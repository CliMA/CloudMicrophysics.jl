"""
    BulkMicrophysicsTendencies

Fused bulk microphysics source terms for atmospheric models.

Provides a dispatch-based API to compute all microphysics tendencies
in a single function call, enabling:
- Simplified integration in atmospheric models
- Point-evaluation suitable for quadrature over subgrid-scale fluctuations
- GPU-friendly pure function design

# Usage

```julia
using CloudMicrophysics.BulkMicrophysicsTendencies

# For 1-moment microphysics
tendencies = bulk_microphysics_tendencies(
    Microphysics1Moment(),
    lcl, icl, rai, sno, ce, aps, vel,
    tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno
)
(; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt) = tendencies
```
"""
module BulkMicrophysicsTendencies

import ..Parameters as CMP
import ..Utilities as UT
import ..Microphysics0M as CM0
import ..Microphysics1M as CM1
import ..Microphysics2M as CM2
import ..MicrophysicsNonEq as CMNonEq
import ..P3Scheme as CMP3
import ..HetIceNucleation as CM_HetIce
import ..AerosolActivation as CMAA
import ...ThermodynamicsInterface as TDI
import ..Common as CO

export MicrophysicsScheme,
    Microphysics0Moment,
    Microphysics1Moment,
    Microphysics2Moment,
    bulk_microphysics_tendencies,
    average_bulk_microphysics_tendencies,
    bulk_microphysics_derivatives,
    bulk_microphysics_cloud_derivatives,
    activation_source,
    repair_ice_state

#####
##### Singleton types for dispatch
#####

"""
    MicrophysicsScheme

Abstract type for microphysics scheme dispatch.
"""
abstract type MicrophysicsScheme end

"""
    Microphysics0Moment <: MicrophysicsScheme

Singleton type for 0-moment microphysics scheme dispatch.
"""
struct Microphysics0Moment <: MicrophysicsScheme end

"""
    Microphysics1Moment <: MicrophysicsScheme

Singleton type for 1-moment microphysics scheme dispatch.
"""
struct Microphysics1Moment <: MicrophysicsScheme end

"""
    Microphysics2Moment <: MicrophysicsScheme

Singleton type for 2-moment microphysics scheme dispatch.

This unified scheme handles:
- Warm rain only (Seifert-Beheng 2006) when ice parameters are not provided
- Warm rain + P3 ice when ice state is provided
"""
struct Microphysics2Moment <: MicrophysicsScheme end

# --- Helper Functions ---

"""
    warm_accretion_melt_factor(tps, sno, T)

Ratio of sensible heat available from warm liquid water to latent heat required
for melting: α = cv_l / L_f × (T - T_freeze).

Returns 0 if T ≤ T_freeze.
"""
@inline function warm_accretion_melt_factor(tps, sno, T)
    L_f = TDI.Lf(tps, T)
    cv_l = TDI.cv_l(tps)
    ΔT = T - sno.T_freeze
    # Branchless: return 0 when cold
    is_cold = (T <= sno.T_freeze)
    return ifelse(is_cold, zero(T), cv_l / L_f * ΔT)
end

# --- Aerosol activation schemes (dispatch on AbstractActivationScheme) ---

"""
    activation_source(scheme::AbstractActivationScheme, tps, ρ, T, q_tot,
                      q_lcl, q_ice, n_lcl, w, p)

Return the cloud-droplet number source `dn_lcl_activation_dt`
[kg⁻¹ air / s] implied by `scheme` at the given local state.

All arguments are passed by value so the function is pure and suitable
for both broadcast and GPU kernels. Inputs that are not needed by a
particular scheme (e.g. `p` for `DiagnosticNc`) are simply ignored.

# Arguments
- `scheme`: activation scheme; dispatch target.
- `tps`: thermodynamics parameters (for supersaturation computations).
- `ρ`: air density [kg/m³].
- `T`: air temperature [K].
- `q_tot`: total water specific content [kg/kg].
- `q_lcl`: cloud liquid specific content [kg/kg].
- `q_ice`: ice specific content [kg/kg] (for supersaturation).
- `n_lcl`: current cloud-droplet number per mass of air [kg⁻¹].
- `w`: vertical velocity [m/s] (0 when unused).
- `p`: air pressure [Pa] (0 when unused).
"""
@inline function activation_source(::CMP.NoActivation, tps, ρ, T, q_tot,
        q_lcl, q_ice, n_lcl, w, p)
    return zero(n_lcl)
end

@inline function activation_source(
        scheme::CMP.DiagnosticNc, tps, ρ, T, q_tot, q_lcl, q_ice, n_lcl, w, p,
)
    FT = typeof(n_lcl)
    target = ifelse(q_lcl > FT(scheme.q_thresh), FT(scheme.N_c), zero(FT))
    return (target - n_lcl) / FT(scheme.τ_relax)
end

@inline function activation_source(
        scheme::CMP.TwomeyActivation, tps, ρ, T, q_tot, q_lcl, q_ice, n_lcl, w, p,
)
    FT = typeof(n_lcl)
    # Twomey: N_CCN = C * S^k whenever w > w_min and cloud is present.
    gate_open = (q_lcl > FT(scheme.q_thresh)) & (w > FT(scheme.w_min))
    S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl, q_ice, ρ, T)
    # Supersaturation can be negative in subsaturated cells; guard before `S^k`
    S_pos = max(zero(FT), S)
    target = ifelse(gate_open, FT(scheme.C) * S_pos^FT(scheme.k), zero(FT))
    return (target - n_lcl) / FT(scheme.τ_relax)
end

@inline function activation_source(
        scheme::CMP.FixedARGActivation, tps, ρ, T, q_tot, q_lcl, q_ice, n_lcl, w, p,
)
    FT = typeof(n_lcl)
    # Early exit for inadmissible ARG inputs (the kernel uses sqrt(α·w/G) and
    # assumes p > 0). A zero source is returned rather than throwing so the
    # BMT call remains pure even in evaporating / descending cells.
    q_gate = q_lcl > FT(scheme.q_thresh)
    if !q_gate || w <= zero(FT) || p <= zero(FT) || !isfinite(T) || !isfinite(p)
        return zero(FT)
    end
    S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl, q_ice, ρ, T)
    if S <= zero(FT)
        return zero(FT)
    end
    # AerosolActivation uses AirProperties from the scheme owner; the closure
    # parameters already carry one. `CMAA.total_N_activated` returns #/m³, so
    # convert to per-mass of air by dividing by ρ before relaxing.
    # The ARG kernel signature expects `N_liq` as a per-volume quantity.
    N_liq = n_lcl * ρ
    # Air properties live on the aerosol-activation parameters only via the
    # separate `air_properties` object; to keep the scheme self-contained we
    # recompute them from `tps` via the existing ARG helper API. The scheme
    # must therefore carry an `air_properties` field in `act_params`'s owner;
    # since `AerosolActivationParameters` already contains the required
    # constants, we pass `tps`'s air block implicitly by calling the
    # 4-argument ARG variant.
    #
    # The 12-arg dispatch in AerosolActivation requires an `AirProperties`
    # parameter. Users that want ARG must therefore provide an `act_params`
    # struct that already embeds it, or call via the high-level helper
    # (not exposed here). For BMT we require the caller to pass an
    # `AerosolActivationParameters` struct that is complete, and fetch
    # `AirProperties` from the scheme's stored `distribution` via the
    # enclosing microphysics parameters.  See `FixedARGActivation` docstring.
    #
    # Pragmatically: we inline just enough of the ARG logic to return N_act.
    # Compute max-supersaturation with a default AirProperties provided by
    # the activation parameters' parent, relayed through a helper method.
    return _fixed_arg_relax(scheme, tps, ρ, T, p, w, q_tot, q_lcl, q_ice, n_lcl)
end

@inline function _fixed_arg_relax(
        scheme::CMP.FixedARGActivation, tps, ρ, T, p, w, q_tot, q_lcl, q_ice, n_lcl,
)
    FT = typeof(n_lcl)
    # Use CMAA.total_N_activated with the scheme's stored `act_params`,
    # `air_properties`, and `distribution`. Any DomainError raised by the
    # kernel (sqrt of negative, etc.) is treated as "kernel cannot apply"
    # → no activation this step.
    try
        N_act = CMAA.total_N_activated(
            scheme.act_params, scheme.distribution,
            scheme.air_properties, tps,
            T, p, w, q_tot, q_lcl, q_ice, n_lcl * ρ, FT(0),
        ) / ρ
        return ifelse(!isfinite(N_act) || N_act <= n_lcl,
            zero(FT),
            (N_act - n_lcl) / FT(scheme.τ_relax))
    catch
        return zero(FT)
    end
end

# --- 1-Moment Microphysics ---

"""
    bulk_microphysics_tendencies(
        ::Microphysics1Moment, mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute all 1-moment microphysics tendencies in one fused call.

Returns a NamedTuple with all source/sink terms for hydrometeor species.
This is a pure function of local thermodynamic state, suitable for:
- Point quadrature over subgrid-scale (T, q_tot) distributions
- GPU kernel evaluation
- Unit testing in isolation

# Arguments
- `mp`: NamedTuple of microphysics parameters with fields:
  - `lcl`: CloudLiquid parameters
  - `icl`: CloudIce parameters
  - `rai`: Rain parameters
  - `sno`: Snow parameters
  - `ce`: CollisionEff parameters
  - `aps`: AirProperties parameters
  - `vel`: Blk1MVelType parameters
  - `var`: VarTimescaleAcnv parameters (optional, for 2M autoconversion)
- `tps`: Thermodynamics parameters
- `ρ`: Air density [kg/m³]
- `T`: Temperature [K]
- `q_tot`: Total water specific content [kg/kg]
- `q_lcl`: Cloud liquid water specific content [kg/kg]
- `q_icl`: Cloud ice specific content [kg/kg]
- `q_rai`: Rain specific content [kg/kg]
- `q_sno`: Snow specific content [kg/kg]
- `N_lcl`: Cloud droplet number concentration [1/m³], optional, default=0
  When N_lcl > 0, uses 2M autoconversion; otherwise uses 1M

# Returns
`NamedTuple` with fields:
- `dq_lcl_dt`: Cloud liquid tendency [kg/kg/s]
- `dq_icl_dt`: Cloud ice tendency [kg/kg/s]
- `dq_rai_dt`: Rain tendency [kg/kg/s]
- `dq_sno_dt`: Snow tendency [kg/kg/s]

# Input Validation
- Negative specific contents are clamped to zero for robustness against numerical errors.

# Notes
- This function does NOT apply timestep-dependent limiters to prevent negative concentrations.
  Limiters depend on timestep `dt` and should be applied by the caller after computing tendencies.
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, N_lcl = zero(ρ),
)
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)

    # Unpack microphysics parameter container
    lcl = mp.cloud.liquid
    icl = mp.cloud.ice
    rai = mp.precip.rain
    sno = mp.precip.snow
    ce = mp.collision
    aps = mp.air_properties
    vel = mp.terminal_velocity
    var = mp.autoconv_2M

    # Initialize tendencies
    dq_lcl_dt = zero(T)
    dq_icl_dt = zero(T)
    dq_rai_dt = zero(T)
    dq_sno_dt = zero(T)

    # --- Cloud condensate formation (non-equilibrium) ---

    # Condensation/evaporation of cloud liquid
    S_lcl_cond = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(lcl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    dq_lcl_dt += S_lcl_cond

    # Deposition/sublimation of cloud ice (INP limiter applied inside conv_q_vap_to_q_lcl_icl_MM2015)
    S_icl_dep = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(icl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    dq_icl_dt += S_icl_dep

    # --- Autoconversion (cloud → precipitation) ---

    # Cloud liquid → rain
    # Use 2M autoconversion when N_lcl > 0 and var provided, otherwise 1M
    S_acnv_1M = CM1.conv_q_lcl_to_q_rai(rai.acnv1M, q_lcl, true)
    S_acnv_2M = isnothing(var) ? S_acnv_1M : CM2.conv_q_lcl_to_q_rai(var, q_lcl, ρ, N_lcl)
    S_acnv_lcl = ifelse(N_lcl > zero(N_lcl), S_acnv_2M, S_acnv_1M)
    dq_lcl_dt -= S_acnv_lcl
    dq_rai_dt += S_acnv_lcl

    # Cloud ice → snow (no supersaturation version for simplicity)
    S_acnv_icl = CM1.conv_q_icl_to_q_sno_no_supersat(sno.acnv1M, q_icl, true)
    dq_icl_dt -= S_acnv_icl
    dq_sno_dt += S_acnv_icl

    # --- Accretion (collisions between species) ---

    # Cloud liquid + rain → rain
    S_accr_lcl_rai = CM1.accretion(lcl, rai, vel.rain, ce, q_lcl, q_rai, ρ)
    dq_lcl_dt -= S_accr_lcl_rai
    dq_rai_dt += S_accr_lcl_rai

    # Cloud liquid + snow → snow (riming, cold) or rain + snow melt (warm)
    # Use branchless ifelse for GPU compatibility
    S_accr_lcl_sno = CM1.accretion(lcl, sno, vel.snow, ce, q_lcl, q_sno, ρ)
    dq_lcl_dt -= S_accr_lcl_sno

    is_warm = T >= sno.T_freeze
    # Cold: riming → snow; Warm: shedding → rain
    dq_sno_dt += ifelse(is_warm, zero(T), S_accr_lcl_sno)
    dq_rai_dt += ifelse(is_warm, S_accr_lcl_sno, zero(T))
    # Thermal melting (α=0 when cold, so these terms vanish)
    # Physical consistency: melt is mass of snow melted, which is α * mass of warm liquid
    α = warm_accretion_melt_factor(tps, sno, T)
    S_accr_melt = α * S_accr_lcl_sno
    dq_sno_dt -= S_accr_melt
    dq_rai_dt += S_accr_melt

    # Cloud ice + rain → snow
    S_accr_icl_rai = CM1.accretion(icl, rai, vel.rain, ce, q_icl, q_rai, ρ)
    dq_icl_dt -= S_accr_icl_rai
    dq_sno_dt += S_accr_icl_rai

    # Cloud ice + snow → snow
    S_accr_icl_sno = CM1.accretion(icl, sno, vel.snow, ce, q_icl, q_sno, ρ)
    dq_icl_dt -= S_accr_icl_sno
    dq_sno_dt += S_accr_icl_sno

    # Rain + cloud ice → sink of rain (forms snow)
    S_accr_rai_icl = CM1.accretion_rain_sink(rai, icl, vel.rain, ce, q_icl, q_rai, ρ)
    dq_rai_dt -= S_accr_rai_icl
    dq_sno_dt += S_accr_rai_icl

    # Rain + snow collisions (temperature dependent)
    # Use branchless ifelse for GPU compatibility
    # Cold: rain + snow → snow (rain freezes); Warm: snow + rain → rain (snow melts)
    # Note: accretion_snow_rain arguments differ by case (collector vs collected swap)
    S_accr_rai_sno = CM1.accretion_snow_rain(sno, rai, vel.snow, vel.rain, ce, q_sno, q_rai, ρ)
    S_accr_sno_rai = CM1.accretion_snow_rain(rai, sno, vel.rain, vel.snow, ce, q_rai, q_sno, ρ)

    is_warm = T >= sno.T_freeze
    # Cold pathway: rain → snow
    dq_rai_dt -= ifelse(is_warm, zero(T), S_accr_rai_sno)
    dq_sno_dt += ifelse(is_warm, zero(T), S_accr_rai_sno)
    # Warm pathway: snow → rain
    dq_sno_dt -= ifelse(is_warm, S_accr_sno_rai, zero(T))
    dq_rai_dt += ifelse(is_warm, S_accr_sno_rai, zero(T))
    # Thermal melting from warm rain (α=0 when cold)
    # Physical consistency: melt is mass of snow melted, which is α * mass of warm rain
    S_accr_melt = α * S_accr_rai_sno
    dq_sno_dt -= ifelse(is_warm, S_accr_melt, zero(T))
    dq_rai_dt += ifelse(is_warm, S_accr_melt, zero(T))

    # --- Evaporation and sublimation ---

    # Rain evaporation
    S_evap_rai = CM1.evaporation_sublimation(rai, vel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    dq_rai_dt += S_evap_rai  # negative tendency (evaporation)

    # Snow sublimation/deposition
    S_subl_sno = CM1.evaporation_sublimation(sno, vel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    dq_sno_dt += S_subl_sno  # can be positive (deposition) or negative (sublimation)

    # --- Melting ---

    # Snow melt → rain
    S_melt_sno = CM1.snow_melt(sno, vel.snow, aps, tps, q_sno, ρ, T)
    dq_sno_dt -= S_melt_sno
    dq_rai_dt += S_melt_sno

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end

"""
Construct a local linear approximation of 1-moment microphysics tendencies:

    dq/dt ≈ M * q + e

using a donor-based linearization:
- donor → receiver transfers are represented as `D * q_donor`
- vapor → condensate sources are treated as constants (`e`)
- condensate sinks are treated as linear sinks (`-D * q`)

All coefficients use `D = S / max(q_min, q_donor)` for robustness.

With this formulation, sink terms behave like exponential decays over a timestep
in the implicit solve. The resulting operator is sparse and only stores the
nonzero entries needed by the specialized solver.

Returns a `NamedTuple` containing the nonzero entries of `M` and `e`.
"""
@inline function _bulk_microphysics_linearized_operator(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, N_lcl = zero(ρ),
)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)

    lcl = mp.cloud.liquid
    icl = mp.cloud.ice
    rai = mp.precip.rain
    sno = mp.precip.snow
    ce = mp.collision
    aps = mp.air_properties
    vel = mp.terminal_velocity
    var = mp.autoconv_2M

    FT = promote_type(
        typeof(ρ), typeof(T), typeof(q_tot), typeof(q_lcl),
        typeof(q_icl), typeof(q_rai), typeof(q_sno),
    )

    M11 = zero(FT)
    M22 = zero(FT)
    M31 = zero(FT)
    M33 = zero(FT)
    M34 = zero(FT)
    M41 = zero(FT)
    M42 = zero(FT)
    M43 = zero(FT)
    M44 = zero(FT)
    e1 = zero(FT)
    e2 = zero(FT)
    e4 = zero(FT)

    # minimum specific humidity threshold for regularizing S/q
    q_min = TDI.TD.Parameters.q_min(tps)

    # ------------------------------------------------------------
    # Vapor -> cloud condensate: constant sources
    # ------------------------------------------------------------
    S_lcl_cond = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(
        lcl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
    )
    D = S_lcl_cond / max(q_min, q_lcl)
    is_source = S_lcl_cond >= zero(FT)
    e1 += ifelse(is_source, S_lcl_cond, zero(FT))
    M11 += ifelse(is_source, zero(FT), D)

    S_icl_dep = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(
        icl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
    )
    D = S_icl_dep / max(q_min, q_icl)
    is_source = S_icl_dep >= zero(FT)
    e2 += ifelse(is_source, S_icl_dep, zero(FT))
    M22 += ifelse(is_source, zero(FT), D)

    # ------------------------------------------------------------
    # Autoconversion: donor-based transfer
    # ------------------------------------------------------------
    S_acnv_1M = CM1.conv_q_lcl_to_q_rai(rai.acnv1M, q_lcl, true)
    S_acnv_2M = isnothing(var) ? S_acnv_1M : CM2.conv_q_lcl_to_q_rai(var, q_lcl, ρ, N_lcl)
    S_acnv_lcl = ifelse(N_lcl > zero(N_lcl), S_acnv_2M, S_acnv_1M)
    D = S_acnv_lcl / max(q_min, q_lcl)
    M11 -= D
    M31 += D

    S_acnv_icl = CM1.conv_q_icl_to_q_sno_no_supersat(sno.acnv1M, q_icl, true)
    D = S_acnv_icl / max(q_min, q_icl)
    M22 -= D
    M42 += D

    # ------------------------------------------------------------
    # Accretion: donor-based transfer
    # ------------------------------------------------------------
    S_accr_lcl_rai = CM1.accretion(lcl, rai, vel.rain, ce, q_lcl, q_rai, ρ)
    D = S_accr_lcl_rai / max(q_min, q_lcl)
    M11 -= D
    M31 += D

    S_accr_lcl_sno = CM1.accretion(lcl, sno, vel.snow, ce, q_lcl, q_sno, ρ)
    D = S_accr_lcl_sno / max(q_min, q_lcl)
    M11 -= D
    is_warm = T >= sno.T_freeze
    M31 += ifelse(is_warm, D, zero(FT))
    M41 += ifelse(is_warm, zero(FT), D)

    α = warm_accretion_melt_factor(tps, sno, T)
    S_accr_melt = α * S_accr_lcl_sno
    # this is snow -> rain melt induced by warm collected liquid
    # donor is snow
    D = S_accr_melt / max(q_min, q_sno)
    M44 -= D
    M34 += D

    S_accr_icl_rai = CM1.accretion(icl, rai, vel.rain, ce, q_icl, q_rai, ρ)
    D = S_accr_icl_rai / max(q_min, q_icl)
    M22 -= D
    M42 += D

    S_accr_icl_sno = CM1.accretion(icl, sno, vel.snow, ce, q_icl, q_sno, ρ)
    D = S_accr_icl_sno / max(q_min, q_icl)
    M22 -= D
    M42 += D

    S_accr_rai_icl = CM1.accretion_rain_sink(rai, icl, vel.rain, ce, q_icl, q_rai, ρ)
    D = S_accr_rai_icl / max(q_min, q_rai)
    M33 -= D
    M43 += D

    S_accr_rai_sno = CM1.accretion_snow_rain(sno, rai, vel.snow, vel.rain, ce, q_sno, q_rai, ρ)
    S_accr_sno_rai = CM1.accretion_snow_rain(rai, sno, vel.rain, vel.snow, ce, q_rai, q_sno, ρ)

    # snow -> rain
    D = S_accr_sno_rai / max(q_min, q_sno)
    M44 -= ifelse(is_warm, D, zero(FT))
    M34 += ifelse(is_warm, D, zero(FT))

    S_accr_melt2 = α * S_accr_rai_sno
    D = S_accr_melt2 / max(q_min, q_sno)
    M44 -= ifelse(is_warm, D, zero(FT))
    M34 += ifelse(is_warm, D, zero(FT))

    # rain -> snow
    D = S_accr_rai_sno / max(q_min, q_rai)
    M33 -= ifelse(is_warm, zero(FT), D)
    M43 += ifelse(is_warm, zero(FT), D)

    # ------------------------------------------------------------
    # Rain evaporation: sink to vapor (always zero or negative)
    # ------------------------------------------------------------
    S_evap_rai = CM1.evaporation_sublimation(
        rai, vel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
    )
    D = (-S_evap_rai) / max(q_min, q_rai)
    M33 -= D

    # ------------------------------------------------------------
    # Snow sublimation/deposition
    # ------------------------------------------------------------
    S_subl_sno = CM1.evaporation_sublimation(
        sno, vel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T,
    )
    D = S_subl_sno / max(q_min, q_sno)
    is_source = S_subl_sno >= zero(FT)
    e4 += ifelse(is_source, S_subl_sno, zero(FT))
    M44 += ifelse(is_source, zero(FT), D)

    # ------------------------------------------------------------
    # Snow melt: snow -> rain
    # ------------------------------------------------------------
    S_melt_sno = CM1.snow_melt(sno, vel.snow, aps, tps, q_sno, ρ, T)
    D = S_melt_sno / max(q_min, q_sno)
    M44 -= D
    M34 += D

    return (
        M11 = M11, M22 = M22,
        M31 = M31, M33 = M33, M34 = M34,
        M41 = M41, M42 = M42, M43 = M43, M44 = M44,
        e1 = e1, e2 = e2, e4 = e4,
    )
end

"""
Compute time-averaged 1-moment microphysics tendencies over a single linearized substep.

Solves the linearized implicit system

    (q* - q⁰) / Δt = M q* + e

and returns the average tendency

    dq/dt = (q* - q⁰) / Δt.

The system uses a sparse structure specific to the 1-moment microphysics model:
`q_lcl` and `q_icl` are solved independently, while `q_rai` and `q_sno` are
solved from a coupled 2×2 system.

Because sinks are linearized as `-D q`, they are effectively integrated as
exponential decays over the substep.
"""
@inline function _average_bulk_microphysics_tendencies(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, N_lcl = zero(ρ),
)

    FT = typeof(q_tot)

    lin = _bulk_microphysics_linearized_operator(
        Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, N_lcl,
    )

    invΔt = one(FT) / Δt

    # A = I/Δt - M
    a11 = invΔt - lin.M11
    a22 = invΔt - lin.M22
    a31 = -lin.M31
    a33 = invΔt - lin.M33
    a34 = -lin.M34
    a41 = -lin.M41
    a42 = -lin.M42
    a43 = -lin.M43
    a44 = invΔt - lin.M44

    # rhs = e + q_0/Δt
    # e3 = 0 by the 1m model
    b1 = lin.e1 + invΔt * q_lcl
    b2 = lin.e2 + invΔt * q_icl
    b3 = invΔt * q_rai
    b4 = lin.e4 + invΔt * q_sno

    q_lcl_new = b1 / a11
    q_icl_new = b2 / a22

    # Reduced 2x2 system for q_rai_new, q_sno_new
    r3 = muladd(-a31, q_lcl_new, b3)
    r4 = muladd(-a41, q_lcl_new, muladd(-a42, q_icl_new, b4))

    det = muladd(-a34, a43, a33 * a44)
    # det is a positive number because a44 and a33 are positive (greater than invΔt) 
    # and a34 and a43 are non-positive so we don't need to safeguard division by det.
    q_rai_new = (r3 * a44 - a34 * r4) / det
    q_sno_new = (a33 * r4 - r3 * a43) / det

    dq_lcl_dt = (q_lcl_new - q_lcl) * invΔt
    dq_icl_dt = (q_icl_new - q_icl) * invΔt
    dq_rai_dt = (q_rai_new - q_rai) * invΔt
    dq_sno_dt = (q_sno_new - q_sno) * invΔt

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end

"""
Compute average 1-moment microphysics tendencies over `Δt` using repeated
linearized implicit substeps.

The interval `Δt` is divided into `nsub` equal substeps. At each substep, a local
linearized microphysics system is rebuilt from the current state and solved
implicitly for cloud liquid, cloud ice, rain, and snow. Temperature is then
updated from the latent heating implied by the substep tendencies.

The returned tendencies are the net change in the hydrometeor species over the
full interval divided by `Δt`.

Increasing `nsub` improves how well the method captures nonlinear changes in the
active microphysical processes, including regime changes near freezing.
"""
@inline function average_bulk_microphysics_tendencies(
    cm::Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    Δt,
    nsub = 1,
    N_lcl = zero(ρ),
)
    FT = typeof(q_tot)

    q_lcl_0 = q_lcl
    q_icl_0 = q_icl
    q_rai_0 = q_rai
    q_sno_0 = q_sno

    Δt_sub = Δt / FT(nsub)

    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    for _ in 1:nsub
        rates = _average_bulk_microphysics_tendencies(
            cm,
            mp,
            tps,
            ρ,
            T,
            q_tot,
            q_lcl,
            q_icl,
            q_rai,
            q_sno,
            Δt_sub,
            N_lcl,
        )

        q_lcl += rates.dq_lcl_dt * Δt_sub
        q_icl += rates.dq_icl_dt * Δt_sub
        q_rai += rates.dq_rai_dt * Δt_sub
        q_sno += rates.dq_sno_dt * Δt_sub

        T +=
            (
                Lv_over_cp * (rates.dq_lcl_dt + rates.dq_rai_dt) +
                Ls_over_cp * (rates.dq_icl_dt + rates.dq_sno_dt)
            ) * Δt_sub
    end

    dq_lcl_dt = (q_lcl - q_lcl_0) / Δt
    dq_icl_dt = (q_icl - q_icl_0) / Δt
    dq_rai_dt = (q_rai - q_rai_0) / Δt
    dq_sno_dt = (q_sno - q_sno_0) / Δt

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end

"""
    bulk_microphysics_derivatives(
        ::Microphysics1Moment,
        mp,
        tps,
        ρ,
        T,
        q_tot,
        q_lcl,
        q_icl,
        q_rai,
        q_sno,
    )

Compute all 1-moment microphysics tendency derivatives in one fused call.

Returns a NamedTuple with the leading-order derivative of each species tendency w.r.t.
its own specific content. Uses the same input clamping and clipping logic
as `bulk_microphysics_tendencies`.

# Returns
`NamedTuple` with fields:
- `∂tendency_∂q_lcl`: Derivative of cloud liquid tendency w.r.t. q_lcl [1/s]
- `∂tendency_∂q_icl`: Derivative of cloud ice tendency w.r.t. q_icl [1/s]
- `∂tendency_∂q_rai`: Derivative of rain tendency w.r.t. q_rai [1/s]
- `∂tendency_∂q_sno`: Derivative of snow tendency w.r.t. q_sno [1/s]
"""
@inline function bulk_microphysics_derivatives(
    ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    N_lcl = zero(ρ),
)
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)

    # Unpack microphysics parameter container
    lcl = mp.cloud.liquid
    icl = mp.cloud.ice
    rai = mp.precip.rain
    sno = mp.precip.snow
    ce = mp.collision
    aps = mp.air_properties
    vel = mp.terminal_velocity

    # --- Cloud condensate derivatives (reuse dedicated function) ---
    (; ∂tendency_∂q_lcl, ∂tendency_∂q_icl) = bulk_microphysics_cloud_derivatives(
        Microphysics1Moment(), mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

    # --- Rain: evaporation + sink self-derivatives ---
    ∂tendency_∂q_rai =
        CM1.∂evaporation_sublimation_∂q_precip(rai, vel.rain, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    # Accretion rain_sink: ice + rain → snow (rate/q_rai approximation)
    S_accr_rai_icl = CM1.accretion_rain_sink(rai, icl, vel.rain, ce, q_icl, q_rai, ρ)
    # Snow-rain accretion (cold): rain → snow
    is_warm = T >= sno.T_freeze
    S_accr_rai_sno = ifelse(is_warm, zero(T),
        CM1.accretion_snow_rain(sno, rai, vel.snow, vel.rain, ce, q_sno, q_rai, ρ))
    ∂tendency_∂q_rai -=
        (S_accr_rai_icl + S_accr_rai_sno) /
        max(q_rai, UT.ϵ_numerics(typeof(q_rai)))

    # --- Snow: sublimation/deposition + melting + sink self-derivatives ---
    dS_subl_sno =
        CM1.∂evaporation_sublimation_∂q_precip(sno, vel.snow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    dS_melt_sno = CM1.∂snow_melt_∂q_sno(sno, vel.snow, aps, tps, q_sno, ρ, T)
    ∂tendency_∂q_sno = dS_subl_sno - dS_melt_sno
    # Snow-rain accretion (warm): snow → rain
    S_accr_sno_rai = ifelse(is_warm,
        CM1.accretion_snow_rain(rai, sno, vel.rain, vel.snow, ce, q_rai, q_sno, ρ), zero(T))
    # Thermal melting from warm rain collision -> snow sink proportional to S_accr_rai_sno
    # (Since this is a diagonal derivative, we approximate ∂/∂q_sno by rate/q_sno, using the rain sink form S_accr_rai_sno)
    S_accr_sno_rai_melt = ifelse(is_warm,
        CM1.accretion_snow_rain(sno, rai, vel.snow, vel.rain, ce, q_sno, q_rai, ρ), zero(T))
    # Accretion melt from warm liquid-snow collision
    α = warm_accretion_melt_factor(tps, sno, T)
    S_accr_lcl_sno_melt = α * CM1.accretion(lcl, sno, vel.snow, ce, q_lcl, q_sno, ρ)
    ∂tendency_∂q_sno -=
        (S_accr_sno_rai + α * S_accr_sno_rai_melt + S_accr_lcl_sno_melt) /
        max(q_sno, UT.ϵ_numerics(typeof(q_sno)))

    return (; ∂tendency_∂q_lcl, ∂tendency_∂q_icl, ∂tendency_∂q_rai, ∂tendency_∂q_sno)
end

"""
    bulk_microphysics_cloud_derivatives(
        ::Microphysics1Moment,
        mp,
        tps,
        ρ,
        T,
        q_tot,
        q_lcl,
        q_icl,
        q_rai,
        q_sno,
    )

Compute 1-moment cloud condensate Jacobian derivatives only.

Returns a NamedTuple with the leading-order derivative of each cloud species
tendency w.r.t. its own specific content. Precipitation derivatives are omitted
(the Jacobian uses quadrature-consistent S/q for precipitation instead).

This is cheaper than `bulk_microphysics_derivatives` because it skips rain
evaporation, snow sublimation/melting, and precipitation accretion derivatives
and avoids quadratures over the derivatives.

# Returns
`NamedTuple` with fields:
- `∂tendency_∂q_lcl`: Derivative of cloud liquid tendency w.r.t. q_lcl [1/s]
- `∂tendency_∂q_icl`: Derivative of cloud ice tendency w.r.t. q_icl [1/s]
"""
@inline function bulk_microphysics_cloud_derivatives(
    ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
)
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)

    # Unpack microphysics parameter container
    lcl = mp.cloud.liquid
    icl = mp.cloud.ice
    rai = mp.precip.rain
    sno = mp.precip.snow
    ce = mp.collision
    vel = mp.terminal_velocity

    # --- Cloud liquid: condensation + sink self-derivatives ---
    ∂tendency_∂q_lcl = -1 / CMNonEq.τ_relax(lcl)
    # Autoconversion q_lcl → q_rai (sink)
    S_acnv_lcl = CM1.conv_q_lcl_to_q_rai(rai.acnv1M, q_lcl, true)
    # Accretion q_lcl + q_rai → q_rai (rate is exactly linear in q_lcl)
    S_accr_lcl_rai = CM1.accretion(lcl, rai, vel.rain, ce, q_lcl, q_rai, ρ)
    # Accretion q_lcl + q_sno → q_sno (rate is exactly linear in q_lcl)
    S_accr_lcl_sno = CM1.accretion(lcl, sno, vel.snow, ce, q_lcl, q_sno, ρ)
    ∂tendency_∂q_lcl -=
        (S_acnv_lcl + S_accr_lcl_rai + S_accr_lcl_sno) /
        max(q_lcl, UT.ϵ_numerics(typeof(q_lcl)))

    # --- Cloud ice: deposition + sink self-derivatives ---
    ∂tendency_∂q_icl = -1 / CMNonEq.τ_relax(icl)
    # Autoconversion q_icl → q_sno
    S_acnv_icl = CM1.conv_q_icl_to_q_sno_no_supersat(sno.acnv1M, q_icl, true)
    # Accretion q_icl + q_rai → q_sno (rate is exactly linear in q_icl)
    S_accr_icl_rai = CM1.accretion(icl, rai, vel.rain, ce, q_icl, q_rai, ρ)
    # Accretion q_icl + q_sno → q_sno (rate is exactly linear in q_icl)
    S_accr_icl_sno = CM1.accretion(icl, sno, vel.snow, ce, q_icl, q_sno, ρ)
    ∂tendency_∂q_icl -=
        (S_acnv_icl + S_accr_icl_rai + S_accr_icl_sno) /
        max(q_icl, UT.ϵ_numerics(typeof(q_icl)))

    return (; ∂tendency_∂q_lcl, ∂tendency_∂q_icl)
end

# --- 2-Moment Microphysics derivatives ---

"""
    bulk_microphysics_derivatives(
        ::Microphysics2Moment,
        mp,
        tps,
        ρ,
        T,
        q_tot,
        q_lcl,
        q_icl,
        q_rai,
        q_sno,
        n_lcl,
        n_rai,
    )

Compute 2-moment microphysics tendency derivatives in one fused call.

Returns leading-order derivatives of species tendencies w.r.t. their own
specific content (q) and number (n). Rain uses the 2-moment rain evaporation
derivatives; snow and cloud formation derivatives are zero for now.

# Arguments
- `mp`: Microphysics2MParams (warm rain; ice optional)
- `n_lcl`, `n_rai`: Cloud and rain number per kg air (1/kg); N_rai = ρ * n_rai used for rain evaporation derivative

# Returns
`NamedTuple` with fields: `∂tendency_∂q_lcl`, `∂tendency_∂q_icl`, `∂tendency_∂q_rai`, `∂tendency_∂q_sno`, `∂tendency_∂n_lcl`, `∂tendency_∂n_rai`
"""
@inline function bulk_microphysics_derivatives(
    ::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE},
    tps,
    ρ,
    T,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    n_lcl,
    n_rai,
) where {WR, ICE}
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)
    n_lcl = UT.clamp_to_nonneg(n_lcl)
    n_rai = UT.clamp_to_nonneg(n_rai)

    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties
    N_rai = ρ * n_rai

    # TODO: Cloud formation — 2M bulk_microphysics_tendencies does not call cloud formation yet;
    # once it does, set these from MM2015 (same as 1M) via -1/τ_relax.
    ∂tendency_∂q_lcl = zero(ρ)
    ∂tendency_∂q_icl = zero(ρ)

    # 2-moment rain evaporation derivatives (return order: ∂n_rai, ∂q_rai, matching rain_evaporation)
    (; ∂N_rai, ∂q_rai) =
        CM2.∂rain_evaporation_∂N_rai_∂q_rai(sb, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)
    ∂tendency_∂q_rai = ∂q_rai
    ∂tendency_∂n_rai = ∂N_rai / ρ

    # TODO: 2M scheme also applies number/mass adjustment terms (e.g. number_increase_for_mass_limit,
    # number_decrease_for_mass_limit) to keep the assumed size distribution stable; consider adding
    # their derivatives here later.

    # Snow derivative zero for now
    ∂tendency_∂q_sno = zero(ρ)

    # Cloud number derivative zero for now
    ∂tendency_∂n_lcl = zero(ρ)

    return (;
        ∂tendency_∂q_lcl,
        ∂tendency_∂q_icl,
        ∂tendency_∂q_rai,
        ∂tendency_∂q_sno,
        ∂tendency_∂n_lcl,
        ∂tendency_∂n_rai,
    )
end

# --- 0-Moment Microphysics ---
"""
    bulk_microphysics_tendencies(::Microphysics0Moment, mp, tps, T, q_lcl, q_icl)
    bulk_microphysics_tendencies(::Microphysics0Moment, mp, tps, T, q_lcl, q_icl, q_vap_sat)

Compute 0-moment microphysics tendencies in one fused call.

Returns the total water tendency `dq_tot_dt` (a scalar, in kg/kg/s) from precipitation removal.

The first form uses the fixed condensate threshold `qc_0`;
the second form uses the supersaturation threshold `S_0 * q_vap_sat`.

# Arguments
- `mp`: Microphysics0MParams (contains τ_precip, qc_0, S_0)
- `tps`: Thermodynamics parameters
- `T`: Temperature [K]
- `q_lcl`: Cloud liquid specific content [kg/kg]
- `q_icl`: Cloud ice specific content [kg/kg]
- `q_vap_sat`: (second method only) Saturation specific humidity [kg/kg]

# Notes
- Does NOT apply limiters (caller applies based on timestep)
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics0Moment, mp::CMP.Microphysics0MParams, tps,
    T, q_lcl, q_icl,
)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    dq_tot_dt = CM0.remove_precipitation(mp.precip, q_lcl, q_icl)
    return dq_tot_dt
end
@inline function bulk_microphysics_tendencies(
    ::Microphysics0Moment,
    mp::CMP.Microphysics0MParams,
    tps,
    T,
    q_lcl,
    q_icl,
    q_vap_sat,
)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    dq_tot_dt = CM0.remove_precipitation(mp.precip, q_lcl, q_icl, q_vap_sat)
    return dq_tot_dt
end

# --- 2-Moment Microphysics Helper Functions ---

"""
    warm_rain_tendencies_2m(sb, q_lcl, q_rai, ρ, n_lcl, n_rai)

Internal helper function that computes core 2M warm rain processes:
autoconversion, self-collection, accretion, and rain breakup.

Used by both warm-only and warm+ice dispatch methods to reduce code duplication.

# Arguments
- `sb`: SB2006 parameters
- `q_lcl`: Cloud liquid specific content (kg/kg)
- `q_rai`: Rain specific content (kg/kg)
- `ρ`: Air density (kg/m³)
- `n_lcl`: Cloud droplet number per kg air (1/kg)
- `n_rai`: Rain number per kg air (1/kg)

# Returns
`NamedTuple` with core warm rain tendencies:
- `dq_lcl_dt`: Cloud liquid mass tendency (kg/kg/s)
- `dq_rai_dt`: Rain mass tendency (kg/kg/s)
- `dn_lcl_dt`: Cloud number tendency (1/kg/s)
- `dn_rai_dt`: Rain number tendency (1/kg/s)
"""
@inline function warm_rain_tendencies_2m(
    warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai,
    w = zero(ρ), p = zero(ρ),
)

    # Unpack parameters
    sb = warm_rain.seifert_beheng
    aps = warm_rain.air_properties
    condevap = warm_rain.condevap

    # Convert to number densities for CM2 functions
    N_lcl = ρ * n_lcl
    N_rai = ρ * n_rai

    # Initialize tendencies
    FT = typeof(ρ)
    dq_lcl_dt = zero(FT)
    dq_rai_dt = zero(FT)
    dn_lcl_dt = zero(FT)
    dn_rai_dt = zero(FT)

    # --- Aerosol activation (scheme is in `warm_rain.activation_scheme`) ---
    # Naturally belongs here: activation is a warm-rain droplet source. The
    # scheme produces only a number tendency; mass is sourced from MM2015
    # condensation in the next block (cloud droplet starter mass is small
    # enough that the rate is dominated by the depositional growth term).
    # `w`, `p` are the per-cell ambient state inputs that w/p-dependent
    # schemes (Twomey, FixedARG) need; positional so the call site can use
    # `@.` broadcast over ClimaCore Fields. Schemes that don't need them
    # (DiagnosticNc, NoActivation) ignore the values, so the `zero(ρ)`
    # defaults are safe for the diagnostic-Nc / no-activation cases.
    dn_lcl_activation_dt = activation_source(
        warm_rain.activation_scheme, tps, ρ, T, q_tot, q_lcl, q_ice, n_lcl, w, p,
    )
    dn_lcl_dt += dn_lcl_activation_dt

    # --- Condensation of vapor / evaporation of cloud liquid water ---
    ∂ₜq_lcl_cond = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(condevap, tps, q_tot, q_lcl, q_ice, q_rai, zero(q_ice), ρ, T)
    dq_lcl_dt += ∂ₜq_lcl_cond
    # dn_lcl_dt += zero(∂ₜq_lcl_cond)  # neglect number change from condensation/evaporation

    # --- Evaporation of rain ---
    evap = CM2.rain_evaporation(sb, aps, tps, q_tot, q_lcl, q_ice, q_rai, zero(q_ice), ρ, N_rai, T)
    dq_rai_dt += evap.∂ₜq_rai
    dn_rai_dt += evap.∂ₜρn_rai / ρ

    # --- Autoconversion ---
    acnv = CM2.autoconversion(sb.acnv, sb.pdf_c, q_lcl, q_rai, ρ, N_lcl)
    dq_lcl_dt += acnv.dq_lcl_dt
    dq_rai_dt += acnv.dq_rai_dt
    dn_lcl_dt += acnv.dN_lcl_dt / ρ
    dn_rai_dt += acnv.dN_rai_dt / ρ

    # --- Cloud liquid self-collection ---
    ∂ₜN_lcl_sc = CM2.cloud_liquid_self_collection(sb.acnv, sb.pdf_c, q_lcl, ρ, acnv.dN_lcl_dt)
    dn_lcl_dt += ∂ₜN_lcl_sc / ρ

    # --- Accretion ---
    accr = CM2.accretion(sb, q_lcl, q_rai, ρ, N_lcl)
    dq_lcl_dt += accr.dq_lcl_dt
    dq_rai_dt += accr.dq_rai_dt
    dn_lcl_dt += accr.dN_lcl_dt / ρ

    # --- Rain self-collection ---
    ∂ₜN_rai_sc = CM2.rain_self_collection(sb.pdf_r, sb.self, q_rai, ρ, N_rai)
    dn_rai_dt += ∂ₜN_rai_sc / ρ

    # --- Rain breakup ---
    ∂ₜN_rai_br = CM2.rain_breakup(sb.pdf_r, sb.brek, q_rai, ρ, N_rai, ∂ₜN_rai_sc)
    dn_rai_dt += ∂ₜN_rai_br / ρ

    # --- Number adjustment for mass limits ---
    # Cloud liquid
    numadj_lcl = (; sb.numadj.τ, x_min = sb.pdf_c.xc_min, x_max = sb.pdf_c.xc_max)
    ∂ₜn_lcl_numadj = CM2.number_tendency_from_mass_limits(numadj_lcl, q_lcl, n_lcl)
    dn_lcl_dt += ∂ₜn_lcl_numadj
    # Rain
    numadj_rai = (; sb.numadj.τ, x_min = sb.pdf_r.xr_min, x_max = sb.pdf_r.xr_max)
    ∂ₜn_rai_numadj = CM2.number_tendency_from_mass_limits(numadj_rai, q_rai, n_rai)
    dn_rai_dt += ∂ₜn_rai_numadj

    return (; dq_lcl_dt, dq_rai_dt, dn_lcl_dt, dn_rai_dt, dn_lcl_activation_dt)
end

# --- 2-Moment Microphysics (Unified Warm + Optional Ice) ---


"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{WR, Nothing},
        ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
    )

Compute 2-moment **warm rain only** microphysics tendencies (Seifert-Beheng 2006).

This method is type-stable and GPU-optimized for warm rain processes only.
For warm rain + P3 ice, see the method that accepts `Microphysics2MParams{FT, WR, <:P3IceParams}`.

# Arguments
- `mp`: Microphysics2MParams with `mp.ice == nothing` (warm rain only)
- `tps`: Thermodynamics parameters
- `ρ`: Air density (kg/m³)
- `T`: Temperature (K)
- `q_tot`: Total water specific content (kg/kg)
- `q_lcl`: Cloud liquid specific content (kg/kg)
- `n_lcl`: Cloud droplet number per kg air (1/kg)
- `q_rai`: Rain specific content (kg/kg)
- `n_rai`: Rain number per kg air (1/kg)

# Returns
`NamedTuple` with warm rain tendency fields:
- `dq_lcl_dt`: Cloud liquid tendency (kg/kg/s)
- `dn_lcl_dt`: Cloud number tendency (1/kg/s)
- `dq_rai_dt`: Rain tendency (kg/kg/s)
- `dn_rai_dt`: Rain number tendency (1/kg/s)
- `dq_ice_dt`: Ice tendency (always zero for warm-only)
- `dq_rim_dt`: Rime mass tendency (always zero for warm-only)
- `db_rim_dt`: Rime volume tendency (always zero for warm-only)
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics2Moment, mp::CMP.Microphysics2MParams{WR, Nothing}, tps,
    ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
    q_ice = zero(ρ), n_ice = zero(ρ), q_rim = zero(ρ), b_rim = zero(ρ), logλ = zero(ρ),
    inpc_log_shift = zero(ρ),
    w = zero(ρ), p = zero(ρ),
) where {WR}
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    n_lcl = UT.clamp_to_nonneg(n_lcl)
    n_rai = UT.clamp_to_nonneg(n_rai)
    q_ice = UT.clamp_to_nonneg(q_ice)
    n_ice = UT.clamp_to_nonneg(n_ice)
    q_rim = UT.clamp_to_nonneg(q_rim)
    b_rim = UT.clamp_to_nonneg(b_rim)

    # Initialize ice-related tendencies (always zero for warm-only)
    dq_ice_dt = zero(ρ)
    dq_rim_dt = zero(ρ)
    db_rim_dt = zero(ρ)

    # --- Core Warm Rain Processes (shared helper, includes activation) ---
    warm = warm_rain_tendencies_2m(mp.warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai,
                                   w, p)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt
    dn_lcl_activation_dt = warm.dn_lcl_activation_dt

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt,
              dq_ice_dt, dq_rim_dt, db_rim_dt,
              dn_lcl_activation_dt)
end

"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{FT, WR, <:P3IceParams},
        ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
        q_ice, n_ice, q_rim, b_rim, logλ,
    )

Compute 2-moment **warm rain + P3 ice** microphysics tendencies.

This method is type-stable and GPU-optimized. The P3 ice parameters are guaranteed
to be non-Nothing, eliminating runtime type checks and dynamic dispatch.

# Arguments
## Required
- `mp`: Microphysics2MParams with P3 ice parameters present
- `tps`: Thermodynamics parameters
- `ρ`: Air density (kg/m³)
- `T`: Temperature (K)
- `q_tot`: Total water specific content (kg/kg)
- `q_lcl`: Cloud liquid specific content (kg/kg)
- `n_lcl`: Cloud droplet number per kg air (1/kg)
- `q_rai`: Rain specific content (kg/kg)
- `n_rai`: Rain number per kg air (1/kg)

## Optional (P3 ice state)
- `q_ice`: Ice specific content (kg/kg), default = 0
- `n_ice`: Ice number per kg air (1/kg), default = 0
- `q_rim`: Rime mass (kg/kg), default = 0
- `b_rim`: Rime volume (m³/kg), default = 0
- `logλ`: Log of P3 distribution slope parameter, log(1/m), default = 0

# Returns
`NamedTuple` with all tendency fields:
- `dq_lcl_dt`: Cloud liquid tendency (kg/kg/s)
- `dn_lcl_dt`: Cloud number tendency (1/kg/s)
- `dq_rai_dt`: Rain tendency (kg/kg/s)
- `dn_rai_dt`: Rain number tendency (1/kg/s)
- `dq_ice_dt`: Ice tendency (kg/kg/s)
- `dq_rim_dt`: Rime mass tendency (kg/kg/s)
- `db_rim_dt`: Rime volume tendency (m³/kg/s)
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics2Moment, mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai,
    q_ice, n_ice, q_rim, b_rim, logλ,
    inpc_log_shift = zero(ρ),
    w = zero(ρ), p = zero(ρ),
    n_INP_used = zero(ρ),
) where {WR, ICE <: CMP.P3IceParams}
    FT = eltype(ρ)
    ϵₘ = UT.ϵ_numerics_2M_M(FT)
    ϵₙ = UT.ϵ_numerics_2M_N(FT)
    ϵB = UT.ϵ_numerics_P3_B(FT)
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    n_lcl = UT.clamp_to_nonneg(n_lcl)
    n_rai = UT.clamp_to_nonneg(n_rai)
    q_ice = UT.clamp_to_nonneg(q_ice)
    n_ice = UT.clamp_to_nonneg(n_ice)
    q_rim = UT.clamp_to_nonneg(q_rim)
    b_rim = UT.clamp_to_nonneg(b_rim)

    # Convert to volumetric quantities for P3 functions
    L_lcl = q_lcl * ρ  # [kg lcl / m³ air]
    L_rai = q_rai * ρ  # [kg rai / m³ air]
    N_lcl = n_lcl * ρ  # [1 / m³ air]
    N_rai = n_rai * ρ  # [1 / m³ air]
    L_ice = q_ice * ρ  # [kg ice / m³ air]
    N_ice = n_ice * ρ  # [1 / m³ air]
    L_rim = q_rim * ρ  # [kg rim / m³ air]
    B_rim = b_rim * ρ  # [m³ rim / m³ air]

    # P3State construction handles the F_rim/ρ_rim regularisation and
    # domain clamps internally — see `get_state_from_prognostic`. The
    # regularised ratios (`UT.rime_mass_fraction`, `UT.rime_density`) used
    # there smoothly blend to zero as their denominators shrink, avoiding
    # the hard discontinuity at `q_ice = ϵ` / `b_rim = ϵ`.
    state = CMP3.get_state_from_prognostic(mp.ice.scheme, L_ice, N_ice, L_rim, B_rim)

    # Unpack warm rain parameters
    aps = mp.warm_rain.air_properties
    subdep = mp.warm_rain.subdep

    # Initialize ice-related tendencies
    dq_ice_dt = zero(ρ)
    dn_ice_dt = zero(ρ)
    dq_rim_dt = zero(ρ)
    db_rim_dt = zero(ρ)

    # --- Core Warm Rain Processes (includes activation, see warm_rain_tendencies_2m) ---
    warm = warm_rain_tendencies_2m(mp.warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai, w, p)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt
    dn_lcl_activation_dt = warm.dn_lcl_activation_dt
    # NOTE on latent-heat coupling: per-process phase-change rates
    # (`S_cond`, `S_dep`, `S_frz_net`) are no longer returned. Hosts that
    # need a combined LH rate can compute it from the prognostic mass
    # tendencies via the identity
    #   ḣ_lh = Lv(T) · (dq_lcl/dt + dq_rai/dt) + Ls(T) · dq_ice/dt
    # which holds because Ls = Lv + Lf (sublimation = vaporisation +
    # fusion) makes the L_f · S_frz term cancel when expressed in
    # species coordinates. This keeps BMT dt-agnostic and the host's LH
    # rate self-consistent with whatever clipping the host applies to
    # the species tendencies.

    # --- P3 Ice Processes ---
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf
    ice_nucleation = mp.ice.ice_nucleation
    inp_depletion_model = mp.ice.inp_depletion_model

    quad = CMP3.ChebyshevGauss(mp.ice.quadrature_order)

    # Only compute ice processes if there is ice mass/number present
    if q_ice > ϵₘ && n_ice > ϵₙ

        # --- Liquid-ice collisions ---
        coll = CMP3.bulk_liquid_ice_collision_sources(
            state, logλ, pdf_c, pdf_r, L_lcl, N_lcl, L_rai, N_rai, aps, tps, vel, ρ, T;
            quad,
        )
        dq_lcl_dt += coll.∂ₜq_c
        dq_rai_dt += coll.∂ₜq_r
        dn_lcl_dt += coll.∂ₜN_c / ρ
        dn_rai_dt += coll.∂ₜN_r / ρ
        dq_ice_dt += coll.∂ₜL_ice / ρ
        dq_rim_dt += coll.∂ₜL_rim / ρ
        db_rim_dt += coll.∂ₜB_rim / ρ

        # --- Ice self-collection (aggregation) ---
        S_ice_agg = CMP3.ice_self_collection(state, logλ, vel, ρ; quad)
        dn_ice_dt -= S_ice_agg.dNdt / ρ

        # Ice melting (above freezing temperature)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        melt = ifelse(T > T_freeze,
            CMP3.ice_melt(vel, aps, tps, T, ρ, state, logλ; quad),
            (; dNdt = zero(ρ), dLdt = zero(ρ))
        )
        # Specific (per-kg-air) ice-mass melt rate.
        ∂ₜq_ice_melt = melt.dLdt / ρ
        ∂ₜn_ice_melt = melt.dNdt / ρ
        # Melting converts ice to rain.
        dq_rai_dt += ∂ₜq_ice_melt
        dn_rai_dt += ∂ₜn_ice_melt  # Melted ice becomes rain drops
        dq_ice_dt -= ∂ₜq_ice_melt
        dn_ice_dt -= ∂ₜn_ice_melt  # Ice particles consumed by melting
        # Rim mass and rim volume drain proportionally to ice mass during
        # melting (uniform-melt assumption — F_rim and ρ_rim are
        # preserved through the melt). Without this, q_rim and b_rim
        # were conserved while q_ice drained, so q_rim/q_ice → ∞ at the
        # 0 °C isoline (analysis showed q_rim/q_ice up to ~470× just
        # below the freezing line in newBMT_tau1d_seed20260423).
        #
        #   ∂ₜq_rim = -F_rim · ∂ₜq_ice
        #   ∂ₜb_rim = -F_rim/ρ_rim · ∂ₜq_ice
        #
        # together preserve both F_rim and ρ_rim. We pull F_rim and
        # ρ_rim from `state` (the regularised P3State built at the top
        # of the function), so trace-mass cells with q_rim near or
        # above q_ice see the clamped F_rim ≤ 1 − ε rather than the raw
        # prognostic ratio. Guard ρ_rim ≈ 0 (no rim volume → no
        # b_rim drain) — `_regularised_ratio` returns 0 for `ρ_rim`
        # when `b_rim` is below the smoothing scale.
        dq_rim_dt -= ∂ₜq_ice_melt * state.F_rim
        db_rim_dt -= ifelse(state.ρ_rim > 0,
            ∂ₜq_ice_melt * state.F_rim / state.ρ_rim, zero(FT),
        )
    end

    # --- ----------------------------- ---
    # --- Ice nucleation (F23 + Bigg)   ---
    # --- ----------------------------- ---
    #
    # Two parallel heterogeneous-freezing channels, MM15-aligned:
    #
    #   Pathway 1 — F23 deposition / condensation-freezing nucleation
    #               (analog of Fortran `qinuc`): vapor → pristine ice on
    #               an INP. F_rim = 0 at genesis; subsequent growth via
    #               Pathway 5 (MM2015 vapor deposition).
    #
    #   Pathway 2 — F23-bounded Bigg immersion freezing of cloud drops
    #               (analog of Fortran `qcheti`): drop + INP → fully-rimed
    #               ice (embryo graupel). F_rim = 1 at genesis; mass and
    #               number drained from q_lcl, n_lcl.
    #
    # Both channels read the same Frostenberg-2023 INPC climatology
    # (with a stochastic OU-driven shift `inpc_log_shift` supplied by
    # the driver). They differ in mass source, gates, and F_rim.
    #
    # Pathway 3 (rain immersion `qrheti`) appears further below; it is
    # MM15's Bigg formula on the rain PSD without an F23 cap. Pathway 4
    # (homogeneous, `qchom`/`qrhom`) is not yet implemented in BMT.
    #
    # See `kid_jouan_ou/F23_REGIME_AWARE_FIX.md` for the full pathway
    # inventory + the MM15 ↔ Fortran reference.

    # F23 activation timescale, taken from the depletion model so the
    # host can override it by passing e.g.
    # `inp_depletion_model = NIceProxyDepletion(τ_act = FT(60))` when
    # constructing the P3IceParams.
    τ_act = inp_depletion_model.τ_act
    # Vapor deposition nucleation size
    D_nuc = FT(10e-6)  # 10 μm nascent crystal - small-D tail of the P3
    m_nuc = p3.ρ_i * CO.volume_sphere_D(D_nuc)

    # F23 INP-activation depletion proxy. With `NIceProxyDepletion`
    # (default) this returns `n_ice` — the legacy always-on proxy.
    # With `PrognosticINPDecay` it returns `n_INP_used` — the host's
    # Phillips-style activation-memory tracer. The host is responsible
    # for advancing `n_INP_used` per the source/sink budget written in
    # `F23_REGIME_AWARE_FIX.md`; BMT just dispatches.
    n_active = CM_HetIce.n_active(inp_depletion_model, n_ice, n_INP_used)

    # ---- Pathway 1: F23 deposition nucleation (vapor → pristine ice) ----
    #
    # Target-relaxation form, structurally identical to Cooper-style qinuc
    # with F23's log-normal target substituted for Cooper's exponential.
    # We route through `CM_HetIce.f23_deposition_rate`, which uses the
    # strict-MM15 doc gates `T < T_freeze − 15 K ∧ S_i > 0.05` by
    # default. Hosts that need the legacy always-on form can override
    # by passing e.g. `T_thresh = FT(2000), S_i_thresh = FT(-2)` to the
    # call below.
    f23_dep = CM_HetIce.f23_deposition_rate(
        ice_nucleation, tps, T, ρ, q_tot, q_lcl + q_rai, q_ice, n_active;
        m_nuc, τ_act, inpc_log_shift,
    )

    dn_ice_dt += f23_dep.∂ₜn_frz
    dq_ice_dt += f23_dep.∂ₜq_frz
    # NO contribution to q_rim, b_rim — pristine deposition crystals have F_rim = 0.

    # ---- Pathway 2: F23-bounded Bigg immersion freezing of cloud drops ----
    #
    # Bigg (1953) volume-weighted kinetics over the SB2006 cloud-drop PSD
    # (vanilla MM15 / qcheti, both rates ≥ 0 by construction):
    #
    #     ∂ₜn_imm^Bigg = J_bigg(T)       · (π/6)  · M_D³(N_lcl, λ_c, ν_c, μ_c)
    #     ∂ₜq_imm^Bigg = J_bigg(T) · ρ_w · (π/6)² · M_D⁶(N_lcl, λ_c, ν_c, μ_c)
    #
    # F23 INPC imposes an upper bound on the activation rate (the new
    # piece for clean-air / OU-SIF regimes; ≥ 0):
    #
    #     ∂ₜn_imm^INPC = INPC / τ_act
    #
    # Realised rate is the smaller of the two; mass tracks the limiting
    # branch via the size-weighted Bigg mass:
    #
    #     ∂ₜn_imm = min(∂ₜn_imm^Bigg, ∂ₜn_imm^INPC)
    #     ∂ₜq_imm = ∂ₜq_imm^Bigg · ∂ₜn_imm / ∂ₜn_imm^Bigg
    #
    # Per MM15 (Fortran `qcheti`) frozen cloud drops are fully rimed
    # (embryo graupel) → mass adds to BOTH qitot and qirim/birim.
    cld_bigg = CM_HetIce.liquid_freezing_rate(
        mp.ice.rain_freezing, pdf_c, tps, q_lcl, ρ, N_lcl, T,
    )
    cld_cap = CM_HetIce.f23_immersion_limit_rate(
        ice_nucleation, T, ρ; τ = τ_act, inpc_log_shift, n_active,
    )
    ∂ₜn_imm = min(cld_bigg.∂ₜn_frz, cld_cap.∂ₜn_frz)
    # Mass = (size-weighted Bigg mass) × (rate ratio). When Bigg is
    # silent (off-gate or zero N_lcl/q_lcl), `cld_bigg.∂ₜn_frz` is exactly
    # zero — and so is `∂ₜn_imm` — so any safe ratio gives ∂ₜq_imm = 0.
    # Use 0 as the safe value to avoid 0/0 NaN.
    ∂ₜq_imm = ifelse(cld_bigg.∂ₜn_frz > zero(FT), cld_bigg.∂ₜq_frz * ∂ₜn_imm / cld_bigg.∂ₜn_frz, zero(FT))

    # Drain liquid:
    dq_lcl_dt -= ∂ₜq_imm
    dn_lcl_dt -= ∂ₜn_imm
    # Add to ice as fully-rimed embryo graupel:
    dq_ice_dt += ∂ₜq_imm
    dn_ice_dt += ∂ₜn_imm
    dq_rim_dt += ∂ₜq_imm           # F_rim = 1 (frozen drop)
    db_rim_dt += ∂ₜq_imm / p3.ρ_i  # solid-ice rime volume

    # --- Ice Sublimation / Deposition ---
    n_per_q_ice = ifelse(q_ice > ϵₘ, n_ice / q_ice, zero(n_ice))
    # Deposition/sublimation of cloud ice
    ∂ₜq_ice_dep = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(subdep, tps, q_tot, q_lcl, q_ice, q_rai, zero(q_ice), ρ, T)
    # No ice deposition above freezing (lack of INPs)
    ∂ₜq_ice_dep = ifelse(T > tps.T_freeze, min(∂ₜq_ice_dep, zero(T)), ∂ₜq_ice_dep)
    # During sublimation, the number of ice particles decreases in proportion to the mean ice mass
    # During deposition, the number of ice particles remain unchanged
    ∂ₜn_ice_dep = ifelse(∂ₜq_ice_dep < 0, n_per_q_ice * ∂ₜq_ice_dep, zero(∂ₜq_ice_dep))
    dq_ice_dt += ∂ₜq_ice_dep
    dn_ice_dt += ∂ₜn_ice_dep
    # Rim mass and rim volume drain proportionally during *sublimation*
    # (∂ₜq_ice_dep < 0). Same uniform-melt-style assumption as the melt
    # branch: F_rim and ρ_rim are preserved through the phase change, so
    # q_rim drains at F_rim · ∂ₜq_ice_sub and b_rim drains at
    # (F_rim/ρ_rim) · ∂ₜq_ice_sub. Per MM15 Eqs. for S_qrim and S_Br,
    # `−(qrim/qi)·QISUB` and `−(qrim/(ρ_r·qi))·QISUB` (matching the
    # `−QIMLT` rim drains we already have). Without this branch, the
    # cirrus-level F_rim ≈ 1 band (where ice nucleates, partially
    # sublimates, and leaves rim mass behind) would persist
    # indefinitely — same root cause as the melt-side bug, just on the
    # other side of the freezing line. Deposition (∂ₜq_ice_dep > 0)
    # adds pristine ice mass with F_rim = 0, so q_rim/b_rim are
    # unaffected; we only drain on the sublimation branch.
    ∂ₜq_ice_sub = min(∂ₜq_ice_dep, 0)   # ≤ 0; zero on the deposition branch
    dq_rim_dt += ∂ₜq_ice_sub * state.F_rim
    db_rim_dt += ifelse(state.ρ_rim > 0, ∂ₜq_ice_sub * state.F_rim / state.ρ_rim, zero(FT))

    # --- Ice number adjustment for mass limits ---
    # Number adjustment for ice mass limits (Horn 2012, DOI: 10.5194/gmd-5-345-2012).
    # Nudges n_ice toward [q_ice / x_max, q_ice / x_min] over timescale τ.
    numadj = (;
        τ = FT(100), 
        x_min = FT(1e-12),  # min mean ice particle mass [kg] (~10 μm crystal)
        x_max = FT(1e-5),   # max mean ice particle mass [kg] (~5 mm aggregate)
    )
    ∂ₜn_ice_numadj = CM2.number_tendency_from_mass_limits(numadj, q_ice, n_ice)
    dn_ice_dt += ∂ₜn_ice_numadj

    # --- Rain Heterogeneous Freezing (Bigg 1953) ---
    # Immersion freezing of rain drops, as in Morrison & Milbrandt (2015).
    rain_frz = CM_HetIce.liquid_freezing_rate(mp.ice.rain_freezing, pdf_r, tps, q_rai, ρ, N_rai, T)

    # Rain → ice (frozen rain is fully rimed, per MM15)
    dq_rai_dt -= rain_frz.∂ₜq_frz
    dn_rai_dt -= rain_frz.∂ₜn_frz
    dq_ice_dt += rain_frz.∂ₜq_frz
    dn_ice_dt += rain_frz.∂ₜn_frz
    dq_rim_dt += rain_frz.∂ₜq_frz
    db_rim_dt += rain_frz.∂ₜq_frz / p3.ρ_i  # ρ_i = 916.7 kg m⁻³, the density of solid bulk ice

    # Aerosol activation is folded into `warm_rain_tendencies_2m` above —
    # `dn_lcl_activation_dt` from `warm` is already included in `dn_lcl_dt`.

    # Source rate for the host's F23 activation-memory tracer
    # (`n_INP_used` in `PrognosticINPDecay` mode). Hosts in
    # `NIceProxyDepletion` mode can ignore this field.
    dn_INP_used_source_dt = f23_dep.∂ₜn_frz + ∂ₜn_imm

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt,
              dq_ice_dt, dn_ice_dt, dq_rim_dt, db_rim_dt,
              dn_lcl_activation_dt, dn_INP_used_source_dt)
end

"""
    repair_ice_state(N_ice, L_ice, ρ; m_crystal_max = 1e-5, q_ice_floor = 1e-14)

Enforce two complementary physical constraints on the mean P3 ice-crystal
mass `L_ice / N_ice`:

1. **Upper-mass bound** (mid-column "bowling ball" pathology): when
   `L_ice > q_ice_floor · ρ` and the implied mean mass would exceed
   `m_crystal_max`, raise `N_ice` to `L_ice / m_crystal_max`.
2. **Trace-ice phantom-number** (top-of-column F23 + sedimentation
   pathology): when `L_ice ≤ q_ice_floor · ρ` (mass effectively zero),
   force `N_ice = 0`. Without this, a column top where mass-weighted
   sedimentation drains all ice mass faster than number can leak F23
   nucleation events into a runaway `N_ice` count with no corresponding
   mass.

Both branches are mass-preserving — only `N_ice` is touched.

# Background (issue 012 + top-of-domain debug 2026-04-24)

In a 1-D column with P3 sedimentation, the mass-weighted fall velocity
differs from the number-weighted fall velocity. Two distinct symptoms
result:

- Mid-column "donor" cells lose mass faster than number → end with
  almost-zero `N_ice` but non-trivial `L_ice` → mean mass ~10^−3 kg
  (heavier than a marble) → runaway riming. Branch (1) above repairs
  this.
- Top-of-column "trace ice" cells — where F23 nucleation continually
  sources new particles but mass-weighted sedimentation drains the
  starter mass DOWN every step — accumulate `N_ice` with no
  corresponding `L_ice`. Branch (2) above repairs this. Without it,
  the top cell can show `N_ice ~ 10^9 m⁻³` (cirrus×10³) at `L_ice = 0`
  after a few hours.

Two design properties make this the right layer:

1. **Mass-preserving.** `L_ice` is left untouched. Only the
   bookkeeping-error variable (`N_ice`) is nudged.
2. **Opt-in.** BMT does NOT call this internally — callers decide when to
   apply it (typically from a `DiscreteCallback` so the correction hits
   exactly once per step). Surreptitiously mutating caller state inside
   BMT would hide the very bug a user is trying to diagnose.

# Units convention

This helper operates on **per-volume** state, using the CliMA casing
convention `lower q/n = per-kg`, `upper L/N = per-m³`, matching how
`ClimaCore` fields store the prognostic variables in both KiD and
ClimaAtmos:

- `N_ice` [1/m³] — number concentration per unit volume
- `L_ice` [kg/m³] — mass per unit volume (callers typically pass
   `Y.ρq_ice` or similar — same quantity, different spelling)
- `ρ` [kg/m³] — air density

Inputs in specific units (`q_ice` [kg/kg], `n_ice` [1/kg]) are NOT
supported: the `q_ice_floor · ρ` threshold assumes the first positional
argument is per-volume.

# Arguments
- `N_ice`: current cloud-ice number concentration [1/m³].
- `L_ice`: current cloud-ice mass per unit volume [kg/m³].
- `ρ`: air density [kg/m³].

# Keyword arguments
- `m_crystal_max = 1e-5 [kg]`: upper bound on an individual P3 ice
   crystal mass. 1e-5 kg corresponds to a ~5 mm aggregate, which the P3
   `numadj.x_max` already uses.
- `q_ice_floor = 1e-14 [kg/kg]`: specific-humidity floor; below this
   (after conversion to per-volume via `× ρ`), the cell is in trace-ice
   territory where the BMT input clamp handles it. Repair is a no-op
   below the floor to avoid creating droplets from numerical noise.

# Return
A revised `N_ice` value of the same type as the input.

# Notes
Callers integrating into host codes (KiD / ClimaAtmos) typically apply
it as a broadcast over the column:
`@. Y.N_ice = BMT.repair_ice_state(Y.N_ice, Y.ρq_ice, ρ)`.
"""
@inline function repair_ice_state(N_ice, L_ice, ρ;
        m_crystal_max = oftype(N_ice, 1e-5),
        q_ice_floor = oftype(N_ice, 1e-14))
    FT = typeof(N_ice)
    # threshold and L_ice both in [kg/m³]; q_ice_floor is in [kg/kg]
    # and × ρ converts it to the per-volume scale.
    threshold = FT(q_ice_floor) * ρ
    # imputed: [kg/m³] / [kg/crystal] = [crystals/m³], matching N_ice.
    imputed = L_ice / FT(m_crystal_max)
    # Branch by trace-ice condition:
    #   L_ice > threshold  → upper-mass bound (raise N_ice if needed)
    #   L_ice ≤ threshold  → trace ice; collapse N_ice to 0 (no phantom number)
    return ifelse(L_ice > threshold, max(N_ice, imputed), zero(N_ice))
end

end # module BulkMicrophysicsTendencies
