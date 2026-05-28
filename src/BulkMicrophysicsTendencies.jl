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
    Microphysics1Moment(), mp, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno
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
import ...ThermodynamicsInterface as TDI

export MicrophysicsScheme,
    Microphysics0Moment,
    Microphysics1Moment,
    Microphysics2Moment,
    TendencyMode,
    Instantaneous,
    InstantaneousVerbose,
    LinearizedAverage,
    bulk_microphysics_tendencies

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

# --- Tendency output mode dispatch ---

"""
    TendencyMode

Abstract type for selecting the output mode of `bulk_microphysics_tendencies`.
"""
abstract type TendencyMode end

"""
    Instantaneous <: TendencyMode

Return raw nonlinear tendencies from a single evaluation of all microphysical
processes (no linearization, no time-averaging).
"""
struct Instantaneous <: TendencyMode end

"""
    InstantaneousVerbose <: TendencyMode

Return all individual source terms alongside aggregated `dq_*_dt` tendencies.
Useful for model diagnostics. Only works with instantaneous (nonlinear) evaluation.
"""
struct InstantaneousVerbose <: TendencyMode end

"""
    LinearizedAverage <: TendencyMode

Return time-averaged tendencies computed via repeated linearized implicit substeps.
This is the mode used operationally by ClimaAtmos.
"""
struct LinearizedAverage <: TendencyMode end

# --- 1-Moment Microphysics ---

# --- Internal helpers ---

"""
Compute all individual 1-moment microphysics source terms in a single pass.

This is the **single source of truth** for which microphysical processes are
called and with what arguments. Both the raw tendency aggregation and the
linearized operator construction consume this output.

Constructs two `NamedTuple`s that are passed to all process functions
(see `Microphysics1M` module docs for the full convention):
- `micro = (; q_tot, q_lcl, q_icl, q_rai, q_sno)` — specific humidities (kg/kg)
- `thermo = (; ρ, T)` — air density (kg/m³) and temperature (K)

Naming convention: `S_process_species1_species2`
 - process: physical mechanism (phase_change, acnv, accr, melt, accr_melt, accr_freeze)
 - species1, species2: interacting pair (not from/to)
 - `_cold` / `_warm` suffix for two-sided collision arms (inactive arm = zero)

Returns a `NamedTuple` of ~19 scalar source terms.  All two-sided collision
processes are pre-routed by temperature, so consumers never need `is_warm`.
"""
@inline function _microphysics_source_terms(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    # Clamp negative inputs to zero (robustness against numerical errors)
    ρ = UT.clamp_to_nonneg(ρ)
    q_tot = UT.clamp_to_nonneg(q_tot)
    q_lcl = UT.clamp_to_nonneg(q_lcl)
    q_icl = UT.clamp_to_nonneg(q_icl)
    q_rai = UT.clamp_to_nonneg(q_rai)
    q_sno = UT.clamp_to_nonneg(q_sno)

    FT = typeof(q_tot)
    opts = mp.options

    # Construct state tuples (reused across all process calls)
    micro = (; q_tot, q_lcl, q_icl, q_rai, q_sno)
    thermo = (; ρ, T)

    # --- Phase change: vapor ↔ cloud condensate (bidirectional, ±) ---
    S_phase_change_vap_lcl = CMNonEq.conv_q_vap_to_q_lcl(opts.cloud_liquid_formation, mp, tps, micro, thermo)
    S_phase_change_vap_icl = CMNonEq.conv_q_vap_to_q_icl(opts.cloud_ice_formation, mp, tps, micro, thermo)

    # --- Autoconversion (cloud → precipitation, ≥ 0) ---
    S_acnv_lcl_rai = CM1.conv_q_lcl_to_q_rai(opts.rain_autoconversion, mp, tps, micro, thermo)
    S_acnv_icl_sno = CM1.conv_q_icl_to_q_sno(opts.snow_autoconversion, mp, tps, micro, thermo)

    # --- Accretion (collisions between species) ---
    is_warm = T >= TDI.T_freeze(tps)

    # Cloud liquid + rain → rain
    S_accr_lcl_rai = CM1.accretion(opts.cloud_liquid_rain_accretion, mp, tps, micro, thermo)

    # Cloud liquid + snow: product goes to sno (cold) or rai (warm), plus thermal melt
    (; S_accr, S_melt) = CM1.accretion(opts.cloud_liquid_snow_accretion, mp, tps, micro, thermo)
    S_accr_lcl_sno_cold = ifelse(is_warm, zero(FT), S_accr)    # lcl → sno (cold)
    S_accr_lcl_sno_warm = ifelse(is_warm, S_accr, zero(FT))    # lcl → rai (warm)
    S_accr_melt_lcl_sno = S_melt                                # thermal melt of sno from warm lcl (already zero when cold)

    # Cloud ice + rain → snow (ice-side sink)
    S_accr_icl_rai = CM1.accretion(opts.cloud_ice_rain_accretion, mp, tps, micro, thermo)

    # Rain frozen in cloud ice + rain collision → snow (rain sink)
    S_accr_freeze_icl_rai = CM1.accretion_rain_sink(opts.cloud_ice_rain_accretion, mp, tps, micro, thermo)

    # Cloud ice + snow → snow
    S_accr_icl_sno = CM1.accretion(opts.cloud_ice_snow_accretion, mp, tps, micro, thermo)

    # Rain-snow collisions: split into cold/warm arms (inactive arm = zero)
    (; S_rai_sno, S_sno_rai, S_melt) = CM1.accretion_snow_rain(opts.rain_snow_accretion, mp, tps, micro, thermo)
    S_accr_rai_sno_cold = ifelse(is_warm, zero(FT), S_rai_sno) # cold arm: rai freezes → sno
    S_accr_rai_sno_warm = ifelse(is_warm, S_sno_rai, zero(FT)) # warm arm: sno melts → rai
    S_accr_melt_rai_sno = ifelse(is_warm, S_melt, zero(FT))    # thermal melt of sno from warm rai

    # --- Phase change: precipitation ↔ vapor ---
    S_phase_change_vap_rai = CM1.conv_q_rai_to_q_vap(opts.rain_condensation_evaporation, mp, tps, micro, thermo)
    S_phase_change_vap_sno = CM1.conv_q_sno_to_q_vap(opts.snow_deposition_sublimation, mp, tps, micro, thermo)

    # --- Melting ---
    S_melt_icl_lcl = CM1.conv_q_icl_to_q_lcl(opts.cloud_ice_melt, mp, tps, micro, thermo)
    S_melt_sno_rai = CM1.conv_q_sno_to_q_rai(opts.snow_melt, mp, tps, micro, thermo)

    return (;
        S_phase_change_vap_lcl, S_phase_change_vap_icl,
        S_acnv_lcl_rai, S_acnv_icl_sno,
        S_accr_lcl_rai, S_accr_lcl_sno_cold, S_accr_lcl_sno_warm, S_accr_melt_lcl_sno,
        S_accr_icl_rai, S_accr_freeze_icl_rai, S_accr_icl_sno,
        S_accr_rai_sno_cold, S_accr_rai_sno_warm, S_accr_melt_rai_sno,
        S_phase_change_vap_rai, S_phase_change_vap_sno,
        S_melt_icl_lcl, S_melt_sno_rai,
    )
end

"""
Aggregate individual source terms into the four hydrometeor tendency totals.

This is the **single location** where the sign convention of source terms
to tendency accumulators is defined.  All temperature-dependent routing is
pre-applied in `_microphysics_source_terms` (cold/warm arms), so every term
here appears with a fixed sign — no `ifelse` branching.
"""
@inline function _aggregate_tendencies(src)
    dq_lcl_dt =
        src.S_phase_change_vap_lcl - src.S_acnv_lcl_rai - src.S_accr_lcl_rai -
        src.S_accr_lcl_sno_cold - src.S_accr_lcl_sno_warm + src.S_melt_icl_lcl

    dq_icl_dt =
        src.S_phase_change_vap_icl - src.S_acnv_icl_sno - src.S_accr_icl_rai -
        src.S_accr_icl_sno - src.S_melt_icl_lcl

    dq_rai_dt =
        src.S_acnv_lcl_rai + src.S_accr_lcl_rai +
        src.S_accr_lcl_sno_warm + src.S_accr_melt_lcl_sno -
        src.S_accr_freeze_icl_rai -
        src.S_accr_rai_sno_cold + src.S_accr_rai_sno_warm + src.S_accr_melt_rai_sno +
        src.S_phase_change_vap_rai + src.S_melt_sno_rai

    dq_sno_dt =
        src.S_acnv_icl_sno +
        src.S_accr_lcl_sno_cold - src.S_accr_melt_lcl_sno +
        src.S_accr_icl_rai + src.S_accr_freeze_icl_rai +
        src.S_accr_icl_sno +
        src.S_accr_rai_sno_cold - src.S_accr_rai_sno_warm - src.S_accr_melt_rai_sno +
        src.S_phase_change_vap_sno - src.S_melt_sno_rai

    return (; dq_lcl_dt, dq_icl_dt, dq_rai_dt, dq_sno_dt)
end

"""
Construct a local linear approximation of 1-moment microphysics tendencies
from pre-computed source terms:

    dq/dt ≈ M * q + e

using a donor-based linearization:
- donor → receiver transfers are represented as `D * q_donor`
- vapor → condensate sources are treated as constants (`e`)
- condensate sinks are treated as linear sinks (`-D * q`)

All coefficients use `D = S / max(q_min, q_donor)` for robustness.

Returns a `NamedTuple` containing the nonzero entries of `M` and `e`.
"""
@inline function _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min)
    FT = typeof(src.S_phase_change_vap_lcl)

    M11 = zero(FT)
    M12 = zero(FT)
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

    # --- Phase change: vapor ↔ cloud condensate ---
    D = src.S_phase_change_vap_lcl / max(q_min, q_lcl)
    is_source = src.S_phase_change_vap_lcl >= zero(FT)
    e1 += ifelse(is_source, src.S_phase_change_vap_lcl, zero(FT))
    M11 += ifelse(is_source, zero(FT), D)

    D = src.S_phase_change_vap_icl / max(q_min, q_icl)
    is_source = src.S_phase_change_vap_icl >= zero(FT)
    e2 += ifelse(is_source, src.S_phase_change_vap_icl, zero(FT))
    M22 += ifelse(is_source, zero(FT), D)

    # --- Melt: ice cloud → liquid cloud ---
    D = src.S_melt_icl_lcl / max(q_min, q_icl)
    M22 -= D
    M12 += D

    # --- Autoconversion: donor-based transfer ---
    D = src.S_acnv_lcl_rai / max(q_min, q_lcl)
    M11 -= D
    M31 += D

    D = src.S_acnv_icl_sno / max(q_min, q_icl)
    M22 -= D
    M42 += D

    # --- Accretion: donor-based transfer ---
    D = src.S_accr_lcl_rai / max(q_min, q_lcl)
    M11 -= D
    M31 += D

    # lcl + sno accretion (cold/warm arms already zeroed)
    D_cold = src.S_accr_lcl_sno_cold / max(q_min, q_lcl)
    D_warm = src.S_accr_lcl_sno_warm / max(q_min, q_lcl)
    M11 -= D_cold + D_warm
    M31 += D_warm           # warm: lcl → rai
    M41 += D_cold           # cold: lcl → sno

    # thermal melt of sno from warm lcl
    D = src.S_accr_melt_lcl_sno / max(q_min, q_sno)
    M44 -= D
    M34 += D

    D = src.S_accr_icl_rai / max(q_min, q_icl)
    M22 -= D
    M42 += D

    D = src.S_accr_icl_sno / max(q_min, q_icl)
    M22 -= D
    M42 += D

    # rain frozen in icl + rai collision
    D = src.S_accr_freeze_icl_rai / max(q_min, q_rai)
    M33 -= D
    M43 += D

    # warm arm: sno melts → rai (already zero when cold)
    D = src.S_accr_rai_sno_warm / max(q_min, q_sno)
    M44 -= D
    M34 += D

    # thermal melt of sno from warm rai (already zero when cold)
    D = src.S_accr_melt_rai_sno / max(q_min, q_sno)
    M44 -= D
    M34 += D

    # cold arm: rai freezes → sno (already zero when warm)
    D = src.S_accr_rai_sno_cold / max(q_min, q_rai)
    M33 -= D
    M43 += D

    # --- Rain phase change: sink to vapor (always zero or negative) ---
    D = (-src.S_phase_change_vap_rai) / max(q_min, q_rai)
    M33 -= D

    # --- Snow phase change: deposition/sublimation ---
    D = src.S_phase_change_vap_sno / max(q_min, q_sno)
    is_source = src.S_phase_change_vap_sno >= zero(FT)
    e4 += ifelse(is_source, src.S_phase_change_vap_sno, zero(FT))
    M44 += ifelse(is_source, zero(FT), D)

    # --- Snow melt: snow → rain ---
    D = src.S_melt_sno_rai / max(q_min, q_sno)
    M44 -= D
    M34 += D

    return (
        M11 = M11, M12 = M12, M22 = M22,
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

The system uses a sparse structure specific to the 1-moment microphysics model.
`q_lcl` and `q_icl` as well as `q_rai` and `q_sno` are solved from a coupled 2×2 system.

Because sinks are linearized as `-D q`, they are effectively integrated as
exponential decays over the substep.
"""
@inline function _linearized_implicit_step(
    ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt,
)

    FT = typeof(q_tot)

    src = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    q_min = TDI.TD.Parameters.q_min(tps)
    lin = _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min)

    invΔt = one(FT) / Δt

    # A = I/Δt - M
    a11 = invΔt - lin.M11
    a12 = -lin.M12
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

    # Solve 2×2 system for q_lcl, q_icl (coupled via ice melt M12)
    det12 = a11 * a22  # a21 = 0
    q_lcl_new = (b1 * a22 - a12 * b2) / det12
    q_icl_new = a11 * b2 / det12

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

# --- Public API: bulk_microphysics_tendencies with TendencyMode dispatch ---

"""
    bulk_microphysics_tendencies(
        ::Instantaneous, ::Microphysics1Moment, mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute all 1-moment microphysics tendencies in one fused call.

Returns a NamedTuple with all source/sink terms for hydrometeor species.
This is a pure function of local thermodynamic state, suitable for:
- Point quadrature over subgrid-scale (T, q_tot) distributions
- GPU kernel evaluation
- Unit testing in isolation

# Arguments
- `mp`: Microphysics1MParams parameter container
- `tps`: Thermodynamics parameters
- `ρ`: Air density [kg/m³]
- `T`: Temperature [K]
- `q_tot`: Total water specific content [kg/kg]
- `q_lcl`: Cloud liquid water specific content [kg/kg]
- `q_icl`: Cloud ice specific content [kg/kg]
- `q_rai`: Rain specific content [kg/kg]
- `q_sno`: Snow specific content [kg/kg]

# Returns
`NamedTuple` with fields:
- `dq_lcl_dt`: Cloud liquid tendency [kg/kg/s]
- `dq_icl_dt`: Cloud ice tendency [kg/kg/s]
- `dq_rai_dt`: Rain tendency [kg/kg/s]
- `dq_sno_dt`: Snow tendency [kg/kg/s]

# Notes
- Negative specific contents are clamped to zero for robustness.
- Does NOT apply timestep-dependent limiters.
"""
@inline function bulk_microphysics_tendencies(
    ::Instantaneous, ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    return _aggregate_tendencies(src)
end

"""
    bulk_microphysics_tendencies(
        ::InstantaneousVerbose, ::Microphysics1Moment, mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )

Compute all 1-moment microphysics tendencies and return both aggregated
tendencies (`dq_*_dt`) and all individual source terms (`S_*`).

Useful for model diagnostics. The `dq_*_dt` fields are identical to those
returned by `Instantaneous()`.

# Returns
`NamedTuple` with all fields from `Instantaneous()` plus individual source
terms: `S_phase_change_vap_lcl`, `S_phase_change_vap_icl`, `S_acnv_lcl_rai`,
`S_acnv_icl_sno`, etc.
"""
@inline function bulk_microphysics_tendencies(
    ::InstantaneousVerbose, ::Microphysics1Moment, mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(
        Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    agg = _aggregate_tendencies(src)
    return merge(agg, src)
end

"""
    bulk_microphysics_tendencies(
        ::LinearizedAverage, ::Microphysics1Moment, mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1,
    )

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

# Returns
`NamedTuple` with fields:
- `dq_lcl_dt`: Cloud liquid tendency [kg/kg/s]
- `dq_icl_dt`: Cloud ice tendency [kg/kg/s]
- `dq_rai_dt`: Rain tendency [kg/kg/s]
- `dq_sno_dt`: Snow tendency [kg/kg/s]
"""
@inline function bulk_microphysics_tendencies(
    ::LinearizedAverage,
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
        rates = _linearized_implicit_step(
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
@inline function warm_rain_tendencies_2m(warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai)

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

    # --- Condensation of vapor / evaporation of cloud liquid water ---
    micro_mock = (; q_tot, q_lcl, q_icl = q_ice, q_rai, q_sno = zero(q_ice))
    thermo_mock = (; ρ, T)
    ∂ₜq_lcl_cond = CMNonEq.conv_q_vap_to_q_lcl(
        CMP.CloudLiquidFormation(condevap.τ_relax), nothing, tps, micro_mock, thermo_mock,
    )
    ∂ₜn_lcl_cond = zero(∂ₜq_lcl_cond)  # neglect number change from condensation/evaporation
    dq_lcl_dt += ∂ₜq_lcl_cond
    dn_lcl_dt += ∂ₜn_lcl_cond

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
    dn_lcl_sc = CM2.cloud_liquid_self_collection(sb.acnv, sb.pdf_c, q_lcl, ρ, acnv.dN_lcl_dt)
    dn_lcl_dt += dn_lcl_sc / ρ

    # --- Accretion ---
    accr = CM2.accretion(sb, q_lcl, q_rai, ρ, N_lcl)
    dq_lcl_dt += accr.dq_lcl_dt
    dq_rai_dt += accr.dq_rai_dt
    dn_lcl_dt += accr.dN_lcl_dt / ρ

    # --- Rain self-collection ---
    dn_rai_sc = CM2.rain_self_collection(sb.pdf_r, sb.self, q_rai, ρ, N_rai)
    dn_rai_dt += dn_rai_sc / ρ

    # --- Rain breakup ---
    dn_rai_br = CM2.rain_breakup(sb.pdf_r, sb.brek, q_rai, ρ, N_rai, dn_rai_sc)
    dn_rai_dt += dn_rai_br / ρ

    # --- Number adjustment for mass limits ---
    # Cloud liquid
    dn_lcl_inc = CM2.number_increase_for_mass_limit(sb.numadj, sb.pdf_c.xc_max, q_lcl, ρ, N_lcl)
    dn_lcl_dec = CM2.number_decrease_for_mass_limit(sb.numadj, sb.pdf_c.xc_min, q_lcl, ρ, N_lcl)
    dn_lcl_dt += (dn_lcl_inc + dn_lcl_dec) / ρ

    # Rain
    dn_rai_inc = CM2.number_increase_for_mass_limit(sb.numadj, sb.pdf_r.xr_max, q_rai, ρ, N_rai)
    dn_rai_dec = CM2.number_decrease_for_mass_limit(sb.numadj, sb.pdf_r.xr_min, q_rai, ρ, N_rai)
    dn_rai_dt += (dn_rai_inc + dn_rai_dec) / ρ

    return (; dq_lcl_dt, dq_rai_dt, dn_lcl_dt, dn_rai_dt)
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

    # --- Core Warm Rain Processes (shared helper) ---
    warm = warm_rain_tendencies_2m(mp.warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt, dq_ice_dt, dq_rim_dt, db_rim_dt)
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
    ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai,
    q_ice = zero(ρ), n_ice = zero(ρ), q_rim = zero(ρ), b_rim = zero(ρ), logλ = zero(ρ),
) where {WR, ICE <: CMP.P3IceParams}
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

    # Unpack warm rain parameters (always present)
    aps = mp.warm_rain.air_properties

    # Initialize ice-related tendencies
    dq_ice_dt = zero(ρ)
    # TODO: When ice number concentration becomes prognostic, add:
    # dn_ice_dt = zero(ρ)  # Ice number tendency (changes due to melting, aggregation)
    dq_rim_dt = zero(ρ)
    db_rim_dt = zero(ρ)

    # --- Core Warm Rain Processes (shared helper) ---
    warm = warm_rain_tendencies_2m(mp.warm_rain, tps, T, q_tot, q_lcl, q_rai, q_ice, ρ, n_lcl, n_rai)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt

    # Convert to number densities for remaining functions
    N_lcl = ρ * n_lcl
    N_rai = ρ * n_rai

    # --- P3 Ice Processes ---
    p3 = mp.ice.scheme
    vel = mp.ice.terminal_velocity
    pdf_c = mp.ice.cloud_pdf
    pdf_r = mp.ice.rain_pdf

    # Only compute ice processes if there is ice mass/number present
    if (q_ice > zero(q_ice) || n_ice > zero(n_ice))
        # Convert to volumetric quantities for P3 functions
        L_ice = q_ice * ρ  # [kg/m³]
        N_ice = n_ice * ρ  # [1/m³]
        L_lcl = q_lcl * ρ  # [kg/m³]
        N_lcl = n_lcl * ρ  # [1/m³]
        L_rai = q_rai * ρ  # [kg/m³]
        N_rai = n_rai * ρ  # [1/m³]

        # Compute rime fraction and density
        F_rim = ifelse(q_ice > zero(q_ice), q_rim / q_ice, zero(q_rim))
        ρ_rim = ifelse(b_rim > zero(b_rim), q_rim * ρ / (b_rim * ρ), ρ)  # [kg/m³]

        # Liquid-ice collision sources (core P3 process)
        # Only compute if there is ice present
        if L_ice > zero(L_ice) && N_ice > zero(N_ice)
            coll = CMP3.bulk_liquid_ice_collision_sources(
                p3,
                logλ,
                L_ice,
                N_ice,
                F_rim,
                ρ_rim,
                pdf_c,
                pdf_r,
                L_lcl,
                N_lcl,
                L_rai,
                N_rai,
                aps,
                tps,
                vel,
                ρ,
                T,
            )

            # Add collision tendencies
            dq_lcl_dt += coll.∂ₜq_c
            dq_rai_dt += coll.∂ₜq_r
            dn_lcl_dt += coll.∂ₜN_c / ρ
            dn_rai_dt += coll.∂ₜN_r / ρ
            dq_ice_dt += coll.∂ₜL_ice / ρ
            # TODO: When P3 collision sources return ∂ₜN_ice (aggregation, etc.), add:
            # dn_ice_dt += coll.∂ₜN_ice / ρ
            dq_rim_dt += coll.∂ₜL_rim / ρ
            db_rim_dt += coll.∂ₜB_rim / ρ
        end

        # Ice melting (above freezing temperature)
        T_freeze = TDI.TD.Parameters.T_freeze(tps)
        if T > T_freeze && L_ice > zero(L_ice)
            state = CMP3.P3State(p3, L_ice, N_ice, F_rim, ρ_rim)
            # TODO: Using a function that takes dt as an argument is not compatible with the current API.
            # We should use a function that doesn't take dt as an argument.
            dt_dummy = 1000 * one(T)  # P3 uses dt for limiting, we'll limit later
            melt = CMP3.ice_melt(vel, aps, tps, T, ρ, dt_dummy, state, logλ)

            # Melting converts ice to rain
            dq_ice_dt -= melt.dLdt / ρ
            dq_rai_dt += melt.dLdt / ρ
            # TODO: When ice number concentration is tracked, add:
            # dn_ice_dt -= melt.dNdt / ρ  # Ice particles consumed by melting
            dn_rai_dt += melt.dNdt / ρ  # Melted ice becomes rain drops
        end
    end

    # TODO: When ice number concentration is tracked, add dn_ice_dt to return tuple:
    # return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt, dq_ice_dt, dn_ice_dt, dq_rim_dt, db_rim_dt)
    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt, dq_ice_dt, dq_rim_dt, db_rim_dt)
end

end # module BulkMicrophysicsTendencies
