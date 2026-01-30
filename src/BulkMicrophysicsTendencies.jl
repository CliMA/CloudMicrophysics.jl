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
import ..Common as CO
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

# --- 1-Moment Microphysics ---

"""
    bulk_microphysics_tendencies(
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

# Notes
- This function does NOT apply limiters to prevent negative concentrations.
  Limiters depend on timestep `dt` and should be applied by the caller.
"""
@inline function bulk_microphysics_tendencies(
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

    # Deposition/sublimation of cloud ice
    S_icl_dep = CMNonEq.conv_q_vap_to_q_lcl_icl_MM2015(icl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    # No ice deposition above freezing (lack of INPs) - use ifelse to avoid GPU branching
    S_icl_dep = ifelse(T > tps.T_freeze, min(S_icl_dep, zero(T)), S_icl_dep)
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
    S_accr_melt = α * S_accr_sno_rai
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

# --- 0-Moment Microphysics ---

"""
    bulk_microphysics_tendencies(
        ::Microphysics0Moment,
        mp,
        tps,
        T,
        q_lcl,
        q_icl,
    )

Compute 0-moment microphysics tendencies in one fused call.

Returns a NamedTuple with:
- `dq_tot_dt`: Total water tendency from precipitation removal [kg/kg/s]
- `e_int_precip`: Internal energy of removed condensate [J/kg]

Caller adds geopotential Φ for energy tendency: `e_tot = dq_tot_dt * (e_int_precip + Φ)`

# Arguments
- `mp`: NamedTuple with microphysics parameters:
  - `params_0M`: Parameters0M struct (contains τ_precip, qc_0 or S_0)
- `tps`: Thermodynamics parameters
- `T`: Temperature [K]
- `q_lcl`: Cloud liquid specific content [kg/kg]
- `q_icl`: Cloud ice specific content [kg/kg]

# Notes
- Does NOT apply limiters (caller applies based on timestep)
- Does NOT include geopotential (caller adds Φ for energy tendency)
"""
@inline function bulk_microphysics_tendencies(
    ::Microphysics0Moment,
    mp::CMP.Microphysics0MParams,
    tps,
    T,
    q_lcl,
    q_icl,
)
    # Unpack microphysics parameters
    params_0M = mp.precip

    # Precipitation removal rate
    dq_tot_dt = CM0.remove_precipitation(params_0M, q_lcl, q_icl)

    # Internal energy of removed condensate (liquid fraction weighted)
    λ = TDI.liquid_fraction(tps, T, q_lcl, q_icl)
    I_liq = TDI.internal_energy_liquid(tps, T)
    I_ice = TDI.internal_energy_ice(tps, T)
    e_int_precip = λ * I_liq + (1 - λ) * I_ice

    return (; dq_tot_dt, e_int_precip)
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
@inline function warm_rain_tendencies_2m(sb, q_lcl, q_rai, ρ, n_lcl, n_rai)
    # Convert to number densities for CM2 functions
    N_lcl = ρ * n_lcl
    N_rai = ρ * n_rai

    # Initialize tendencies
    FT = typeof(ρ)
    dq_lcl_dt = zero(FT)
    dq_rai_dt = zero(FT)
    dn_lcl_dt = zero(FT)
    dn_rai_dt = zero(FT)

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

    return (; dq_lcl_dt, dq_rai_dt, dn_lcl_dt, dn_rai_dt)
end

# --- 2-Moment Microphysics (Unified Warm + Optional Ice) ---


"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{FT, WR, Nothing},
        ...
    )

Compute 2-moment **warm rain only** microphysics tendencies (Seifert-Beheng 2006).

This method is type-stable and GPU-optimized for warm rain processes only.
For warm rain + P3 ice, see the method that accepts `Microphysics2MParams{FT, WR, <:P3IceParams}`.

# Arguments
- `mp`: Microphysics2MParams with `mp.ice == nothing` (warm rain only)
- `tps`: Thermodynamics parameters
- `ρ`: Air density (kg/m³)
- `T`: Temperature (K)
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
    ::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{FT, WR, Nothing},
    tps,
    ρ,
    T,
    q_lcl,
    n_lcl,
    q_rai,
    n_rai,
    q_ice = zero(ρ),
    n_ice = zero(ρ),
    q_rim = zero(ρ),
    b_rim = zero(ρ),
    logλ = zero(ρ),
) where {FT, WR}
    # Unpack warm rain parameters
    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties

    # Initialize ice-related tendencies (always zero for warm-only)
    dq_ice_dt = zero(typeof(ρ))
    dq_rim_dt = zero(typeof(ρ))
    db_rim_dt = zero(typeof(ρ))

    # --- Core Warm Rain Processes (shared helper) ---
    warm = warm_rain_tendencies_2m(sb, q_lcl, q_rai, ρ, n_lcl, n_rai)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt

    # Convert to number densities for remaining functions
    N_lcl = ρ * n_lcl
    N_rai = ρ * n_rai


    # --- Rain evaporation ---
    evap = CM2.rain_evaporation(sb, aps, tps, zero(ρ), q_lcl, zero(ρ), q_rai, zero(ρ), ρ, N_rai, T)
    dq_rai_dt += evap.evap_rate_1
    dn_rai_dt += evap.evap_rate_0 / ρ

    # --- Number adjustment for mass limits ---
    # Cloud liquid
    dn_lcl_inc = CM2.number_increase_for_mass_limit(sb.numadj, sb.pdf_c.xc_max, q_lcl, ρ, N_lcl)
    dn_lcl_dec = CM2.number_decrease_for_mass_limit(sb.numadj, sb.pdf_c.xc_min, q_lcl, ρ, N_lcl)
    dn_lcl_dt += (dn_lcl_inc + dn_lcl_dec) / ρ

    # Rain
    dn_rai_inc = CM2.number_increase_for_mass_limit(sb.numadj, sb.pdf_r.xr_max, q_rai, ρ, N_rai)
    dn_rai_dec = CM2.number_decrease_for_mass_limit(sb.numadj, sb.pdf_r.xr_min, q_rai, ρ, N_rai)
    dn_rai_dt += (dn_rai_inc + dn_rai_dec) / ρ

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt, dq_ice_dt, dq_rim_dt, db_rim_dt)
end

"""
    bulk_microphysics_tendencies(
        ::Microphysics2Moment,
        mp::Microphysics2MParams{FT, WR, <:P3IceParams},
        ...
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
    ::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{FT, WR, ICE},
    tps,
    ρ,
    T,
    q_lcl,
    n_lcl,
    q_rai,
    n_rai,
    q_ice = zero(ρ),
    n_ice = zero(ρ),
    q_rim = zero(ρ),
    b_rim = zero(ρ),
    logλ = zero(ρ),
) where {FT, WR, ICE <: CMP.P3IceParams}
    # Unpack warm rain parameters (always present)
    sb = mp.warm_rain.seifert_beheng
    aps = mp.warm_rain.air_properties

    # Initialize ice-related tendencies
    dq_ice_dt = zero(typeof(ρ))
    # TODO: When ice number concentration becomes prognostic, add:
    # dn_ice_dt = zero(typeof(ρ))  # Ice number tendency (changes due to melting, aggregation)
    dq_rim_dt = zero(typeof(ρ))
    db_rim_dt = zero(typeof(ρ))

    # --- Core Warm Rain Processes (shared helper) ---
    warm = warm_rain_tendencies_2m(sb, q_lcl, q_rai, ρ, n_lcl, n_rai)
    dq_lcl_dt = warm.dq_lcl_dt
    dn_lcl_dt = warm.dn_lcl_dt
    dq_rai_dt = warm.dq_rai_dt
    dn_rai_dt = warm.dn_rai_dt

    # Convert to number densities for remaining functions
    N_lcl = ρ * n_lcl
    N_rai = ρ * n_rai

    # --- Rain evaporation ---
    evap = CM2.rain_evaporation(sb, aps, tps, zero(ρ), q_lcl, zero(ρ), q_rai, zero(ρ), ρ, N_rai, T)
    dq_rai_dt += evap.evap_rate_1
    dn_rai_dt += evap.evap_rate_0 / ρ

    # --- Number adjustment for mass limits ---
    # Cloud liquid
    dn_lcl_inc = CM2.number_increase_for_mass_limit(sb.numadj, sb.pdf_c.xc_max, q_lcl, ρ, N_lcl)
    dn_lcl_dec = CM2.number_decrease_for_mass_limit(sb.numadj, sb.pdf_c.xc_min, q_lcl, ρ, N_lcl)
    dn_lcl_dt += (dn_lcl_inc + dn_lcl_dec) / ρ

    # Rain
    dn_rai_inc = CM2.number_increase_for_mass_limit(sb.numadj, sb.pdf_r.xr_max, q_rai, ρ, N_rai)
    dn_rai_dec = CM2.number_decrease_for_mass_limit(sb.numadj, sb.pdf_r.xr_min, q_rai, ρ, N_rai)
    dn_rai_dt += (dn_rai_inc + dn_rai_dec) / ρ

    # --- P3 Ice Processes ---
    # Ice parameters are guaranteed to be present (enforced by type signature)
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
