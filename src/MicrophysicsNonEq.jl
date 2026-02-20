"""
Non-equilibrium bulk microphysics scheme for cloud condensate formation.

Implements relaxation-to-equilibrium approach for:
- Condensation and evaporation of cloud liquid water
- Deposition and sublimation of cloud ice

See also: `Microphysics1M` for precipitating hydrometeor processes.
"""
module MicrophysicsNonEq

import ..Parameters as CMP
import ..ThermodynamicsInterface as TDI
import ..Common as CO
import ..Utilities as UT

export τ_relax
export conv_q_vap_to_q_lcl_icl
export conv_q_vap_to_q_lcl_icl_MM2015
export ∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld
export INP_limiter
export limit_MM2015_sinks
export dqcld_dT
export gamma_helper
export d2qcld_dT2

"""
    τ_relax(liquid::CloudLiquid)
    τ_relax(ice::CloudIce)

Returns the relaxation timescale for phase change processes.

# Arguments
- `liquid::CloudLiquid` or `ice::CloudIce` - cloud condensate parameters struct

# Returns
- Relaxation timescale [s] for condensation/evaporation (liquid) or
  deposition/sublimation (ice)
"""
@inline τ_relax(p::CMP.CloudLiquid) = p.τ_relax
@inline τ_relax(p::CMP.CloudIce) = p.τ_relax

"""
    conv_q_vap_to_q_lcl_icl(params::CloudLiquid, q_sat_liq, q_lcl)
    conv_q_vap_to_q_lcl_icl(params::CloudIce, q_sat_ice, q_icl)

Computes the tendency of cloud condensate specific content using a
simple relaxation-to-equilibrium formulation with a constant timescale.

# Arguments
- `params` - cloud liquid or ice parameters struct containing `τ_relax`
- `q_sat_liq` or `q_sat_ice` - saturation specific humidity [kg/kg]
- `q_lcl` or `q_icl` - cloud liquid water or ice specific content [kg/kg]

# Returns
- Cloud condensate tendency in kg/kg/s, positive for condensation/deposition,
  negative for evaporation/sublimation

The tendency is computed as:
```math
\\frac{dq}{dt} = \\frac{q_{sat} - q}{\\tau_{relax}}
```
"""
@inline function conv_q_vap_to_q_lcl_icl((; τ_relax)::CMP.CloudLiquid, q_sat_liq, q_lcl)
    return (q_sat_liq - q_lcl) / τ_relax
end

@inline function conv_q_vap_to_q_lcl_icl((; τ_relax)::CMP.CloudIce, q_sat_ice, q_icl)
    return (q_sat_ice - q_icl) / τ_relax
end

"""
    limit_MM2015_sinks(tendency, q_cld)

Checks if the cloud evaporation or sublimation tendency need to be
limited to prevent negative cloud condensate specific content.
The tendency is limited if it is negative and the cloud condensate specific
content is zero or negative.

# Arguments
- `tendency` - cloud condensate tendency [kg/kg/s]
- `q_cld` - cloud liquid water or ice specific content [kg/kg]

# Returns
- `true` if the tendency needs to be limited, `false` otherwise
"""
@inline function limit_MM2015_sinks(tendency, q_cld)
    return tendency <= zero(tendency) && q_cld <= zero(q_cld)
end

"""
    INP_limiter(tendency, tps, T)

Returns `true` when ice deposition should be suppressed:
positive tendency (deposition) at T > T_freeze (no INPs available).
"""
@inline function INP_limiter(tendency, tps, T)
    return T > TDI.T_freeze(tps) && tendency > zero(tendency)
end

"""
    dqcld_dT(qᵥ_sat, L, Rᵥ, T)

Computes the derivative of the saturation specific humidity with respect to
temperature for a given phase of water.

# Arguments
- `qᵥ_sat` - saturation specific humidity [kg/kg]
- `L` - latent heat [J/kg]
- `Rᵥ` - gas constant for water vapor [J/kg/K]
- `T` - temperature [K]
"""
@inline function dqcld_dT(qᵥ_sat, L, Rᵥ, T)
    return qᵥ_sat * (L / (Rᵥ * T^2) - 1 / T)
end

"""
    gamma_helper(L, cₚ_air, dqcld_dT)

Computes the thermodynamic adjustment factor Γ.

# Arguments
- `L` - latent heat [J/kg]
- `cₚ_air` - specific heat capacity of air [J/kg/K]
- `dqcld_dT` - derivative of saturation specific humidity with respect to temperature [kg/kg/K]
"""
@inline function gamma_helper(L, cₚ_air, dqcld_dT)
    return 1 + (L / cₚ_air) * dqcld_dT
end

"""
    d2qcld_dT2(qᵥ_sat, L, Rᵥ, T)

Computes the second derivative of the saturation specific humidity with respect to
temperature for a given phase of water.

# Arguments
- `qᵥ_sat` - saturation specific humidity [kg/kg]
- `L` - latent heat [J/kg]
- `Rᵥ` - gas constant for water vapor [J/kg/K]
- `T` - temperature [K]
"""
@inline function d2qcld_dT2(qᵥ_sat, L, Rᵥ, T)
    return qᵥ_sat * ((L / Rᵥ / T^2 - 1 / T)^2 + (1 / T^2 - 2 * L / Rᵥ / T^3))
end

"""
    conv_q_vap_to_q_lcl_icl_MM2015(params, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)

Computes cloud condensate tendency using the formulation from
Morrison & Grabowski (2008), https://doi.org/10.1175/2007JAS2374.1, and
Morrison & Milbrandt (2015), https://doi.org/10.1175/JAS-D-14-0065.1.

This formulation includes a thermodynamic adjustment factor Γ that
accounts for latent heat release modifying the saturation state.

# Arguments
- `params` - cloud liquid or ice parameters struct containing `τ_relax`
- `tps` - thermodynamics parameters struct
- `q_tot` - total water specific content [kg/kg]
- `q_lcl` - cloud liquid water specific content [kg/kg]
- `q_icl` - cloud ice specific content [kg/kg]
- `q_rai` - rain specific content [kg/kg]
- `q_sno` - snow specific content [kg/kg]
- `ρ` - air density [kg/m³]
- `T` - air temperature [K]

# Returns
- Cloud condensate tendency [kg/kg/s]

# Notes
This function does NOT apply limiters for small or negative specific humidities.
Users should apply appropriate bounds checking when integrating in a model.
"""
function conv_q_vap_to_q_lcl_icl_MM2015(
    (; τ_relax)::CMP.CloudLiquid,
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T,
)
    Rᵥ = TDI.Rᵥ(tps)
    Lᵥ = TDI.Lᵥ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)
    qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

    dqsl_dT = dqcld_dT(qᵥ_sat_liq, Lᵥ, Rᵥ, T)
    Γₗ = gamma_helper(Lᵥ, cₚ_air, dqsl_dT)

    # compute the tendency
    tendency = (qᵥ - qᵥ_sat_liq) / (τ_relax * Γₗ)

    return ifelse(limit_MM2015_sinks(tendency, q_lcl), zero(tendency), tendency)
end
function conv_q_vap_to_q_lcl_icl_MM2015(
    (; τ_relax)::CMP.CloudIce,
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T,
)
    Rᵥ = TDI.Rᵥ(tps)
    Lₛ = TDI.Lₛ(tps, T)
    cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)

    qᵥ = TDI.q_vap(q_tot, q_lcl + q_rai, q_icl + q_sno)
    qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)

    dqsi_dT = dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
    Γᵢ = gamma_helper(Lₛ, cₚ_air, dqsi_dT)

    # compute the tendency
    tendency = (qᵥ - qᵥ_sat_ice) / (τ_relax * Γᵢ)

    limiter = limit_MM2015_sinks(tendency, q_icl) || INP_limiter(tendency, tps, T)
    return ifelse(limiter, zero(tendency), tendency)
end

"""
    ∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(params::CloudLiquid, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T; simplified = true)
    ∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(params::CloudIce,    tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T; simplified = true)

Returns the derivative of the cloud condensate tendency with respect to the
corresponding cloud condensate species:
 - `∂(tendency)/∂q_lcl` for the `CloudLiquid` dispatch
 - `∂(tendency)/∂q_icl` for the `CloudIce` dispatch

# Keyword Arguments
- `simplified` — if `true` (default), returns the leading-order approximation (`-1/τ_relax`);
  if `false`, returns the total derivative accounting for implicit temperature feedback.
"""
function ∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
    lcl::CMP.CloudLiquid,
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T;
    simplified = true,
)
    tendency = conv_q_vap_to_q_lcl_icl_MM2015(lcl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    ∂tendency_∂q_lcl = -1 / lcl.τ_relax
    if !simplified
        Rᵥ = TDI.Rᵥ(tps)
        Lᵥ = TDI.Lᵥ(tps, T)
        cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
        qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)

        dqsl_dT = dqcld_dT(qᵥ_sat_liq, Lᵥ, Rᵥ, T)
        Γₗ = gamma_helper(Lᵥ, cₚ_air, dqsl_dT)
        d²qᵥ_sat_liq_dT² = d2qcld_dT2(qᵥ_sat_liq, Lᵥ, Rᵥ, T)

        ∂tendency_∂q_lcl -= tendency / Γₗ * (Lᵥ / cₚ_air)^2 * d²qᵥ_sat_liq_dT²
    end

    return ifelse(limit_MM2015_sinks(tendency, q_lcl), zero(∂tendency_∂q_lcl), ∂tendency_∂q_lcl)
end
function ∂conv_q_vap_to_q_lcl_icl_MM2015_∂q_cld(
    icl::CMP.CloudIce,
    tps::TDI.PS,
    q_tot,
    q_lcl,
    q_icl,
    q_rai,
    q_sno,
    ρ,
    T;
    simplified = true,
)
    tendency = conv_q_vap_to_q_lcl_icl_MM2015(icl, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    ∂tendency_∂q_icl = -1 / icl.τ_relax
    if !simplified
        Rᵥ = TDI.Rᵥ(tps)
        Lₛ = TDI.Lₛ(tps, T)
        cₚ_air = TDI.cpₘ(tps, q_tot, q_lcl + q_rai, q_icl + q_sno)
        qᵥ_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)

        dqsi_dT = dqcld_dT(qᵥ_sat_ice, Lₛ, Rᵥ, T)
        Γᵢ = gamma_helper(Lₛ, cₚ_air, dqsi_dT)
        d²qᵥ_sat_ice_dT² = d2qcld_dT2(qᵥ_sat_ice, Lₛ, Rᵥ, T)

        ∂tendency_∂q_icl -= tendency / Γᵢ * (Lₛ / cₚ_air)^2 * d²qᵥ_sat_ice_dT²
    end

    limiter = limit_MM2015_sinks(tendency, q_icl) || INP_limiter(tendency, tps, T)
    return ifelse(limiter, zero(∂tendency_∂q_icl), ∂tendency_∂q_icl)
end

"""
    terminal_velocity(sediment::CloudLiquid, vel::StokesRegimeVelType, ρₐ, q)
    terminal_velocity(sediment::CloudIce, vel::Chen2022VelTypeSmallIce, ρₐ, q)

Computes mass-weighted average terminal velocity for cloud droplets or ice crystals
assuming a monodisperse size distribution.

- **Cloud Liquid**: Uses Stokes Law (v ∝ D²), valid for small droplets (Re < 1).
- **Cloud Ice**: Uses Chen et al. (2022) parameterization,
  [DOI: 10.1016/j.atmosres.2022.106171](https://doi.org/10.1016/j.atmosres.2022.106171)

# Arguments
- `sediment` - cloud liquid or ice parameters struct (provides density and N_0)
- `vel` - velocity parameters (StokesRegimeVelType for liquid, Chen2022VelTypeSmallIce for ice)
- `ρₐ` - air density [kg/m³]
- `q` - cloud liquid water or ice specific content [kg/kg]

# Returns
- Mass-weighted terminal velocity [m/s]
"""
function terminal_velocity(
    (; ρw, N_0)::CMP.CloudLiquid{FT},
    vel::CMP.StokesRegimeVelType{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > UT.ϵ_numerics(FT)
        # Stokes law: v(D) = C * D^2, valid for D < ~80 μm (Re < 1)
        v_term = CO.particle_terminal_velocity(vel, ρₐ)
        # Mean volume diameter from assumed number concentration
        D = cbrt(FT(6 / π) * ρₐ * q / N_0 / ρw)
        fall_w = v_term(D)
    end
    return fall_w
end

function terminal_velocity(
    (; ρᵢ, N_0)::CMP.CloudIce{FT},
    vel::CMP.Chen2022VelTypeSmallIce{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > UT.ϵ_numerics(FT)
        v_term = CO.particle_terminal_velocity(vel, ρₐ, ρᵢ)
        # Mean volume diameter from assumed number concentration
        D = cbrt(FT(6 / π) * ρₐ * q / N_0 / ρᵢ)
        fall_w = max(FT(0), v_term(D))
    end
    return fall_w
end

end #module MicrophysicsNonEq.jl
