"""
Non-equilibrium bulk microphysics scheme for cloud condensate formation.

Implements relaxation-to-equilibrium approach for:
- Condensation and evaporation of cloud liquid water
- Deposition and sublimation of cloud ice

See also: [`Microphysics1M`](@ref) for precipitating hydrometeor processes.
"""
module MicrophysicsNonEq

import ..Parameters as CMP
import ..ThermodynamicsInterface as TDI
import ...Common as CO

export τ_relax
export conv_q_vap_to_q_lcl_icl
export conv_q_vap_to_q_lcl_icl_MM2015

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
- Cloud condensate tendency [kg/kg/s] (positive for condensation/deposition,
  negative for evaporation/sublimation)

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

    pᵥ_sat_liq = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    qᵥ_sat_liq = TDI.p2q(tps, T, ρ, pᵥ_sat_liq)

    dqsldT = qᵥ_sat_liq * (Lᵥ / (Rᵥ * T^2) - 1 / T)
    Γₗ = 1 + (Lᵥ / cₚ_air) * dqsldT

    return (qᵥ - qᵥ_sat_liq) / (τ_relax * Γₗ)
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

    pᵥ_sat_ice = TDI.saturation_vapor_pressure_over_ice(tps, T)
    qᵥ_sat_ice = TDI.p2q(tps, T, ρ, pᵥ_sat_ice)

    dqsidT = qᵥ_sat_ice * (Lₛ / (Rᵥ * T^2) - 1 / T)
    Γᵢ = 1 + (Lₛ / cₚ_air) * dqsidT

    return (qᵥ - qᵥ_sat_ice) / (τ_relax * Γᵢ)
end

"""
    terminal_velocity(sediment::CloudLiquid, vel::Chen2022VelTypeRain, ρₐ, q)
    terminal_velocity(sediment::CloudIce, vel::Chen2022VelTypeSmallIce, ρₐ, q)

Computes mass-weighted average terminal velocity for cloud droplets or ice crystals
assuming a monodisperse size distribution with prescribed number concentration.

Uses Chen et al. (2022) terminal velocity parameterization.
See [DOI: 10.1016/j.atmosres.2022.106171](https://doi.org/10.1016/j.atmosres.2022.106171)

# Arguments
- `sediment` - cloud liquid or ice parameters struct (provides density)
- `vel` - Chen2022 terminal velocity parameters struct
- `ρₐ` - air density [kg/m³]
- `q` - cloud liquid water or ice specific content [kg/kg]

# Returns
- Mass-weighted terminal velocity [m/s]

# Notes
The Chen et al. (2022) coefficients are calibrated for D > 100 μm and may
not be accurate for smaller cloud droplets. A correction factor is applied
for cloud liquid. Number concentration is assumed fixed at 500 × 10⁶ m⁻³.
"""
function terminal_velocity(
    (; ρw)::CMP.CloudLiquid{FT},
    vel::CMP.Chen2022VelTypeRain{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > CO.ϵ_numerics(FT)
        # TODO: Coefficients from Table B1 from Chen et. al. 2022 are only valid
        # for D > 100mm. We should look for a different parameterization
        # that is more suited for cloud droplets. For now I'm just multiplying
        # by an arbitrary correction factor.
        v_term = CO.particle_terminal_velocity(vel, ρₐ)
        # The 1M scheme does not assume any cloud droplet size distribution.
        # TODO - For now I compute a mean volume radius assuming a fixed value
        # for the total number concentration of droplets.
        N = FT(500 * 1e6)
        D = cbrt(ρₐ * q / N / ρw)
        corr = FT(0.1)
        # assuming spherical particles (ϕ = 1)
        fall_w = max(FT(0), corr * v_term(D))
    end
    return fall_w
end

function terminal_velocity(
    (; ρᵢ)::CMP.CloudIce{FT},
    vel::CMP.Chen2022VelTypeSmallIce{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > CO.ϵ_numerics(FT)
        v_term = CO.particle_terminal_velocity(vel, ρₐ, ρᵢ)
        # Assume fixed ice crystal number concentration (see comment for liquid above)
        N = FT(500 * 1e6)
        D = cbrt(ρₐ * q / N / ρᵢ)
        # assuming spherical particles (ϕ = 1)
        fall_w = max(FT(0), v_term(D))
    end
    return fall_w
end

end #module MicrophysicsNonEq.jl

