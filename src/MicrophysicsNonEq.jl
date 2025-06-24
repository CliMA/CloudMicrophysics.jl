"""
    Non-equilibrium bulk microphysics scheme, which includes:

  - condensation and evaporation of cloud liquid water and
    deposition and sublimation of cloud ice (relaxation to equilibrium)
"""
module MicrophysicsNonEq

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import CloudMicrophysics.Common as CO

import ..Parameters as CMP

export τ_relax
export conv_q_vap_to_q_liq_ice
export conv_q_vap_to_q_liq_ice_MM2015

"""
    τ_relax(liquid)
    τ_relax(ice)

 - `liquid` or `ice` - a type for cloud liquid water or ice

Returns the relaxation timescale for condensation and evaporation of
cloud liquid water or the relaxation timescale for sublimation and
deposition of cloud ice.
"""
τ_relax(p::CMP.CloudLiquid) = p.τ_relax
τ_relax(p::CMP.CloudIce) = p.τ_relax

"""
    conv_q_vap_to_q_liq_ice(params, q_sat, q_)

 - `params` - a struct with cloud water or ice free parameters
 - `q_sat` - liquid or ice equilibrium specific contents
 - `q_`` - current liquid or ice specific contants

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale and is based on the difference between
specific contents in equilibrium at the current temperature
and the current cloud condensate.
"""
function conv_q_vap_to_q_liq_ice((; τ_relax)::CMP.CloudLiquid, q_sat_liq, q_liq)
    return (q_sat_liq - q_liq) / τ_relax
end
function conv_q_vap_to_q_liq_ice((; τ_relax)::CMP.CloudIce, q_sat_ice, q_ice)
    return (q_sat_ice - q_ice) / τ_relax
end

"""
    conv_q_vap_to_q_liq_ice_MM2015(params, tps, q_tot, q_liq, q_ice, q_rai, q_sno, ρ, T)
    conv_q_vap_to_q_liq_ice_MM2015(params, tps, q::PhasePartition, q_rai, q_sno, ρ, T)

- `params` - a struct with cloud water or ice free parameters
- `tps` - thermodynamics parameters struct
- `q_tot`, `q_liq`, `q_ice`, `q_rai`, `q_sno` - specific contents of total water, cloud liquid and ice, rain and snow,
- `q` - current PhasePartition
- `ρ` - air density [kg/m3]
- `T` - air temperature [K]

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The formulation is based on Morrison and Grabowski 2008 and
Morrison and Milbrandt 2015. It does NOT screen for small or
negative values for humidities, so we suggest applying a limiter
of choice to this function when running it in a model.
"""
function conv_q_vap_to_q_liq_ice_MM2015(
    (; τ_relax)::CMP.CloudLiquid,
    tps::TDP.ThermodynamicsParameters,
    q_tot,
    q_liq,
    q_ice,
    q_rai,
    q_sno,
    ρ,
    T,
)
    FT = eltype(tps)

    Rᵥ = TD.Parameters.R_v(tps)
    Lᵥ = TD.latent_heat_vapor(tps, T)

    qᵥ = q_tot - q_liq - q_ice - q_rai - q_sno
    cₚ_air = TD.cp_m(tps, q_tot, q_liq + q_rai, q_ice + q_sno)

    pᵥ_sat_liq = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    qᵥ_sat_liq = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sat_liq)

    dqsldT = qᵥ_sat_liq * (Lᵥ / (Rᵥ * T^2) - 1 / T)
    Γₗ = 1 + (Lᵥ / cₚ_air) * dqsldT

    return (qᵥ - qᵥ_sat_liq) / (τ_relax * Γₗ)
end
function conv_q_vap_to_q_liq_ice_MM2015(
    (; τ_relax)::CMP.CloudIce,
    tps::TDP.ThermodynamicsParameters,
    q_tot,
    q_liq,
    q_ice,
    q_rai,
    q_sno,
    ρ,
    T,
)
    FT = eltype(tps)
    Rᵥ = TD.Parameters.R_v(tps)
    Lₛ = TD.latent_heat_sublim(tps, T)

    qᵥ = q_tot - q_liq - q_ice - q_rai - q_sno
    cₚ_air = TD.cp_m(tps, q_tot, q_liq + q_rai, q_ice + q_sno)

    pᵥ_sat_ice = TD.saturation_vapor_pressure(tps, T, TD.Ice())
    qᵥ_sat_ice = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sat_ice)

    dqsidT = qᵥ_sat_ice * (Lₛ / (Rᵥ * T^2) - 1 / T)
    Γᵢ = 1 + (Lₛ / cₚ_air) * dqsidT

    return (qᵥ - qᵥ_sat_ice) / (τ_relax * Γᵢ)
end

conv_q_vap_to_q_liq_ice_MM2015(
    params::CMP.CloudLiquid,
    tps::TDP.ThermodynamicsParameters,
    q::TD.PhasePartition,
    q_rai,
    q_sno,
    ρ,
    T,
) = conv_q_vap_to_q_liq_ice_MM2015(
    params::CMP.CloudLiquid,
    tps::TDP.ThermodynamicsParameters,
    q.tot,
    q.liq - q_rai,
    q.ice - q_sno,
    q_rai,
    q_sno,
    ρ,
    T,
)

conv_q_vap_to_q_liq_ice_MM2015(
    params::CMP.CloudIce,
    tps::TDP.ThermodynamicsParameters,
    q::TD.PhasePartition,
    q_rai,
    q_sno,
    ρ,
    T,
) = conv_q_vap_to_q_liq_ice_MM2015(
    params::CMP.CloudIce,
    tps::TDP.ThermodynamicsParameters,
    q.tot,
    q.liq - q_rai,
    q.ice - q_sno,
    q_rai,
    q_sno,
    ρ,
    T,
)

"""
    terminal_velocity(sediment, vel, ρ, q)

 - `sediment` - a struct with sedimentation type (cloud liquid or ice)
 - `vel` - a struct with terminal velocity parameters
 - `ρₐ` - air density
 - `q` - cloud liquid or ice specific content

Returns the mass weighted average terminal velocity assuming a
monodisperse size distribution with prescribed number concentration.
The fall velocity of individual particles is parameterized following
Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171
"""
function terminal_velocity(
    (; ρw)::CMP.CloudLiquid{FT},
    vel::CMP.Chen2022VelTypeRain{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > FT(0)
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
        # assuming ϕ = 1 (spherical)
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
    if q > FT(0)
        v_term = CO.particle_terminal_velocity(vel, ρₐ, ρᵢ)
        # See the comment for liquid droplets above
        N = FT(500 * 1e6)
        D = cbrt(ρₐ * q / N / ρᵢ)
        # assuming ϕ = 1 (spherical)
        fall_w = max(FT(0), v_term(D))
    end
    return fall_w
end

end #module MicrophysicsNonEq.jl
