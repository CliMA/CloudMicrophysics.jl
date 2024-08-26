"""
    Non-equilibrium bulk microphysics scheme, which includes:

  - condensation and evaporation of cloud liquid water and
    deposition and sublimation of cloud ice (relaxation to equilibrium)
"""
module MicrophysicsNonEq

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

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
    conv_q_vap_to_q_liq_ice(liquid, q_sat, q)
    conv_q_vap_to_q_liq_ice(ice, q_sat, q)

 - `liquid` or `ice` - a struct with cloud water or ice free parameters
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale and is based on the difference between
specific humidity in equilibrium at the current temperature
and the current cloud condensate.
"""
function conv_q_vap_to_q_liq_ice(
    (; τ_relax)::CMP.CloudLiquid{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT}
    return (q_sat.liq - q.liq) / τ_relax
end
function conv_q_vap_to_q_liq_ice(
    (; τ_relax)::CMP.CloudIce{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
) where {FT}
    return (q_sat.ice - q.ice) / τ_relax
end

"""
    conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, q_sat, q, T)
    conv_q_vap_to_q_liq_ice_MM2015(ice, tps, q_sat, q, T)

- `liquid` or `ice` - a struct with cloud water or ice free parameters
- `tps` - thermodynamics parameters struct
- `q_sat` - PhasePartition of the q_vap values at saturation for liquid and ice
- `q` - current PhasePartition
- `ρ` - air density [kg/m3]
- `T` - air temperature [K]

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The formulation is based on Morrison and Grabowski 2008 and
Morrison and Milbrandt 2015
"""
function conv_q_vap_to_q_liq_ice_MM2015(
    (; τ_relax)::CMP.CloudLiquid{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
    T::FT,
) where {FT}
    Rᵥ = TD.Parameters.R_v(tps)
    cₚ_air = TD.cp_m(tps, q)
    Lᵥ = TD.latent_heat_vapor(tps, T)
    qᵥ = TD.vapor_specific_humidity(q)

    dqsldT = q_sat.liq * (Lᵥ / (Rᵥ * T^2) - 1 / T)
    Γₗ = FT(1) + (Lᵥ / cₚ_air) * dqsldT

    return (qᵥ - q_sat.liq) / (τ_relax * Γₗ)
end
function conv_q_vap_to_q_liq_ice_MM2015(
    (; τ_relax)::CMP.CloudIce{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
    T::FT,
) where {FT}
    Rᵥ = TD.Parameters.R_v(tps)
    cₚ_air = TD.cp_m(tps, q)
    Lₛ = TD.latent_heat_sublim(tps, T)
    qᵥ = TD.vapor_specific_humidity(q)

    dqsidT = q_sat.ice * (Lₛ / (Rᵥ * T^2) - 1 / T)
    Γᵢ = FT(1) + (Lₛ / cₚ_air) * dqsidT

    return (qᵥ - q_sat.ice) / (τ_relax * Γᵢ)
end

end #module MicrophysicsNonEq.jl
