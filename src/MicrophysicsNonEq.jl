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

 - `liquid` or `ice` - a type for cloud water or ice
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.
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

function conv_q_vap_to_q_liq_ice(
    tps::TDP.ThermodynamicsParameters{FT},
    liquid::CMP.CloudLiquid{FT},
    ice::CMP.CloudIce{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
    T:: FT, # temperature
    Sₗ::FT, # liquid saturation ratio
    const_dt:: FT,
) where {FT}

    # (might want to change the name of this function at some point)
    # implementing this -- first the simplest form that
    # uses the assumption that dqs/dT = 1 and sets A_c = 1
    # and currently this is in specific humidity, but it technically should
    # be converted to mixing ratio

    dqsdT = FT(1)
    A_c = FT(1) # i need to actually make this into something

    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_v = TD.latent_heat_vapor(tps,T)

    Γₗ = 1 + (L_v/cp_air)*dqsdT
    Γᵢ = 1 + (L_subl/cp_air)*dqsdT

    τ = (liquid.τ_relax^(-1) + (1 + (L_subl/cp_air))*ice.τ_relax^(-1) / Γᵢ)^(-1)

    # this does need to be mixing ratio instead
    delta_0 =  (Sₗ-1)*q_sat.liq

    # solving for new Sl after timestep delta t:
    Sₗ = (1/q_sat.liq)(A_c * τ/(ice.τ_relax * Γₗ) + (delta_0 - A_c*τ)
        *τ/(const_dt*ice.τ_relax*Γₗ)*(FT(1) - exp( - const_dt/τ))) + 1

    return Sₗ
end

end #module MicrophysicsNonEq.jl