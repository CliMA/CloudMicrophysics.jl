"""
    Non-equilibrium bulk microphysics scheme, which includes:

  - condensation and evaporation of cloud liquid water and
    deposition and sublimation of cloud ice (relaxation to equilibrium)
"""
module MicrophysicsNonEq

import Thermodynamics as TD

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
    (; τₗ)::CMP.CloudLiquid{FT},
    (; τᵢ)::CMP.CloudIce{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
    T:: FT, # temperature
    const_dt:: FT,
) where {FT}

    # implementing this -- first the simplest form that drops the
    # derivatives and sets A_c = 1

    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)

    Γₗ = 1 + (L_subl/cp_air)
    Γᵢ = Γₗ

    τ = τₗ + (1 + (L_subl/cp_air))*τᵢ^(-1) / Γᵢ

    # replace this w QCCONN -- see what Jordan did in his milbrandt code

    A_c = 1
    delta_0 =  q_sat # ??? this is maybe the part that scares me most
    # and yeah also no clue what tau c is

    # something like this -- need to get everything into the right format etc etc.
    result = A_c * τ/(τᵢ * Γₗ) + (delta_0 - A_c*τ)*τ/(const_dt*τᵢ*Γₗ)(1-np.exp(-const_dt/τ))

    # beware that this version is mixing ratio not specific humidity
    # so may need to convert back to that at some point

    return result
end

end #module MicrophysicsNonEq.jl