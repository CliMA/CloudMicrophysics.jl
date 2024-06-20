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
    T::FT, # temperature
    Sₗ::FT, # liquid saturation ratio
    w::FT, # upwards velocity
    p_air::FT, # air pressure
    e::FT, # vapor pressure
    ρ_air::FT, # air density
    const_dt:: FT,
) where {FT}

    # (might want to change the name of this function at some point)
    # might also want to make turn these parameters into a struct or something
    # implementing this -- first the simplest form that
    # uses the assumption that dqs/dT = 1 and sets A_c = 1
    # and currently this is in specific humidity, but it technically should
    # be converted to mixing ratio

    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_v = TD.latent_heat_vapor(tps, T)
    g = TD.Parameters.grav(tps)

    # need to calculate the internal energy now...
    nonequilibrium_phase = TD.PhaseNonEquil(tps, e_int, ρ_air, q) # is this q correct?
    λ = TD.liquid_fraction(tps, T, TD.Liquid(), q) # the phase type here is currently incorrect
    L = TD.weighted_latent_heat(tps, T, λ)

    dqsldt = TD.∂q_vap_sat_∂T(tps,λ,T,q_sat,L)

    dqsidT = FT(1) # i dont see anything in thermo to calculate this? might need to do myself

    #A_c = FT(1) # i need to actually make this into something

    A_c_WBF = - (q_sat.liq - q_sat.ice)/(ice.τ_relax*Γᵢ)*(1+(L_subl/cp_air)*dqsidT)
    e_s = e/Sₗ # assuming this is ok but would like to double check
    A_c_uplift = -(q_sat.liq * g * w * ρ_air)/(p-e_s) + dqsldT*w*g/c_p

    A_c =  A_c_uplift + A_c_WBF

    Γₗ = FT(1) + (L_v/cp_air)*dqsldT
    Γᵢ = FT(1) + (L_subl/cp_air)*dqsidT

    τ = (liquid.τ_relax^(-1) + (1 + (L_subl/cp_air))*ice.τ_relax^(-1) / Γᵢ)^(-1)
       
    # this does need to be mixing ratio instead
    δ_0 = (Sₗ-1)*q_sat.liq

    # solving for new Sl after timestep delta t:
    Sₗ = (1/q_sat.liq)*(A_c * τ/(ice.τ_relax * Γₗ) + (δ_0 - A_c*τ)*τ/(const_dt*ice.τ_relax*Γₗ)*(FT(1) - exp(- const_dt/τ))) + 1

    return Sₗ
end

end #module MicrophysicsNonEq.jl