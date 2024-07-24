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

    # going to need to either change this function to have ice functionality or do something else?


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


    # e_int = TD.internal_energy(tps,T)

    nonequil_phase = TD.PhaseNonEquil # 
    λ = TD.liquid_fraction(tps, T, nonequil_phase, q)
    L = TD.weighted_latent_heat(tps, T, λ)

    dqsldT = TD.∂q_vap_sat_∂T(tps,λ,T,q_sat.liq,L)

    dqsidT = TD.∂q_vap_sat_∂T(tps,λ,T,q_sat.ice,L) # just changing the phase? is that enough?
    #FT(1) # i dont see anything in thermo to calculate this? might need to do myself
    # ask Amy

    Γₗ = FT(1) + (L_v/cp_air)*dqsldT
    Γᵢ = FT(1) + (L_subl/cp_air)*dqsidT

    #A_c = FT(1) # i need to actually make this into something

    A_c_WBF = - (q_sat.liq - q_sat.ice)/(ice.τ_relax*Γᵢ)*(1+(L_subl/cp_air)*dqsidT)
    e_s = e/Sₗ # assuming this is ok but would like to double check
    A_c_uplift = -(q_sat.liq * g * w * ρ_air)/(p_air-e_s) + dqsldT*w*g/cp_air

    A_c =  A_c_uplift + A_c_WBF

    τ = (liquid.τ_relax^(-1) + (1 + (L_subl/cp_air))*ice.τ_relax^(-1) / Γᵢ)^(-1)
       
    # this does need to be mixing ratio instead -- double check this calc
    δ_0ₗ = (Sₗ-1)*q_sat.liq
    # just doing a quick hacky version of calculating S_i this way -- not entirely sure
    # taken from parcel common
    Sᵢ = TD.saturation_vapor_pressure(tps, T, TD.Liquid()) / TD.saturation_vapor_pressure(tps, T, TD.Ice()) * Sₗ 
    
    δ_0ᵢ = (Sᵢ-1)*q_sat.ice

    # solving for new Sl after timestep delta t:
    Sₗ = (1/q_sat.liq)*(A_c * τ/(liquid.τ_relax * Γₗ) + (δ_0ₗ - A_c*τ)*τ/(const_dt*liquid.τ_relax*Γₗ)*(FT(1) - exp(- const_dt/τ))) + 1

    Sᵢ = (1/q_sat.ice)*(A_c * τ/(ice.τ_relax * Γᵢ) + (δ_0ᵢ - A_c*τ)*τ/(const_dt*ice.τ_relax*Γᵢ)*(FT(1) - exp(- const_dt/τ)) + (q_sat.liq + q_sat.ice)/(ice.τ_relax*Γᵢ)) + 1

    # a question -- is this Sl or total S?

    return Sₗ, Sᵢ
end

end #module MicrophysicsNonEq.jl