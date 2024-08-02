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
    T::FT
) where {FT}
    # THE MORRISON & MILBRANDT VERSION WITHOUT THE THE INTEGRATOR

    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_v = TD.latent_heat_vapor(tps, T)

    nonequil_phase = TD.PhaseNonEquil
    λ = TD.liquid_fraction(tps, T, nonequil_phase, q)
    L = TD.weighted_latent_heat(tps, T, λ)

    # honestly both of these i dont like
    dqsldT = TD.∂q_vap_sat_∂T(tps,λ,T,q_sat.liq,L)
    # oh lol this is probably still wrong -- probably need to think more abt it
    dqsidT = TD.∂q_vap_sat_∂T(tps,λ,T,q_sat.ice,L)

    Γₗ = FT(1) + (L_v/cp_air)*dqsldT
    Γᵢ = FT(1) + (L_subl/cp_air)*dqsidT

    τ = (liquid.τ_relax^(-1) + (1 + (L_subl/cp_air))*ice.τ_relax^(-1) / Γᵢ)^(-1)
    
    q_v = q.tot - q.liq - q.ice
    cond_rate = (q_v - q_sat.liq) / τ*Γₗ
    dep_rate = (q_v - q_sat.ice) / τ*Γᵢ

    return cond_rate, dep_rate
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
    type::String,
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

    A_c_WBF = (q_sat.liq - q_sat.ice)/(ice.τ_relax*Γᵢ)*(1+(L_subl/cp_air)*dqsldT)
    #A_c_WBF = 0
    e_sl = e/Sₗ # assuming this is ok but would like to double check
    Sᵢ = TD.saturation_vapor_pressure(tps, T, TD.Liquid()) / TD.saturation_vapor_pressure(tps, T, TD.Ice()) * Sₗ 
    e_si = e/Sᵢ

    A_c_uplift_l = -(q_sat.liq * g * w * ρ_air)/(p_air-e_sl) + dqsldT*w*g/cp_air
    A_c_uplift_i = -(q_sat.ice * g * w * ρ_air)/(p_air-e_si) + dqsidT*w*g/cp_air

    A_c_l =  A_c_uplift_l - A_c_WBF
    A_c_i = A_c_uplift_i + A_c_WBF
    
    τ = (liquid.τ_relax^(-1) + (1 + (L_subl/cp_air))*ice.τ_relax^(-1) / Γᵢ)^(-1)
    
    # this does need to be mixing ratio instead -- double check this calc
    #q.tot - q.liq -q.ice - q_sat.liq
    #δ_0 = (Sₗ-1)*q_sat.liq
    q_v = q.tot - q.liq - q.ice
    δ_0_l = q_v - q_sat.liq #(Sₗ-1)*q_sat.liq
    δ_0_i = q_v - q_sat.ice #(Sₗ-1)*q_sat.liq

    @info("", q_v, q_sat.liq, q_sat.ice, δ_0_l, δ_0_i)
    # just doing a quick hacky version of calculating S_i this way -- not entirely sure
    # taken from parcel common


    #@info("",A_c_uplift, A_c_WBF, A_c, e_s, e_si, e)
    #@info("", ice.τ_relax, liquid.τ_relax)


    # DONT DO THIS -- CALCULATE DQDT INSTEAD.

    if type == "condensation"
        cond_rate = A_c_l * τ/(liquid.τ_relax * Γₗ) + (δ_0_l - A_c_l*τ)*τ/(const_dt*liquid.τ_relax*Γₗ)*(FT(1) - exp(- const_dt/τ))
        #term1 = A_c * τ/(liquid.τ_relax * Γₗ)
        #term2 = (δ_0 - A_c*τ)*τ/(const_dt*liquid.τ_relax*Γₗ)*(FT(1) - exp(- const_dt/τ))
        @info("",cond_rate)#, term1, term2, q.tot-q_sat.liq)
        return cond_rate
    elseif type == "deposition"
        dep_rate = A_c_i * τ/(ice.τ_relax * Γᵢ) + (δ_0_i - A_c_i*τ)*τ/(const_dt*ice.τ_relax*Γᵢ)*(FT(1) - exp(- const_dt/τ)) #+ (q_sat.liq - q_sat.ice)/(ice.τ_relax*Γᵢ)
        #term1 = A_c * τ/(ice.τ_relax * Γᵢ)
        #term2 = (δ_0 - A_c*τ)*τ/(const_dt*ice.τ_relax*Γᵢ)*(FT(1) - exp(- const_dt/τ))
        #term3 = (q_sat.liq - q_sat.ice)/(ice.τ_relax*Γᵢ)
        #avg_δ = A_c * τ/ + (δ_0 - A_c*τ)*τ/(const_dt)*(FT(1) - exp(- const_dt/τ))
        #@info("",dep_rate, term1, term2, term3, avg_δ, q.tot-q_sat.ice, q_sat.ice)
        @info("",dep_rate)# - A_c*τ, τ/(const_dt*ice.τ_relax*Γᵢ), (FT(1) - exp(- const_dt/τ)))
        #@info("", τ)
        return dep_rate
    end

    # solving for new Sl after timestep delta t:
    #Sₗ = (1/q_sat.liq)*() + 1

    #Sᵢ = (1/q_sat.ice)*() + 1

    # a question -- is this Sl or total S?

    #return Sₗ, Sᵢ
end

end #module MicrophysicsNonEq.jl