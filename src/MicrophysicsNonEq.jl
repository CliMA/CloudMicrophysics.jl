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
first simple formulation:

    conv_q_vap_to_q_liq_ice(liquid, q_sat, q)
    conv_q_vap_to_q_liq_ice(ice, q_sat, q)

 - `liquid` or `ice` - a type for cloud water or ice
 - `q_sat` - PhasePartition at equilibrium
 - `q` - current PhasePartition

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The tendency is obtained assuming a relaxation to equilibrium with
a constant timescale.

second simple formulation:

    conv_q_vap_to_q_liq_ice(tps, liquid, q_sat, q, T)
    conv_q_vap_to_q_liq_ice(tps, ice, q_sat, q, T)

- `tps` - thermodynamics Parameters
- `liquid` or `ice` - a type for cloud water or ice
- `q_sat` - PhasePartition of water vapor specific humidity at saturation (different than above!!)
- `q` - current PhasePartition
- `T` - temperature in Kelvin

third complex formulation

    conv_q_vap_to_q_liq_ice(tps, liquid, ice, q_sat, q, T, Sₗ, w, p_air, e, ρ_air, const_dt, "condensation")
    conv_q_vap_to_q_liq_ice(tps, liquid, ice, q_sat, q, T, Sₗ, w, p_air, e, ρ_air, const_dt, "deposition")

- `tps` - thermodynamics Parameters
- `liquid` - a type for cloud water
- `ice` - a type for cloud ice
- `q_sat` - PhasePartition of water vapor specific humidity at saturation
- `q` - current PhasePartition
- `T` - temperature in Kelvin
- `w` - vertical uplift velocity
- `p_air` - air pressure
- `e` - water vapor pressure
- `ρ_air` - air density
- `const_dt` - length of time step
- `"condensation"` or `"deposition"` - type of process to calculate and output

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
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
    T::FT,
) where {FT}
    # condensation version

    cp_air = TD.cp_m(tps, q)
    L_v = TD.latent_heat_vapor(tps, T)

    nonequil_phase = TD.PhaseNonEquil
    λ = TD.liquid_fraction(tps, T, nonequil_phase, q)
    L = TD.weighted_latent_heat(tps, T, λ)

    # honestly both of these i dont like
    dqsldT = TD.∂q_vap_sat_∂T(tps, λ, T, q_sat.liq, L)

    Γₗ = FT(1) + (L_v / cp_air) * dqsldT

    q_v = q.tot - q.liq - q.ice
    cond_rate = (q_v - q_sat.liq) / liquid.τ_relax * Γₗ

    return cond_rate
end

function conv_q_vap_to_q_liq_ice(
    tps::TDP.ThermodynamicsParameters{FT},
    ice::CMP.CloudIce{FT},
    q_sat::TD.PhasePartition{FT},
    q::TD.PhasePartition{FT},
    T::FT,
) where {FT}
    # deposition version

    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)

    nonequil_phase = TD.PhaseNonEquil
    λ = TD.liquid_fraction(tps, T, nonequil_phase, q)
    L = TD.weighted_latent_heat(tps, T, λ)

    # oh lol this is probably still wrong -- probably need to think more abt it
    dqsidT = TD.∂q_vap_sat_∂T(tps, λ, T, q_sat.ice, L)

    Γᵢ = FT(1) + (L_subl / cp_air) * dqsidT

    q_v = q.tot - q.liq - q.ice
    dep_rate = (q_v - q_sat.ice) / ice.τ_relax * Γᵢ

    return dep_rate
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
    const_dt::FT,
    type::String,
) where {FT}

    cp_air = TD.cp_m(tps, q)
    L_subl = TD.latent_heat_sublim(tps, T)
    L_v = TD.latent_heat_vapor(tps, T)
    g = TD.Parameters.grav(tps)

    nonequil_phase = TD.PhaseNonEquil # 
    λ = TD.liquid_fraction(tps, T, nonequil_phase, q)
    L = TD.weighted_latent_heat(tps, T, λ)

    # a bit unsure that this is the correct way to calculate this
    dqsldT = TD.∂q_vap_sat_∂T(tps, λ, T, q_sat.liq, L)
    dqsidT = TD.∂q_vap_sat_∂T(tps, λ, T, q_sat.ice, L)

    Γₗ = FT(1) + (L_v / cp_air) * dqsldT
    Γᵢ = FT(1) + (L_subl / cp_air) * dqsidT

    A_c_WBF =
        (q_sat.liq - q_sat.ice) / (ice.τ_relax * Γᵢ) * (1 + (L_subl / cp_air) * dqsldT)
    #A_c_WBF = 0
    e_sl = e / Sₗ # assuming this is ok but would like to double check
    Sᵢ =
        TD.saturation_vapor_pressure(tps, T, TD.Liquid()) /
        TD.saturation_vapor_pressure(tps, T, TD.Ice()) * Sₗ
    e_si = e / Sᵢ

    A_c_uplift_l = -(q_sat.liq * g * w * ρ_air) / (p_air - e_sl) + dqsldT * w * g / cp_air
    A_c_uplift_i = -(q_sat.ice * g * w * ρ_air) / (p_air - e_si) + dqsidT * w * g / cp_air

    A_c_l = A_c_uplift_l - A_c_WBF
    A_c_i = A_c_uplift_i + A_c_WBF

    τ = (liquid.τ_relax^(-1) + (1 + (L_subl / cp_air)) * ice.τ_relax^(-1) / Γᵢ)^(-1)

    q_v = q.tot - q.liq - q.ice
    δ_0_l = q_v - q_sat.liq #(Sₗ-1)*q_sat.liq
    δ_0_i = q_v - q_sat.ice #(Sₗ-1)*q_sat.liq

    if type == "condensation"
        cond_rate =
            A_c_l * τ / (liquid.τ_relax * Γₗ) +
            (δ_0_l - A_c_l * τ) * τ / (const_dt * liquid.τ_relax * Γₗ) *
            (FT(1) - exp(-const_dt / τ))
        return cond_rate
    elseif type == "deposition"
        dep_rate =
            A_c_i * τ / (ice.τ_relax * Γᵢ) +
            (δ_0_i - A_c_i * τ) * τ / (const_dt * ice.τ_relax * Γᵢ) *
            (FT(1) - exp(-const_dt / τ)) #+ (q_sat.liq - q_sat.ice)/(ice.τ_relax*Γᵢ)
        return dep_rate
    end

end

end #module MicrophysicsNonEq.jl
