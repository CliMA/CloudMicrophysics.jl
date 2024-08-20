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
    conv_q_vap_to_q_liq_ice_MM2015(liquid, tps, q, ρ, T)
    conv_q_vap_to_q_liq_ice_MM2015(ice, tps, q, ρ, T)

- `liquid` or `ice` - a struct with cloud water or ice free parameters
- `tps` - thermodynamics parameters struct
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
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT}
    Rᵥ = TD.Parameters.R_v(tps)
    cₚ_air = TD.cp_m(tps, q)
    Lᵥ = TD.latent_heat_vapor(tps, T)
    qᵥ = TD.vapor_specific_humidity(q)

    pᵥ_sat_liq = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    qᵥ_sat_liq = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sat_liq)

    dqsldT = qᵥ_sat_liq * (Lᵥ / (Rᵥ * T^2) - 1 / T)
    Γₗ = FT(1) + (Lᵥ / cₚ_air) * dqsldT

    return (qᵥ - qᵥ_sat_liq) / τ_relax * Γₗ
end
function conv_q_vap_to_q_liq_ice_MM2015(
    (; τ_relax)::CMP.CloudIce{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT}
    Rᵥ = TD.Parameters.R_v(tps)
    cₚ_air = TD.cp_m(tps, q)
    Lₛ = TD.latent_heat_sublim(tps, T)
    qᵥ = TD.vapor_specific_humidity(q)

    pᵥ_sat_ice = TD.saturation_vapor_pressure(tps, T, TD.Ice())
    qᵥ_sat_ice = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sat_ice)

    dqsidT = qᵥ_sat_ice * (Lₛ / (Rᵥ * T^2) - 1 / T)
    Γᵢ = FT(1) + (Lₛ / cₚ_air) * dqsidT

    return (qᵥ - qᵥ_sat_ice) / τ_relax * Γᵢ
end

"""
    conv_q_vap_to_q_liq_ice_MM2015_timeintegrator(liquid, ice, tps, q, ρ, T, w, p_air, const_dt, type)

- `liquid` - a struct with cloud water free parameters
- `ice` - a struct with cloud ice free parameters
- `tps` - thermodynamics parameters struct
- `q` - current PhasePartition
- `ρ` - air density [kg/m3]
- `T` - air temperature [K]
- `w` - vertical velocity [m/s]
- `p_air` - air pressure 
- `const_dt`-
- `type` -

Returns the cloud water tendency due to condensation and evaporation
or cloud ice tendency due to sublimation and vapor deposition.
The formulation is based on Morrison and Grabowski 2008 and
Morrison and Milbrandt 2015
"""
function conv_q_vap_to_q_liq_ice_MM2015_timeintegrator(
    liquid::CMP.CloudLiquid{FT},
    ice::CMP.CloudIce{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
    w::FT,
    p_air::FT,
    const_dt::FT,
    type::String,
) where {FT}

    cₚ_air = TD.cp_m(tps, q)
    Lₛ = TD.latent_heat_sublim(tps, T)
    Lᵥ = TD.latent_heat_vapor(tps, T)
    Rᵥ = TD.Parameters.R_v(tps)
    g = TD.Parameters.grav(tps)
    qᵥ = TD.vapor_specific_humidity(q)

    pᵥ_sat_liq = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    qᵥ_sat_liq = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sat_liq)

    pᵥ_sat_ice = TD.saturation_vapor_pressure(tps, T, TD.Ice())
    qᵥ_sat_ice = TD.q_vap_saturation_from_density(tps, T, ρ, pᵥ_sat_ice)

    dqsldT = qᵥ_sat_liq * (Lᵥ/(Rᵥ * T^2) - 1 / T)
    dqsidT = qᵥ_sat_ice * (Lₛ/(Rᵥ * T^2) - 1 / T)

    Γₗ = FT(1) + (Lᵥ / cₚ_air) * dqsldT
    Γᵢ = FT(1) + (Lₛ / cₚ_air) * dqsidT

    A_c_WBF =
        (qᵥ_sat_liq - qᵥ_sat_ice) / (ice.τ_relax * Γᵢ) * (1 + (Lₛ / cₚ_air) * dqsldT)

    A_c_uplift = -(qᵥ_sat_liq * g * w * ρ) / (p_air - pᵥ_sat_liq) + dqsldT * w * g / cₚ_air

    A_c = A_c_uplift - A_c_WBF

    τ = (liquid.τ_relax^(-1) + (1 + (Lₛ / cₚ_air)*dqsldT) * ice.τ_relax^(-1) / Γᵢ)^(-1)
    
    δ_0 = qᵥ - qᵥ_sat_liq

    if type == "condensation"
        cond_rate =
            A_c * τ / (liquid.τ_relax * Γₗ) +
            (δ_0 - A_c * τ) * τ / (const_dt * liquid.τ_relax * Γₗ) *
            (FT(1) - exp(-const_dt / τ))
        return cond_rate
    elseif type == "deposition"
        dep_rate =
            A_c * τ / (ice.τ_relax * Γᵢ) +
            (δ_0 - A_c * τ) * τ / (const_dt * ice.τ_relax * Γᵢ) *
            (FT(1) - exp(-const_dt / τ)) + (qᵥ_sat_liq - qᵥ_sat_ice)/(ice.τ_relax*Γᵢ)
        return dep_rate
    end

end

end
