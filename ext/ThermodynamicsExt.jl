module ThermodynamicsExt

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP
const TPS = TD.Parameters.ThermodynamicsParameters

#! format: off

###
### Gas constants
###
Rᵥ(tps::TPS) = TD.Parameters.R_v(tps)

Lᵥ(tps::TPS, T) = TD.latent_heat_vapor(tps, T)
Lₛ(tps::TPS, T) = TD.latent_heat_sublim(tps, T)
Lf(tps::TPS, T) = TD.latent_heat_fusion(tps, T)

cpₘ(tps::TPS, qₜ, qₗ, qᵢ) = TD.cp_m(tps, qₜ, qₗ, qᵢ)

###
### Supersaturations
###
saturation_vapor_pressure_over_liquid(tps::TPS, T) = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
saturation_vapor_pressure_over_ice(tps::TPS, T) = TD.saturation_vapor_pressure(tps, T, TD.Ice())

supersaturation_over_liquid(tps::TPS, q_tot, q_liq, q_ice, q_rai, q_sno, ρ, T) =
    TD.supersaturation(tps, TD.PhasePartition(q_tot, q_liq + q_rai, q_ice + q_sno), ρ, T, TD.Liquid())
supersaturation_over_ice(tps::TPS, q_tot, q_liq, q_ice, q_rai, q_sno, ρ, T) =
    TD.supersaturation(tps, TD.PhasePartition(q_tot, q_liq + q_rai, q_ice + q_sno), ρ, T, TD.Ice())

###
### Utility functions
###
q_vap(q_tot, q_liq, q_ice) = q_tot - q_liq - q_ice
p2q(tps::TPS, T, ρ, p) = TD.q_vap_saturation_from_density(tps, T, ρ, p)

###
### Wrappers for using TD.PhasePartition in different CloudMicrophysics modules
###

# 0-moment

remove_precipitation(params::CMP.Parameters0M, q::TD.PhasePartition) = remove_precipitation(params, q.liq, q.ice)
remove_precipitation(params::CMP.Parameters0M, q::TD.PhasePartition, q_vap_sat) = remove_precipitation(params, q.liq, q.ice, q_vap_sat)

# Assuming q = (q_tot, q_liq + q_rai, q_ice + q_sno)
q_(q::TD.PhasePartition, q_rai, q_sno) = (q.tot, q.liq - q_rai, q.ice - q_sno)

# Non-equilibrium
conv_q_vap_to_q_liq_ice_MM2015(prs::CMP.CloudLiquid, tps::TPS, q::TD.PhasePartition, q_rai, q_sno, ρ, T) =
conv_q_vap_to_q_liq_ice_MM2015(prs::CMP.CloudLiquid, tps::TPS, q_(q::TD.PhasePartition, q_rai, q_sno)..., ρ, T)

conv_q_vap_to_q_liq_ice_MM2015(prs::CMP.CloudIce, tps::TPS, q::TD.PhasePartition, q_rai, q_sno, ρ, T) =
conv_q_vap_to_q_liq_ice_MM2015(prs::CMP.CloudIce, tps::TPS, q_(q::TD.PhasePartition, q_rai, q_sno)..., ρ, T)

# 1-moment
conv_q_ice_to_q_sno(
    ice_params::CMP.CloudIce, aps::CMP.AirProperties, tps::TPS, q::TD.PhasePartition, q_rai, q_sno, ρ, T) =
conv_q_ice_to_q_sno(
    ice_params::CMP.CloudIce, aps::CMP.AirProperties, tps::TPS, q_(q::TD.PhasePartition, q_rai, q_sno)..., ρ, T)

evaporation_sublimation(
    rain_params::CMP.Rain, vel::CMP.Blk1MVelTypeRain, aps::CMP.AirProperties, tps::TPS, q::TD.PhasePartition, q_rai, q_sno, ρ, T) =
evaporation_sublimation(
    rain_params::CMP.Rain, vel::CMP.Blk1MVelTypeRain, aps::CMP.AirProperties, tps::TPS, q_(q::TD.PhasePartition, q_rai, q_sno)..., ρ, T)

evaporation_sublimation(
    snow_params::CMP.Snow, vel::CMP.Blk1MVelTypeSnow, aps::CMP.AirProperties, tps::TPS, q::TD.PhasePartition, q_rai, q_sno, ρ, T) =
evaporation_sublimation(
    snow_params::CMP.Snow, vel::CMP.Blk1MVelTypeSnow, aps::CMP.AirProperties, tps::TPS, q_(q::TD.PhasePartition, q_rai, q_sno)..., ρ, T)

#! format: on
end
