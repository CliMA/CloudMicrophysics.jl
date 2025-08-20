module ThermodynamicsInterface

import Thermodynamics as TD
const PS = TD.Parameters.ThermodynamicsParameters
const PP = TD.PhasePartition  # Still used in Atmos EquilMoist model rain formation

###
### Constants
###
grav(tps::PS) = TD.Parameters.grav(tps) # Needed in parcel model

Rᵥ(tps::PS) = TD.Parameters.R_v(tps)
Rd(tps::PS) = TD.Parameters.R_d(tps)
Rd_over_Rv(tps::PS) = 1 / TD.Parameters.Rv_over_Rd(tps)

# Gas constant for moist air
Rₘ(tps::PS, qₜ, qₗ, qᵢ) = TD.gas_constant_air(tps, qₜ, qₗ, qᵢ)

Lᵥ(tps::PS, T) = TD.latent_heat_vapor(tps, T)
Lₛ(tps::PS, T) = TD.latent_heat_sublim(tps, T)
Lf(tps::PS, T) = TD.latent_heat_fusion(tps, T)

cpₘ(tps::PS, qₜ, qₗ, qᵢ) = TD.cp_m(tps, qₜ, qₗ, qᵢ)

###
### Utility functions
###

# Get vapor specific content from total water, total liquid water and total ice
# water specific contents.
q_vap(q_tot, q_liq, q_ice) = q_tot - q_liq - q_ice
# Get vapor specific content from total water, cloud liquid water, cloud ice,
# rain and snow specific contents.
q_vap(q_tot, q_liq, q_ice, q_rai, q_sno) = q_tot - q_liq - q_ice - q_rai - q_sno

# Get specific content from partial pressure
p2q(tps::PS, T, ρ, pᵥ) = TD.q_vap_from_p_vap(tps, T, ρ, pᵥ)

# Get partial pressure from specific content
q2p(tps::PS, T, ρ, qᵥ) = qᵥ * ρ * Rᵥ(tps) * T

# Get air density from temperature, pressure and water content
# (only used in tests)
air_density(tps::PS, T, p, q_tot, q_liq, q_ice) =
    TD.air_density(tps, T, p, PP(q_tot, q_liq, q_ice))

###
### Supersaturations
###
saturation_vapor_pressure_over_liquid(tps::PS, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Liquid())
saturation_vapor_pressure_over_ice(tps::PS, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Ice())

# (only used in tests)
saturation_vapor_specific_content_over_liquid(tps::PS, T, ρ) =
    TD.q_vap_saturation_generic(tps, T, ρ, TD.Liquid())
saturation_vapor_specific_content_over_ice(tps::PS, T, ρ) =
    TD.q_vap_saturation_generic(tps, T, ρ, TD.Ice())

function supersaturation_over_liquid(tps::PS, qₜ, qₗ, qᵢ, ρ, T)
    pᵥ_sat = saturation_vapor_pressure_over_liquid(tps, T)
    qᵥ = q_vap(qₜ, qₗ, qᵢ)
    pᵥ = q2p(tps, T, ρ, qᵥ)
    return pᵥ / pᵥ_sat - 1
end
function supersaturation_over_ice(tps::PS, qₜ, qₗ, qᵢ, ρ, T)
    pᵥ_sat = saturation_vapor_pressure_over_ice(tps, T)
    qᵥ = q_vap(qₜ, qₗ, qᵢ)
    pᵥ = q2p(tps, T, ρ, qᵥ)
    return pᵥ / pᵥ_sat - 1
end

end
