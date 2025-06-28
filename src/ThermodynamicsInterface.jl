module ThermodynamicsInterface

import Thermodynamics as TD
const PS = TD.Parameters.ThermodynamicsParameters
const PP = TD.PhasePartition


# Constants

grav(tps::PS) = TD.Parameters.grav(tps) # Needed in parcel model


# Constants for moist air

Rᵥ(tps::PS) = TD.Parameters.R_v(tps)
Rd(tps::PS) = TD.Parameters.R_d(tps)
Rd_over_Rv(tps::PS) = 1 / TD.Parameters.molmass_ratio(tps)
# Gas constant for mois air
Rₘ(tps::PS, qₜ, qₗ, qᵢ) = TD.gas_constant_air(tps, PP(qₜ, qₗ, qᵢ)) #TODO - PP

Lᵥ(tps::PS, T) = TD.latent_heat_vapor(tps, T)
Lₛ(tps::PS, T) = TD.latent_heat_sublim(tps, T)
Lf(tps::PS, T) = TD.latent_heat_fusion(tps, T)

cpₘ(tps::PS, qₜ, qₗ, qᵢ) = TD.cp_m(tps, qₜ, qₗ, qᵢ)

# Supersaturations

saturation_vapor_pressure_over_liquid(tps::PS, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Liquid())
saturation_vapor_pressure_over_ice(tps::PS, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Ice())

# (only used in tests)
saturation_vapor_specific_content_over_liquid(tps::PS, T, ρ) =
    TD.q_vap_saturation_generic(tps, T, ρ, TD.Liquid())
saturation_vapor_specific_content_over_ice(tps::PS, T, ρ) =
    TD.q_vap_saturation_generic(tps, T, ρ, TD.Ice())

supersaturation_over_liquid(tps::PS, qₜ, qₗ, qᵢ, ρ, T) =
    TD.supersaturation(tps, PP(qₜ, qₗ, qᵢ), ρ, T, TD.Liquid()) # TODO - PP
supersaturation_over_ice(tps::PS, qₜ, qₗ, qᵢ, ρ, T) =
    TD.supersaturation(tps, PP(qₜ, qₗ, qᵢ), ρ, T, TD.Ice()) # TODO - PP


# Utility functions

# Get vapor specific content from total water, total liquid water and total ice
# water specific contents.
q_vap(q_tot, q_liq, q_ice) = q_tot - q_liq - q_ice
# Get vapor specific content from total water, cloud liquid water, cloud ice,
# rain and snow specific contents.
q_vap(q_tot, q_liq, q_ice, q_rai, q_sno) = q_tot - q_liq - q_ice - q_rai - q_sno

# Get specific content from partial pressure
p2q(tps::PS, T, ρ, p) = TD.q_vap_saturation_from_density(tps, T, ρ, p)


# Get air density from temperature, pressure and water content
# (only used in tests)
air_density(tps::PS, T, p, q_tot, q_liq, q_ice) =
    TD.air_density(tps, T, p, PP(q_tot, q_liq, q_ice))

# Get a tuple containing total water, cloud liquid water and cloud ice specific
# contents from TD.PhasePartition and rain and snow specific contents.
# Assumes that q::PhasePartition = (q_tot, q_cloud_liq + q_rai, q_cloud_ice + q_sno)
q_(q::PP, q_rai, q_sno) = (q.tot, q.liq - q_rai, q.ice - q_sno)

end
