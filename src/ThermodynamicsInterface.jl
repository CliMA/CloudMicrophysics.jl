module ThermodynamicsInterface

import Thermodynamics as TD
const PS = TD.Parameters.ThermodynamicsParameters

###
### Constants
###
grav(tps::PS) = TD.Parameters.grav(tps) # Needed in parcel model
T_freeze(tps::PS) = TD.Parameters.T_freeze(tps)

Rᵥ(tps::PS) = TD.Parameters.R_v(tps)
Rd(tps::PS) = TD.Parameters.R_d(tps)
Rd_over_Rv(tps::PS) = 1 / TD.Parameters.Rv_over_Rd(tps)

# Gas constant for moist air
Rₘ(tps::PS, qₜ, qₗ, qᵢ) = TD.gas_constant_air(tps, qₜ, qₗ, qᵢ)

Lᵥ(tps::PS, T) = TD.latent_heat_vapor(tps, T)
Lₛ(tps::PS, T) = TD.latent_heat_sublim(tps, T)
Lf(tps::PS, T) = TD.latent_heat_fusion(tps, T)

cpₘ(tps::PS, qₜ, qₗ, qᵢ) = TD.cp_m(tps, qₜ, qₗ, qᵢ)
cv_l(tps::PS) = TD.Parameters.cv_l(tps)
cp_l(tps::PS) = TD.Parameters.cp_l(tps)

###
### Internal energies and liquid fraction (needed for energy sources from 0M)
###
liquid_fraction(tps::PS, T, q_lcl, q_icl) =
    TD.liquid_fraction(tps, T, q_lcl, q_icl)
internal_energy_liquid(tps::PS, T) = TD.internal_energy_liquid(tps, T)
internal_energy_ice(tps::PS, T) = TD.internal_energy_ice(tps, T)

###
### Utility functions
###

import ..Utilities as UT

"""
    q_vap(q_tot, q_liq, q_ice)
    q_vap(q_tot, q_lcl, q_icl, q_rai, q_sno)

Compute vapor specific content from total water and condensed phase specific contents.

# Arguments
- `q_tot`: Total water specific content [kg/kg]
- `q_liq` or `q_lcl`: Liquid water specific content [kg/kg]
- `q_ice` or `q_icl`: Ice specific content [kg/kg]
- `q_rai`: Rain specific content [kg/kg] (5-argument version only)
- `q_sno`: Snow specific content [kg/kg] (5-argument version only)

# Returns
- Vapor specific content [kg/kg], clamped to be non-negative

# Notes
- Negative values from `q_tot - Σq_condensed` are clamped to zero using AD-compatible operations
"""
q_vap(q_tot, q_liq, q_ice) = UT.clamp_to_nonneg(q_tot - q_liq - q_ice)
q_vap(q_tot, q_lcl, q_icl, q_rai, q_sno) = UT.clamp_to_nonneg(q_tot - q_lcl - q_icl - q_rai - q_sno)

# Get specific content from partial pressure
p2q(tps::PS, T, ρ, pᵥ) = TD.q_vap_from_p_vap(tps, T, ρ, pᵥ)

# Get partial pressure from specific content
q2p(tps::PS, T, ρ, qᵥ) = qᵥ * ρ * Rᵥ(tps) * T

# Get air density from temperature, pressure and water content
# (only used in tests)
air_density(tps::PS, T, p, q_tot, q_liq, q_ice) =
    TD.air_density(tps, T, p, q_tot, q_liq, q_ice)

# Get water vapor specific content from relative humidity over liquid
# (used in documentation plots for P3 scheme)
q_vap_from_RH_over_liquid(tps::PS, p, T, RH) =
    TD.q_vap_from_RH(tps, p, T, RH, TD.Liquid())

###
### Supersaturations
###
saturation_vapor_pressure_over_liquid(tps::PS, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Liquid())
saturation_vapor_pressure_over_ice(tps::PS, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Ice())

# (only used in tests)
saturation_vapor_specific_content_over_liquid(tps::PS, T, ρ) =
    TD.q_vap_saturation(tps, T, ρ, TD.Liquid())
saturation_vapor_specific_content_over_ice(tps::PS, T, ρ) =
    TD.q_vap_saturation(tps, T, ρ, TD.Ice())

"""
    supersaturation_over_liquid(tps, qₜ, qₗ, qᵢ, ρ, T)
    supersaturation_over_ice(tps, qₜ, qₗ, qᵢ, ρ, T)

Compute supersaturation with respect to liquid water or ice.

Uses clamped vapor specific content calculation to ensure robustness against
negative humidity inputs while preserving AD compatibility.

# Arguments
- `tps`: Thermodynamics parameters
- `qₜ`: Total water specific content [kg/kg]
- `qₗ`: Total liquid water specific content [kg/kg]
- `qᵢ`: Total ice specific content [kg/kg]
- `ρ`: Air density [kg/m³]
- `T`: Temperature [K]

# Returns
- Supersaturation: `S = (pᵥ / pᵥ_sat) - 1` (dimensionless)
  - `S > 0`: supersaturated (condensation/deposition occurs)
  - `S < 0`: subsaturated (evaporation/sublimation occurs)

# Notes
- Internally uses `q_vap` which clamps negative vapor values to zero
- This prevents unphysical supersaturation values from negative humidity inputs
"""
function supersaturation_over_liquid(tps::PS, qₜ, qₗ, qᵢ, ρ, T)
    qᵥ = q_vap(qₜ, qₗ, qᵢ)
    return TD.supersaturation(tps, qᵥ, ρ, T, TD.Liquid())
end
function supersaturation_over_ice(tps::PS, qₜ, qₗ, qᵢ, ρ, T)
    qᵥ = q_vap(qₜ, qₗ, qᵢ)
    return TD.supersaturation(tps, qᵥ, ρ, T, TD.Ice())
end

end
