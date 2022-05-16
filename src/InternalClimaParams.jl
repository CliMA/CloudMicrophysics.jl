module InternalClimaParams

import CLIMAParameters
const CP = CLIMAParameters
const CP_planet = CLIMAParameters.Planet
const CP_micro = CLIMAParameters.Atmos.Microphysics
const CP_0M = CLIMAParameters.Atmos.Microphysics_0M

const APS = CP.AbstractParameterSet

a0_rai(param_set::APS) = CP_micro.a0_rai(param_set)
ae_rai(param_set::APS) = CP_micro.ae_rai(param_set)
a0_sno(param_set::APS) = CP_micro.a0_sno(param_set)
ae_sno(param_set::APS) = CP_micro.ae_sno(param_set)

a_vent_rai(param_set::APS) = CP_micro.a_vent_rai(param_set)
a_vent_sno(param_set::APS) = CP_micro.a_vent_sno(param_set)
b_vent_rai(param_set::APS) = CP_micro.b_vent_rai(param_set)
b_vent_sno(param_set::APS) = CP_micro.b_vent_sno(param_set)

C_drag(param_set::APS) = CP_micro.C_drag(param_set)

χm_ice(param_set::APS) = CP_micro.χm_ice(param_set)
χm_rai(param_set::APS) = CP_micro.χm_rai(param_set)
χa_rai(param_set::APS) = CP_micro.χa_rai(param_set)
χv_rai(param_set::APS) = CP_micro.χv_rai(param_set)
χm_sno(param_set::APS) = CP_micro.χm_sno(param_set)
χa_sno(param_set::APS) = CP_micro.χa_sno(param_set)
χv_sno(param_set::APS) = CP_micro.χv_sno(param_set)

D_vapor(param_set::APS) = CP_micro.D_vapor(param_set)

Δm_ice(param_set::APS) = CP_micro.Δm_ice(param_set)
Δm_rai(param_set::APS) = CP_micro.Δm_rai(param_set)
Δa_rai(param_set::APS) = CP_micro.Δa_rai(param_set)
Δv_rai(param_set::APS) = CP_micro.Δv_rai(param_set)
Δm_sno(param_set::APS) = CP_micro.Δm_sno(param_set)
Δa_sno(param_set::APS) = CP_micro.Δa_sno(param_set)
Δv_sno(param_set::APS) = CP_micro.Δv_sno(param_set)

E_liq_rai(param_set::APS) = CP_micro.E_liq_rai(param_set)
E_liq_sno(param_set::APS) = CP_micro.E_liq_sno(param_set)
E_ice_rai(param_set::APS) = CP_micro.E_ice_rai(param_set)
E_ice_sno(param_set::APS) = CP_micro.E_ice_sno(param_set)
E_rai_sno(param_set::APS) = CP_micro.E_rai_sno(param_set)

gas_constant() = CP.gas_constant()
grav(param_set::APS) = CP_planet.grav(param_set)

K_therm(param_set::APS) = CP_micro.K_therm(param_set)

molmass_ratio(param_set::APS) = CP_planet.molmass_ratio(param_set)
molmass_water(param_set::APS) = CP_planet.molmass_water(param_set)

m0_ice(param_set::APS) = CP_micro.m0_ice(param_set)
m0_rai(param_set::APS) = CP_micro.m0_rai(param_set)
me_rai(param_set::APS) = CP_micro.me_rai(param_set)
m0_sno(param_set::APS) = CP_micro.m0_sno(param_set)
me_ice(param_set::APS) = CP_micro.me_ice(param_set)
me_sno(param_set::APS) = CP_micro.me_sno(param_set)

n0_ice(param_set::APS) = CP_micro.n0_ice(param_set)
n0_rai(param_set::APS) = CP_micro.n0_rai(param_set)

μ_sno(param_set::APS) = CP_micro.μ_sno(param_set)
ν_sno(param_set::APS) = CP_micro.ν_sno(param_set)

ν_air(param_set::APS) = CP_micro.ν_air(param_set)

qc_0(param_set::APS) = CP_0M.qc_0(param_set)
q_liq_threshold(param_set::APS) = CP_micro.q_liq_threshold(param_set)
q_ice_threshold(param_set::APS) = CP_micro.q_ice_threshold(param_set)

r0_rai(param_set::APS) = CP_micro.r0_rai(param_set)
r0_ice(param_set::APS) = CP_micro.r0_ice(param_set)
r0_sno(param_set::APS) = CP_micro.r0_sno(param_set)

r_ice_snow(param_set::APS) = CP_micro.r_ice_snow(param_set)

R_v(param_set::APS) = CP_planet.R_v(param_set)

ρ_cloud_liq(param_set::APS) = CP_planet.ρ_cloud_liq(param_set)

S_0(param_set::APS) = CP_0M.S_0(param_set)
surface_tension_coeff(param_set::APS) =
    CP_planet.surface_tension_coeff(param_set)

T_freeze(param_set::APS) = CP_planet.T_freeze(param_set)

τ_cond_evap(param_set::APS) = CP_micro.τ_cond_evap(param_set)
τ_sub_dep(param_set::APS) = CP_micro.τ_sub_dep(param_set)
τ_acnv_rai(param_set::APS) = CP_micro.τ_acnv_rai(param_set)
τ_acnv_sno(param_set::APS) = CP_micro.τ_acnv_sno(param_set)
τ_precip(param_set::APS) = CP_0M.τ_precip(param_set)

ve_rai(param_set::APS) = CP_micro.ve_rai(param_set)
v0_sno(param_set::APS) = CP_micro.v0_sno(param_set)
ve_sno(param_set::APS) = CP_micro.ve_sno(param_set)

end
