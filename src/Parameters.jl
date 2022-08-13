module Parameters

import Thermodynamics
const TD = Thermodynamics
const TDPS = TD.Parameters.ThermodynamicsParameters

abstract type AbstractCloudMicrophysicsParameters end

# TODO: add doc strings
# Cloud microphysics parameters
"""
    CloudMicrophysicsParameters

"""
Base.@kwdef struct CloudMicrophysicsParameters{FT, TP} <:
                   AbstractCloudMicrophysicsParameters
    K_therm::FT
    D_vapor::FT
    ρ_cloud_liq::FT
    surface_tension_coeff::FT
    τ_cond_evap::FT
    τ_sub_dep::FT
    τ_precip::FT
    qc_0::FT
    S_0::FT
    C_drag::FT
    ν_air::FT
    r_ice_snow::FT
    n0_ice::FT
    r0_ice::FT
    me_ice::FT
    χm_ice::FT
    Δm_ice::FT
    a_vent_rai::FT
    b_vent_rai::FT
    n0_rai::FT
    r0_rai::FT
    me_rai::FT
    ae_rai::FT
    ve_rai::FT
    χm_rai::FT
    Δm_rai::FT
    χa_rai::FT
    Δa_rai::FT
    χv_rai::FT
    Δv_rai::FT
    τ_acnv_rai::FT
    q_liq_threshold::FT
    q_ice_threshold::FT
    τ_acnv_sno::FT
    a_vent_sno::FT
    b_vent_sno::FT
    ν_sno::FT
    μ_sno::FT
    r0_sno::FT
    me_sno::FT
    ae_sno::FT
    ve_sno::FT
    χm_sno::FT
    Δm_sno::FT
    χa_sno::FT
    Δa_sno::FT
    χv_sno::FT
    Δv_sno::FT
    E_liq_rai::FT
    E_liq_sno::FT
    E_ice_rai::FT
    E_ice_sno::FT
    E_rai_sno::FT
    ρ_cloud_ice::FT
    gas_constant::FT
    thermo_params::TP
    D_acnv_TC1980::FT
    a_acnv_TC1980::FT
    b_acnv_TC1980::FT
    r_0_acnv_TC1980::FT
    A_acc_TC1980::FT
    me_liq_TC1980::FT
    C_acnv_B1994::FT
    a_acnv_B1994::FT
    b_acnv_B1994::FT
    c_acnv_B1994::FT
    d_low_acnv_B1994::FT
    d_high_acnv_B1994::FT
    N_0_B1994::FT
    A_acc_B1994::FT
    A_acnv_KK2000::FT
    a_acnv_KK2000::FT
    b_acnv_KK2000::FT
    c_acnv_KK2000::FT
    A_acc_KK2000::FT
    a_acc_KK2000::FT
    b_acc_KK2000::FT
    R_6C_coeff_LD2004::FT
    E_0_LD2004::FT
    k_thrshld_stpnss::FT
end

const CMPS = CloudMicrophysicsParameters

thermodynamics_params(cmp::CMPS) = cmp.thermo_params

# Derived parameters
N_Sc(ps::CMPS) = ps.ν_air / ps.D_vapor
m0_ice(ps::CMPS{FT}) where {FT} =
    FT(4 / 3) * π * ps.ρ_cloud_ice * ps.r0_ice^ps.me_ice
m0_rai(ps::CMPS{FT}) where {FT} =
    FT(4 / 3) * π * ps.ρ_cloud_liq * ps.r0_rai^ps.me_rai
a0_rai(ps::CMPS) = π * ps.r0_rai^ps.ae_rai
m0_sno(ps::CMPS{FT}) where {FT} = FT(1e-1) * ps.r0_sno^ps.me_sno
a0_sno(ps::CMPS{FT}) where {FT} = FT(0.3) * π * ps.r0_sno^ps.ae_sno
v0_sno(ps::CMPS{FT}) where {FT} = FT(2^(9 / 4)) * ps.r0_sno^ps.ve_sno

m0_liq_coeff_TC1980(ps::CMPS{FT}) where {FT} = FT(4 / 3) * π * ps.ρ_cloud_liq

# For example: ρ_cloud_ice(ps::CMPS) = ps.ρ_cloud_ice
for var in filter(x -> x ≠ :thermo_params, fieldnames(CMPS))
    @eval $var(ps::CMPS) = ps.$var
end

# Parameters forwarded to Thermodynamics
for var in fieldnames(TDPS)
    @eval $var(ps::CMPS) = TD.Parameters.$var(thermodynamics_params(ps))
end

# Thermodynamics derived parameters
R_d(ps::CMPS) = TD.Parameters.R_d(thermodynamics_params(ps))
R_v(ps::CMPS) = TD.Parameters.R_v(thermodynamics_params(ps))
molmass_ratio(ps::CMPS) = TD.Parameters.molmass_ratio(thermodynamics_params(ps))
LH_f0(ps::CMPS) = TD.Parameters.LH_f0(thermodynamics_params(ps))
e_int_v0(ps::CMPS) = TD.Parameters.e_int_v0(thermodynamics_params(ps))
e_int_i0(ps::CMPS) = TD.Parameters.e_int_i0(thermodynamics_params(ps))
cp_d(ps::CMPS) = TD.Parameters.cp_d(thermodynamics_params(ps))
cv_d(ps::CMPS) = TD.Parameters.cv_d(thermodynamics_params(ps))
cv_v(ps::CMPS) = TD.Parameters.cv_v(thermodynamics_params(ps))
cv_l(ps::CMPS) = TD.Parameters.cv_l(thermodynamics_params(ps))
cv_i(ps::CMPS) = TD.Parameters.cv_i(thermodynamics_params(ps))


end # module
