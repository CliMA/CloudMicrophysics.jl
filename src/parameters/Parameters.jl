module Parameters

using DocStringExtensions
import CLIMAParameters
import Thermodynamics
const CP = CLIMAParameters
const TD = Thermodynamics
const TDPS = TD.Parameters.ThermodynamicsParameters

abstract type AbstractCloudMicrophysicsParameters end
Base.broadcastable(ps::AbstractCloudMicrophysicsParameters) = tuple(ps)

export CloudMicrophysicsParameters0M
export CloudMicrophysicsParametersP3

include("cmp_0m.jl")
include("cmp_p3.jl")
include("modal_nucleation.jl")

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
    α_var_time_scale_acnv::FT
    k_thrshld_stpnss::FT
    kcc_SB2006::FT
    kcr_SB2006::FT
    krr_SB2006::FT
    κrr_SB2006::FT
    xr_min_SB2006::FT
    xr_max_SB2006::FT
    νc_SB2006::FT
    ρ0_SB2006::FT
    A_phi_au_SB2006::FT
    a_phi_au_SB2006::FT
    b_phi_au_SB2006::FT
    τ0_phi_ac_SB2006::FT
    c_phi_ac_SB2006::FT
    d_sc_SB2006::FT
    Deq_br_SB2006::FT
    Dr_th_br_SB2006::FT
    kbr_SB2006::FT
    κbr_SB2006::FT
    aR_tv_SB2006::FT
    bR_tv_SB2006::FT
    cR_tv_SB2006::FT
    av_evap_SB2006::FT
    bv_evap_SB2006::FT
    α_evap_SB2006::FT
    β_evap_SB2006::FT
    N0_min_SB2006::FT
    N0_max_SB2006::FT
    λ_min_SB2006::FT
    λ_max_SB2006::FT
    q_coeff_rain_Ch2022::FT
    a1_coeff_rain_Ch2022::FT
    a2_coeff_rain_Ch2022::FT
    a3_coeff_rain_Ch2022::FT
    a3_pow_coeff_rain_Ch2022::FT
    b1_coeff_rain_Ch2022::FT
    b2_coeff_rain_Ch2022::FT
    b3_coeff_rain_Ch2022::FT
    b_rho_coeff_rain_Ch2022::FT
    c1_coeff_rain_Ch2022::FT
    c2_coeff_rain_Ch2022::FT
    c3_coeff_rain_Ch2022::FT
    As_coeff_1_Ch2022::FT
    As_coeff_2_Ch2022::FT
    As_coeff_3_Ch2022::FT
    Bs_coeff_1_Ch2022::FT
    Bs_coeff_2_Ch2022::FT
    Bs_coeff_3_Ch2022::FT
    Cs_coeff_1_Ch2022::FT
    Cs_coeff_2_Ch2022::FT
    Cs_coeff_3_Ch2022::FT
    Cs_coeff_4_Ch2022::FT
    Es_coeff_1_Ch2022::FT
    Es_coeff_2_Ch2022::FT
    Es_coeff_3_Ch2022::FT
    Fs_coeff_1_Ch2022::FT
    Fs_coeff_2_Ch2022::FT
    Fs_coeff_3_Ch2022::FT
    Gs_coeff_1_Ch2022::FT
    Gs_coeff_2_Ch2022::FT
    Gs_coeff_3_Ch2022::FT
    Si_max_Mohler2006::FT
    T_thr_Mohler2006::FT
    S0_warm_ATD_Mohler2006::FT
    S0_cold_ATD_Mohler2006::FT
    a_warm_ATD_Mohler2006::FT
    a_cold_ATD_Mohler2006::FT
    S0_warm_DD_Mohler2006::FT
    S0_cold_DD_Mohler2006::FT
    a_warm_DD_Mohler2006::FT
    a_cold_DD_Mohler2006::FT
    J_ABIFM_m_KA2013_DesertDust::FT
    J_ABIFM_c_KA2013_DesertDust::FT
    J_ABIFM_m_KA2013_Kaolinite::FT
    J_ABIFM_c_KA2013_Kaolinite::FT
    J_ABIFM_m_KA2013_Illite::FT
    J_ABIFM_c_KA2013_Illite::FT
    Koop2000_min_delta_aw::FT
    Koop2000_max_delta_aw::FT
    Koop2000_J_hom_c1::FT
    Koop2000_J_hom_c2::FT
    Koop2000_J_hom_c3::FT
    Koop2000_J_hom_c4::FT
    H2SO4_sol_T_max::FT
    H2SO4_sol_T_min::FT
    H2SO4_sol_w_2::FT
    H2SO4_sol_c1::FT
    H2SO4_sol_c2::FT
    H2SO4_sol_c3::FT
    H2SO4_sol_c4::FT
    H2SO4_sol_c5::FT
    H2SO4_sol_c6::FT
    H2SO4_sol_c7::FT
    molmass_seasalt::FT
    rho_seasalt::FT
    osm_coeff_seasalt::FT
    N_ion_seasalt::FT
    water_soluble_mass_frac_seasalt::FT
    kappa_seasalt::FT
    molmass_sulfate::FT
    rho_sulfate::FT
    osm_coeff_sulfate::FT
    N_ion_sulfate::FT
    water_soluble_mass_frac_sulfate::FT
    kappa_sulfate::FT
    α_va_BF1995::FT
    β_va_BF1995::FT
    σ_M1996::FT
    γ_M1996::FT
    p_b_n::FT
    p_b_i::FT
    u_b_n::FT
    u_b_i::FT
    v_b_n::FT
    v_b_i::FT
    w_b_n::FT
    w_b_i::FT
    p_t_n::FT
    p_t_i::FT
    u_t_n::FT
    u_t_i::FT
    v_t_n::FT
    v_t_i::FT
    w_t_n::FT
    w_t_i::FT
    p_A_n::FT
    p_A_i::FT
    a_n::FT
    a_i::FT
    a_1::FT
    a_2::FT
    a_3::FT
    a_4::FT
    a_5::FT
    Y_MTO3::FT
    Y_MTOH::FT
    k_MTO3::FT
    k_MTOH::FT
    exp_MTO3::FT
    exp_MTOH::FT
    k_H2SO4org::FT
end

Base.eltype(::CloudMicrophysicsParameters{FT}) where {FT} = FT

const CMPS = CloudMicrophysicsParameters

thermodynamics_params(cmp::CMPS) = cmp.thermo_params

# Derived parameters
N_Sc(ps::CMPS) = ps.ν_air / ps.D_vapor
m0_ice(ps::CMPS{FT}) where {FT} =
    FT(4 / 3) * π * ps.ρ_cloud_ice * ps.r0_ice^ps.me_ice
m0_rai(ps::CMPS{FT}) where {FT} =
    FT(4 / 3) * π * ps.ρ_cloud_liq * ps.r0_rai^ps.me_rai
a0_rai(ps::CMPS{FT}) where {FT} = FT(π) * ps.r0_rai^ps.ae_rai
m0_sno(ps::CMPS{FT}) where {FT} = FT(1e-1) * ps.r0_sno^ps.me_sno
a0_sno(ps::CMPS{FT}) where {FT} = FT(0.3) * π * ps.r0_sno^ps.ae_sno
v0_sno(ps::CMPS{FT}) where {FT} = FT(2^(9 / 4)) * ps.r0_sno^ps.ve_sno

m0_liq_coeff_TC1980(ps::CMPS{FT}) where {FT} = FT(4 / 3) * π * ps.ρ_cloud_liq

# For example: ρ_cloud_ice(ps::CMPS) = ps.ρ_cloud_ice
for var in filter(x -> !(x in [:thermo_params]), fieldnames(CMPS))
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
