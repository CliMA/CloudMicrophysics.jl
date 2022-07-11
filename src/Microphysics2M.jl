"""
    Two-moment bulk microphysics autoconversion and accretion rates

"""
module Microphysics2M

import SpecialFunctions
const SF = SpecialFunctions

import ..Parameters
const CMP = Parameters
const APS = CMP.AbstractCloudMicrophysicsParameters

export conv_q_liq_to_q_rai_KK2000
export conv_q_liq_to_q_rai_B1994
export conv_q_liq_to_q_rai_TC1980
export conv_q_liq_to_q_rai_LD2004
export accretion_KK2000
export accretion_B1994
export accretion_TC198

# autoconversion rates
function conv_q_liq_to_q_rai_KK2000(param_set::APS, q_liq::FT, ρ::FT; N_d::FT = 1e8) where {FT <: Real}

    q_liq = max(0.0, q_liq)

    A::FT = CMP.A_acnv_KK2000(param_set)
    a::FT = CMP.a_acnv_KK2000(param_set)
    b::FT = CMP.b_acnv_KK2000(param_set)
    c::FT = CMP.c_acnv_KK2000(param_set)

    return A * (q_liq * ρ)^a * N_d^b * ρ^c / ρ
end

function conv_q_liq_to_q_rai_B1994(param_set::APS, q_liq::FT, ρ::FT; N_d::FT = 1e8) where {FT <: Real}

    q_liq = max(0.0, q_liq)

    C::FT = CMP.C_acnv_B1994(param_set)
    a::FT = CMP.a_acnv_B1994(param_set)
    b::FT = CMP.b_acnv_B1994(param_set)
    c::FT = CMP.c_acnv_B1994(param_set)
    N_0::FT = CMP.N_0_B1994(param_set)
    d::FT = N_d::FT >= N_0 ? CMP.d_low_acnv_B1994(param_set) : CMP.d_high_acnv_B1994(param_set)

    return C * d^a * (q_liq * ρ)^b * N_d^c / ρ
end

function conv_q_liq_to_q_rai_TC1980(param_set::APS, q_liq::FT, ρ::FT; N_d::FT = 1e8) where {FT <: Real}

    q_liq = max(0.0, q_liq)

    m0_liq_coeff::FT = CMP.m0_liq_coeff_TC1980(param_set)
    me_liq::FT = CMP.me_liq_TC1980(param_set)

    D::FT = CMP.D_acnv_TC1980(param_set)
    a::FT = CMP.a_acnv_TC1980(param_set)
    b::FT = CMP.b_acnv_TC1980(param_set)
    r_0::FT = CMP.r_0_acnv_TC1980(param_set)

    #q_liq_threshold::FT = 4/3 * pi * CMP.ρ_cloud_liq(param_set) * N_d * (7 * 1e-6)^3
    #return 3268.0 * (q_liq * 1.2)^(7.0/3) / N_d^(1.0/3) * max(0.0, (q_liq * 1.2) - q_liq_threshold) * 1e3
    q_liq_threshold::FT = m0_liq_coeff * N_d * r_0^me_liq
    return D * (q_liq * ρ)^a * N_d^b * max(0.0, (q_liq * ρ) - q_liq_threshold) * 1e3 / ρ
end

function conv_q_liq_to_q_rai_LD2004(param_set::APS, q_liq::FT, ρ::FT; N_d::FT = 1e8) where {FT <: Real}

    q_liq = max(0.0, q_liq)

    ρ_w::FT = CMP.ρ_cloud_liq(param_set)
    R_6C_0::FT = CMP.R_6C_coeff_LD2004(param_set)
    E_0::FT = CMP.E_0_LD2004(param_set)

    # Mean volume radius in microns (assuming spherical cloud droplets)
    r_vol = (3 * (q_liq * ρ) / 4 / π / ρ_w / N_d)^(1/3) * 1e6

    # Assumed size distribution: modified gamma distribution
    β_6::FT = ((r_vol + 3) / r_vol)^(1/3)
    E::FT = 1.08e10 * β_6^6
    R_6::FT = β_6 * r_vol
    R_6C::FT = R_6C_0 / (q_liq * ρ)^(1/6) / R_6^(1/2)

    return E * (q_liq * ρ)^3 / N_d * max(0.0, R_6 - R_6C) / ρ
end

# accretion rates
function accretion_KK2000(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    q_liq = max(0.0, q_liq)
    q_rai = max(0.0, q_rai)

    A::FT = CMP.A_acc_KK2000(param_set)
    a::FT = CMP.a_acc_KK2000(param_set)
    b::FT = CMP.b_acc_KK2000(param_set)

    return A * (q_liq * ρ * q_rai * ρ)^a * ρ^b / ρ
end

function accretion_B1994(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    q_liq = max(0.0, q_liq)
    q_rai = max(0.0, q_rai)

    A::FT = CMP.A_acc_B1994(param_set)

    return A * q_liq * ρ * q_rai
end

function accretion_TC1980(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    q_liq = max(0.0, q_liq)
    q_rai = max(0.0, q_rai)

    A::FT = CMP.A_acc_TC1980(param_set)

    return A * q_liq * ρ * q_rai
end

end
