"""
    Two-moment bulk microphysics autoconversion and accretion rates from:
      - Khairoutdinov and Kogan 2000,
      - Beheng 1994,
      - Tripoli and Cotton 1980,
      - Liu and Daum 2004.
"""
module Microphysics2M

import ..Common
const CO = Common

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

"""
    A Heaviside step function
"""
function heaviside(x::FT) where {FT <: Real}
    return FT(x > 0)
end

# autoconversion rates
"""
    conv_q_liq_to_q_rai_KK2000(param_set, q_liq, ρ, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water specific humidity
 - `ρ` - air density
 - `N_d` - cloud droplet number concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Khairoutdinov and Kogan (2000).
"""
function conv_q_liq_to_q_rai_KK2000(
    param_set::APS,
    q_liq::FT,
    ρ::FT;
    N_d::FT = 1e8,
) where {FT <: Real}

    q_liq = max(0, q_liq)

    A::FT = CMP.A_acnv_KK2000(param_set)
    a::FT = CMP.a_acnv_KK2000(param_set)
    b::FT = CMP.b_acnv_KK2000(param_set)
    c::FT = CMP.c_acnv_KK2000(param_set)

    return A * q_liq^a * N_d^b * ρ^c
end

"""
    conv_q_liq_to_q_rai_B1994(param_set, q_liq, ρ, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water specific humidity
 - `ρ` - air density
 - `N_d` - cloud droplet number concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Beheng (1994).
"""
function conv_q_liq_to_q_rai_B1994(
    param_set::APS,
    q_liq::FT,
    ρ::FT;
    N_d::FT = 1e8,
    smooth_transition::Bool = false,
) where {FT <: Real}

    q_liq = max(0, q_liq)

    C::FT = CMP.C_acnv_B1994(param_set)
    a::FT = CMP.a_acnv_B1994(param_set)
    b::FT = CMP.b_acnv_B1994(param_set)
    c::FT = CMP.c_acnv_B1994(param_set)
    N_0::FT = CMP.N_0_B1994(param_set)
    d::FT = FT(0)
    if smooth_transition
        _k::FT = CMP.k_thrshld_stpnss(param_set)
        _d_low_acnv_fraction::FT = CO.logistic_function(N_d, N_0, _k)
        _d_high_acnv_fraction::FT = 1 - _d_low_acnv_fraction
        d =
            _d_low_acnv_fraction * CMP.d_low_acnv_B1994(param_set) +
            _d_high_acnv_fraction * CMP.d_high_acnv_B1994(param_set)
    else
        d =
            N_d >= N_0 ? CMP.d_low_acnv_B1994(param_set) :
            CMP.d_high_acnv_B1994(param_set)
    end

    return C * d^a * (q_liq * ρ)^b * N_d^c / ρ
end

"""
    conv_q_liq_to_q_rai_TC1980(param_set, q_liq, ρ, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water water specific humidity
 - `ρ` - air density
 - `N_d` - cloud droplet number concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Tripoli and Cotton (1980).
"""
function conv_q_liq_to_q_rai_TC1980(
    param_set::APS,
    q_liq::FT,
    ρ::FT;
    N_d::FT = 1e8,
    smooth_transition::Bool = false,
) where {FT <: Real}
    #TODO - The original paper is actually formulated for mixing ratios, not specific humidities

    q_liq = max(0, q_liq)

    m0_liq_coeff::FT = CMP.m0_liq_coeff_TC1980(param_set)
    me_liq::FT = CMP.me_liq_TC1980(param_set)

    D::FT = CMP.D_acnv_TC1980(param_set)
    a::FT = CMP.a_acnv_TC1980(param_set)
    b::FT = CMP.b_acnv_TC1980(param_set)
    r_0::FT = CMP.r_0_acnv_TC1980(param_set)

    q_liq_threshold::FT = m0_liq_coeff * N_d / ρ * r_0^me_liq

    _output::FT = FT(0)
    if smooth_transition
        _k::FT = CMP.k_thrshld_stpnss(param_set)
        _output = CO.logistic_function(q_liq, q_liq_threshold, _k)
    else
        _output = heaviside(q_liq - q_liq_threshold)
    end
    return D * q_liq^a * N_d^b * _output
end

"""
    conv_q_liq_to_q_rai_LD2004(param_set, q_liq, ρ, N_d)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water specific humidity
 - `ρ` - air density
 - `N_d` - cloud droplet number concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Liu and Daum (2004).
"""
function conv_q_liq_to_q_rai_LD2004(
    param_set::APS,
    q_liq::FT,
    ρ::FT;
    N_d::FT = 1e8,
    smooth_transition::Bool = false,
) where {FT <: Real}

    if q_liq <= eps(FT)
        return FT(0)
    else
        ρ_w::FT = CMP.ρ_cloud_liq(param_set)
        R_6C_0::FT = CMP.R_6C_coeff_LD2004(param_set)
        E_0::FT = CMP.E_0_LD2004(param_set)

        # Mean volume radius in microns (assuming spherical cloud droplets)
        r_vol = (3 * (q_liq * ρ) / 4 / π / ρ_w / N_d)^(1 / 3) * 1e6

        # Assumed size distribution: modified gamma distribution
        β_6::FT = ((r_vol + 3) / r_vol)^(1 / 3)
        E::FT = E_0 * β_6^6
        R_6::FT = β_6 * r_vol
        R_6C::FT = R_6C_0 / (q_liq * ρ)^(1 / 6) / R_6^(1 / 2)

        _output::FT = FT(0)
        if smooth_transition
            _k::FT = CMP.k_thrshld_stpnss(param_set)
            _output = CO.logistic_function(R_6, R_6C, _k)
        else
            _output = heaviside(R_6 - R_6C)
        end
        return E * (q_liq * ρ)^3 / N_d / ρ * _output
    end
end

# accretion rates

"""
    accretion_KK2000(param_set, q_liq, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

 Returns the accretion rate of rain, parametrized
 following Khairoutdinov and Kogan (2000).
"""
function accretion_KK2000(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)

    A::FT = CMP.A_acc_KK2000(param_set)
    a::FT = CMP.a_acc_KK2000(param_set)
    b::FT = CMP.b_acc_KK2000(param_set)

    return A * (q_liq * q_rai)^a * ρ^b
end

"""
    accretion_B1994(param_set, q_liq, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

 Returns the accretion rate of rain, parametrized
 following Beheng (1994).
"""
function accretion_B1994(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)

    A::FT = CMP.A_acc_B1994(param_set)

    return A * q_liq * ρ * q_rai
end

"""
    accretion_TC1980(param_set, q_liq, q_rai)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

 Returns the accretion rate of rain, parametrized
 following Tripoli and Cotton (1980).
"""
function accretion_TC1980(
    param_set::APS,
    q_liq::FT,
    q_rai::FT,
) where {FT <: Real}
    #TODO - The original paper is actually formulated for mixing ratios, not specific humidities

    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)

    A::FT = CMP.A_acc_TC1980(param_set)

    return A * q_liq * q_rai
end

end
