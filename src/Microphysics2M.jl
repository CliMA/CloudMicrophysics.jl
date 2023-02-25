"""
Double-moment bulk microphysics parametrizations including:
 - autoconversion, accretion, and self-collection rates from Seifert and Beheng 2001,
 - additional double-moment bulk microphysics autoconversion and accretion rates
   from: Khairoutdinov and Kogan 2000, Beheng 1994, Tripoli and Cotton 1980, and 
   Liu and Daum 2004.
"""
module Microphysics2M

import ..Common
const CO = Common

import ..CommonTypes
const CT = CommonTypes

import SpecialFunctions
const SF = SpecialFunctions

import Thermodynamics
const TD = Thermodynamics

import ..Parameters
const CMP = Parameters
const APS = CMP.AbstractCloudMicrophysicsParameters

export autoconversion
export accretion
export liquid_self_collection
export rain_self_collection
export conv_q_liq_to_q_rai

"""
A structure containing the rates of change of the specific humidities and number
densities of liquid and rain water.
"""
Base.@kwdef struct LiqRaiRates{FT <: Real}
    "Rate of change of the liquid water specific humidity"
    dq_liq_dt::FT
    "Rate of change of the liquid water number density"
    dN_liq_dt::FT
    "Rate of change of the rain water specific humidity"
    dq_rai_dt::FT
    "Rate of change of the rain water number density"
    dN_rai_dt::FT
end

# Double-moment bulk microphysics autoconversion, accretion, and self-collection
# rates from Seifert and Beheng 2001

"""
    autoconversion(param_set, scheme, q_liq, q_rai, ρ, N_liq)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between cloud droplets (autoconversion) for `scheme == SB2001Type`
"""
function autoconversion(
    param_set::APS,
    scheme::CT.SB2001Type,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT <: Real}

    if q_liq < eps(FT)
        return LiqRaiRates(
            dq_liq_dt = FT(0),
            dN_liq_dt = FT(0),
            dq_rai_dt = FT(0),
            dN_rai_dt = FT(0),
        )
    end

    kc::FT = CMP.kc_SB2001(param_set)
    ν::FT = CMP.ν_SB2001(param_set)
    xstar::FT = CMP.xstar_SB2001(param_set)
    A::FT = CMP.A_phi_au_SB2001(param_set)
    a::FT = CMP.a_phi_au_SB2001(param_set)
    b::FT = CMP.b_phi_au_SB2001(param_set)

    L_liq = ρ * q_liq
    x_liq = L_liq / N_liq
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_au = A * τ^a * (FT(1) - τ^a)^b

    dL_rai_dt =
        kc / 20 / xstar * (ν + 2) * (ν + 4) / (ν + 1)^2 *
        L_liq^2 *
        x_liq^2 *
        (1 + ϕ_au / (1 - τ)^2)
    dN_rai_dt = dL_rai_dt / xstar
    dL_liq_dt = -dL_rai_dt
    dN_liq_dt = -2 * dN_rai_dt

    return LiqRaiRates(
        dq_liq_dt = dL_liq_dt / ρ,
        dN_liq_dt = dN_liq_dt,
        dq_rai_dt = dL_rai_dt / ρ,
        dN_rai_dt = dN_rai_dt,
    )
end

"""
    accretion(param_set, scheme, q_liq, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment accretion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between raindrops and cloud droplets (accretion) for `scheme == SB2001Type`
"""
function accretion(
    param_set::APS,
    scheme::CT.SB2001Type,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT <: Real}

    if q_liq < eps(FT) || q_rai < eps(FT)
        return LiqRaiRates(
            dq_liq_dt = FT(0),
            dN_liq_dt = FT(0),
            dq_rai_dt = FT(0),
            dN_rai_dt = FT(0),
        )
    end

    kr::FT = CMP.kr_SB2001(param_set)
    τ0::FT = CMP.τ_0_phi_ac_SB2001(param_set)
    c::FT = CMP.c_phi_ac_SB2001(param_set)

    L_liq = ρ * q_liq
    L_rai = ρ * q_rai
    x_liq = L_liq / N_liq
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_ac = (τ / (τ + τ0))^c

    dL_rai_dt = kr * L_liq * L_rai * ϕ_ac
    dN_rai_dt = FT(0)
    dL_liq_dt = -dL_rai_dt
    dN_liq_dt = dL_liq_dt / x_liq

    return LiqRaiRates(
        dq_liq_dt = dL_liq_dt / ρ,
        dN_liq_dt = dN_liq_dt,
        dq_rai_dt = dL_rai_dt / ρ,
        dN_rai_dt = dN_rai_dt,
    )
end

"""
    liquid_self_collection(param_set, scheme, q_liq, ρ, dN_liq_dt_au)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_liq` - cloud water specific humidity
 - `ρ` - air density
 - `dN_liq_dt_au` - rate of change of cloud droplets number density due to autoconversion

Returns the cloud droplets number density tendency due to collisions of cloud droplets
that produce larger cloud droplets (self-collection) for `scheme == SB2001Type`
"""
function liquid_self_collection(
    param_set::APS,
    scheme::CT.SB2001Type,
    q_liq::FT,
    ρ::FT,
    dN_liq_dt_au::FT,
) where {FT <: Real}

    if q_liq < eps(FT)
        return FT(0)
    end

    kc::FT = CMP.kc_SB2001(param_set)
    ν::FT = CMP.ν_SB2001(param_set)

    L_liq = ρ * q_liq

    dN_liq_dt_sc = -kc * (ν + 2) / (ν + 1) * L_liq^2 - dN_liq_dt_au

    return dN_liq_dt_sc
end

"""
    rain_self_collection(param_set, scheme, q_rai, ρ, N_rai)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain self-collection parameterization
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns the raindrops number density tendency due to collisions of raindrops
that produce larger raindrops (self-collection) for `scheme == SB2001Type`
"""
function rain_self_collection(
    param_set::APS,
    scheme::CT.SB2001Type,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT <: Real}

    if q_rai < eps(FT)
        return FT(0)
    end

    kr::FT = CMP.kr_SB2001(param_set)

    L_rai = ρ * q_rai

    dN_rai_dt_sc = -kr * N_rai * L_rai

    return dN_rai_dt_sc
end

# Additional double moment autoconversion and accretion parametrizations:
# - Khairoutdinov and Kogan (2000)
# - Beheng (1994)
# - Tripoli and Cotton (1980)
# - Liu and Daum (2004)

"""
    conv_q_liq_to_q_rai(param_set, scheme, q_liq, ρ; N_d, smooth_transition)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific humidity
 - `ρ` - air density
 - `N_d` - prescribed cloud droplet number concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following:
 - Khairoutdinov and Kogan (2000) for `scheme == KK2000Type`
 - Beheng (1994) for `scheme == B1994Type`
 - Tripoli and Cotton (1980) for `scheme == TC1980Type`
 - Liu and Daum (2004) for `scheme ==LD2004Type`

`N_d` is an optional argument with the default value of 100 cm-3

The `Beheng1994Type`, `TC1980Type` and `LD2004Type` of schemes
additionally accept `smooth_transition` flag that
smoothes their thershold behaviour if set to `true`.
The default value is `false`.
"""
function conv_q_liq_to_q_rai(
    param_set::APS,
    scheme::CT.KK2000Type,
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
function conv_q_liq_to_q_rai(
    param_set::APS,
    scheme::CT.B1994Type,
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
function conv_q_liq_to_q_rai(
    param_set::APS,
    scheme::CT.TC1980Type,
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
        _output = CO.heaviside(q_liq - q_liq_threshold)
    end
    return D * q_liq^a * N_d^b * _output
end
function conv_q_liq_to_q_rai(
    param_set::APS,
    scheme::CT.LD2004Type,
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
            _output = CO.heaviside(R_6 - R_6C)
        end
        return E * (q_liq * ρ)^3 / N_d / ρ * _output
    end
end

"""
    accretion(param_set, scheme, q_liq, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain accretion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density (for `KK2000Type` and `Beheng1994Type`)

 Returns the accretion rate of rain, parametrized following
 - Khairoutdinov and Kogan (2000) for `scheme == KK2000Type`
 - Beheng (1994) for `scheme == B1994Type`
 - Tripoli and Cotton (1980) for `scheme == TC1980Type`
"""
function accretion(
    param_set::APS,
    scheme::CT.KK2000Type,
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
function accretion(
    param_set::APS,
    scheme::CT.B1994Type,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)

    A::FT = CMP.A_acc_B1994(param_set)

    return A * q_liq * ρ * q_rai
end
function accretion(
    param_set::APS,
    scheme::CT.TC1980Type,
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
