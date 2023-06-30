"""
Double-moment bulk microphysics parametrizations including:
 - autoconversion, accretion, self-collection, breakup, mean terminal velocity of raindrops and rain
    evaporation rates from Seifert and Beheng 2006,
 - additional double-moment bulk microphysics autoconversion and accretion rates
   from: Khairoutdinov and Kogan 2000, Beheng 1994, Tripoli and Cotton 1980, and
   Liu and Daum 2004.
"""
module Microphysics2M

import ..Common as CO

import ..CommonTypes as CT

import SpecialFunctions as SF

import Thermodynamics as TD

import ..Parameters as CMP
const APS = CMP.AbstractCloudMicrophysicsParameters

export autoconversion
export accretion
export liquid_self_collection
export autoconversion_and_liquid_self_collection
export rain_self_collection
export rain_breakup
export rain_self_collection_and_breakup
export rain_terminal_velocity
export rain_evaporation
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

# Double-moment bulk microphysics autoconversion, accretion, self-collection, breakup,
# mean terminal velocity of raindrops, and rain evaporation rates from Seifert and Beheng 2001

"""
    raindrops_limited_vars(param_set, q_rai, ρ, N_rai)

 - `param_set` - abstract set with Earth parameters
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` raindrops number density

Returns a named tupple containing the mean mass of raindrops, xr, and the rate parameter of the assumed
size distribution of raindrops (based on drops diameter), λr, limited within prescribed ranges
"""
function raindrops_limited_vars(
    param_set::APS,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT <: Real}

    xr_min::FT = CMP.xr_min_SB2006(param_set)
    xr_max::FT = CMP.xr_max_SB2006(param_set)
    N0_min::FT = CMP.N0_min_SB2006(param_set)
    N0_max::FT = CMP.N0_max_SB2006(param_set)
    λ_min::FT = CMP.λ_min_SB2006(param_set)
    λ_max::FT = CMP.λ_max_SB2006(param_set)
    ρw::FT = CMP.ρ_cloud_liq(param_set)

    L_rai = ρ * q_rai
    xr_0 = L_rai / N_rai
    xr_hat = max(xr_min, min(xr_max, xr_0))
    N0 = max(N0_min, min(N0_max, N_rai * (FT(π) * ρw / xr_hat)^FT(1 / 3)))
    λr = max(λ_min, min(λ_max, (FT(π) * ρw * N0 / L_rai)^FT(1 / 4)))
    xr = max(xr_min, min(xr_max, L_rai * λr / N0))

    return (; λr, xr)
end

"""
    autoconversion(param_set, scheme, q_liq, q_rai, ρ, N_liq)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between cloud droplets (autoconversion) for `scheme == SB2006Type`
"""
function autoconversion(
    param_set::APS,
    scheme::CT.SB2006Type,
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

    kcc::FT = CMP.kcc_SB2006(param_set)
    νc::FT = CMP.νc_SB2006(param_set)
    x_star::FT = CMP.xr_min_SB2006(param_set)
    ρ0::FT = CMP.ρ0_SB2006(param_set)
    A::FT = CMP.A_phi_au_SB2006(param_set)
    a::FT = CMP.a_phi_au_SB2006(param_set)
    b::FT = CMP.b_phi_au_SB2006(param_set)

    L_liq = ρ * q_liq
    x_liq = min(x_star, L_liq / N_liq)
    q_rai = max(FT(0), q_rai)
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_au = A * τ^a * (FT(1) - τ^a)^b

    dL_rai_dt =
        kcc / 20 / x_star * (νc + 2) * (νc + 4) / (νc + 1)^2 *
        L_liq^2 *
        x_liq^2 *
        (1 + ϕ_au / (1 - τ)^2) *
        ρ0 / ρ
    dN_rai_dt = dL_rai_dt / x_star
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
    accretion(param_set, scheme, q_liq, q_rai, ρ, N_liq)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment accretion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between raindrops and cloud droplets (accretion) for `scheme == SB2006Type`
"""
function accretion(
    param_set::APS,
    scheme::CT.SB2006Type,
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

    kcr::FT = CMP.kcr_SB2006(param_set)
    τ0::FT = CMP.τ0_phi_ac_SB2006(param_set)
    ρ0::FT = CMP.ρ0_SB2006(param_set)
    c::FT = CMP.c_phi_ac_SB2006(param_set)

    L_liq = ρ * q_liq
    L_rai = ρ * q_rai
    x_liq = L_liq / N_liq
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_ac = (τ / (τ + τ0))^c

    dL_rai_dt = kcr * L_liq * L_rai * ϕ_ac * sqrt(ρ0 / ρ)
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
that produce larger cloud droplets (self-collection) for `scheme == SB2006Type`
"""
function liquid_self_collection(
    param_set::APS,
    scheme::CT.SB2006Type,
    q_liq::FT,
    ρ::FT,
    dN_liq_dt_au::FT,
) where {FT <: Real}

    if q_liq < eps(FT)
        return FT(0)
    end

    kcc::FT = CMP.kcc_SB2006(param_set)
    ρ0::FT = CMP.ρ0_SB2006(param_set)
    νc::FT = CMP.νc_SB2006(param_set)

    L_liq = ρ * q_liq

    dN_liq_dt_sc =
        -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * L_liq^2 - dN_liq_dt_au

    return dN_liq_dt_sc
end

"""
    autoconversion_and_liquid_self_collection(param_set, scheme, q_liq, q_rai, ρ, N_liq)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a named tupple containing a LiqRaiRates object for the autoconversion rate and
the liquid self-collection rate for `scheme == SB2006Type`
"""
function autoconversion_and_liquid_self_collection(
    param_set::APS,
    scheme::CT.SB2006Type,
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT <: Real}

    au = autoconversion(param_set, scheme, q_liq, q_rai, ρ, N_liq)
    sc = liquid_self_collection(param_set, scheme, q_liq, ρ, au.dN_liq_dt)

    return (; au, sc)
end

"""
    rain_self_collection(param_set, scheme, q_rai, ρ, N_rai)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment rain self-collection parameterization
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns the raindrops number density tendency due to collisions of raindrops
that produce larger raindrops (self-collection) for `scheme == SB2006Type`
"""
function rain_self_collection(
    param_set::APS,
    scheme::CT.SB2006Type,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT <: Real}

    if q_rai < eps(FT)
        return FT(0)
    end

    krr::FT = CMP.krr_SB2006(param_set)
    κrr::FT = CMP.κrr_SB2006(param_set)
    d::FT = CMP.d_sc_SB2006(param_set)
    ρw::FT = CMP.ρ_cloud_liq(param_set)
    ρ0::FT = CMP.ρ0_SB2006(param_set)

    L_rai = ρ * q_rai
    λr =
        raindrops_limited_vars(param_set, q_rai, ρ, N_rai).λr *
        (SF.gamma(FT(4)) / FT(π) / ρw)^FT(1 / 3)

    dN_rai_dt_sc = -krr * N_rai * L_rai * sqrt(ρ0 / ρ) * (1 + κrr / λr)^d

    return dN_rai_dt_sc
end

"""
    rain_breakup(param_set, scheme, q_rai, ρ, dN_rai_dt_sc)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density
 - `dN_rai_dt_sc` - rate of change of raindrops number density due to self-collection

Returns the raindrops number density tendency due to breakup of raindrops
that produce smaller raindrops for `scheme == SB2006Type`
"""
function rain_breakup(
    param_set::APS,
    scheme::CT.SB2006Type,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
    dN_rai_dt_sc::FT,
) where {FT <: Real}

    if q_rai < eps(FT)
        return FT(0)
    end

    ρw::FT = CMP.ρ_cloud_liq(param_set)
    Deq::FT = CMP.Deq_br_SB2006(param_set)
    Dr_th::FT = CMP.Dr_th_br_SB2006(param_set)
    kbr::FT = CMP.kbr_SB2006(param_set)
    κbr::FT = CMP.κbr_SB2006(param_set)

    xr = raindrops_limited_vars(param_set, q_rai, ρ, N_rai).xr
    Dr = (xr * 6 / FT(π) / ρw)^FT(1 / 3)
    ΔD = Dr - Deq
    phi_br =
        (Dr < Dr_th) ? FT(-1) : ((ΔD <= 0) ? kbr * ΔD : 2 * (exp(κbr * ΔD) - 1))
    dN_rai_dt_br = -(phi_br + 1) * dN_rai_dt_sc

    return dN_rai_dt_br
end

"""
    rain_sef_collection_and_breakup(param_set, scheme, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns a named tupple containing the raindrops self-collection and breakup rates
for `scheme == SB2006Type`
"""
function rain_self_collection_and_breakup(
    param_set::APS,
    scheme::CT.SB2006Type,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT <: Real}

    sc = rain_self_collection(param_set, scheme, q_rai, ρ, N_rai)
    br = rain_breakup(param_set, scheme, q_rai, ρ, N_rai, sc)

    return (; sc, br)
end

"""
    rain_terminal_velocity(param_set, scheme, velo_scheme, q_rai, ρ, N_rai)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment parameterization
 - `velo_scheme` - type for terminal velocity parameterization
 - `q_rai` - rain water specific humidity [kg/kg]
 - `ρ` - air density [kg/m^3]
 - `N_rai` - raindrops number density [1/m^3]

Returns a tuple containing the number and mass weigthed mean fall velocities of raindrops in [m/s].
Assuming an exponential size distribution from Seifert and Beheng 2006 for `scheme == SB2006Type`
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical relation similar to Rogers (1993) for `velo_scheme == SB2006VelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171 for `velo_scheme == Chen2022Type`
"""
function rain_terminal_velocity(
    param_set::APS,
    scheme::CT.SB2006Type,
    velo_scheme::CT.SB2006VelType,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT <: Real}
    if q_rai < eps(FT)
        return (FT(0), FT(0))
    end

    ρ0::FT = CMP.ρ0_SB2006(param_set)
    aR::FT = CMP.aR_tv_SB2006(param_set)
    bR::FT = CMP.bR_tv_SB2006(param_set)
    cR::FT = CMP.cR_tv_SB2006(param_set)

    λr = raindrops_limited_vars(param_set, q_rai, ρ, N_rai).λr
    vt0 = max(FT(0), sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)))
    vt1 = max(FT(0), sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)^FT(4)))

    return (vt0, vt1)
end
function rain_terminal_velocity(
    param_set::APS,
    scheme::CT.SB2006Type,
    velo_scheme::CT.Chen2022Type,
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT <: Real}
    if q_rai < eps(FT)
        return (FT(0), FT(0))
    end
    # coefficients from Table B1 from Chen et. al. 2022
    aiu, bi, ciu = CO.Chen2022_vel_coeffs(param_set, CT.RainType(), ρ)
    # size distribution parameter
    λ = raindrops_limited_vars(param_set, q_rai, ρ, N_rai).λr

    # eq 20 from Chen et al 2022
    vt0 = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, λ, 0))
    vt3 = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, λ, 3))

    vt0 = max(FT(0), vt0)
    vt3 = max(FT(0), vt3)
    # It should be (ϕ^κ * vt0, ϕ^κ * vt3), but for rain drops ϕ = 1 and κ = 0
    return (vt0, vt3)
end

"""
    rain_evaporation(param_set, scheme, q, q_rai, ρ, N_rai, T)

 - `param_set` - abstract set with Earth parameters
 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q` - phase partition
 - `q_rai` - rain specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density
 - `T` - air temperature

Returns a tupple containing the tendency of raindrops number density and rain water
specific humidity due to rain rain_evaporation, assuming a power law velocity relation for
fall velocity of individual drops and an exponential size distribution, for `scheme == SB2006Type`
"""
function rain_evaporation(
    param_set::APS,
    scheme::CT.SB2006Type,
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
    T::FT,
) where {FT <: Real}

    evap_rate_0 = FT(0)
    evap_rate_1 = FT(0)
    thermo_params = CMP.thermodynamics_params(param_set)
    S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && S < FT(0))

        ν_air::FT = CMP.ν_air(param_set)
        D_vapor::FT = CMP.D_vapor(param_set)
        ρw::FT = CMP.ρ_cloud_liq(param_set)
        ρ0::FT = CMP.ρ0_SB2006(param_set)
        av::FT = CMP.av_evap_SB2006(param_set)
        bv::FT = CMP.bv_evap_SB2006(param_set)
        α::FT = CMP.α_evap_SB2006(param_set)
        β::FT = CMP.β_evap_SB2006(param_set)
        x_star::FT = CMP.xr_min_SB2006(param_set)

        G::FT = CO.G_func(param_set, T, TD.Liquid())

        xr = raindrops_limited_vars(param_set, q_rai, ρ, N_rai).xr
        Dr = (FT(6) / FT(π) / ρw)^FT(1 / 3) * xr^FT(1 / 3)
        t_star = (FT(6) * x_star / xr)^FT(1 / 3)
        a_vent_0::FT = av * SF.gamma(FT(-1), t_star) / FT(6)^FT(-2 / 3)
        b_vent_0::FT =
            bv * SF.gamma(FT(-1 / 2) + FT(3 / 2) * β, t_star) /
            FT(6)^FT(β / 2 - 1 / 2)
        a_vent_1::FT = av * SF.gamma(FT(2)) / FT(6)^FT(1 / 3)
        b_vent_1::FT =
            bv * SF.gamma(FT(5 / 2) + FT(3 / 2) * β) / FT(6)^FT(β / 2 + 1 / 2)

        N_Re = α * xr^β * sqrt(ρ0 / ρ) * Dr / ν_air
        Fv0::FT = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
        Fv1::FT = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

        evap_rate_0 = min(0, 2 * FT(π) * G * S * N_rai * Dr * Fv0 / xr)
        evap_rate_1 = min(0, 2 * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)
    end

    return (evap_rate_0, evap_rate_1)
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
    N_d::FT = FT(1e8),
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
    N_d::FT = FT(1e8),
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
    N_d::FT = FT(1e8),
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
    N_d::FT = FT(1e8),
    smooth_transition::Bool = false,
) where {FT <: Real}

    if q_liq <= eps(FT)
        return FT(0)
    else
        ρ_w::FT = CMP.ρ_cloud_liq(param_set)
        R_6C_0::FT = CMP.R_6C_coeff_LD2004(param_set)
        E_0::FT = CMP.E_0_LD2004(param_set)

        # Mean volume radius in microns (assuming spherical cloud droplets)
        r_vol = (3 * (q_liq * ρ) / 4 / FT(π) / ρ_w / N_d)^FT(1 / 3) * 1e6

        # Assumed size distribution: modified gamma distribution
        β_6::FT = ((r_vol + 3) / r_vol)^FT(1 / 3)
        E::FT = E_0 * β_6^6
        R_6::FT = β_6 * r_vol
        R_6C::FT = R_6C_0 / (q_liq * ρ)^FT(1 / 6) / R_6^FT(1 / 2)

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
