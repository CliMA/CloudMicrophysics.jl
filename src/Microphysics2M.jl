"""
Double-moment bulk microphysics parametrizations including:
 - autoconversion, accretion, self-collection, breakup, mean terminal velocity of raindrops and rain
    evaporation rates from Seifert and Beheng 2006,
 - additional double-moment bulk microphysics autoconversion and accretion rates
   from: Khairoutdinov and Kogan 2000, Beheng 1994, Tripoli and Cotton 1980, and
   Liu and Daum 2004.
"""
module Microphysics2M

import SpecialFunctions as SF

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import ..Common as CO
import ..Parameters as CMP

export autoconversion,
    accretion,
    liquid_self_collection,
    autoconversion_and_liquid_self_collection,
    rain_terminal_velocity,
    conv_q_liq_to_q_rai,
    rain_evaporation,
    rain_self_collection,
    rain_breakup,
    rain_self_collection_and_breakup

"""
A structure containing the rates of change of the specific humidities and number
densities of liquid and rain water.
"""
@kwdef struct LiqRaiRates{FT}
    "Rate of change of the liquid water specific humidity"
    dq_liq_dt::FT = FT(0)
    "Rate of change of the liquid water number density"
    dN_liq_dt::FT = FT(0)
    "Rate of change of the rain water specific humidity"
    dq_rai_dt::FT = FT(0)
    "Rate of change of the rain water number density"
    dN_rai_dt::FT = FT(0)
end

# Double-moment bulk microphysics autoconversion, accretion, self-collection, breakup,
# mean terminal velocity of raindrops, and rain evaporation rates from Seifert and Beheng 2001

"""
    raindrops_limited_vars(pdf, q_rai, ρ, N_rai)

 - `pdf` - a struct with SB2006 size distribution parameters
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` raindrops number density

Returns a named tupple containing the mean mass of raindrops, xr, and the rate parameter of the assumed
size distribution of raindrops (based on drops diameter), λr, limited within prescribed ranges
"""
function raindrops_limited_vars(
    pdf::CMP.ParticlePDF_SB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}
    (; xr_min, xr_max, N0_min, N0_max, λ_min, λ_max, ρw) = pdf

    L_rai = ρ * q_rai
    xr_0 = L_rai / N_rai
    xr_hat = max(xr_min, min(xr_max, xr_0))
    N0 = max(N0_min, min(N0_max, N_rai * (FT(π) * ρw / xr_hat)^FT(1 / 3)))
    λr = max(λ_min, min(λ_max, (FT(π) * ρw * N0 / L_rai)^FT(1 / 4)))
    xr = max(xr_min, min(xr_max, L_rai * λr / N0))

    return (; λr, xr)
end

"""
    autoconversion(scheme, q_liq, q_rai, ρ, N_liq)

 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between cloud droplets (autoconversion) for `scheme == SB2006Type`
"""
function autoconversion(
    acnv::CMP.AcnvSB2006{FT},
    q_liq,
    q_rai,
    ρ,
    N_liq,
) where {FT}

    if q_liq < eps(FT)
        return LiqRaiRates{FT}()
    end

    (; kcc, νc, x_star, ρ0, A, a, b) = acnv

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
    accretion(scheme, q_liq, q_rai, ρ, N_liq)

 - `scheme` - type for 2-moment accretion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between raindrops and cloud droplets (accretion) for `scheme == SB2006Type`
"""
function accretion((; accr)::CMP.SB2006{FT}, q_liq, q_rai, ρ, N_liq) where {FT}

    if q_liq < eps(FT) || q_rai < eps(FT)
        return LiqRaiRates{FT}()
    end

    (; kcr, τ0, ρ0, c) = accr
    L_liq = ρ * q_liq
    L_rai = ρ * q_rai
    x_liq = L_liq / N_liq
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_ac = (τ / (τ + τ0))^c

    dL_rai_dt = kcr * L_liq * L_rai * ϕ_ac * sqrt(ρ0 / ρ)
    dN_rai_dt = zero(N_liq)
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
    liquid_self_collection(scheme, q_liq, ρ, dN_liq_dt_au)

 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_liq` - cloud water specific humidity
 - `ρ` - air density
 - `dN_liq_dt_au` - rate of change of cloud droplets number density due to autoconversion

Returns the cloud droplets number density tendency due to collisions of cloud droplets
that produce larger cloud droplets (self-collection) for `scheme == SB2006Type`
"""
function liquid_self_collection(
    acnv::CMP.AcnvSB2006{FT},
    q_liq::FT,
    ρ::FT,
    dN_liq_dt_au::FT,
) where {FT}

    if q_liq < eps(FT)
        return FT(0)
    end
    (; kcc, ρ0, νc) = acnv

    L_liq = ρ * q_liq

    dN_liq_dt_sc =
        -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * L_liq^2 - dN_liq_dt_au

    return dN_liq_dt_sc
end

"""
    autoconversion_and_liquid_self_collection(scheme, q_liq, q_rai, ρ, N_liq)

 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a named tupple containing a LiqRaiRates object for the autoconversion rate and
the liquid self-collection rate for `scheme == SB2006Type`
"""
function autoconversion_and_liquid_self_collection(
    (; acnv)::CMP.SB2006{FT},
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT}

    au = autoconversion(acnv, q_liq, q_rai, ρ, N_liq)
    sc = liquid_self_collection(acnv, q_liq, ρ, au.dN_liq_dt)

    return (; au, sc)
end

"""
    rain_self_collection(scheme, q_rai, ρ, N_rai)

 - `scheme` - type for 2-moment rain self-collection parameterization
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns the raindrops number density tendency due to collisions of raindrops
that produce larger raindrops (self-collection) for `scheme == SB2006Type`
"""
function rain_self_collection(
    pdf::CMP.ParticlePDF_SB2006{FT},
    self::CMP.SelfColSB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}

    if q_rai < eps(FT)
        return FT(0)
    end

    (; krr, κrr, d) = self
    (; ρ0, ρw) = pdf

    L_rai = ρ * q_rai
    λr =
        raindrops_limited_vars(pdf, q_rai, ρ, N_rai).λr *
        (SF.gamma(FT(4)) / FT(π) / ρw)^FT(1 / 3)

    dN_rai_dt_sc = -krr * N_rai * L_rai * sqrt(ρ0 / ρ) * (1 + κrr / λr)^d

    return dN_rai_dt_sc
end

"""
    rain_breakup(scheme, q_rai, ρ, dN_rai_dt_sc)

 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density
 - `dN_rai_dt_sc` - rate of change of raindrops number density due to self-collection

Returns the raindrops number density tendency due to breakup of raindrops
that produce smaller raindrops for `scheme == SB2006Type`
"""
function rain_breakup(
    pdf::CMP.ParticlePDF_SB2006{FT},
    brek::CMP.BreakupSB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
    dN_rai_dt_sc::FT,
) where {FT}

    if q_rai < eps(FT)
        return FT(0)
    end
    (; Deq, Dr_th, kbr, κbr) = brek
    ρw = pdf.ρw
    xr = raindrops_limited_vars(pdf, q_rai, ρ, N_rai).xr
    Dr = (xr * 6 / FT(π) / ρw)^FT(1 / 3)
    ΔD = Dr - Deq
    phi_br =
        (Dr < Dr_th) ? FT(-1) : ((ΔD <= 0) ? kbr * ΔD : 2 * (exp(κbr * ΔD) - 1))
    dN_rai_dt_br = -(phi_br + 1) * dN_rai_dt_sc

    return dN_rai_dt_br
end

"""
    rain_self_collection_and_breakup(SB2006, q_rai, ρ, N_rai)

 - `SB2006` - a struct with SB2006 parameters for size distribution
    self collection and breakup
 - `q_rai` - rain water specific humidity
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns a named tupple containing the raindrops self-collection and breakup rates
for `scheme == SB2006Type`
"""
function rain_self_collection_and_breakup(
    (; pdf, self, brek)::CMP.SB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}

    sc = rain_self_collection(pdf, self, q_rai, ρ, N_rai)
    br = rain_breakup(pdf, brek, q_rai, ρ, N_rai, sc)

    return (; sc, br)
end

"""
    rain_terminal_velocity(SB2006, vel, q_rai, ρ, N_rai)

 - `SB2006` - a struct with SB2006 rain size distribution parameters
 - `vel` - a struct with terminal velocity parameters
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
    (; pdf)::CMP.SB2006{FT},
    (; ρ0, aR, bR, cR)::CMP.SB2006VelType{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}
    # TODO: Input argument list needs to be redesigned
    if q_rai < eps(FT)
        return (FT(0), FT(0))
    end
    λr = raindrops_limited_vars(pdf, q_rai, ρ, N_rai).λr
    vt0 = max(FT(0), sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)))
    vt1 = max(FT(0), sqrt(ρ0 / ρ) * (aR - bR / (1 + cR / λr)^FT(4)))
    return (vt0, vt1)
end
function rain_terminal_velocity(
    (; pdf)::CMP.SB2006{FT},
    vel::CMP.Chen2022VelTypeRain{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}
    if q_rai < eps(FT)
        return (FT(0), FT(0))
    end
    # coefficients from Table B1 from Chen et. al. 2022
    aiu, bi, ciu = CO.Chen2022_vel_coeffs(vel, ρ)
    # size distribution parameter
    λ = raindrops_limited_vars(pdf, q_rai, ρ, N_rai).λr

    # eq 20 from Chen et al 2022
    vt0 = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, λ, 0))
    vt3 = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, λ, 3))

    vt0 = max(FT(0), vt0)
    vt3 = max(FT(0), vt3)
    # It should be (ϕ^κ * vt0, ϕ^κ * vt3), but for rain drops ϕ = 1 and κ = 0
    return (vt0, vt3)
end

"""
    Returns the approximation of an incomplete gamma function
    for a ∈ {-1.0, -0.101}, and x in [0.067 1.82]
"""
function Γ_incl(a::FT, x::FT) where {FT}
    #return exp(-x) / ((FT(1.5) - FT(0.54) * a) * x^(FT(0.46) - FT(0.75) * a))
    return exp(-x) / (
        (FT(0.33) - FT(0.7) * a) * x^(FT(0.08) - FT(0.93) * a) +
        (FT(1.34) - FT(0.1) * a) * x^(FT(0.8) - a)
    )
end

"""
    rain_evaporation(evap, aps, tps, q, q_rai, ρ, N_rai, T)

 - `evap` - evaporation parameterization scheme
 - `aps` - air properties
 - `tps` - thermodynamics parameters
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
    (; pdf, evap)::CMP.SB2006{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
    T::FT,
) where {FT}

    evap_rate_0 = FT(0)
    evap_rate_1 = FT(0)
    S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && S < FT(0))

        (; ν_air, D_vapor) = aps
        (; av, bv, α, β, ρ0) = evap
        x_star = pdf.xr_min
        ρw = pdf.ρw
        G = CO.G_func(aps, tps, T, TD.Liquid())

        xr = raindrops_limited_vars(pdf, q_rai, ρ, N_rai).xr
        Dr = (FT(6) / FT(π) / ρw)^FT(1 / 3) * xr^FT(1 / 3)

        t_star = (FT(6) * x_star / xr)^FT(1 / 3)
        a_vent_0 = av * Γ_incl(FT(-1), t_star) / FT(6)^FT(-2 / 3)
        b_vent_0 =
            bv * Γ_incl(-FT(0.5) + FT(1.5) * β, t_star) /
            FT(6)^FT(β / 2 - FT(0.5))

        a_vent_1 = av * SF.gamma(FT(2)) / FT(6)^FT(1 / 3)
        b_vent_1 =
            bv * SF.gamma(FT(5 / 2) + FT(3 / 2) * β) / FT(6)^FT(β / 2 + 1 / 2)

        N_Re = α * xr^β * sqrt(ρ0 / ρ) * Dr / ν_air
        Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
        Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

        evap_rate_0 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv0 / xr)
        evap_rate_1 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)
    end

    return (evap_rate_0, evap_rate_1)
end

# Additional double moment autoconversion and accretion parametrizations:
# - Khairoutdinov and Kogan (2000)
# - Beheng (1994)
# - Tripoli and Cotton (1980)
# - Liu and Daum (2004)
# - variable timescale autoconversion Azimi (2023)

"""
    conv_q_liq_to_q_rai(acnv, q_liq, ρ; N_d, smooth_transition)

 - `acnv` - 2-moment rain autoconversion parameterization
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
    (; acnv)::CMP.KK2000{FT},
    q_liq,
    ρ;
    N_d = FT(1e8),
) where {FT}
    q_liq = max(0, q_liq)
    (; A, a, b, c) = acnv
    return A * q_liq^a * N_d^b * ρ^c
end
function conv_q_liq_to_q_rai(
    (; acnv)::CMP.B1994{FT},
    q_liq,
    ρ;
    N_d = FT(1e8),
    smooth_transition = false,
) where {FT}
    q_liq = max(0, q_liq)
    (; C, a, b, c, N_0, k, d_low, d_high) = acnv
    d = FT(0)
    if smooth_transition
        d_low_acnv_fraction = CO.logistic_function(N_d, N_0, k)
        d_high_acnv_fraction = FT(1) - d_low_acnv_fraction
        d = d_low_acnv_fraction * d_low + d_high_acnv_fraction * d_high
    else
        d = N_d >= N_0 ? d_low : d_high
    end
    return C * d^a * (q_liq * ρ)^b * N_d^c / ρ
end
function conv_q_liq_to_q_rai(
    (; acnv)::CMP.TC1980{FT},
    q_liq,
    ρ;
    N_d = FT(1e8),
    smooth_transition = false,
) where {FT}
    #TODO - The original paper is actually formulated for mixing ratios, not specific humidities
    q_liq = max(0, q_liq)
    (; m0_liq_coeff, me_liq, D, a, b, r_0, k) = acnv
    q_liq_threshold::FT = m0_liq_coeff * N_d / ρ * r_0^me_liq
    output =
        smooth_transition ? CO.logistic_function(q_liq, q_liq_threshold, k) :
        CO.heaviside(q_liq - q_liq_threshold)
    return D * q_liq^a * N_d^b * output
end
function conv_q_liq_to_q_rai(
    (; ρ_w, R_6C_0, E_0, k)::CMP.LD2004{FT},
    q_liq,
    ρ;
    N_d = FT(1e8),
    smooth_transition = false,
) where {FT}
    if q_liq <= eps(FT)
        return FT(0)
    else
        # Mean volume radius in microns (assuming spherical cloud droplets)
        r_vol =
            (FT(3) * (q_liq * ρ) / FT(4) / FT(π) / ρ_w / N_d)^FT(1 / 3) *
            FT(1e6)

        # Assumed size distribution: modified gamma distribution
        β_6 = ((r_vol + FT(3)) / r_vol)^FT(1 / 3)
        E = E_0 * β_6^6
        R_6 = β_6 * r_vol
        R_6C = R_6C_0 / (q_liq * ρ)^FT(1 / 6) / R_6^FT(1 / 2)

        output =
            smooth_transition ? CO.logistic_function(R_6, R_6C, k) :
            CO.heaviside(R_6 - R_6C)
        return E * (q_liq * ρ)^3 / N_d / ρ * output
    end
end
function conv_q_liq_to_q_rai(
    (; τ, α)::CMP.VarTimescaleAcnv{FT},
    q_liq::FT,
    ρ::FT;
    N_d::FT = FT(1e8),
) where {FT}
    return max(0, q_liq) / (τ * (N_d / 1e8)^α)
end

"""
    accretion(accretion_scheme, q_liq, q_rai, ρ)

 - `accretion_scheme` - type for 2-moment rain accretion parameterization
 - `q_liq` - cloud water specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density (for `KK2000Type` and `Beheng1994Type`)

 Returns the accretion rate of rain, parametrized following
 - Khairoutdinov and Kogan (2000) for `scheme == KK2000Type`
 - Beheng (1994) for `scheme == B1994Type`
 - Tripoli and Cotton (1980) for `scheme == TC1980Type`
"""
function accretion((; accr)::CMP.KK2000{FT}, q_liq, q_rai, ρ) where {FT}
    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)
    (; A, a, b) = accr
    return A * (q_liq * q_rai)^a * ρ^b
end

function accretion((; accr)::CMP.B1994{FT}, q_liq, q_rai, ρ) where {FT}
    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_liq * ρ * q_rai
end

function accretion((; accr)::CMP.TC1980{FT}, q_liq, q_rai) where {FT}
    #TODO - The original paper is actually formulated for mixing ratios, not specific humidities
    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_liq * q_rai
end

end
