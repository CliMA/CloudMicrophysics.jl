"""
    One-moment bulk microphysics scheme, which includes:

  - terminal velocity of precipitation
  - autoconversion of cloud liquid water into rain and of cloud ice into snow
  - accretion due to collisions between categories of condensed species
  - evaporation and sublimation of hydrometeors
  - melting of snow into rain
"""
module Microphysics1M

import SpecialFunctions
const SF = SpecialFunctions

import Thermodynamics
import CloudMicrophysics

const SF = SpecialFunctions
const TD = Thermodynamics
const CO = CloudMicrophysics.Common

import CloudMicrophysics.Microphysics_1M_Parameters

import ..CommonTypes
const CT = CommonTypes

import ..Common
const CO = Common

E(param_set::Microphysics_1M_Parameters, ::CT.LiquidType, ::CT.RainType) =
    param_set.E_liq_rai
E(param_set::Microphysics_1M_Parameters, ::CT.LiquidType, ::CT.SnowType) =
    param_set.E_liq_sno
E(param_set::Microphysics_1M_Parameters, ::CT.IceType, ::CT.RainType) =
    param_set.E_ice_rai
E(param_set::Microphysics_1M_Parameters, ::CT.IceType, ::CT.SnowType) =
    param_set.E_ice_sno
E(param_set::Microphysics_1M_Parameters, ::CT.RainType, ::CT.SnowType) =
    param_set.E_rai_sno
E(param_set::Microphysics_1M_Parameters, ::CT.SnowType, ::CT.RainType) =
    param_set.E_rai_sno

export terminal_velocity

export conv_q_liq_to_q_rai
export conv_q_ice_to_q_sno

export accretion
export accretion_rain_sink
export accretion_snow_rain

export evaporation_sublimation
export snow_melt

"""
    v0_rai(param_set, ρ)

 - `param_set` - abstract set with Earth parameters
 - `ρ` air density

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
function v0_rai(param_set::Microphysics_1M_Parameters, ρ::FT) where {FT <: Real}

    ρ_cloud_liq = param_set.ρ_cloud_liq
    C_drag = param_set.C_drag
    grav = param_set.grav
    r0_rai = param_set.r0_rai

    return sqrt(FT(8 / 3) / C_drag * (ρ_cloud_liq / ρ - FT(1)) * grav * r0_rai)
end

"""
    n0_sno(param_set, q_sno, ρ)

 - `param_set` - abstract set with Earth parameters
 - `q_sno` -  snow specific humidity
 - `ρ` - air density

Returns the intercept parameter of the assumed Marshall-Palmer distribution of
snow particles.
"""
function n0_sno(
    param_set::Microphysics_1M_Parameters,
    q_sno::FT,
    ρ::FT,
) where {FT <: Real}

    ν_sno = param_set.ν_sno
    μ_sno = param_set.μ_sno

    # TODO               this max should be replaced by
    #                    limiting inside a PhasePartition struct for
    #                    precipitation (once it is implemented)
    return μ_sno * (ρ * max(0, q_sno))^ν_sno
end

"""
    unpack_params(param_set, micro, ρ, q_)

 - `param_set` - abstract set with Earth parameters
 - `micro` - type for cloud ice, rain or snow
 - `q_` - specific humidity
 - `ρ` - air density

Utility function that unpacks microphysics parameters.
"""
function unpack_params(
    param_set::Microphysics_1M_Parameters,
    ice::CT.IceType,
    ρ::FT,
    q_ice::FT,
) where {FT <: Real}
    #TODO - make ρ and q_ice optional
    n0_ice = param_set.n0_ice
    r0_ice = param_set.r0_ice

    m0_ice = param_set.m0_ice
    me_ice = param_set.me_ice

    χm_ice = param_set.χm_ice
    Δm_ice = param_set.Δm_ice

    return (n0_ice, r0_ice, m0_ice, me_ice, χm_ice, Δm_ice)
end
function unpack_params(
    param_set::Microphysics_1M_Parameters,
    rain::CT.RainType,
    ρ::FT,
    q_rai::FT,
) where {FT <: Real}
    #TODO - make q_rai optional
    n0_rai = param_set.n0_rai
    r0_rai = param_set.r0_rai

    m0_rai = param_set.m0_rai
    me_rai = param_set.me_rai
    a0_rai = param_set.a0_rai
    ae_rai = param_set.ae_rai
    _v0_rai = v0_rai(param_set, ρ)
    ve_rai = param_set.ve_rai

    χm_rai = param_set.χm_rai
    Δm_rai = param_set.Δm_rai
    χa_rai = param_set.χa_rai
    Δa_rai = param_set.Δa_rai
    χv_rai = param_set.χv_rai
    Δv_rai = param_set.Δv_rai

    return (
        n0_rai,
        r0_rai,
        m0_rai,
        me_rai,
        χm_rai,
        Δm_rai,
        a0_rai,
        ae_rai,
        χa_rai,
        Δa_rai,
        _v0_rai,
        ve_rai,
        χv_rai,
        Δv_rai,
    )
end
function unpack_params(
    param_set::Microphysics_1M_Parameters,
    snow::CT.SnowType,
    ρ::FT,
    q_sno::FT,
) where {FT <: Real}

    _n0_sno = n0_sno(param_set, q_sno, ρ)
    r0_sno = param_set.r0_sno

    m0_sno = param_set.m0_sno
    me_sno = param_set.me_sno
    a0_sno = param_set.a0_sno
    ae_sno = param_set.ae_sno
    v0_sno = param_set.v0_sno
    ve_sno = param_set.ve_sno

    χm_sno = param_set.χm_sno
    Δm_sno = param_set.Δm_sno
    χa_sno = param_set.χa_sno
    Δa_sno = param_set.Δa_sno
    χv_sno = param_set.χv_sno
    Δv_sno = param_set.Δv_sno

    return (
        _n0_sno,
        r0_sno,
        m0_sno,
        me_sno,
        χm_sno,
        Δm_sno,
        a0_sno,
        ae_sno,
        χa_sno,
        Δa_sno,
        v0_sno,
        ve_sno,
        χv_sno,
        Δv_sno,
    )
end


"""
    lambda(q, ρ, n0, m0, me, r0, χm, Δm)

 - `q` - specific humidity of rain, ice or snow
 - `ρ` - air density
 - `n0` - size distribution parameter
 - `m0`, `me`, `χm`, `Δm`, `r0` - mass(radius) parameters

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(
    q::FT,
    ρ::FT,
    n0::FT,
    m0::FT,
    me::FT,
    r0::FT,
    χm::FT,
    Δm::FT,
) where {FT <: Real}

    λ::FT = FT(0)

    if q > FT(0)
        λ =
            (
                χm * m0 * n0 * SF.gamma(me + Δm + FT(1)) / ρ / q / r0^(me + Δm)
            )^FT(1 / (me + Δm + 1))
    end
    return λ
end

"""

    terminal_velocity(param_set, precip, ρ, q_)

 - `param_set` - abstract set with Earth parameters
 - `precip` - a type for rain or snow
 - `ρ` - air density
 - `q_` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of rain drops and snow crystals.
"""
function terminal_velocity(
    param_set::Microphysics_1M_Parameters,
    precip::CT.AbstractPrecipType,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        (n0, r0, m0, me, χm, Δm, a0, ae, χa, Δa, v0, ve, χv, Δv) =
            unpack_params(param_set, precip, ρ, q_)

        λ::FT = lambda(q_, ρ, n0, m0, me, r0, χm, Δm)

        fall_w =
            χv *
            v0 *
            (λ * r0)^(-ve - Δv) *
            SF.gamma(me + ve + Δm + Δv + FT(1)) / SF.gamma(me + Δm + FT(1))
    end

    return fall_w
end

"""

    conv_q_liq_to_q_rai(param_set, q_liq)

 - `param_set` - abstract set with Earth parameters
 - `q_liq` - liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Kessler (1995).
"""
function conv_q_liq_to_q_rai(
    param_set::Microphysics_1M_Parameters,
    q_liq::FT,
) where {FT <: Real}

    τ_acnv_rai = param_set.τ_acnv_rai
    q_liq_threshold = param_set.q_liq_threshold

    return max(0, q_liq - q_liq_threshold) / τ_acnv_rai
end

"""
    conv_q_ice_to_q_sno_no_supersat(param_set, q_ice)

 - `param_set` - abstract set with Earth parameters
 - `q_ice` -  cloud ice specific humidity

Returns the q_sno tendency due to autoconversion from ice.
This is a simplified version of a snow autoconversion rate that can be used in
simulations where there is no supersaturation
(for example in TC.jl when using saturation adjustment).
"""
function conv_q_ice_to_q_sno_no_supersat(
    param_set::Microphysics_1M_Parameters,
    q_ice::FT,
) where {FT <: Real}

    τ_acnv_sno = param_set.τ_acnv_sno
    q_ice_threshold = param_set.q_ice_threshold

    return max(0, q_ice - q_ice_threshold) / τ_acnv_sno
end

"""
    conv_q_ice_to_q_sno(param_set, q, ρ, T)

 - `param_set` - abstract set with Earth parameters
 - `q` - phase partition
 - `ρ` - air density
 - `T` - air temperature

Returns the q_sno tendency due to autoconversion from ice.
Parameterized following Harrington et al. (1996) and Kaul et al. (2015).
"""
function conv_q_ice_to_q_sno(
    param_set::Microphysics_1M_Parameters,
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT <: Real}
    acnv_rate = FT(0)
    S::FT = TD.supersaturation(param_set.TPS, q, ρ, T, TD.Ice())

    if (q.ice > FT(0) && S > FT(0))

        G::FT = CO.G_func(param_set, T, TD.Ice())

        r_ice_snow::FT = param_set.r_ice_snow

        (n0, r0, m0, me, χm, Δm) =
            unpack_params(param_set, CT.IceType(), ρ, q.ice)

        λ::FT = lambda(q.ice, ρ, n0, m0, me, r0, χm, Δm)

        acnv_rate =
            4 * FT(π) * S * G * n0 / ρ *
            exp(-λ * r_ice_snow) *
            (r_ice_snow^FT(2) / (me + Δm) + (r_ice_snow * λ + FT(1)) / λ^FT(2))
    end
    return acnv_rate
end

"""
    accretion(param_set, cloud, precip, q_clo, q_pre, ρ)

 - `param_set` - abstract set with Earth parameters
 - `cloud` - type for cloud water or cloud ice
 - `precip` - type for rain or snow
 - `q_clo` - cloud water or cloud ice specific humidity
 - `q_pre` - rain water or snow specific humidity
 - `ρ` - rain water or snow specific humidity

Returns the source of precipitating water (rain or snow)
due to collisions with cloud water (liquid or ice).
"""
function accretion(
    param_set::Microphysics_1M_Parameters,
    cloud::CT.AbstractCloudType,
    precip::CT.AbstractPrecipType,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        (n0, r0, m0, me, χm, Δm, a0, ae, χa, Δa, v0, ve, χv, Δv) =
            unpack_params(param_set, precip, ρ, q_pre)

        λ::FT = lambda(q_pre, ρ, n0, m0, me, r0, χm, Δm)
        _E::FT = E(param_set, cloud, precip)

        accr_rate =
            q_clo * _E * n0 * a0 * v0 * χa * χv / λ *
            SF.gamma(ae + ve + Δa + Δv + FT(1)) / (λ * r0)^(ae + ve + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_rain_sink(param_set, q_ice, q_rai, ρ)

 - `param_set` - abstract set with Earth parameters
 - `q_ice` - cloud ice specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

Returns the sink of rain water (partial source of snow) due to collisions
with cloud ice.
"""
function accretion_rain_sink(
    param_set::Microphysics_1M_Parameters,
    q_ice::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_ice > FT(0) && q_rai > FT(0))

        (n0_ice, r0_ice, m0_ice, me_ice, χm_ice, Δm_ice) =
            unpack_params(param_set, CT.IceType(), ρ, q_ice)

        (
            n0_rai,
            r0_rai,
            m0_rai,
            me_rai,
            χm_rai,
            Δm_rai,
            a0_rai,
            ae_rai,
            χa_rai,
            Δa_rai,
            v0_rai,
            ve_rai,
            χv_rai,
            Δv_rai,
        ) = unpack_params(param_set, CT.RainType(), ρ, q_rai)

        _E::FT = E(param_set, CT.IceType(), CT.RainType())

        λ_rai::FT =
            lambda(q_rai, ρ, n0_rai, m0_rai, me_rai, r0_rai, χm_rai, Δm_rai)
        λ_ice::FT =
            lambda(q_ice, ρ, n0_ice, m0_ice, me_ice, r0_ice, χm_ice, Δm_ice)

        accr_rate =
            _E / ρ *
            n0_rai *
            n0_ice *
            m0_rai *
            a0_rai *
            v0_rai *
            χm_rai *
            χa_rai *
            χv_rai / λ_ice / λ_rai * SF.gamma(
                me_rai + ae_rai + ve_rai + Δm_rai + Δa_rai + Δv_rai + FT(1),
            ) /
            (
                r0_rai * λ_rai
            )^(me_rai + ae_rai + ve_rai + Δm_rai + Δa_rai + Δv_rai)
    end
    return accr_rate
end

"""
    accretion_snow_rain(param_set, type_i, type_j, q_i, q_j, ρ)

 - `i` - snow for temperatures below freezing
         or rain for temperatures above freezing
 - `j` - rain for temperatures below freezing
         or rain for temperatures above freezing
 - `param_set` - abstract set with Earth parameters
 - `type_i`, `type_j` - a type for snow or rain
 - `q_` - specific humidity of snow or rain
 - `ρ` - air density

Returns the accretion rate between rain and snow.
Collisions between rain and snow result in
snow at temperatures below freezing and in rain at temperatures above freezing.
"""
function accretion_snow_rain(
    param_set::Microphysics_1M_Parameters,
    type_i::CT.AbstractPrecipType,
    type_j::CT.AbstractPrecipType,
    q_i::FT,
    q_j::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        (
            n0_i,
            r0_i,
            m0_i,
            me_i,
            χm_i,
            Δm_i,
            a0_i,
            ae_i,
            χa_i,
            Δa_i,
            v0_i,
            ve_i,
            χv_i,
            Δv_i,
        ) = unpack_params(param_set, type_i, ρ, q_i)
        (
            n0_j,
            r0_j,
            m0_j,
            me_j,
            χm_j,
            Δm_j,
            a0_j,
            ae_j,
            χa_j,
            Δa_j,
            v0_j,
            ve_j,
            χv_j,
            Δv_j,
        ) = unpack_params(param_set, type_j, ρ, q_j)

        E_ij::FT = E(param_set, type_i, type_j)

        λ_i::FT = lambda(q_i, ρ, n0_i, m0_i, me_i, r0_i, χm_i, Δm_i)
        λ_j::FT = lambda(q_j, ρ, n0_j, m0_j, me_j, r0_j, χm_j, Δm_j)

        v_ti = terminal_velocity(param_set, type_i, ρ, q_i)
        v_tj = terminal_velocity(param_set, type_j, ρ, q_j)

        accr_rate =
            FT(π) / ρ * n0_i * n0_j * m0_j * χm_j * E_ij * abs(v_ti - v_tj) /
            r0_j^(me_j + Δm_j) * (
                FT(2) * SF.gamma(me_j + Δm_j + FT(1)) / λ_i^FT(3) /
                λ_j^(me_j + Δm_j + FT(1)) +
                FT(2) * SF.gamma(me_j + Δm_j + FT(2)) / λ_i^FT(2) /
                λ_j^(me_j + Δm_j + FT(2)) +
                SF.gamma(me_j + Δm_j + FT(3)) / λ_i / λ_j^(me_j + Δm_j + FT(3))
            )
    end
    return accr_rate
end

"""
    evaporation_sublimation(param_set, rain, q, q_rai, ρ, T)
    evaporation_sublimation(param_set, snow, q, q_sno, ρ, T)

 - `param_set` - abstract set with Earth parameters
 - `rain` - a type for rain
 - `snow` - a type for snow
 - `q` - phase partition
 - `q_rai` - rain specific humidity
 - `q_sno` - snow specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to rain evaporation or snow sublimation.
"""
function evaporation_sublimation(
    param_set::Microphysics_1M_Parameters,
    rain::CT.RainType,
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    S::FT = TD.supersaturation(param_set.TPS, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && S < FT(0))

        a_vent = param_set.a_vent_rai
        b_vent = param_set.b_vent_rai
        ν_air = param_set.ν_air
        D_vapor = param_set.D_vapor

        G::FT = CO.G_func(param_set, T, TD.Liquid())

        (n0, r0, m0, me, χm, Δm, a0, ae, χa, Δa, v0, ve, χv, Δv) =
            unpack_params(param_set, rain, ρ, q_rai)

        λ::FT = lambda(q_rai, ρ, n0, m0, me, r0, χm, Δm)

        evap_subl_rate =
            4 * FT(π) * n0 / ρ * S * G / λ^FT(2) * (
                a_vent +
                b_vent * (ν_air / D_vapor)^FT(1 / 3) /
                (r0 * λ)^((ve + Δv) / FT(2)) *
                (FT(2) * v0 * χv / ν_air / λ)^FT(1 / 2) *
                SF.gamma((ve + Δv + FT(5)) / FT(2))
            )
    end
    # only evaporation is considered for rain
    return min(0, evap_subl_rate)
end
function evaporation_sublimation(
    param_set::Microphysics_1M_Parameters,
    snow::CT.SnowType,
    q::TD.PhasePartition{FT},
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    if q_sno > FT(0)

        a_vent = param_set.a_vent_sno
        b_vent = param_set.b_vent_sno
        ν_air = param_set.ν_air
        D_vapor = param_set.D_vapor

        S::FT = TD.supersaturation(param_set.TPS, q, ρ, T, TD.Ice())
        G::FT = CO.G_func(param_set, T, TD.Ice())

        (n0, r0, m0, me, χm, Δm, a0, ae, χa, Δa, v0, ve, χv, Δv) =
            unpack_params(param_set, snow, ρ, q_sno)
        λ::FT = lambda(q_sno, ρ, n0, m0, me, r0, χm, Δm)

        evap_subl_rate =
            4 * FT(π) * n0 / ρ * S * G / λ^FT(2) * (
                a_vent +
                b_vent * (ν_air / D_vapor)^FT(1 / 3) /
                (r0 * λ)^((ve + Δv) / FT(2)) *
                (FT(2) * v0 * χv / ν_air / λ)^FT(1 / 2) *
                SF.gamma((ve + Δv + FT(5)) / FT(2))
            )
    end
    return evap_subl_rate
end

"""
    snow_melt(param_set, q_sno, ρ, T)

 - `param_set` - abstract set with Earth parameters
 - `q_sno` - snow water specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to snow melt.
"""
function snow_melt(
    param_set::Microphysics_1M_Parameters,
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}

    snow_melt_rate = FT(0)
    T_freeze = param_set.T_freeze

    if (q_sno > FT(0) && T > T_freeze)

        a_vent::FT = param_set.a_vent_sno
        b_vent::FT = param_set.b_vent_sno
        ν_air::FT = param_set.ν_air
        D_vapor::FT = param_set.D_vapor
        K_therm::FT = param_set.K_therm

        L = TD.latent_heat_fusion(param_set.TPS, T)

        (n0, r0, m0, me, χm, Δm, a0, ae, χa, Δa, v0, ve, χv, Δv) =
            unpack_params(param_set, CT.SnowType(), ρ, q_sno)
        λ::FT = lambda(q_sno, ρ, n0, m0, me, r0, χm, Δm)

        snow_melt_rate =
            4 * FT(π) * n0 / ρ * K_therm / L * (T - T_freeze) / λ^FT(2) * (
                a_vent +
                b_vent * (ν_air / D_vapor)^FT(1 / 3) /
                (r0 * λ)^((ve + Δv) / FT(2)) *
                (FT(2) * v0 * χv / ν_air / λ)^FT(1 / 2) *
                SF.gamma((ve + Δv + FT(5)) / FT(2))
            )
    end
    return snow_melt_rate
end

end #module Microphysics1M.jl
