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
const TD = Thermodynamics

import CLIMAParameters
const CP = CLIMAParameters
const APS = CP.AbstractParameterSet

import ..CommonTypes
const CT = CommonTypes

import ..Common
const CO = Common

import ..InternalClimaParams
const ICP = InternalClimaParams

E(param_set::APS, ::CT.LiquidType, ::CT.RainType) =
    ICP.E_liq_rai(param_set::APS)
E(param_set::APS, ::CT.LiquidType, ::CT.SnowType) =
    ICP.E_liq_sno(param_set::APS)
E(param_set::APS, ::CT.IceType, ::CT.RainType) = ICP.E_ice_rai(param_set::APS)
E(param_set::APS, ::CT.IceType, ::CT.SnowType) = ICP.E_ice_sno(param_set::APS)
E(param_set::APS, ::CT.RainType, ::CT.SnowType) = ICP.E_rai_sno(param_set::APS)
E(param_set::APS, ::CT.SnowType, ::CT.RainType) = ICP.E_rai_sno(param_set::APS)

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
function v0_rai(param_set::APS, ρ::FT) where {FT <: Real}

    _ρ_cloud_liq::FT = ICP.ρ_cloud_liq(param_set)
    _C_drag::FT = ICP.C_drag(param_set)
    _grav::FT = ICP.grav(param_set)
    _r0_rai::FT = ICP.r0_rai(param_set)

    return sqrt(
        FT(8 / 3) / _C_drag * (_ρ_cloud_liq / ρ - FT(1)) * _grav * _r0_rai,
    )
end

"""
    n0_sno(param_set, q_sno, ρ)

 - `param_set` - abstract set with Earth parameters
 - `q_sno` -  snow specific humidity
 - `ρ` - air density

Returns the intercept parameter of the assumed Marshall-Palmer distribution of
snow particles.
"""
function n0_sno(param_set::APS, q_sno::FT, ρ::FT) where {FT <: Real}

    _ν_sno::FT = ICP.ν_sno(param_set)
    _μ_sno::FT = ICP.μ_sno(param_set)

    # TODO               this max should be replaced by
    #                    limiting inside a PhasePartition struct for
    #                    precipitation (once it is implemented)
    return _μ_sno * (ρ * max(0, q_sno))^_ν_sno
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
    param_set::APS,
    ice::CT.IceType,
    ρ::FT,
    q_ice::FT,
) where {FT <: Real}
    #TODO - make ρ and q_ice optional
    _n0_ice::FT = ICP.n0_ice(param_set)
    _r0_ice::FT = ICP.r0_ice(param_set)

    _m0_ice::FT = ICP.m0_ice(param_set)
    _me_ice::FT = ICP.me_ice(param_set)

    _χm_ice::FT = ICP.χm_ice(param_set)
    _Δm_ice::FT = ICP.Δm_ice(param_set)

    return (_n0_ice, _r0_ice, _m0_ice, _me_ice, _χm_ice, _Δm_ice)
end
function unpack_params(
    param_set::APS,
    rain::CT.RainType,
    ρ::FT,
    q_rai::FT,
) where {FT <: Real}
    #TODO - make q_rai optional
    _n0_rai::FT = ICP.n0_rai(param_set)
    _r0_rai::FT = ICP.r0_rai(param_set)

    _m0_rai::FT = ICP.m0_rai(param_set)
    _me_rai::FT = ICP.me_rai(param_set)
    _a0_rai::FT = ICP.a0_rai(param_set)
    _ae_rai::FT = ICP.ae_rai(param_set)
    _v0_rai::FT = v0_rai(param_set, ρ)
    _ve_rai::FT = ICP.ve_rai(param_set)

    _χm_rai::FT = ICP.χm_rai(param_set)
    _Δm_rai::FT = ICP.Δm_rai(param_set)
    _χa_rai::FT = ICP.χa_rai(param_set)
    _Δa_rai::FT = ICP.Δa_rai(param_set)
    _χv_rai::FT = ICP.χv_rai(param_set)
    _Δv_rai::FT = ICP.Δv_rai(param_set)

    return (
        _n0_rai,
        _r0_rai,
        _m0_rai,
        _me_rai,
        _χm_rai,
        _Δm_rai,
        _a0_rai,
        _ae_rai,
        _χa_rai,
        _Δa_rai,
        _v0_rai,
        _ve_rai,
        _χv_rai,
        _Δv_rai,
    )
end
function unpack_params(
    param_set::APS,
    snow::CT.SnowType,
    ρ::FT,
    q_sno::FT,
) where {FT <: Real}

    _n0_sno::FT = n0_sno(param_set, q_sno, ρ)
    _r0_sno::FT = ICP.r0_sno(param_set)

    _m0_sno::FT = ICP.m0_sno(param_set)
    _me_sno::FT = ICP.me_sno(param_set)
    _a0_sno::FT = ICP.a0_sno(param_set)
    _ae_sno::FT = ICP.ae_sno(param_set)
    _v0_sno::FT = ICP.v0_sno(param_set)
    _ve_sno::FT = ICP.ve_sno(param_set)

    _χm_sno::FT = ICP.χm_sno(param_set)
    _Δm_sno::FT = ICP.Δm_sno(param_set)
    _χa_sno::FT = ICP.χa_sno(param_set)
    _Δa_sno::FT = ICP.Δa_sno(param_set)
    _χv_sno::FT = ICP.χv_sno(param_set)
    _Δv_sno::FT = ICP.Δv_sno(param_set)

    return (
        _n0_sno,
        _r0_sno,
        _m0_sno,
        _me_sno,
        _χm_sno,
        _Δm_sno,
        _a0_sno,
        _ae_sno,
        _χa_sno,
        _Δa_sno,
        _v0_sno,
        _ve_sno,
        _χv_sno,
        _Δv_sno,
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
    param_set::APS,
    precip::CT.AbstractPrecipType,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, precip, ρ, q_)

        _λ::FT = lambda(q_, ρ, _n0, _m0, _me, _r0, _χm, _Δm)

        fall_w =
            _χv *
            _v0 *
            (_λ * _r0)^(-_ve - _Δv) *
            SF.gamma(_me + _ve + _Δm + _Δv + FT(1)) /
            SF.gamma(_me + _Δm + FT(1))
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
function conv_q_liq_to_q_rai(param_set::APS, q_liq::FT) where {FT <: Real}

    _τ_acnv_rai::FT = ICP.τ_acnv_rai(param_set)
    _q_liq_threshold::FT = ICP.q_liq_threshold(param_set)

    return max(0, q_liq - _q_liq_threshold) / _τ_acnv_rai
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
    param_set::APS,
    q_ice::FT,
) where {FT <: Real}

    _τ_acnv_sno::FT = ICP.τ_acnv_sno(param_set)
    _q_ice_threshold::FT = ICP.q_ice_threshold(param_set)

    return max(0, q_ice - _q_ice_threshold) / _τ_acnv_sno
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
    param_set::APS,
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT <: Real}
    acnv_rate = FT(0)
    _S::FT = TD.supersaturation(param_set, q, ρ, T, TD.Ice())

    if (q.ice > FT(0) && _S > FT(0))

        _G::FT = CO.G_func(param_set, T, TD.Ice())

        _r_ice_snow::FT = ICP.r_ice_snow(param_set)

        (_n0, _r0, _m0, _me, _χm, _Δm) =
            unpack_params(param_set, CT.IceType(), ρ, q.ice)

        _λ::FT = lambda(q.ice, ρ, _n0, _m0, _me, _r0, _χm, _Δm)

        acnv_rate =
            4 * FT(π) * _S * _G * _n0 / ρ *
            exp(-_λ * _r_ice_snow) *
            (
                _r_ice_snow^FT(2) / (_me + _Δm) +
                (_r_ice_snow * _λ + FT(1)) / _λ^FT(2)
            )
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
    param_set::APS,
    cloud::CT.AbstractCloudType,
    precip::CT.AbstractPrecipType,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, precip, ρ, q_pre)

        _λ::FT = lambda(q_pre, ρ, _n0, _m0, _me, _r0, _χm, _Δm)
        _E::FT = E(param_set, cloud, precip)

        accr_rate =
            q_clo * _E * _n0 * _a0 * _v0 * _χa * _χv / _λ *
            SF.gamma(_ae + _ve + _Δa + _Δv + FT(1)) /
            (_λ * _r0)^(_ae + _ve + _Δa + _Δv)
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
    param_set::APS,
    q_ice::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_ice > FT(0) && q_rai > FT(0))

        (_n0_ice, _r0_ice, _m0_ice, _me_ice, _χm_ice, _Δm_ice) =
            unpack_params(param_set, CT.IceType(), ρ, q_ice)

        (
            _n0_rai,
            _r0_rai,
            _m0_rai,
            _me_rai,
            _χm_rai,
            _Δm_rai,
            _a0_rai,
            _ae_rai,
            _χa_rai,
            _Δa_rai,
            _v0_rai,
            _ve_rai,
            _χv_rai,
            _Δv_rai,
        ) = unpack_params(param_set, CT.RainType(), ρ, q_rai)

        _E::FT = E(param_set, CT.IceType(), CT.RainType())

        _λ_rai::FT = lambda(
            q_rai,
            ρ,
            _n0_rai,
            _m0_rai,
            _me_rai,
            _r0_rai,
            _χm_rai,
            _Δm_rai,
        )
        _λ_ice::FT = lambda(
            q_ice,
            ρ,
            _n0_ice,
            _m0_ice,
            _me_ice,
            _r0_ice,
            _χm_ice,
            _Δm_ice,
        )

        accr_rate =
            _E / ρ *
            _n0_rai *
            _n0_ice *
            _m0_rai *
            _a0_rai *
            _v0_rai *
            _χm_rai *
            _χa_rai *
            _χv_rai / _λ_ice / _λ_rai * SF.gamma(
                _me_rai +
                _ae_rai +
                _ve_rai +
                _Δm_rai +
                _Δa_rai +
                _Δv_rai +
                FT(1),
            ) /
            (
                _r0_rai * _λ_rai
            )^(_me_rai + _ae_rai + _ve_rai + _Δm_rai + _Δa_rai + _Δv_rai)
    end
    return accr_rate
end

"""
    accretion_snow_rain(param_set, type_i, type_j, q_i, q_j, ρ)

 - `i` - snow for temperatures below freezing
         or rain for temperatures above freezing
 - `j` - rain for temperatures below freezing
         or snow for temperatures above freezing
 - `param_set` - abstract set with Earth parameters
 - `type_i`, `type_j` - a type for snow or rain
 - `q_` - specific humidity of snow or rain
 - `ρ` - air density

Returns the accretion rate between rain and snow.
Collisions between rain and snow result in
snow at temperatures below freezing and in rain at temperatures above freezing.
"""
function accretion_snow_rain(
    param_set::APS,
    type_i::CT.AbstractPrecipType,
    type_j::CT.AbstractPrecipType,
    q_i::FT,
    q_j::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        (
            _n0_i,
            _r0_i,
            _m0_i,
            _me_i,
            _χm_i,
            _Δm_i,
            _a0_i,
            _ae_i,
            _χa_i,
            _Δa_i,
            _v0_i,
            _ve_i,
            _χv_i,
            _Δv_i,
        ) = unpack_params(param_set, type_i, ρ, q_i)
        (
            _n0_j,
            _r0_j,
            _m0_j,
            _me_j,
            _χm_j,
            _Δm_j,
            _a0_j,
            _ae_j,
            _χa_j,
            _Δa_j,
            _v0_j,
            _ve_j,
            _χv_j,
            _Δv_j,
        ) = unpack_params(param_set, type_j, ρ, q_j)

        _E_ij::FT = E(param_set, type_i, type_j)

        _λ_i::FT = lambda(q_i, ρ, _n0_i, _m0_i, _me_i, _r0_i, _χm_i, _Δm_i)
        _λ_j::FT = lambda(q_j, ρ, _n0_j, _m0_j, _me_j, _r0_j, _χm_j, _Δm_j)

        _v_ti = terminal_velocity(param_set, type_i, ρ, q_i)
        _v_tj = terminal_velocity(param_set, type_j, ρ, q_j)

        accr_rate =
            FT(π) / ρ *
            _n0_i *
            _n0_j *
            _m0_j *
            _χm_j *
            _E_ij *
            abs(_v_ti - _v_tj) / _r0_j^(_me_j + _Δm_j) * (
                FT(2) * SF.gamma(_me_j + _Δm_j + FT(1)) / _λ_i^FT(3) /
                _λ_j^(_me_j + _Δm_j + FT(1)) +
                FT(2) * SF.gamma(_me_j + _Δm_j + FT(2)) / _λ_i^FT(2) /
                _λ_j^(_me_j + _Δm_j + FT(2)) +
                SF.gamma(_me_j + _Δm_j + FT(3)) / _λ_i /
                _λ_j^(_me_j + _Δm_j + FT(3))
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
    param_set::APS,
    rain::CT.RainType,
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    _S::FT = TD.supersaturation(param_set, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && _S < FT(0))

        _a_vent::FT = ICP.a_vent_rai(param_set)
        _b_vent::FT = ICP.b_vent_rai(param_set)
        _ν_air::FT = ICP.ν_air(param_set)
        _D_vapor::FT = ICP.D_vapor(param_set)

        _G::FT = CO.G_func(param_set, T, TD.Liquid())

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, rain, ρ, q_rai)

        _λ::FT = lambda(q_rai, ρ, _n0, _m0, _me, _r0, _χm, _Δm)

        evap_subl_rate =
            4 * FT(π) * _n0 / ρ * _S * _G / _λ^FT(2) * (
                _a_vent +
                _b_vent * (_ν_air / _D_vapor)^FT(1 / 3) /
                (_r0 * _λ)^((_ve + _Δv) / FT(2)) *
                (FT(2) * _v0 * _χv / _ν_air / _λ)^FT(1 / 2) *
                SF.gamma((_ve + _Δv + FT(5)) / FT(2))
            )
    end
    # only evaporation is considered for rain
    return min(0, evap_subl_rate)
end
function evaporation_sublimation(
    param_set::APS,
    snow::CT.SnowType,
    q::TD.PhasePartition{FT},
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    if q_sno > FT(0)

        _a_vent::FT = ICP.a_vent_sno(param_set)
        _b_vent::FT = ICP.b_vent_sno(param_set)
        _ν_air::FT = ICP.ν_air(param_set)
        _D_vapor::FT = ICP.D_vapor(param_set)

        _S::FT = TD.supersaturation(param_set, q, ρ, T, TD.Ice())
        _G::FT = CO.G_func(param_set, T, TD.Ice())

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, snow, ρ, q_sno)
        _λ::FT = lambda(q_sno, ρ, _n0, _m0, _me, _r0, _χm, _Δm)

        evap_subl_rate =
            4 * FT(π) * _n0 / ρ * _S * _G / _λ^FT(2) * (
                _a_vent +
                _b_vent * (_ν_air / _D_vapor)^FT(1 / 3) /
                (_r0 * _λ)^((_ve + _Δv) / FT(2)) *
                (FT(2) * _v0 * _χv / _ν_air / _λ)^FT(1 / 2) *
                SF.gamma((_ve + _Δv + FT(5)) / FT(2))
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
function snow_melt(param_set::APS, q_sno::FT, ρ::FT, T::FT) where {FT <: Real}

    snow_melt_rate = FT(0)
    _T_freeze::FT = ICP.T_freeze(param_set)

    if (q_sno > FT(0) && T > _T_freeze)

        _a_vent::FT = ICP.a_vent_sno(param_set)
        _b_vent::FT = ICP.b_vent_sno(param_set)
        _ν_air::FT = ICP.ν_air(param_set)
        _D_vapor::FT = ICP.D_vapor(param_set)
        _K_therm::FT = ICP.K_therm(param_set)

        L = TD.latent_heat_fusion(param_set, T)

        (_n0, _r0, _m0, _me, _χm, _Δm, _a0, _ae, _χa, _Δa, _v0, _ve, _χv, _Δv) =
            unpack_params(param_set, CT.SnowType(), ρ, q_sno)
        _λ::FT = lambda(q_sno, ρ, _n0, _m0, _me, _r0, _χm, _Δm)

        snow_melt_rate =
            4 * FT(π) * _n0 / ρ * _K_therm / L * (T - _T_freeze) / _λ^FT(2) * (
                _a_vent +
                _b_vent * (_ν_air / _D_vapor)^FT(1 / 3) /
                (_r0 * _λ)^((_ve + _Δv) / FT(2)) *
                (FT(2) * _v0 * _χv / _ν_air / _λ)^FT(1 / 2) *
                SF.gamma((_ve + _Δv + FT(5)) / FT(2))
            )
    end
    return snow_melt_rate
end

end #module Microphysics1M.jl
