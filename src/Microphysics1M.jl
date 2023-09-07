"""
    One-moment bulk microphysics scheme, which includes:

  - terminal velocity of precipitation
  - autoconversion of cloud liquid water into rain and of cloud ice into snow
  - accretion due to collisions between categories of condensed species
  - evaporation and sublimation of hydrometeors
  - melting of snow into rain
"""
module Microphysics1M

import SpecialFunctions as SF

import Thermodynamics as TD

import ..CommonTypes as CT
import ..Common as CO
import ..Parameters as CMP

const APS = CMP.AbstractCloudMicrophysicsParameters

export terminal_velocity

export conv_q_liq_to_q_rai
export conv_q_ice_to_q_sno

export accretion
export accretion_rain_sink
export accretion_snow_rain

export evaporation_sublimation
export snow_melt

"""
    E(prs, type_1, type_2)

 - `prs` - set with model parameters
 - `type_1` and `type_2` - types of colliding species

Returns collision efficiency for two colliding species
"""
E(prs::APS, ::CT.LiquidType, ::CT.RainType) = CMP.E_liq_rai(prs::APS)
E(prs::APS, ::CT.LiquidType, ::CT.SnowType) = CMP.E_liq_sno(prs::APS)
E(prs::APS, ::CT.IceType, ::CT.RainType) = CMP.E_ice_rai(prs::APS)
E(prs::APS, ::CT.IceType, ::CT.SnowType) = CMP.E_ice_sno(prs::APS)
E(prs::APS, ::CT.RainType, ::CT.SnowType) = CMP.E_rai_sno(prs::APS)
E(prs::APS, ::CT.SnowType, ::CT.RainType) = CMP.E_rai_sno(prs::APS)

"""
    n0(prs, q_sno, ρ, snow_type)

 - `prs` - abstract set with Earth parameters
 - `q_sno` -  snow specific humidity
 - `ρ` - air density
 - `type` - type for dispatch

Returns the intercept parameter of the assumed Marshall-Palmer distribution
"""
function n0(prs::APS, q_sno::FT, ρ::FT, ::CT.SnowType) where {FT <: Real}

    _ν_sno::FT = CMP.ν_sno(prs)
    _μ_sno::FT = CMP.μ_sno(prs)

    # TODO               this max should be replaced by
    #                    limiting inside a PhasePartition struct for
    #                    precipitation (once it is implemented)
    return _μ_sno * (ρ * max(0, q_sno))^_ν_sno
end
n0(prs::APS, ::Any, ::Any, ::CT.IceType) = CMP.n0_ice(prs)
n0(prs::APS, ::Any, ::Any, ::CT.RainType) = CMP.n0_rai(prs)

"""
    v0(prs, ρ, rain_type)

 - `prs` - abstract set with Earth parameters
 - `ρ` - air density
 - `type` - type for dispatch

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
function v0(prs::APS, ρ::FT, ::CT.RainType) where {FT <: Real}

    _ρ_cloud_liq::FT = CMP.ρ_cloud_liq(prs)
    _C_drag::FT = CMP.C_drag(prs)
    _grav::FT = CMP.grav(prs)
    _r0_rai::FT = CMP.r0_rai(prs)

    return sqrt(
        FT(8 / 3) / _C_drag * (_ρ_cloud_liq / ρ - FT(1)) * _grav * _r0_rai,
    )
end
v0(prs::APS, ::Any, ::CT.SnowType) = CMP.v0_sno(prs)

# Other ice/rain/snow parameters to dispatch over
a_vent(prs::APS, ::CT.RainType) = CMP.a_vent_rai(prs)
b_vent(prs::APS, ::CT.RainType) = CMP.b_vent_rai(prs)
a_vent(prs::APS, ::CT.SnowType) = CMP.a_vent_sno(prs)
b_vent(prs::APS, ::CT.SnowType) = CMP.b_vent_sno(prs)

r0(prs::APS, ::CT.IceType) = CMP.r0_ice(prs)
m0(prs::APS, ::CT.IceType) = CMP.m0_ice(prs)
me(prs::APS, ::CT.IceType) = CMP.me_ice(prs)
χm(prs::APS, ::CT.IceType) = CMP.χm_ice(prs)
Δm(prs::APS, ::CT.IceType) = CMP.Δm_ice(prs)

r0(prs::APS, ::CT.RainType) = CMP.r0_rai(prs)
m0(prs::APS, ::CT.RainType) = CMP.m0_rai(prs)
me(prs::APS, ::CT.RainType) = CMP.me_rai(prs)
a0(prs::APS, ::CT.RainType) = CMP.a0_rai(prs)
ae(prs::APS, ::CT.RainType) = CMP.ae_rai(prs)
ve(prs::APS, ::CT.RainType) = CMP.ve_rai(prs)
χm(prs::APS, ::CT.RainType) = CMP.χm_rai(prs)
Δm(prs::APS, ::CT.RainType) = CMP.Δm_rai(prs)
χa(prs::APS, ::CT.RainType) = CMP.χa_rai(prs)
Δa(prs::APS, ::CT.RainType) = CMP.Δa_rai(prs)
χv(prs::APS, ::CT.RainType) = CMP.χv_rai(prs)
Δv(prs::APS, ::CT.RainType) = CMP.Δv_rai(prs)

r0(prs::APS, ::CT.SnowType) = CMP.r0_sno(prs)
m0(prs::APS, ::CT.SnowType) = CMP.m0_sno(prs)
me(prs::APS, ::CT.SnowType) = CMP.me_sno(prs)
a0(prs::APS, ::CT.SnowType) = CMP.a0_sno(prs)
ae(prs::APS, ::CT.SnowType) = CMP.ae_sno(prs)
ve(prs::APS, ::CT.SnowType) = CMP.ve_sno(prs)
χm(prs::APS, ::CT.SnowType) = CMP.χm_sno(prs)
Δm(prs::APS, ::CT.SnowType) = CMP.Δm_sno(prs)
χa(prs::APS, ::CT.SnowType) = CMP.χa_sno(prs)
Δa(prs::APS, ::CT.SnowType) = CMP.Δa_sno(prs)
χv(prs::APS, ::CT.SnowType) = CMP.χv_sno(prs)
Δv(prs::APS, ::CT.SnowType) = CMP.Δv_sno(prs)

"""
    lambda(q, ρ, n0, m0, me, r0, χm, Δm)

 - `prs` - set with free parameters
 - `precip` - a type for cloud ice, rain or snow
 - `q` - specific humidity of rain, ice or snow
 - `ρ` - air density

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(
    prs::APS,
    precip::Union{CT.IceType, CT.RainType, CT.SnowType},
    q::FT,
    ρ::FT,
) where {FT <: Real}

    _n0::FT = n0(prs, q, ρ, precip)
    _r0::FT = r0(prs, precip)
    _m0::FT = m0(prs, precip)
    _me::FT = me(prs, precip)
    _Δm::FT = Δm(prs, precip)
    _χm::FT = χm(prs, precip)

    λ::FT = FT(0)

    if q > FT(0)
        λ =
            (
                _χm * _m0 * _n0 * SF.gamma(_me + _Δm + FT(1)) / ρ / q /
                _r0^(_me + _Δm)
            )^FT(1 / (_me + _Δm + 1))
    end
    return λ
end

"""
    terminal_velocity(prs, precip, velo_scheme, ρ, q_)

 - `prs` - abstract set with Earth parameters
 - `precip` - a type for ice, rain or snow
 - `velo_scheme` - type for terminal velocity parameterization
 - `ρ` - air density
 - `q_` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of particles.
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical power-law relations for `velo_scheme == Blk1MVelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171, for `velo_scheme == Chen2022Type`
"""
function terminal_velocity(
    prs::APS,
    precip::CT.AbstractPrecipType,
    velo_scheme::CT.Blk1MVelType,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        _r0::FT = r0(prs, precip)
        _me::FT = me(prs, precip)
        _Δm::FT = Δm(prs, precip)
        _χm::FT = χm(prs, precip)
        _χv::FT = χv(prs, precip)
        _v0::FT = v0(prs, ρ, precip)
        _ve::FT = ve(prs, precip)
        _Δv::FT = Δv(prs, precip)
        _λ::FT = lambda(prs, precip, q_, ρ)

        fall_w =
            _χv *
            _v0 *
            (_λ * _r0)^(-_ve - _Δv) *
            SF.gamma(_me + _ve + _Δm + _Δv + FT(1)) /
            SF.gamma(_me + _Δm + FT(1))
    end

    return fall_w
end
function terminal_velocity(
    prs::APS,
    precip::CT.RainType,
    velo_scheme::CT.Chen2022Type,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(prs, precip, ρ)
        # size distribution parameter
        _λ::FT = lambda(prs, precip, q_, ρ)

        # eq 20 from Chen et al 2022
        fall_w = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, _λ, 3))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end
function terminal_velocity(
    prs::APS,
    precip::CT.IceType,
    velo_scheme::CT.Chen2022Type,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        ρ_i::FT = CMP.ρ_cloud_ice(prs)
        _λ::FT = lambda(prs, precip, q_, ρ)

        # coefficients from Appendix B from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(prs, precip, ρ)

        # eq 20 from Chen et al 2022
        fall_w = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, _λ, 3))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end
function terminal_velocity(
    prs::APS,
    precip::CT.SnowType,
    velo_scheme::CT.Chen2022Type,
    ρ::FT,
    q_::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q_ > FT(0)

        _r0::FT = r0(prs, precip)
        _λ::FT = lambda(prs, precip, q_, ρ)

        m0c::FT = m0(prs, precip) * χm(prs, precip)
        a0c::FT = a0(prs, precip) * χa(prs, precip)
        mec::FT = me(prs, precip) + Δm(prs, precip)
        aec::FT = ae(prs, precip) + Δa(prs, precip)

        ρ_i::FT = CMP.ρ_cloud_ice(prs)

        # coefficients from Appendix B from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(prs, precip, ρ)

        κ = FT(-1 / 3) #oblate
        k = 3 # mass weighted

        tmp =
            _λ^(k + 1) *
            ((16 * a0c^3 * ρ_i^2) / (9 * π * m0c^2 * _r0^(3 * aec - 2 * mec)))^κ
        ci_pow =
            (2 .* ciu .+ _λ) .^
            (.-(3 .* aec .* κ .- 2 .* mec .* κ .+ bi .+ k .+ 1))

        ti = tmp .* aiu .* FT(2) .^ bi .* ci_pow

        Chen2022_vel_add_sno(t, b, aec, mec, κ, k) =
            t * SF.gamma(3 * κ * aec - 2 * κ * mec + b + k + 1) /
            SF.gamma(k + 1)

        fall_w = sum(Chen2022_vel_add_sno.(ti, bi, aec, mec, κ, k))
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end

"""
    conv_q_liq_to_q_rai(prs, q_liq)

 - `prs` - abstract set with Earth parameters
 - `q_liq` - liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Kessler (1995).
"""
function conv_q_liq_to_q_rai(
    prs::APS,
    q_liq::FT;
    smooth_transition::Bool = false,
) where {FT <: Real}

    _τ_acnv_rai::FT = CMP.τ_acnv_rai(prs)
    _q_liq_threshold::FT = CMP.q_liq_threshold(prs)

    _output::FT = FT(0)
    if smooth_transition
        _k::FT = CMP.k_thrshld_stpnss(prs)
        _output = CO.logistic_function_integral(q_liq, _q_liq_threshold, _k)
    else
        _output = max(0, q_liq - _q_liq_threshold)
    end
    return _output / _τ_acnv_rai

end

"""
    conv_q_ice_to_q_sno_no_supersat(prs, q_ice)

 - `prs` - abstract set with Earth parameters
 - `q_ice` -  cloud ice specific humidity

Returns the q_sno tendency due to autoconversion from ice.
This is a simplified version of a snow autoconversion rate that can be used in
simulations where there is no supersaturation
(for example in TC.jl when using saturation adjustment).
"""
function conv_q_ice_to_q_sno_no_supersat(
    prs::APS,
    q_ice::FT;
    smooth_transition::Bool = false,
) where {FT <: Real}

    _τ_acnv_sno::FT = CMP.τ_acnv_sno(prs)
    _q_ice_threshold::FT = CMP.q_ice_threshold(prs)

    _output::FT = FT(0)
    if smooth_transition
        _k::FT = CMP.k_thrshld_stpnss(prs)
        _output = CO.logistic_function_integral(q_ice, _q_ice_threshold, _k)
    else
        _output = max(0, q_ice - _q_ice_threshold)
    end
    return _output / _τ_acnv_sno

end

"""
    conv_q_ice_to_q_sno(prs, q, ρ, T)

 - `prs` - abstract set with Earth parameters
 - `q` - phase partition
 - `ρ` - air density
 - `T` - air temperature

Returns the q_sno tendency due to autoconversion from ice.
Parameterized following Harrington et al. (1996) and Kaul et al. (2015).
"""
function conv_q_ice_to_q_sno(
    prs::APS,
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT <: Real}
    acnv_rate = FT(0)
    thermo_params = CMP.thermodynamics_params(prs)
    _S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())

    if (q.ice > FT(0) && _S > FT(0))

        _G::FT = CO.G_func(prs, T, TD.Ice())

        _r_ice_snow::FT = CMP.r_ice_snow(prs)
        _n0::FT = n0(prs, FT(0), ρ, CT.IceType())
        _me::FT = me(prs, CT.IceType())
        _Δm::FT = Δm(prs, CT.IceType())
        _λ::FT = lambda(prs, CT.IceType(), q.ice, ρ)

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
    accretion(prs, cloud, precip, q_clo, q_pre, ρ)

 - `prs` - abstract set with Earth parameters
 - `cloud` - type for cloud water or cloud ice
 - `precip` - type for rain or snow
 - `q_clo` - cloud water or cloud ice specific humidity
 - `q_pre` - rain water or snow specific humidity
 - `ρ` - rain water or snow specific humidity

Returns the source of precipitating water (rain or snow)
due to collisions with cloud water (liquid or ice).
"""
function accretion(
    prs::APS,
    cloud::CT.AbstractCloudType,
    precip::CT.AbstractPrecipType,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        _n0::FT = n0(prs, q_pre, ρ, precip)
        _r0::FT = r0(prs, precip)
        _χv::FT = χv(prs, precip)
        _v0::FT = v0(prs, ρ, precip)
        _ve::FT = ve(prs, precip)
        _Δv::FT = Δv(prs, precip)
        _a0::FT = a0(prs, precip)
        _ae::FT = ae(prs, precip)
        _χa::FT = χa(prs, precip)
        _Δa::FT = Δa(prs, precip)
        _λ::FT = lambda(prs, precip, q_pre, ρ)
        _E::FT = E(prs, cloud, precip)

        accr_rate =
            q_clo * _E * _n0 * _a0 * _v0 * _χa * _χv / _λ *
            SF.gamma(_ae + _ve + _Δa + _Δv + FT(1)) /
            (_λ * _r0)^(_ae + _ve + _Δa + _Δv)
    end
    return accr_rate
end

"""
    accretion_rain_sink(prs, q_ice, q_rai, ρ)

 - `prs` - abstract set with Earth parameters
 - `q_ice` - cloud ice specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

Returns the sink of rain water (partial source of snow) due to collisions
with cloud ice.
"""
function accretion_rain_sink(
    prs::APS,
    q_ice::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_ice > FT(0) && q_rai > FT(0))

        _n0_ice::FT = n0(prs, FT(0), ρ, CT.IceType())
        _n0_rai::FT = n0(prs, q_rai, ρ, CT.RainType())
        _r0_rai::FT = r0(prs, CT.RainType())
        _m0_rai::FT = m0(prs, CT.RainType())
        _me_rai::FT = me(prs, CT.RainType())
        _Δm_rai::FT = Δm(prs, CT.RainType())
        _χm_rai::FT = χm(prs, CT.RainType())
        _χv_rai::FT = χv(prs, CT.RainType())
        _v0_rai::FT = v0(prs, ρ, CT.RainType())
        _ve_rai::FT = ve(prs, CT.RainType())
        _Δv_rai::FT = Δv(prs, CT.RainType())
        _a0_rai::FT = a0(prs, CT.RainType())
        _ae_rai::FT = ae(prs, CT.RainType())
        _χa_rai::FT = χa(prs, CT.RainType())
        _Δa_rai::FT = Δa(prs, CT.RainType())
        _E::FT = E(prs, CT.IceType(), CT.RainType())

        _λ_rai::FT = lambda(prs, CT.RainType(), q_rai, ρ)
        _λ_ice::FT = lambda(prs, CT.IceType(), q_ice, ρ)

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
            )^FT(_me_rai + _ae_rai + _ve_rai + _Δm_rai + _Δa_rai + _Δv_rai)
    end
    return accr_rate
end

"""
    accretion_snow_rain(prs, type_i, type_j, q_i, q_j, ρ)

 - `i` - snow for temperatures below freezing
         or rain for temperatures above freezing
 - `j` - rain for temperatures below freezing
         or snow for temperatures above freezing
 - `prs` - abstract set with Earth parameters
 - `type_i`, `type_j` - a type for snow or rain
 - `q_` - specific humidity of snow or rain
 - `ρ` - air density

Returns the accretion rate between rain and snow.
Collisions between rain and snow result in
snow at temperatures below freezing and in rain at temperatures above freezing.
"""
function accretion_snow_rain(
    prs::APS,
    type_i::CT.AbstractPrecipType,
    type_j::CT.AbstractPrecipType,
    q_i::FT,
    q_j::FT,
    ρ::FT,
) where {FT <: Real}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        _n0_i::FT = n0(prs, q_i, ρ, type_i)
        _n0_j::FT = n0(prs, q_j, ρ, type_j)

        _r0_j::FT = r0(prs, type_j)
        _m0_j::FT = m0(prs, type_j)
        _me_j::FT = me(prs, type_j)
        _Δm_j::FT = Δm(prs, type_j)
        _χm_j::FT = χm(prs, type_j)

        _E_ij::FT = E(prs, type_i, type_j)

        _λ_i::FT = lambda(prs, type_i, q_i, ρ)
        _λ_j::FT = lambda(prs, type_j, q_j, ρ)

        _v_ti = terminal_velocity(prs, type_i, CT.Blk1MVelType(), ρ, q_i)
        _v_tj = terminal_velocity(prs, type_j, CT.Blk1MVelType(), ρ, q_j)

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
    evaporation_sublimation(prs, rain, q, q_rai, ρ, T)
    evaporation_sublimation(prs, snow, q, q_sno, ρ, T)

 - `prs` - abstract set with Earth parameters
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
    prs::APS,
    rain::CT.RainType,
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    thermo_params = CMP.thermodynamics_params(prs)
    _S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && _S < FT(0))

        _ν_air::FT = CMP.ν_air(prs)
        _D_vapor::FT = CMP.D_vapor(prs)

        _G::FT = CO.G_func(prs, T, TD.Liquid())

        _n0::FT = n0(prs, q_rai, ρ, rain)
        _r0::FT = r0(prs, rain)
        _χv::FT = χv(prs, rain)
        _v0::FT = v0(prs, ρ, rain)
        _ve::FT = ve(prs, rain)
        _Δv::FT = Δv(prs, rain)

        _a_vent::FT = a_vent(prs, rain)
        _b_vent::FT = b_vent(prs, rain)

        _λ::FT = lambda(prs, rain, q_rai, ρ)

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
    prs::APS,
    snow::CT.SnowType,
    q::TD.PhasePartition{FT},
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    if q_sno > FT(0)
        _ν_air::FT = CMP.ν_air(prs)
        _D_vapor::FT = CMP.D_vapor(prs)

        thermo_params = CMP.thermodynamics_params(prs)
        _S::FT = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        _G::FT = CO.G_func(prs, T, TD.Ice())

        _n0::FT = n0(prs, q_sno, ρ, snow)
        _r0::FT = r0(prs, snow)
        _χv::FT = χv(prs, snow)
        _v0::FT = v0(prs, ρ, snow)
        _ve::FT = ve(prs, snow)
        _Δv::FT = Δv(prs, snow)

        _a_vent::FT = a_vent(prs, snow)
        _b_vent::FT = b_vent(prs, snow)

        _λ::FT = lambda(prs, snow, q_sno, ρ)

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
    snow_melt(prs, q_sno, ρ, T)

 - `prs` - abstract set with Earth parameters
 - `q_sno` - snow water specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to snow melt.
"""
function snow_melt(prs::APS, q_sno::FT, ρ::FT, T::FT) where {FT <: Real}

    snow_melt_rate = FT(0)
    _T_freeze::FT = CMP.T_freeze(prs)

    if (q_sno > FT(0) && T > _T_freeze)

        _ν_air::FT = CMP.ν_air(prs)
        _D_vapor::FT = CMP.D_vapor(prs)
        _K_therm::FT = CMP.K_therm(prs)

        thermo_params = CMP.thermodynamics_params(prs)
        L = TD.latent_heat_fusion(thermo_params, T)

        snow = CT.SnowType()

        _n0::FT = n0(prs, q_sno, ρ, snow)
        _r0::FT = r0(prs, snow)
        _χv::FT = χv(prs, snow)
        _v0::FT = v0(prs, ρ, snow)
        _ve::FT = ve(prs, snow)
        _Δv::FT = Δv(prs, snow)

        _a_vent::FT = a_vent(prs, snow)
        _b_vent::FT = b_vent(prs, snow)

        _λ::FT = lambda(prs, snow, q_sno, ρ)

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
