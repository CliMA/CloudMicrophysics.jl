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

export terminal_velocity,
    conv_q_liq_to_q_rai,
    conv_q_ice_to_q_sno,
    accretion,
    accretion_rain_sink,
    accretion_snow_rain,
    evaporation_sublimation,
    snow_melt

"""
    collision_efficiency(type_1, type_2, ce)

 - `type_1` and `type_2` - types of colliding species
 - `ce` - CollisionEfficiency model parameters

Returns collision efficiency for two colliding species
"""
collision_efficiency(
    ::CT.LiquidType,
    ::CT.RainType,
    (; e_liq_rai)::CT.CollisionEfficiency,
) = e_liq_rai
collision_efficiency(
    ::CT.LiquidType,
    ::CT.SnowType,
    (; e_liq_sno)::CT.CollisionEfficiency,
) = e_liq_sno
collision_efficiency(
    ::CT.IceType,
    ::CT.RainType,
    (; e_ice_rai)::CT.CollisionEfficiency,
) = e_ice_rai
collision_efficiency(
    ::CT.IceType,
    ::CT.SnowType,
    (; e_ice_sno)::CT.CollisionEfficiency,
) = e_ice_sno
collision_efficiency(
    ::CT.RainType,
    ::CT.SnowType,
    (; e_rai_sno)::CT.CollisionEfficiency,
) = e_rai_sno
collision_efficiency(
    ::CT.SnowType,
    ::CT.RainType,
    (; e_rai_sno)::CT.CollisionEfficiency,
) = e_rai_sno

"""
    get_n0(type, q_sno, ρ)

 - `type` - supports `SnowType`, `IceType`, and `RainType`
 - `q_sno` -  snow specific humidity (only for `SnowType`)
 - `ρ` - air density (only for `SnowType`)

Returns the intercept parameter of the assumed Marshall-Palmer distribution
"""
get_n0(
    (; ν, μ)::CT.SnowType{FT},
    q_sno::FT,
    ρ::FT,
) where {FT <: AbstractFloat} = μ * (ρ * max(0, q_sno))^ν
# TODO               this max should be replaced by
#                    limiting inside a PhasePartition struct for
#                    precipitation (once it is implemented)
get_n0((; n0)::Union{CT.IceType, CT.RainType}, args...) = n0

"""
    get_v0(type, ρ)

 - `type` - supports `RainType`and `SnowType`
 - `ρ` - air density (only for `RainType`)

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
get_v0(
    (; ρ_cloud_liq, C_drag, grav, r0)::CT.RainType,
    ρ::FT,
) where {FT <: AbstractFloat} =
    sqrt(FT(8 / 3) / C_drag * (ρ_cloud_liq / ρ - FT(1)) * grav * r0)
get_v0((; v0)::CT.SnowType, args...) = v0

"""
    lambda(precip, q, ρ)

 - `precip` - a type for cloud ice, rain or snow
 - `q` - specific humidity of rain, ice or snow
 - `ρ` - air density

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(
    precip::Union{CT.IceType, CT.RainType, CT.SnowType},
    q::FT,
    ρ::FT,
) where {FT <: AbstractFloat}

    n0::FT = get_n0(precip, q, ρ)
    (; r0, m0, me, Δm, χm) = precip

    return q > FT(0) ?
           (
        χm * m0 * n0 * SF.gamma(me + Δm + FT(1)) / ρ / q / r0^(me + Δm)
    )^FT(1 / (me + Δm + 1)) : FT(0)
end

"""
    terminal_velocity(precip, velo_scheme, ρ, q)

 - `precip` - a type for ice, rain or snow
 - `velo_scheme` - type for terminal velocity parameterization
 - `ρ` - air density
 - `q` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of particles.
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical power-law relations for `velo_scheme == Blk1MVelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171, for `velo_scheme == Chen2022Type`
"""
function terminal_velocity(
    precip::CT.AbstractPrecipType,
    velo_scheme::CT.Blk1MVelType,
    ρ::FT,
    q::FT,
) where {FT <: AbstractFloat}
    if q > FT(0)
        (; r0, me, Δm, χm) = precip
        (; χv, ve, Δv) = velo_scheme
        v0 = get_v0(precip, ρ)
        λ = lambda(precip, q, ρ)

        return χv *
               v0 *
               (λ * r0)^(-ve - Δv) *
               SF.gamma(me + ve + Δm + Δv + FT(1)) / SF.gamma(me + Δm + FT(1))
    else
        return FT(0)
    end
end
function terminal_velocity(
    precip::CT.RainType,
    velo_scheme::CT.Chen2022Type,
    ρ::FT,
    q::FT,
) where {FT <: AbstractFloat}
    fall_w = FT(0)
    if q > FT(0)

        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(precip, velo_scheme, ρ)
        # size distribution parameter
        _λ::FT = lambda(precip, q, ρ)

        # eq 20 from Chen et al 2022
        fall_w = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, _λ, 3))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end
function terminal_velocity(
    precip::CT.IceType,
    velo_scheme::CT.Chen2022Type,
    ρ::FT,
    q::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q > FT(0)

        ρ_i = precip.ρ #::FT = CMP.ρ_cloud_ice(prs)
        _λ::FT = lambda(precip, q, ρ)

        # coefficients from Appendix B from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(precip, velo_scheme, ρ)

        # eq 20 from Chen et al 2022
        fall_w = sum(CO.Chen2022_vel_add.(aiu, bi, ciu, _λ, 3))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end
function terminal_velocity(
    precip::CT.SnowType,
    velo_scheme::CT.Chen2022Type,
    ρ::FT,
    q::FT,
) where {FT <: Real}
    fall_w = FT(0)
    if q > FT(0)

        (; r0, m0, me, Δm, χm, a0, ae, Δa, χa) = precip
        λ::FT = lambda(precip, q, ρ)

        m0c = m0 * χm
        a0c = a0 * χa
        mec = me + Δm
        aec = ae + Δa

        ρ_i = precip.ρ #::FT = CMP.ρ_cloud_ice(prs)

        # coefficients from Appendix B from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(precip, velo_scheme, ρ)

        κ = FT(-1 / 3) #oblate
        k = 3 # mass weighted

        tmp =
            λ^(k + 1) *
            ((16 * a0c^3 * ρ_i^2) / (9 * π * m0c^2 * r0^(3 * aec - 2 * mec)))^κ
        ci_pow =
            (2 .* ciu .+ λ) .^
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
    conv_q_liq_to_q_rai(acnv, q_liq)

 - `acnv` - 1M autoconversion parameters
 - `q_liq` - liquid water specific humidity

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Kessler (1995).
"""
conv_q_liq_to_q_rai(
    (; τ, q_threshold, k)::CT.Autoconversion1M{FT},
    q_liq::FT;
    smooth_transition::Bool = false,
) where {FT <: AbstractFloat} =
    smooth_transition ?
    CO.logistic_function_integral(q_liq, q_threshold, k) / τ :
    max(0, q_liq - q_threshold) / τ

"""
    conv_q_ice_to_q_sno_no_supersat(acnv, q_ice)

 - `acnv` - 1M autoconversion parameters
 - `q_ice` -  cloud ice specific humidity

Returns the q_sno tendency due to autoconversion from ice.
This is a simplified version of a snow autoconversion rate that can be used in
simulations where there is no supersaturation
(for example in TC.jl when using saturation adjustment).
"""
conv_q_ice_to_q_sno_no_supersat(
    (; τ, q_threshold, k)::CT.Autoconversion1M{FT},
    q_ice::FT;
    smooth_transition::Bool = false,
) where {FT <: AbstractFloat} =
    smooth_transition ?
    CO.logistic_function_integral(q_ice, q_threshold, k) / τ :
    max(0, q_ice - q_threshold) / τ


"""
    conv_q_ice_to_q_sno(ice_type, air_props, thermo_params, q, ρ, T)

 - `ice_type` - IceType
 - `air_props` - air properties
 - `thermo_params` - thermodynamics parameters
 - `q` - phase partition
 - `ρ` - air density
 - `T` - air temperature

Returns the q_sno tendency due to autoconversion from ice.
Parameterized following Harrington et al. (1996) and Kaul et al. (2015).
"""
function conv_q_ice_to_q_sno(
    ice_type::CT.IceType{FT},
    air_props,
    thermo_params,
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT <: AbstractFloat}
    acnv_rate = FT(0)
    S = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())

    if (q.ice > FT(0) && S > FT(0))
        (; r_ice_snow, me, Δm) = ice_type
        G = CO.G_func(air_props, thermo_params, T, TD.Ice())
        n0 = get_n0(ice_type, FT(0), ρ)
        λ = lambda(ice_type, q.ice, ρ)

        acnv_rate =
            4 * FT(π) * S * G * n0 / ρ *
            exp(-λ * r_ice_snow) *
            (r_ice_snow^FT(2) / (me + Δm) + (r_ice_snow * λ + FT(1)) / λ^FT(2))
    end
    return acnv_rate
end

"""
    accretion(ce, cloud, precip, q_clo, q_pre, ρ)

 - `ce` - collision efficiency parameters
 - `cloud` - type for cloud water or cloud ice
 - `precip` - type for rain or snow
 - `q_clo` - cloud water or cloud ice specific humidity
 - `q_pre` - rain water or snow specific humidity
 - `ρ` - rain water or snow specific humidity

Returns the source of precipitating water (rain or snow)
due to collisions with cloud water (liquid or ice).
"""
function accretion(
    ce::CT.CollisionEfficiency,
    cloud::CT.AbstractCloudType,
    precip::CT.AbstractPrecipType,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT <: AbstractFloat}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        n0::FT = get_n0(precip, q_pre, ρ)
        v0::FT = get_v0(precip, ρ)
        (; r0, χv, ve, Δv, a0, ae, χa, Δa) = precip
        λ = lambda(precip, q_pre, ρ)
        E = collision_efficiency(cloud, precip, ce)

        accr_rate =
            q_clo * E * n0 * a0 * v0 * χa * χv / λ *
            SF.gamma(ae + ve + Δa + Δv + FT(1)) / (λ * r0)^(ae + ve + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_rain_sink(rain_type, ice_type, ce, q_ice, q_rai, ρ)

 - `rain_type` - rain type parameters
 - `ice_type` - ice type parameters
 - `ce` - collision efficiency parameters
 - `q_ice` - cloud ice specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

Returns the sink of rain water (partial source of snow) due to collisions
with cloud ice.
"""
function accretion_rain_sink(
    rain_type::CT.RainType,
    ice_type::CT.IceType,
    ce::CT.CollisionEfficiency,
    q_ice::FT,
    q_rai::FT,
    ρ::FT,
) where {FT <: AbstractFloat}
    accr_rate = FT(0)
    if (q_ice > FT(0) && q_rai > FT(0))

        n0_ice = get_n0(ice_type, FT(0), ρ)
        λ_ice = lambda(ice_type, q_ice, ρ)

        n0 = get_n0(rain_type, q_rai, ρ)
        v0 = get_v0(rain_type, ρ)
        (; r0, m0, me, Δm, χm, χv, ve, Δv) = rain_type
        (; a0, ae, χa, Δa) = rain_type

        E = collision_efficiency(ice_type, rain_type, ce)

        λ = lambda(rain_type, q_rai, ρ)

        accr_rate =
            E / ρ * n0 * n0_ice * m0 * a0 * v0 * χm * χa * χv / λ_ice / λ *
            SF.gamma(me + ae + ve + Δm + Δa + Δv + FT(1)) /
            (r0 * λ)^FT(me + ae + ve + Δm + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_snow_rain(ce, type_i, type_j, blk1m_type_i, blk1m_type_j, q_i, q_j, ρ)

 - `i` - snow for temperatures below freezing
         or rain for temperatures above freezing
 - `j` - rain for temperatures below freezing
         or snow for temperatures above freezing
 - `ce` - collision efficiency parameters
 - `type_i`, `type_j` - a type for snow or rain
 - `blk1mveltype_ti`, `blk1mveltype_tj` - 1M terminal velocity model
 - `q_` - specific humidity of snow or rain
 - `ρ` - air density

Returns the accretion rate between rain and snow.
Collisions between rain and snow result in
snow at temperatures below freezing and in rain at temperatures above freezing.
"""
function accretion_snow_rain(
    ce::CT.CollisionEfficiency,
    type_i::CT.AbstractPrecipType,
    type_j::CT.AbstractPrecipType,
    blk1mveltype_ti::CT.Blk1MVelType,
    blk1mveltype_tj::CT.Blk1MVelType,
    q_i::FT,
    q_j::FT,
    ρ::FT,
) where {FT <: AbstractFloat}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        n0_i = get_n0(type_i, q_i, ρ)
        n0_j = get_n0(type_j, q_j, ρ)

        r0_j = type_j.r0
        m0_j = type_j.m0
        me_j = type_j.me
        Δm_j = type_j.Δm
        χm_j = type_j.χm

        E_ij = collision_efficiency(type_i, type_j, ce)

        λ_i = lambda(type_i, q_i, ρ)
        λ_j = lambda(type_j, q_j, ρ)

        v_ti = terminal_velocity(type_i, blk1mveltype_ti, ρ, q_i)
        v_tj = terminal_velocity(type_j, blk1mveltype_tj, ρ, q_j)

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
    evaporation_sublimation(rain, air_props, thermo_params, q, q_rai, ρ, T)
    evaporation_sublimation(snow, air_props, thermo_params, q, q_sno, ρ, T)

 - `prs` - abstract set with Earth parameters
 - `rain` - a type for rain
 - `snow` - a type for snow
 - `air_props` - air properties
 - `thermo_params` - thermodynamics parameters
 - `q` - phase partition
 - `q_rai` - rain specific humidity
 - `q_sno` - snow specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to rain evaporation or snow sublimation.
"""
function evaporation_sublimation(
    rain::CT.RainType,
    air_props,
    thermo_params,
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    T::FT,
) where {FT <: AbstractFloat}
    evap_subl_rate = FT(0)
    S = TD.supersaturation(thermo_params, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && S < FT(0))

        (; ν_air, D_vapor) = air_props

        G = CO.G_func(air_props, thermo_params, T, TD.Liquid())

        n0 = get_n0(rain, q_rai, ρ)
        v0 = get_v0(rain, ρ)
        (; r0, χv, ve, Δv, a_vent, b_vent) = rain

        λ = lambda(rain, q_rai, ρ)

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
    snow::CT.SnowType,
    air_props,
    thermo_params,
    q::TD.PhasePartition{FT},
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT <: Real}
    evap_subl_rate = FT(0)
    if q_sno > FT(0)
        (; ν_air, D_vapor) = air_props

        S = TD.supersaturation(thermo_params, q, ρ, T, TD.Ice())
        G = CO.G_func(air_props, thermo_params, T, TD.Ice())

        n0 = get_n0(snow, q_sno, ρ)
        v0 = get_v0(snow, ρ)
        (; r0, χv, ve, Δv, a_vent, b_vent) = snow

        λ = lambda(snow, q_sno, ρ)

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
    snow_melt(snow, air_props, thermo_params, q_sno, ρ, T)

 - `snow` - snow parameters
 - `air_props` - air properties
 - `thermo_params` - thermodynamics parameters
 - `q_sno` - snow water specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to snow melt.
"""
function snow_melt(
    snow::CT.SnowType{FT},
    air_props::CT.AirProperties{FT},
    thermo_params,
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT <: AbstractFloat}
    snow_melt_rate = FT(0)
    (; T_freeze) = snow

    if (q_sno > FT(0) && T > T_freeze)
        (; ν_air, D_vapor, K_therm) = air_props

        L = TD.latent_heat_fusion(thermo_params, T)


        n0 = get_n0(snow, q_sno, ρ)
        v0 = get_v0(snow, ρ)
        (; r0, χv, ve, Δv, a_vent, b_vent) = snow
        λ = lambda(snow, q_sno, ρ)

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
