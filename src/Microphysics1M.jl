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
import Thermodynamics.Parameters as TDP

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

abstract type AbstractSnowShape end
struct Oblate <: AbstractSnowShape end
struct Prolate <: AbstractSnowShape end

"""
    Ec(type_1, type_2, ce)

 - `type_1` and `type_2` - types of colliding species
 - `ce` - a struct with collision efficiency parameters

Returns collision efficiency for two colliding species
"""
Ec(::CMP.CloudLiquid, ::CMP.Rain, (; e_liq_rai)::CMP.CollisionEff) = e_liq_rai
Ec(::CMP.CloudLiquid, ::CMP.Snow, (; e_liq_sno)::CMP.CollisionEff) = e_liq_sno
Ec(::CMP.CloudIce, ::CMP.Rain, (; e_ice_rai)::CMP.CollisionEff) = e_ice_rai
Ec(::CMP.CloudIce, ::CMP.Snow, (; e_ice_sno)::CMP.CollisionEff) = e_ice_sno
Ec(::CMP.Rain, ::CMP.Snow, (; e_rai_sno)::CMP.CollisionEff) = e_rai_sno
Ec(::CMP.Snow, ::CMP.Rain, (; e_rai_sno)::CMP.CollisionEff) = e_rai_sno

"""
    get_n0(pdf, q_sno, ρ)

 - `pdf` -  a struct with parameters for snow, ice, and rain size distribution
 - `q_sno` -  snow specific humidity (only for `Snow`)
 - `ρ` - air density (only for `Snow`)

Returns the intercept parameter of the assumed Marshall-Palmer distribution
"""
get_n0((; ν, μ)::CMP.ParticlePDFSnow{FT}, q_sno::FT, ρ::FT) where {FT} =
    μ * (ρ * max(0, q_sno))^ν
# TODO               this max should be replaced by
#                    limiting inside a PhasePartition struct for
#                    precipitation (once it is implemented)
get_n0((; n0)::CMP.ParticlePDFIceRain{FT}, args...) where {FT} = n0

"""
    get_v0(v, ρ)

 - `v` - a struct with bulk 1-moment terminal velocity parameters
 - `ρ` - air density (only for `Rain`)

Returns the proportionality coefficient in terminal velocity(r/r0).
"""
get_v0((; C_drag, ρw, grav, r0)::CMP.Blk1MVelTypeRain{FT}, ρ::FT) where {FT} =
    sqrt(FT(8 / 3) / C_drag * (ρw / ρ - FT(1)) * grav * r0)
get_v0((; v0)::CMP.Blk1MVelTypeSnow{FT}, args...) where {FT} = v0

"""
    lambda(pdf, mass, q, ρ)

 - `pdf`, `mass` - structs with particle size distribution and mass parameters
 - `q` - specific humidity of rain, ice or snow
 - `ρ` - air density

Returns the rate parameter of the assumed size distribution of
particles (rain drops, ice crystals, snow crystals).
"""
function lambda(
    #(; pdf, mass)::Union{CMP.Snow{FT}, CMP.Rain{FT}, CMP.CloudIce{FT}},
    pdf::Union{CMP.ParticlePDFIceRain{FT}, CMP.ParticlePDFSnow{FT}},
    mass::CMP.ParticleMass{FT},
    q::FT,
    ρ::FT,
) where {FT}
    # size distribution
    n0::FT = get_n0(pdf, q, ρ)
    # mass(size)
    (; r0, m0, me, Δm, χm) = mass

    return q > FT(0) ?
           (
        χm * m0 * n0 * CO.Γ(me + Δm + FT(1)) / ρ / q / r0^(me + Δm)
    )^FT(1 / (me + Δm + 1)) : FT(0)
end

"""
    radar_reflectivity(precip, q, ρ)

    - `precip` - struct with rain free parameters
    - `q` - specific humidity of rain
    - `ρ` - air density

Returns logarithmic radar reflectivity from the assumed rain particle size distribution
normalized by the reflectivty of 1 millimiter drop in a volume of one meter cube
"""
function radar_reflectivity(
    (; pdf, mass)::CMP.Rain{FT},
    q::FT,
    ρ::FT,
) where {FT}

    # change units for accuracy
    n0 = get_n0(pdf) * FT(1e-12)
    λ = lambda(pdf, mass, q, ρ) * FT(1e-3)

    Z = (λ == FT(0)) ? FT(0) : (720 * n0 / λ^7)
    log_10_Z₀ = FT(-18)
    log_Z = FT(10) * (log10(Z) - log_10_Z₀ - FT(9))

    return max(FT(-430), log_Z)
end

"""
    aspect_ratio_coeffs(snow_shape, mass, area, density)

 - `snow_shape` - a struct specifying assumed snow particle shape (Oblate or Prolate)
 - `mass` - a struct with assumed m(r) power law relation parameters
 - `area` - a struct with assumed a(r) power law relation parameters
 - `ρᵢ` particle density

Returns coefficients of the implied power law relationship between aspect ratio
and particle diameter ϕ(D) = ϕ₀ D^α
Also returns the coefficient for the aspect ratio in Chen 2022 terminal velocity
parameterization (κ=1/3 for oblate and κ=-1/6 for prolate).
"""
function aspect_ratio_coeffs(
    snow_shape::Oblate,
    (; r0, m0, me, Δm, χm)::CMP.ParticleMass{FT},
    (; a0, ae, Δa, χa)::CMP.ParticleArea{FT},
    ρᵢ::FT,
) where {FT}
    # ϕ(r) = 3 * sqrt(FT(π)) * mᵢ(r) / (4 * ρᵢ * aᵢ(r)^(3/2)
    α = me + Δm - 3/2 * (ae + Δa)
    ϕ₀ = 3 * sqrt(FT(π)) / 4 / ρᵢ * χm * m0 / (χa * a0)^(3/2) / (2 * r0)^α
    κ = FT(1 / 3)
    return (; ϕ₀, α, κ)
end
function aspect_ratio_coeffs(
    snow_shape::Prolate,
    (; r0, m0, me, Δm, χm)::CMP.ParticleMass{FT},
    (; a0, ae, Δa, χa)::CMP.ParticleArea{FT},
    ρᵢ::FT,
) where {FT}
    # ϕ(r) = 16 * ρᵢ^2 * aᵢ(r)^3 / (9 * π * mᵢ(r)^2)
    α = 3 * (ae + Δa) - 2 * (me + Δm)
    ϕ₀ = 16 * ρᵢ^2 / 9 / FT(π) * (χa * a0)^3 / (χm * m0)^2 / (2 *r0)^α
    κ = FT(-1 / 6)
    return (; ϕ₀, α, κ)
end

"""
    terminal_velocity(precip, vel, ρ, q)

 - `precip` - a struct with precipitation type (rain or snow)
 - `vel` - a struct with terminal velocity parameters
 - `ρ` - air density
 - `q` - rain or snow specific humidity

Returns the mass weighted average terminal velocity assuming
a Marshall-Palmer (1948) distribution of particles.
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical power-law relations for `velocity == Blk1MVelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171, for `velocity == Chen2022VelType`
"""
function terminal_velocity(
    (; pdf, mass)::Union{CMP.Rain{FT}, CMP.Snow{FT}},
    vel::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    ρ::FT,
    q::FT,
) where {FT}
    if q > FT(0)
        # terminal_velocity(size)
        (; χv, ve, Δv) = vel
        v0 = get_v0(vel, ρ)
        # mass(size)
        (; r0, me, Δm, χm) = mass
        # size distrbution
        λ = lambda(pdf, mass, q, ρ)

        return χv * v0 * (λ * r0)^(-ve - Δv) * CO.Γ(me + ve + Δm + Δv + FT(1)) /
               CO.Γ(me + Δm + FT(1))
    else
        return FT(0)
    end
end
function terminal_velocity(
    (; pdf, mass)::CMP.Rain{FT},
    vel::CMP.Chen2022VelTypeRain{FT},
    ρ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > FT(0)
        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs_B1(vel, ρ)
        # size distribution parameter
        λ::FT = lambda(pdf, mass, q, ρ)
        # eq 20 from Chen et al 2022
        fall_w = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, λ, 3))
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end
function terminal_velocity(
    (; pdf, mass, aspr)::CMP.Snow{FT},
    vel::CMP.Chen2022VelTypeSnowIce{FT},
    ρ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    # We assume the B4 table coeffs for snow and B2 table coeffs for cloud ice.
    # aiu, bi, ciu = CO.Chen2022_vel_coeffs_B2(vel, ρ)
    # Instead we should do partial integrals
    # from D=125um to D=625um using B2 and D=625um to inf using B4.
    if q > FT(0)
        # coefficients from Table B4 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs_B4(vel, ρ)
        # size distribution parameter
        λ::FT = lambda(pdf, mass, q, ρ)

        # assume oblate shape and aspect ratio
        (; ϕ, κ) = aspr

        # eq 20 from Chen 2022
        fall_w = ϕ^κ * sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, λ, 3))
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end
function terminal_velocity(
    (; pdf, mass, area)::CMP.Snow{FT},
    vel::CMP.Chen2022VelTypeSnowIce{FT},
    ρ::FT,
    q::FT,
    snow_shape::AbstractSnowShape,
) where {FT}
    fall_w = FT(0)
    # see comments above about B2 vs B4 coefficients
    if q > FT(0)
        # coefficients from Table B4 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs_B4(vel, ρ)
        # size distribution parameter
        λ::FT = lambda(pdf, mass, q, ρ)
        # Compute the mass weighted average aspect ratio ϕ_av
        # As a next step, we could keep ϕ(r) under the integrals
        (ϕ₀, α, κ) = aspect_ratio_coeffs(snow_shape, mass, area, vel.ρᵢ)
        ϕ_av = ϕ₀ / λ^α * CO.Γ(α + 3 + 1) / CO.Γ(3 + 1)
        # eq 20 from Chen 2022
        fall_w = ϕ_av^κ * sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, λ, 3))
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end

"""
    conv_q_liq_to_q_rai(acnv, q_liq, smooth_transition)

 - `acnv` - 1M autoconversion parameters
 - `q_liq` - liquid water specific humidity
 - `smooth_transition` - a flag to switch on smoothing

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following Kessler (1995).
"""
conv_q_liq_to_q_rai(
    (; τ, q_threshold, k)::CMP.Acnv1M{FT},
    q_liq::FT,
    smooth_transition::Bool = false,
) where {FT} =
    smooth_transition ?
    CO.logistic_function_integral(q_liq, q_threshold, k) / τ :
    max(0, q_liq - q_threshold) / τ

"""
    conv_q_ice_to_q_sno_no_supersat(acnv, q_ice, smooth_transition)

 - `acnv` - 1M autoconversion parameters
 - `q_ice` -  cloud ice specific humidity
 - `smooth_transition` - a flag to switch on smoothing

Returns the q_sno tendency due to autoconversion from ice.
This is a simplified version of a snow autoconversion rate that can be used in
simulations where there is no supersaturation
(for example in TC.jl when using saturation adjustment).
"""
conv_q_ice_to_q_sno_no_supersat(
    (; τ, q_threshold, k)::CMP.Acnv1M{FT},
    q_ice::FT,
    smooth_transition::Bool = false,
) where {FT} =
    smooth_transition ?
    CO.logistic_function_integral(q_ice, q_threshold, k) / τ :
    max(0, q_ice - q_threshold) / τ

"""
    conv_q_ice_to_q_sno(ice, aps, tps, q, ρ, T)

 - `ice` - a struct with ice parameters
 - `aps` - a struct with air properties
 - `tps` - a struct with thermodynamics parameters
 - `q` - phase partition
 - `ρ` - air density
 - `T` - air temperature

Returns the q_sno tendency due to autoconversion from ice.
Parameterized following Harrington et al. (1996) and Kaul et al. (2015).
"""
function conv_q_ice_to_q_sno(
    (; r_ice_snow, pdf, mass)::CMP.CloudIce{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    ρ::FT,
    T::FT,
) where {FT}
    acnv_rate = FT(0)
    S = TD.supersaturation(tps, q, ρ, T, TD.Ice())

    if (q.ice > FT(0) && S > FT(0))
        (; me, Δm) = mass
        G = CO.G_func(aps, tps, T, TD.Ice())
        n0 = get_n0(pdf)
        λ = lambda(pdf, mass, q.ice, ρ)

        acnv_rate =
            4 * FT(π) * S * G * n0 / ρ *
            exp(-λ * r_ice_snow) *
            (r_ice_snow^FT(2) / (me + Δm) + (r_ice_snow * λ + FT(1)) / λ^FT(2))
    end
    return acnv_rate
end

"""
    accretion(cloud, precip, vel, ce, q_clo, q_pre, ρ)

 - `cloud` - type for cloud water or cloud ice
 - `precip` - type for rain or snow
 - `vel` - a struct with terminal velocity parameters
 - `ce` - collision efficiency parameters
 - `q_clo` - cloud water or cloud ice specific humidity
 - `q_pre` - rain water or snow specific humidity
 - `ρ` - rain water or snow specific humidity

Returns the source of precipitating water (rain or snow)
due to collisions with cloud water (liquid or ice).
"""
function accretion(
    cloud::CMP.CloudCondensateType{FT},
    precip::CMP.PrecipitationType{FT},
    vel::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    ce::CMP.CollisionEff,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT}

    accr_rate = FT(0)
    if (q_clo > FT(0) && q_pre > FT(0))

        n0::FT = get_n0(precip.pdf, q_pre, ρ)
        v0::FT = get_v0(vel, ρ)

        (; r0) = precip.mass
        (; χv, ve, Δv) = vel
        (; a0, ae, χa, Δa) = precip.area

        λ = lambda(precip.pdf, precip.mass, q_pre, ρ)
        E = Ec(cloud, precip, ce)

        accr_rate =
            q_clo * E * n0 * a0 * v0 * χa * χv / λ *
            CO.Γ(ae + ve + Δa + Δv + FT(1)) / (λ * r0)^(ae + ve + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_rain_sink(rain, ice, vel, ce, q_ice, q_rai, ρ)

 - `rain` - rain type parameters
 - `ice` - ice type parameters
 - `vel` - terminal velocity parameters for rain
 - `ce` - collision efficiency parameters
 - `q_ice` - cloud ice specific humidity
 - `q_rai` - rain water specific humidity
 - `ρ` - air density

Returns the sink of rain water (partial source of snow) due to collisions
with cloud ice.
"""
function accretion_rain_sink(
    rain::CMP.Rain{FT},
    ice::CMP.CloudIce{FT},
    vel::CMP.Blk1MVelTypeRain{FT},
    ce::CMP.CollisionEff,
    q_ice::FT,
    q_rai::FT,
    ρ::FT,
) where {FT}
    accr_rate = FT(0)
    if (q_ice > FT(0) && q_rai > FT(0))

        n0_ice = get_n0(ice.pdf)
        λ_ice = lambda(ice.pdf, ice.mass, q_ice, ρ)

        n0 = get_n0(rain.pdf, q_rai, ρ)
        v0 = get_v0(vel, ρ)
        (; r0, m0, me, Δm, χm) = rain.mass
        (; χv, ve, Δv) = vel
        (; a0, ae, χa, Δa) = rain.area

        E = Ec(ice, rain, ce)

        λ = lambda(rain.pdf, rain.mass, q_rai, ρ)

        accr_rate =
            E / ρ * n0 * n0_ice * m0 * a0 * v0 * χm * χa * χv / λ_ice / λ *
            CO.Γ(me + ae + ve + Δm + Δa + Δv + FT(1)) /
            (r0 * λ)^FT(me + ae + ve + Δm + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_snow_rain(ce, type_i, type_j, blk1m_type_i, blk1m_type_j, q_i, q_j, ρ)

 - `ce` - collision efficiency parameters
 - `i` - snow for temperatures below freezing
         or rain for temperatures above freezing
 - `j` - rain for temperatures below freezing
         or snow for temperatures above freezing
 - `type_i`, `type_j` - a type for snow or rain
 - `blk1mveltype_ti`, `blk1mveltype_tj` - 1M terminal velocity parameters
 - `q_` - specific humidity of snow or rain
 - `ρ` - air density

Returns the accretion rate between rain and snow.
Collisions between rain and snow result in
snow at temperatures below freezing and in rain at temperatures above freezing.
"""
function accretion_snow_rain(
    type_i::CMP.PrecipitationType{FT},
    type_j::CMP.PrecipitationType{FT},
    blk1mveltype_ti::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    blk1mveltype_tj::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    ce::CMP.CollisionEff,
    q_i::FT,
    q_j::FT,
    ρ::FT,
) where {FT}

    accr_rate = FT(0)
    if (q_i > FT(0) && q_j > FT(0))

        n0_i = get_n0(type_i.pdf, q_i, ρ)
        n0_j = get_n0(type_j.pdf, q_j, ρ)

        r0_j = type_j.mass.r0
        m0_j = type_j.mass.m0
        me_j = type_j.mass.me
        Δm_j = type_j.mass.Δm
        χm_j = type_j.mass.χm

        E_ij = Ec(type_i, type_j, ce)

        λ_i = lambda(type_i.pdf, type_i.mass, q_i, ρ)
        λ_j = lambda(type_j.pdf, type_j.mass, q_j, ρ)

        v_ti = terminal_velocity(type_i, blk1mveltype_ti, ρ, q_i)
        v_tj = terminal_velocity(type_j, blk1mveltype_tj, ρ, q_j)

        accr_rate =
            FT(π) / ρ * n0_i * n0_j * m0_j * χm_j * E_ij * abs(v_ti - v_tj) /
            r0_j^(me_j + Δm_j) * (
                FT(2) * CO.Γ(me_j + Δm_j + FT(1)) / λ_i^FT(3) /
                λ_j^(me_j + Δm_j + FT(1)) +
                FT(2) * CO.Γ(me_j + Δm_j + FT(2)) / λ_i^FT(2) /
                λ_j^(me_j + Δm_j + FT(2)) +
                CO.Γ(me_j + Δm_j + FT(3)) / λ_i / λ_j^(me_j + Δm_j + FT(3))
            )
    end
    return accr_rate
end

"""
    evaporation_sublimation(rain, vel, aps, tps, q, q_rai, ρ, T)
    evaporation_sublimation(snow, vel, aps, tps, q, q_sno, ρ, T)

 - `rain` - a struct with rain parameters
 - `snow` - a struct with snow parameters
 - `vel` - a struct with terminal velocity parameters
 - `aps` - a struct with air parameters
 - `tps` - a struct with thermodynamics parameters
 - `q` - phase partition
 - `q_rai` - rain specific humidity
 - `q_sno` - snow specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to rain evaporation or snow sublimation.
"""
function evaporation_sublimation(
    (; pdf, mass, vent)::CMP.Rain{FT},
    vel::CMP.Blk1MVelTypeRain{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    T::FT,
) where {FT}
    evap_subl_rate = FT(0)
    S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())

    if (q_rai > FT(0) && S < FT(0))

        (; ν_air, D_vapor) = aps
        G = CO.G_func(aps, tps, T, TD.Liquid())
        n0 = get_n0(pdf, q_rai, ρ)
        v0 = get_v0(vel, ρ)
        (; χv, ve, Δv) = vel
        (; r0) = mass
        a_vent = vent.a
        b_vent = vent.b

        λ = lambda(pdf, mass, q_rai, ρ)

        evap_subl_rate =
            4 * FT(π) * n0 / ρ * S * G / λ^FT(2) * (
                a_vent +
                b_vent * (ν_air / D_vapor)^FT(1 / 3) /
                (r0 * λ)^((ve + Δv) / FT(2)) *
                (FT(2) * v0 * χv / ν_air / λ)^FT(1 / 2) *
                CO.Γ((ve + Δv + FT(5)) / FT(2))
            )
    end
    # only evaporation is considered for rain
    return min(0, evap_subl_rate)
end
function evaporation_sublimation(
    (; pdf, mass, vent)::CMP.Snow{FT},
    vel::CMP.Blk1MVelTypeSnow{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT}
    evap_subl_rate = FT(0)
    if q_sno > FT(0)
        (; ν_air, D_vapor) = aps

        S = TD.supersaturation(tps, q, ρ, T, TD.Ice())
        G = CO.G_func(aps, tps, T, TD.Ice())

        n0 = get_n0(pdf, q_sno, ρ)
        v0 = get_v0(vel, ρ)
        (; r0) = mass
        (; χv, ve, Δv) = vel

        a_vent = vent.a
        b_vent = vent.b

        λ = lambda(pdf, mass, q_sno, ρ)

        evap_subl_rate =
            4 * FT(π) * n0 / ρ * S * G / λ^FT(2) * (
                a_vent +
                b_vent * (ν_air / D_vapor)^FT(1 / 3) /
                (r0 * λ)^((ve + Δv) / FT(2)) *
                (FT(2) * v0 * χv / ν_air / λ)^FT(1 / 2) *
                CO.Γ((ve + Δv + FT(5)) / FT(2))
            )
    end
    return evap_subl_rate
end

"""
    snow_melt(snow, vel, aps, tps, q_sno, ρ, T)

 - `snow` - snow parameters
 - `vel` - terminal velocity parameters
 - `aps` - air properties
 - `tps` - thermodynamics parameters
 - `q_sno` - snow water specific humidity
 - `ρ` - air density
 - `T` - air temperature

Returns the tendency due to snow melt.
"""
function snow_melt(
    (; T_freeze, pdf, mass, vent)::CMP.Snow{FT},
    vel::CMP.Blk1MVelTypeSnow{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT}
    snow_melt_rate = FT(0)

    if (q_sno > FT(0) && T > T_freeze)
        (; ν_air, D_vapor, K_therm) = aps

        L = TD.latent_heat_fusion(tps, T)

        n0 = get_n0(pdf, q_sno, ρ)
        v0 = get_v0(vel, ρ)
        (; r0) = mass
        (; χv, ve, Δv) = vel

        a_vent = vent.a
        b_vent = vent.b

        λ = lambda(pdf, mass, q_sno, ρ)

        snow_melt_rate =
            4 * FT(π) * n0 / ρ * K_therm / L * (T - T_freeze) / λ^FT(2) * (
                a_vent +
                b_vent * (ν_air / D_vapor)^FT(1 / 3) /
                (r0 * λ)^((ve + Δv) / FT(2)) *
                (FT(2) * v0 * χv / ν_air / λ)^FT(1 / 2) *
                CO.Γ((ve + Δv + FT(5)) / FT(2))
            )
    end
    return snow_melt_rate
end

end #module Microphysics1M.jl
