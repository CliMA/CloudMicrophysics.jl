"""
    Microphysics1M

One-moment bulk microphysics scheme, which includes:
  - terminal velocity of precipitation
  - autoconversion of cloud liquid water into rain and of cloud ice into snow
  - accretion due to collisions between categories of condensed species
  - evaporation and sublimation of hydrometeors
  - melting of snow into rain

# Mathematical Approach

Particle size distributions are assumed to follow the Marshall-Palmer exponential
distribution. Microphysical process rates are derived by integrating particle-scale
physics (terminal velocities, collision kernels, mass transfer) over these assumed
distributions.

# Common argument conventions

All option-dispatched process functions share the signature
`process(opt, mp, tps, micro, thermo)`, where:

- `opt`: process option type (dispatches which parameterization to use)
- `mp`: `Microphysics1MParams` — unified parameter container
- `tps`: thermodynamics parameters (from Thermodynamics.jl)
- `micro`: `NamedTuple` of specific humidities (kg/kg):
  - `q_tot` — total water specific humidity
  - `q_lcl` — cloud liquid water specific humidity
  - `q_icl` — cloud ice specific humidity
  - `q_rai` — rain specific humidity
  - `q_sno` — snow specific humidity
- `thermo`: `NamedTuple` of thermodynamic state:
  - `ρ` — air density (kg/m³)
  - `T` — temperature (K)

Not all processes use every field; each extracts only the fields it needs.
When a process is disabled (option = `nothing`), it returns zero.

Based on:
  - Kessler (1995) - https://doi.org/10.1016/0169-8095(94)00090-Z
  - Grabowski (1998) - https://doi.org/10.1175/1520-0469(1998)055<3283:TCRMOL>2.0.CO;2
  - Kaul et al. (2015) - https://doi.org/10.1175/MWR-D-14-00319.1
"""
module Microphysics1M

import SpecialFunctions as SF

import ..ThermodynamicsInterface as TDI
import ..Common as CO
import ..Parameters as CMP
import ..Utilities as UT


export terminal_velocity,
    conv_q_lcl_to_q_rai,
    conv_q_icl_to_q_sno,
    accretion,
    accretion_rain_sink,
    accretion_snow_rain,
    conv_q_rai_to_q_vap,
    conv_q_sno_to_q_vap,
    conv_q_icl_to_q_lcl,
    conv_q_sno_to_q_rai,
    lambda_inverse

abstract type AbstractSnowShape end
struct Oblate <: AbstractSnowShape end
struct Prolate <: AbstractSnowShape end


"""
    get_n0(pdf::ParticlePDFSnow, q_sno, ρ)
    get_n0(pdf::ParticlePDFIceRain, args...)

Returns the intercept parameter of the assumed Marshall-Palmer distribution
(Marshall and Palmer, 1948).

# Arguments
- `pdf`: size distribution parameters (contains `ν`, `μ` for snow; `n0` for rain/ice)
- `q_sno`: snow specific content (snow only)
- `ρ`: air density (snow only)
"""
@inline get_n0((; ν, μ)::CMP.ParticlePDFSnow{FT}, q_sno::FT, ρ::FT) where {FT} =
    q_sno > UT.ϵ_numerics(FT) ? μ * (ρ * q_sno)^ν : zero(FT)
@inline get_n0((; n0)::CMP.ParticlePDFIceRain{FT}, args...) where {FT} = n0

"""
    get_v0(vel::Blk1MVelTypeRain, ρ)
    get_v0(vel::Blk1MVelTypeSnow, args...)

Returns the proportionality coefficient in terminal velocity(r/r0).

Guards against unphysical density ratios (ρ > ρw) that would cause sqrt of negative.

# Arguments
- `vel`: terminal velocity parameters (contains `C_drag`, `ρw`, `grav`, `r0`, `gamma_term` for rain; `v0`, `gamma_term` for snow)
- `ρ`: air density (rain only)
"""
@inline function get_v0((; C_drag, ρw, grav, r0)::CMP.Blk1MVelTypeRain{FT}, ρ::FT) where {FT}
    # Guard against ρ > ρw (unphysical but could occur from numerical errors)
    density_factor = max(ρw / ρ - 1, zero(FT))
    return sqrt(FT(8 / 3) / C_drag * density_factor * grav * r0)
end
@inline get_v0((; v0)::CMP.Blk1MVelTypeSnow{FT}, args...) where {FT} = v0

"""
    lambda_inverse(pdf, mass::ParticleMass, q, ρ)

Returns the inverse of the rate parameter (λ⁻¹) of the assumed Marshall-Palmer
size distribution of particles (rain drops, ice crystals, snow crystals).

The rate parameter λ is related to the mean particle size: larger λ means smaller
average particles. The value is clipped at `r0 * 1e-5` to prevent numerical issues.

# Arguments
- `pdf`: size distribution parameters (ParticlePDFIceRain or ParticlePDFSnow)
- `mass`: mass(radius) parameters (contains `r0`, `m0`, `me`, `Δm`, `χm`, `gamma_coeff`)
- `q`: specific content of rain, cloud ice, or snow [kg/kg]
- `ρ`: air density [kg/m³]

# Returns
- `λ⁻¹`: inverse rate parameter [m]
"""
@inline function lambda_inverse(
    #(; pdf, mass)::Union{CMP.Snow{FT}, CMP.Rain{FT}, CMP.CloudIce{FT}},
    pdf::Union{CMP.ParticlePDFIceRain{FT}, CMP.ParticlePDFSnow{FT}},
    mass::CMP.ParticleMass{FT},
    q::FT,
    ρ::FT,
) where {FT}
    # size distribution
    n0::FT = get_n0(pdf, q, ρ)
    # mass(size)
    (; r0, m0, me, Δm, χm, gamma_coeff) = mass

    λ_inv = FT(0)
    if q > UT.ϵ_numerics(FT) && ρ > UT.ϵ_numerics(FT)
        # Note: Julia compiles x^y to exp(y * log(x))
        # gamma_coeff is pre-computed in ParticleMass constructor for GPU performance
        λ_inv = (ρ * q * r0^(me + Δm) / (χm * m0 * n0 * gamma_coeff))^(1 / (me + Δm + 1))
    end
    return max(r0 * FT(1e-5), λ_inv)
end

"""
    aspect_ratio_coeffs(snow_shape::Oblate, mass::ParticleMass, area::ParticleArea, ρᵢ)
    aspect_ratio_coeffs(snow_shape::Prolate, mass::ParticleMass, area::ParticleArea, ρᵢ)

Returns coefficients of the implied power law relationship between aspect ratio
and particle diameter φ(D) = φ₀ D^α.
Also returns the coefficient κ for the aspect ratio in Chen et al. (2022)
terminal velocity parameterization (κ=1/3 for oblate, κ=-1/6 for prolate).

# Arguments
- `snow_shape`: assumed snow particle shape (Oblate or Prolate)
- `mass`: mass(radius) parameters (contains `r0`, `m0`, `me`, `Δm`, `χm`, `gamma_coeff`)
- `area`: area(radius) parameters (contains `a0`, `ae`, `Δa`, `χa`)
- `ρᵢ`: particle density
"""
@inline function aspect_ratio_coeffs(
    snow_shape::Oblate,
    (; r0, m0, me, Δm, χm)::CMP.ParticleMass{FT},
    (; a0, ae, Δa, χa)::CMP.ParticleArea{FT},
    ρᵢ::FT,
) where {FT}
    # ϕ(r) = 3 * sqrt(FT(π)) * mᵢ(r) / (4 * ρᵢ * aᵢ(r)^(3/2)
    α = me + Δm - 3 / 2 * (ae + Δa)
    ϕ₀ = 3 * sqrt(FT(π)) / 4 / ρᵢ * χm * m0 / (χa * a0)^FT(3 / 2) / (2 * r0)^α
    κ = FT(1 / 3)
    return (; ϕ₀, α, κ)
end

@inline function aspect_ratio_coeffs(
    snow_shape::Prolate,
    (; r0, m0, me, Δm, χm)::CMP.ParticleMass{FT},
    (; a0, ae, Δa, χa)::CMP.ParticleArea{FT},
    ρᵢ::FT,
) where {FT}
    # ϕ(r) = 16 * ρᵢ^2 * aᵢ(r)^3 / (9 * π * mᵢ(r)^2)
    α = 3 * (ae + Δa) - 2 * (me + Δm)
    ϕ₀ = 16 * ρᵢ^2 / 9 / FT(π) * (χa * a0)^3 / (χm * m0)^2 / (2 * r0)^α
    κ = FT(-1 / 6)
    return (; ϕ₀, α, κ)
end

"""
    terminal_velocity(precip::Rain, vel::Blk1MVelTypeRain, ρ, q)
    terminal_velocity(precip::Snow, vel::Blk1MVelTypeSnow, ρ, q)
    terminal_velocity(precip::Rain, vel::Chen2022VelTypeRain, ρ, q)
    terminal_velocity(precip::Snow, vel::Chen2022VelTypeLargeIce, ρ, q)
    terminal_velocity(precip::Snow, vel::Chen2022VelTypeLargeIce, ρ, q, snow_shape)

Returns the mass-weighted average terminal velocity assuming a Marshall-Palmer
distribution of particles (Ogura and Takahashi, 1971).

The mass-weighted average velocity is computed by integrating the product of
particle mass, terminal velocity, and size distribution, then dividing by the
total mass. This represents the sedimentation velocity of the bulk hydrometeor field.

Fall velocity of individual particles is parameterized:
  - using empirical power-law relations for `Blk1MVelType`
  - following Chen et al. (2022), https://doi.org/10.1016/j.atmosres.2022.106171, for `Chen2022VelType`

# Arguments
- `precip`: precipitation parameters (Rain or Snow, contains `pdf`, `mass`, and for snow: `area`, `ρᵢ`, `aspr`)
- `vel`: terminal velocity parameterization parameters
- `ρ`: air density [kg/m³]
- `q`: rain or snow specific content [kg/kg]
- `snow_shape`: (optional) assumed snow shape (Oblate or Prolate)

# Returns
- Mass-weighted terminal velocity [m/s]
"""
@inline function terminal_velocity(
    (; pdf, mass)::Union{CMP.Rain, CMP.Snow},
    vel::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    ρ::FT,
    q::FT,
) where {FT}
    if q > UT.ϵ_numerics(FT)
        # terminal_velocity(size)
        (; χv, ve, Δv, gamma_term) = vel
        v0 = get_v0(vel, ρ)
        # mass(size)
        (; r0, me, Δm, χm, gamma_coeff) = mass
        # size distribution
        λ_inv = lambda_inverse(pdf, mass, q, ρ)

        # gamma_term = SF.gamma(me + ve + Δm + Δv + 1) (pre-computed in vel)
        # gamma_coeff = SF.gamma(me + Δm + 1) (pre-computed in mass)
        return χv * v0 * (λ_inv / r0)^(ve + Δv) * gamma_term / gamma_coeff
    else
        return FT(0)
    end
end

@inline function terminal_velocity(
    (; pdf, mass)::CMP.Rain,
    vel::CMP.Chen2022VelTypeRain{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > UT.ϵ_numerics(FT)
        # coefficients from Table B1 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(vel, ρₐ)
        # size distribution parameter
        λ_inv_radius::FT = lambda_inverse(pdf, mass, q, ρₐ)
        λ_inv_diameter = 2 * λ_inv_radius
        # eq 20 from Chen et al 2022 (loop unrolled for GPU performance)
        fall_w =
            CO.Chen2022_exponential_pdf(aiu[1], bi[1], ciu[1], λ_inv_diameter, 3) +
            CO.Chen2022_exponential_pdf(aiu[2], bi[2], ciu[2], λ_inv_diameter, 3) +
            CO.Chen2022_exponential_pdf(aiu[3], bi[3], ciu[3], λ_inv_diameter, 3)
        # It should be ϕ^κ * fall_w, but for rain drops ϕ = 1 and κ = 0
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end

@inline function terminal_velocity(
    (; pdf, mass, area, ρᵢ, aspr)::CMP.Snow,
    vel::CMP.Chen2022VelTypeLargeIce{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    # We assume the B4 table coeffs for snow and B2 table coeffs for cloud ice.
    # Instead we should do partial integrals
    # from D=125um to D=625um using B2 and D=625um to inf using B4.
    if q > UT.ϵ_numerics(FT)
        # coefficients from Table B4 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(vel, ρₐ, ρᵢ)
        # size distribution parameter
        λ_inv_radius::FT = lambda_inverse(pdf, mass, q, ρₐ)
        λ_inv_diameter = 2 * λ_inv_radius

        # As a next step, we could keep ϕ(r) under the integrals
        # assume oblate shape and aspect ratio
        (; ϕ, κ) = aspr

        # eq 20 from Chen 2022 (loop unrolled for GPU performance)
        fall_w =
            ϕ^κ * CO.Chen2022_exponential_pdf(aiu[1], bi[1], ciu[1], λ_inv_diameter, 3) +
            ϕ^κ * CO.Chen2022_exponential_pdf(aiu[2], bi[2], ciu[2], λ_inv_diameter, 3)
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end

@inline function terminal_velocity(
    (; pdf, mass, area, ρᵢ, gamma_aspect_oblate, gamma_aspect_prolate)::CMP.Snow,
    vel::CMP.Chen2022VelTypeLargeIce{FT},
    ρₐ::FT,
    q::FT,
    snow_shape::AbstractSnowShape,
) where {FT}
    fall_w = FT(0)
    # see comments above about B2 vs B4 coefficients
    if q > UT.ϵ_numerics(FT)
        # coefficients from Table B4 from Chen et. al. 2022
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(vel, ρₐ, ρᵢ)
        # size distribution parameter
        λ_inv_radius::FT = lambda_inverse(pdf, mass, q, ρₐ)
        λ_inv_diameter = 2 * λ_inv_radius
        # Compute the mass weighted average aspect ratio ϕ_av
        # As a next step, we could keep ϕ(r) under the integrals
        (ϕ₀, α, κ) = aspect_ratio_coeffs(snow_shape, mass, area, ρᵢ)
        # Use pre-computed gamma_aspect from Snow struct
        gamma_aspect = snow_shape isa Oblate ? gamma_aspect_oblate : gamma_aspect_prolate
        ϕ_av = ϕ₀ * λ_inv_radius^α * gamma_aspect
        # eq 20 from Chen 2022 (loop unrolled for GPU performance)
        fall_w =
            ϕ_av^κ * CO.Chen2022_exponential_pdf(aiu[1], bi[1], ciu[1], λ_inv_diameter, 3) +
            ϕ_av^κ * CO.Chen2022_exponential_pdf(aiu[2], bi[2], ciu[2], λ_inv_diameter, 3)
        fall_w = max(FT(0), fall_w)
    end
    return fall_w
end

"""
    conv_q_lcl_to_q_rai(::Nothing, mp, tps, micro, thermo)
    conv_q_lcl_to_q_rai(::Kessler1M, mp, tps, micro, thermo)
    conv_q_lcl_to_q_rai(::PrescribedNd, mp, tps, micro, thermo)

Returns the rain tendency due to autoconversion of cloud liquid, dispatching on
the option stored in `Microphysics1MOptions`.

**Nothing**: returns zero (autoconversion disabled).

**Kessler1M**: Kessler (1995) 1-moment threshold autoconversion
(smooth logistic transition), https://doi.org/10.1016/0169-8095(94)00090-Z.

**PrescribedNd**: Variable-timescale autoconversion following Azimi (2023),
using the prescribed cloud droplet number concentration.

# Arguments
- `option`: `nothing`, `Kessler1M(...)`, or `PrescribedNd(...)`
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters (unused, kept for uniform interface)
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)` (unused for 1M, kept for uniform interface)

# Returns
- Rain autoconversion rate [kg/kg/s]
"""
@inline conv_q_lcl_to_q_rai(::Nothing, mp, tps, micro, thermo) = zero(micro.q_lcl)

@inline function conv_q_lcl_to_q_rai(opt::CMP.Kessler1M, mp, tps, micro, thermo)
    q_lcl = micro.q_lcl
    (; τ, q_threshold, k) = opt.acnv1M
    return CO.logistic_function_integral(q_lcl, q_threshold, k) / τ
end

@inline function conv_q_lcl_to_q_rai(opt::CMP.PrescribedNd, mp, tps, micro, thermo)
    q_lcl = micro.q_lcl
    (; τ, α, Nc) = opt.autoconv
    return max(0, q_lcl) / (τ * (Nc / 100_000_000)^α)
end

"""
    conv_q_icl_to_q_sno(::Nothing, mp, tps, micro, thermo)
    conv_q_icl_to_q_sno(opt::NoSupersaturation, mp, tps, micro, thermo)
    conv_q_icl_to_q_sno(opt::WithSupersaturation, mp, tps, micro, thermo)

Returns the snow tendency due to autoconversion from cloud ice.

**Nothing**: returns zero (snow autoconversion disabled).

**NoSupersaturation**: Simplified Kessler-type threshold autoconversion,
for use in simulations without supersaturation (e.g., with saturation adjustment).

**WithSupersaturation**: Supersaturation-dependent autoconversion following
Harrington et al. (1995) and Kaul et al. (2015).

# Arguments
- `opt`: `nothing`, `NoSupersaturation(...)`, or `WithSupersaturation(...)`
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`
"""
@inline conv_q_icl_to_q_sno(::Nothing, mp, tps, micro, thermo) = zero(micro.q_icl)

@inline function conv_q_icl_to_q_sno(opt::CMP.NoSupersaturation, mp, tps, micro, thermo)
    (; τ, q_threshold, k) = opt.acnv1M
    q_icl = micro.q_icl
    return CO.logistic_function_integral(q_icl, q_threshold, k) / τ
end

@inline function conv_q_icl_to_q_sno(opt::CMP.WithSupersaturation, mp, tps, micro, thermo)
    (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro
    (; ρ, T) = thermo
    r_ice_snow = opt.r_ice_snow
    (; pdf, mass) = mp.cloud.ice
    aps = mp.air_properties
    FT = eltype(ρ)
    acnv_rate = FT(0)
    S = TDI.supersaturation_over_ice(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)

    # Only allow ice autoconversion below freezing with positive supersaturation
    if (q_icl > UT.ϵ_numerics(FT) && S > FT(0) && T < TDI.T_freeze(tps))
        (; me, Δm) = mass
        G = CO.G_func_ice(aps, tps, T)
        n0 = get_n0(pdf)
        λ_inv = lambda_inverse(pdf, mass, q_icl, ρ)

        acnv_rate =
            4 * FT(π) * S * G * n0 / ρ *
            exp(-r_ice_snow / λ_inv) *
            (r_ice_snow^2 / (me + Δm) + (r_ice_snow / λ_inv + 1) * λ_inv^2)
    end
    return acnv_rate
end

"""
    warm_accretion_melt_factor(tps, T)

Ratio of sensible heat from warm collected liquid to latent heat for melting:
`α = cv_l / L_f × (T - T_freeze)`, returning 0 when `T ≤ T_freeze`.

Used by `accretion(::CloudLiquidSnowAccretion, ...)` and
`accretion_snow_rain(::RainSnowAccretion, ...)` to compute the thermal melt
contribution of warm liquid on snow.
"""
@inline function warm_accretion_melt_factor(tps, T)
    L_f = TDI.Lf(tps, T)
    cv_l = TDI.cv_l(tps)
    T_freeze = TDI.T_freeze(tps)
    ΔT = T - T_freeze
    is_cold = (T <= T_freeze)
    return ifelse(is_cold, zero(T), cv_l / L_f * ΔT)
end

"""
    accretion(cloud, precip, vel, E, q_clo, q_pre, ρ)

Returns the source of precipitating water (rain or snow) due to collisions
with cloud water (liquid or ice).

Internal low-level kernel. Prefer the option-dispatched API.

# Arguments
- `cloud`: type for cloud water or cloud ice
- `precip`: type for rain or snow
- `vel`: terminal velocity parameters
- `E`: collision efficiency
- `q_clo`: cloud liquid water or cloud ice specific content
- `q_pre`: rain or snow specific content
- `ρ`: air density
"""
@inline function accretion(
    cloud::CMP.CloudCondensateType,
    precip::CMP.PrecipitationType,
    vel::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    E::FT,
    q_clo::FT,
    q_pre::FT,
    ρ::FT,
) where {FT}

    accr_rate = FT(0)
    if (q_clo > UT.ϵ_numerics(FT) && q_pre > UT.ϵ_numerics(FT))

        n0::FT = get_n0(precip.pdf, q_pre, ρ)
        v0::FT = get_v0(vel, ρ)

        (; r0) = precip.mass
        (; χv, ve, Δv, gamma_accr) = vel
        (; a0, ae, χa, Δa) = precip.area

        λ_inv = lambda_inverse(precip.pdf, precip.mass, q_pre, ρ)

        # gamma_accr = SF.gamma(ae + ve + Δa + Δv + 1) (pre-computed in vel)
        accr_rate =
            q_clo * E * n0 * a0 * v0 * χa * χv * λ_inv *
            gamma_accr / (r0 / λ_inv)^(ae + ve + Δa + Δv)
    end
    return accr_rate
end

# accretion_rain_sink(rain, ice, vel, E, q_icl, q_rai, ρ)
#
# Returns the sink of rain water (partial source of snow) due to collisions
# with cloud ice.
@inline function accretion_rain_sink(
    rain::CMP.Rain,
    ice::CMP.CloudIce,
    vel::CMP.Blk1MVelTypeRain{FT},
    E::FT,
    q_icl::FT,
    q_rai::FT,
    ρ::FT,
) where {FT}
    accr_rate = FT(0)
    if (q_icl > UT.ϵ_numerics(FT) && q_rai > UT.ϵ_numerics(FT))

        n0_ice = get_n0(ice.pdf)
        λ_ice_inv = lambda_inverse(ice.pdf, ice.mass, q_icl, ρ)

        n0 = get_n0(rain.pdf, q_rai, ρ)
        v0 = get_v0(vel, ρ)
        (; r0, m0, me, Δm, χm) = rain.mass
        (; χv, ve, Δv, gamma_accr_rain_sink) = vel
        (; a0, ae, χa, Δa) = rain.area

        λ_inv = lambda_inverse(rain.pdf, rain.mass, q_rai, ρ)

        # gamma_accr_rain_sink = SF.gamma(me + ae + ve + Δm + Δa + Δv + 1) (pre-computed in vel)
        accr_rate =
            E / ρ * n0 * n0_ice * m0 * a0 * v0 * χm * χa * χv * λ_ice_inv * λ_inv *
            gamma_accr_rain_sink /
            (r0 / λ_inv)^FT(me + ae + ve + Δm + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_snow_rain(type_i::PrecipitationType, type_j::PrecipitationType, blk1mveltype_ti, blk1mveltype_tj, ce, q_i, q_j, ρ)

Returns the accretion rate when rain and snow collide.
Collisions result in snow for T < T_freeze and rain for T > T_freeze.

Uses geometric collision kernel assumption: a(r_i, r_j) = π(r_i + r_j)², with
a velocity dispersion correction that assumes that fall velocity standard 
deviations are proportional to the mean fall velocities, with coefficient
`ce.coeff_disp`.

!!! note "Internal use only"
    This low-level positional-argument kernel is kept for internal dispatch.
    Prefer `accretion_snow_rain(::RainSnowAccretion, mp, tps, micro, thermo)`,
    which calls both the cold and warm arms and returns `(; S_rai_sno, S_sno_rai, S_melt)`.

# Arguments
- `type_i`: snow (T < T_freeze) or rain (T > T_freeze)
- `type_j`: rain (T < T_freeze) or snow (T > T_freeze)  
- `blk1mveltype_ti`, `blk1mveltype_tj`: 1M terminal velocity parameters
- `ce`: collision efficiency parameters (contains `e_rai_sno`, `coeff_disp`)
- `q_i`, `q_j`: specific contents of snow or rain [kg/kg]
- `ρ`: air density [kg/m³]
"""
@inline function accretion_snow_rain(
    type_i::CMP.PrecipitationType,
    type_j::CMP.PrecipitationType,
    blk1mveltype_ti,
    blk1mveltype_tj,
    E_ij::FT,
    coeff_disp::FT,
    q_i::FT,
    q_j::FT,
    ρ::FT,
) where {FT}

    accr_rate = FT(0)
    if (q_i > UT.ϵ_numerics(FT) && q_j > UT.ϵ_numerics(FT))

        n0_i = get_n0(type_i.pdf, q_i, ρ)
        n0_j = get_n0(type_j.pdf, q_j, ρ)

        (; r0, m0, me, Δm, χm, gamma_coeff) = type_j.mass
        δ = me + Δm

        λ_i_inv = lambda_inverse(type_i.pdf, type_i.mass, q_i, ρ)
        λ_j_inv = lambda_inverse(type_j.pdf, type_j.mass, q_j, ρ)

        v_ti = terminal_velocity(type_i, blk1mveltype_ti, ρ, q_i)
        v_tj = terminal_velocity(type_j, blk1mveltype_tj, ρ, q_j)

        # Add simple parameterization for velocity dispersion, assuming that fall velocity 
        # standard deviations are proportional to the mean fall velocities, with coefficient 
        # coeff_disp
        Δv_eff = sqrt((v_ti - v_tj)^2 + coeff_disp * (v_ti^2 + v_tj^2))

        # We use the recurrence relation Γ(x+1) = xΓ(x) to simplify gamma terms.
        # gamma_coeff = Γ(δ + 1) is pre-computed.
        accr_rate =
            FT(π) / ρ * n0_i * n0_j * m0 * χm * E_ij * Δv_eff * gamma_coeff /
            r0^δ * (
                2 * λ_i_inv^3 * λ_j_inv^(δ + 1) +
                2 * (δ + 1) * λ_i_inv^2 * λ_j_inv^(δ + 2) +
                (δ + 2) * (δ + 1) * λ_i_inv * λ_j_inv^(δ + 3)
            )
    end
    return accr_rate
end

"""
    accretion(::Nothing, mp, tps, micro, thermo)
    accretion(opt::CloudLiquidRainAccretion, mp, tps, micro, thermo)
    accretion(opt::CloudLiquidSnowAccretion, mp, tps, micro, thermo)
    accretion(opt::CloudIceRainAccretion, mp, tps, micro, thermo)
    accretion(opt::CloudIceSnowAccretion, mp, tps, micro, thermo)
    accretion_snow_rain(::Nothing, mp, tps, micro, thermo)
    accretion_snow_rain(opt::RainSnowAccretion, mp, tps, micro, thermo)

Option-dispatched accretion wrappers. All extract parameters from `mp` and
delegate to the corresponding low-level Marshall-Palmer kernels.
`nothing` variants return zero without computing anything.

# Arguments
- `opt`: accretion option type (carries collision efficiency) or `nothing`
- `mp`: `Microphysics1MParams`
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`
"""
# Nothing dispatch: scalar zero for simple accretion, NamedTuple for processes
# that return NamedTuple shapes. Since `nothing` is a singleton, we can't dispatch
# differently on it by process. BMT calls accretion() for scalar-result processes
# and directly destructures for NamedTuple processes. We provide the scalar zero
# here; NamedTuple-returning processes are handled in BMT by checking `nothing` directly.
@inline accretion(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)

@inline function accretion(opt::CMP.CloudLiquidRainAccretion, mp, tps, micro, thermo)
    q_lcl = micro.q_lcl
    q_rai = micro.q_rai
    ρ = thermo.ρ
    return accretion(mp.cloud.liquid, mp.precip.rain, mp.terminal_velocity.rain, opt.e, q_lcl, q_rai, ρ)
end

@inline function accretion(opt::CMP.CloudLiquidSnowAccretion, mp, tps, micro, thermo)
    q_lcl = micro.q_lcl
    q_sno = micro.q_sno
    ρ = thermo.ρ
    T = thermo.T
    S = accretion(mp.cloud.liquid, mp.precip.snow, mp.terminal_velocity.snow, opt.e, q_lcl, q_sno, ρ)
    α = warm_accretion_melt_factor(tps, T)
    return (; S_accr = S, S_melt = α * S)
end

@inline function accretion(opt::CMP.CloudIceRainAccretion, mp, tps, micro, thermo)
    q_icl = micro.q_icl
    q_rai = micro.q_rai
    ρ = thermo.ρ
    return accretion(mp.cloud.ice, mp.precip.rain, mp.terminal_velocity.rain, opt.e, q_icl, q_rai, ρ)
end

@inline function accretion(opt::CMP.CloudIceSnowAccretion, mp, tps, micro, thermo)
    q_icl = micro.q_icl
    q_sno = micro.q_sno
    ρ = thermo.ρ
    return accretion(mp.cloud.ice, mp.precip.snow, mp.terminal_velocity.snow, opt.e, q_icl, q_sno, ρ)
end

@inline accretion_snow_rain(::Nothing, mp, tps, micro, thermo) =
    (; S_rai_sno = zero(thermo.T), S_sno_rai = zero(thermo.T), S_melt = zero(thermo.T))

@inline function accretion_snow_rain(opt::CMP.RainSnowAccretion, mp, tps, micro, thermo)
    q_rai = micro.q_rai
    q_sno = micro.q_sno
    ρ = thermo.ρ
    T = thermo.T
    vel = mp.terminal_velocity
    sno = mp.precip.snow
    rai = mp.precip.rain
    # cold arm: snow is collector, rain freezes → snow
    S_rai_sno = accretion_snow_rain(sno, rai, vel.snow, vel.rain, opt.e, opt.coeff_disp, q_sno, q_rai, ρ)
    # warm arm: rain is collector, snow melts → rain
    S_sno_rai = accretion_snow_rain(rai, sno, vel.rain, vel.snow, opt.e, opt.coeff_disp, q_rai, q_sno, ρ)
    α = warm_accretion_melt_factor(tps, T)
    return (; S_rai_sno, S_sno_rai, S_melt = α * S_rai_sno)
end

# Rain sink arm of cloud ice + rain accretion
@inline accretion_rain_sink(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)

@inline function accretion_rain_sink(opt::CMP.CloudIceRainAccretion, mp, tps, micro, thermo)
    q_icl = micro.q_icl
    q_rai = micro.q_rai
    ρ = thermo.ρ
    return accretion_rain_sink(mp.precip.rain, mp.cloud.ice, mp.terminal_velocity.rain, opt.e, q_icl, q_rai, ρ)
end

"""
Returns the tendency due to rain evaporation.
Ventilation factor parameterization follows Seifert and Beheng (2006).

Returns the tendency due to rain evaporation.
Ventilation factor parameterization follows Seifert and Beheng (2006).

Only evaporation is considered (sub-saturated over liquid); result is clamped ≤ 0.

# Arguments
- `opt`: `RainEvaporation()` or `nothing`
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`
"""
@inline conv_q_rai_to_q_vap(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)
@inline conv_q_rai_to_q_vap(::Nothing, mp, tps, micro, thermo, p_vs_liq) = zero(thermo.T)

# Convenience method: compute the saturation vapor pressure over liquid, then
# delegate to the variant below. Callers that already hold `p_vs_liq` (e.g. the
# bulk tendency driver, which reuses it across several processes) should call
# the 6-argument method directly to avoid recomputing it.
@inline conv_q_rai_to_q_vap(opt::CMP.RainEvaporation, mp, tps, micro, thermo) =
    conv_q_rai_to_q_vap(
        opt, mp, tps, micro, thermo,
        TDI.saturation_vapor_pressure_over_liquid(tps, thermo.T),
    )

@inline function conv_q_rai_to_q_vap(::CMP.RainEvaporation, mp, tps, micro, thermo, p_vs_liq)
    (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro
    (; ρ, T) = thermo
    (; pdf, mass, vent) = mp.precip.rain
    vel = mp.terminal_velocity.rain
    aps = mp.air_properties
    FT = eltype(ρ)
    evap_rate = FT(0)

    if q_rai > UT.ϵ_numerics(FT)
        S = TDI.supersaturation(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T, p_vs_liq)

        if S < FT(0)
            (; ν_air, D_vapor) = aps
            G = CO.G_func_liquid(aps, tps, T, p_vs_liq)
            n0 = get_n0(pdf, q_rai, ρ)
            v0 = get_v0(vel, ρ)
            (; χv, ve, Δv, gamma_vent) = vel
            (; r0) = mass
            a_vent = vent.a
            b_vent = vent.b

            λ_inv = lambda_inverse(pdf, mass, q_rai, ρ)
            Sc = ν_air / max(D_vapor, UT.ϵ_numerics(FT))

            evap_rate =
                4 * FT(π) * n0 / ρ * S * G * λ_inv^2 *
                (
                    a_vent +
                    b_vent * cbrt(Sc) /
                    (r0 / λ_inv)^((ve + Δv) / 2) *
                    sqrt(2 * v0 * χv / ν_air * λ_inv) *
                    gamma_vent
                )
        end
    end
    return min(0, evap_rate)
end

"""
    conv_q_sno_to_q_vap(::Nothing, mp, tps, micro, thermo)
    conv_q_sno_to_q_vap(::SublimationOnly, mp, tps, micro, thermo)
    conv_q_sno_to_q_vap(::DepositionAndSublimation, mp, tps, micro, thermo)

Returns the tendency due to snow sublimation or sublimation+deposition.
Ventilation factor parameterization follows Seifert and Beheng (2006).

# Arguments
- `opt`: `nothing`, `SublimationOnly()`, or `DepositionAndSublimation()`
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`
"""
@inline conv_q_sno_to_q_vap(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)
@inline conv_q_sno_to_q_vap(::Nothing, mp, tps, micro, thermo, p_vs_ice) = zero(thermo.T)

@inline function conv_q_sno_to_q_vap(::CMP.SublimationOnly, mp, tps, micro, thermo, p_vs_ice)
    return min(0, _snow_subl_dep_rate(mp, tps, micro, thermo, p_vs_ice))
end
@inline conv_q_sno_to_q_vap(opt::CMP.SublimationOnly, mp, tps, micro, thermo) =
    conv_q_sno_to_q_vap(
        opt, mp, tps, micro, thermo,
        TDI.saturation_vapor_pressure_over_ice(tps, thermo.T),
    )

@inline function conv_q_sno_to_q_vap(
    ::CMP.DepositionAndSublimation, mp, tps, micro, thermo, p_vs_ice,
)
    return _snow_subl_dep_rate(mp, tps, micro, thermo, p_vs_ice)
end
@inline conv_q_sno_to_q_vap(opt::CMP.DepositionAndSublimation, mp, tps, micro, thermo) =
    conv_q_sno_to_q_vap(
        opt, mp, tps, micro, thermo,
        TDI.saturation_vapor_pressure_over_ice(tps, thermo.T),
    )

"""Internal helper: snow sublimation/deposition physics kernel."""
@inline _snow_subl_dep_rate(mp, tps, micro, thermo) = _snow_subl_dep_rate(
    mp, tps, micro, thermo, TDI.saturation_vapor_pressure_over_ice(tps, thermo.T),
)
@inline function _snow_subl_dep_rate(mp, tps, micro, thermo, p_vs_ice)
    (; q_tot, q_lcl, q_icl, q_rai, q_sno) = micro
    (; ρ, T) = thermo
    (; pdf, mass, vent) = mp.precip.snow
    vel = mp.terminal_velocity.snow
    aps = mp.air_properties
    FT = eltype(ρ)
    subl_rate = FT(0)

    if q_sno > UT.ϵ_numerics(FT)
        (; ν_air, D_vapor) = aps
        S = TDI.supersaturation(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T, p_vs_ice)
        G = CO.G_func_ice(aps, tps, T, p_vs_ice)
        n0 = get_n0(pdf, q_sno, ρ)
        v0 = get_v0(vel, ρ)
        (; r0) = mass
        (; χv, ve, Δv, gamma_vent) = vel
        a_vent = vent.a
        b_vent = vent.b

        λ_inv = lambda_inverse(pdf, mass, q_sno, ρ)
        Sc = ν_air / max(D_vapor, UT.ϵ_numerics(FT))

        subl_rate =
            4 * FT(π) * n0 / ρ * S * G * λ_inv^2 *
            (
                a_vent +
                b_vent * cbrt(Sc) /
                (r0 / λ_inv)^((ve + Δv) / 2) *
                sqrt(2 * v0 * χv / ν_air * λ_inv) *
                gamma_vent
            )
    end
    return subl_rate
end


"""
    conv_q_icl_to_q_lcl(::Nothing, mp, tps, micro, thermo)
    conv_q_icl_to_q_lcl(::CloudIceMelt, mp, tps, micro, thermo)

Returns the tendency due to cloud ice melt.

# Arguments
- `opt`: `CloudIceMelt()` or `nothing`
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`
"""
@inline conv_q_icl_to_q_lcl(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)

@inline function conv_q_icl_to_q_lcl(::CMP.CloudIceMelt, mp, tps, micro, thermo)
    q_icl = micro.q_icl
    (; ρ, T) = thermo
    (; pdf, mass) = mp.cloud.ice
    (; K_therm) = mp.air_properties
    FT = eltype(ρ)
    cloud_ice_melt_rate = FT(0)
    T_freeze = TDI.T_freeze(tps)

    if (q_icl > UT.ϵ_numerics(FT) && T > T_freeze)
        L = TDI.Lf(tps, T)
        (; n0) = pdf
        λ_inv = lambda_inverse(pdf, mass, q_icl, ρ)
        cloud_ice_melt_rate = 4 * FT(π) * n0 / ρ * K_therm / L * (T - T_freeze) * λ_inv^2
    end
    return cloud_ice_melt_rate
end

"""
    conv_q_sno_to_q_rai(::Nothing, mp, tps, micro, thermo)
    conv_q_sno_to_q_rai(::SnowMelt, mp, tps, micro, thermo)

Returns the tendency due to snow melt.

# Arguments
- `opt`: `SnowMelt()` or `nothing`
- `mp`: 1-moment microphysics parameters
- `tps`: thermodynamics parameters
- `micro`: microphysics state `(; q_tot, q_lcl, q_icl, q_rai, q_sno)`
- `thermo`: thermodynamic state `(; ρ, T)`
"""
@inline conv_q_sno_to_q_rai(::Nothing, mp, tps, micro, thermo) = zero(thermo.T)

@inline function conv_q_sno_to_q_rai(::CMP.SnowMelt, mp, tps, micro, thermo)
    q_sno = micro.q_sno
    (; ρ, T) = thermo
    (; pdf, mass, vent) = mp.precip.snow
    vel = mp.terminal_velocity.snow
    aps = mp.air_properties
    FT = eltype(ρ)
    snow_melt_rate = FT(0)
    T_freeze = TDI.T_freeze(tps)

    if (q_sno > UT.ϵ_numerics(FT) && T > T_freeze)
        (; ν_air, D_vapor, K_therm) = aps

        L = TDI.Lf(tps, T)

        n0 = get_n0(pdf, q_sno, ρ)
        v0 = get_v0(vel, ρ)
        (; r0) = mass
        (; χv, ve, Δv, gamma_vent) = vel

        a_vent = vent.a
        b_vent = vent.b

        λ_inv = lambda_inverse(pdf, mass, q_sno, ρ)

        # Schmidt number (guard against division by near-zero D_vapor)
        Sc = ν_air / max(D_vapor, UT.ϵ_numerics(FT))

        snow_melt_rate =
            4 * FT(π) * n0 / ρ * K_therm / L * (T - T_freeze) * λ_inv^2 *
            (
                a_vent +
                b_vent * cbrt(Sc) /
                (r0 / λ_inv)^((ve + Δv) / 2) *
                sqrt(2 * v0 * χv / ν_air * λ_inv) *
                gamma_vent
            )
    end
    return snow_melt_rate
end

end #module Microphysics1M.jl
