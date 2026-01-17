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


export terminal_velocity,
    conv_q_lcl_to_q_rai,
    conv_q_icl_to_q_sno,
    accretion,
    accretion_rain_sink,
    accretion_snow_rain,
    evaporation_sublimation,
    snow_melt,
    lambda_inverse

abstract type AbstractSnowShape end
struct Oblate <: AbstractSnowShape end
struct Prolate <: AbstractSnowShape end

"""
    Ec(type_1, type_2, ce)

Returns the collision efficiency (Ec) for two colliding species.

Collision efficiency represents the fraction of geometric collisions that result
in coalescence (0 ≤ Ec ≤ 1). A value of 1 means all collisions lead to merging,
while lower values account for particles bouncing apart.

# Arguments
- `type_1`, `type_2`: types of colliding species (CloudLiquid, CloudIce, Rain, Snow)
- `ce`: collision efficiency parameters struct
"""
@inline Ec(::CMP.CloudLiquid, ::CMP.Rain, (; e_lcl_rai)::CMP.CollisionEff) = e_lcl_rai
@inline Ec(::CMP.CloudLiquid, ::CMP.Snow, (; e_lcl_sno)::CMP.CollisionEff) = e_lcl_sno
@inline Ec(::CMP.CloudIce, ::CMP.Rain, (; e_icl_rai)::CMP.CollisionEff) = e_icl_rai
@inline Ec(::CMP.CloudIce, ::CMP.Snow, (; e_icl_sno)::CMP.CollisionEff) = e_icl_sno
@inline Ec(::CMP.Rain, ::CMP.Snow, (; e_rai_sno)::CMP.CollisionEff) = e_rai_sno
@inline Ec(::CMP.Snow, ::CMP.Rain, (; e_rai_sno)::CMP.CollisionEff) = e_rai_sno

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
    μ * (ρ * max(0, q_sno))^ν  # TODO - think about limiters
@inline get_n0((; n0)::CMP.ParticlePDFIceRain{FT}, args...) where {FT} = n0

"""
    get_v0(vel::Blk1MVelTypeRain, ρ)
    get_v0(vel::Blk1MVelTypeSnow, args...)

Returns the proportionality coefficient in terminal velocity(r/r0).

# Arguments
- `vel`: terminal velocity parameters (contains `C_drag`, `ρw`, `grav`, `r0`, `gamma_term` for rain; `v0`, `gamma_term` for snow)
- `ρ`: air density (rain only)
"""
@inline get_v0((; C_drag, ρw, grav, r0)::CMP.Blk1MVelTypeRain{FT}, ρ::FT) where {FT} =
    sqrt(FT(8 / 3) / C_drag * (ρw / ρ - 1) * grav * r0)
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
    if q > CO.ϵ_numerics(FT) && ρ > CO.ϵ_numerics(FT)
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
function aspect_ratio_coeffs(
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

function aspect_ratio_coeffs(
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
function terminal_velocity(
    (; pdf, mass)::Union{CMP.Rain{FT}, CMP.Snow{FT}},
    vel::Union{CMP.Blk1MVelTypeRain{FT}, CMP.Blk1MVelTypeSnow{FT}},
    ρ::FT,
    q::FT,
) where {FT}
    if q > CO.ϵ_numerics(FT)
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

function terminal_velocity(
    (; pdf, mass)::CMP.Rain{FT},
    vel::CMP.Chen2022VelTypeRain{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    if q > CO.ϵ_numerics(FT)
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

function terminal_velocity(
    (; pdf, mass, area, ρᵢ, aspr)::CMP.Snow{FT},
    vel::CMP.Chen2022VelTypeLargeIce{FT},
    ρₐ::FT,
    q::FT,
) where {FT}
    fall_w = FT(0)
    # We assume the B4 table coeffs for snow and B2 table coeffs for cloud ice.
    # Instead we should do partial integrals
    # from D=125um to D=625um using B2 and D=625um to inf using B4.
    if q > CO.ϵ_numerics(FT)
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

function terminal_velocity(
    (; pdf, mass, area, ρᵢ, gamma_aspect_oblate, gamma_aspect_prolate)::CMP.Snow{FT},
    vel::CMP.Chen2022VelTypeLargeIce{FT},
    ρₐ::FT,
    q::FT,
    snow_shape::AbstractSnowShape,
) where {FT}
    fall_w = FT(0)
    # see comments above about B2 vs B4 coefficients
    if q > CO.ϵ_numerics(FT)
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
    conv_q_lcl_to_q_rai(acnv::Acnv1M, q_lcl, smooth_transition)

Returns the rain tendency due to collisions between cloud droplets (autoconversion),
parameterized following Kessler (1995), https://doi.org/10.1016/0169-8095(94)00090-Z.

When `smooth_transition = false`, uses a step function at the threshold.
When `smooth_transition = true`, uses a logistic function to smooth the transition
over the threshold, avoiding discontinuities in the tendency.

# Arguments
- `acnv`: autoconversion parameters (contains `τ`, `q_threshold`, `k`)
- `q_lcl`: cloud liquid water specific content [kg/kg]
- `smooth_transition`: flag to switch on smoothing

# Returns
- Rain autoconversion rate [kg/kg/s]
"""
conv_q_lcl_to_q_rai(
    (; τ, q_threshold, k)::CMP.Acnv1M{FT},
    q_lcl::FT,
    smooth_transition::Bool = false,
) where {FT} =
    smooth_transition ?
    CO.logistic_function_integral(q_lcl, q_threshold, k) / τ :
    max(0, q_lcl - q_threshold) / τ

"""
    conv_q_icl_to_q_sno_no_supersat(acnv::Acnv1M, q_icl, smooth_transition)

Returns the snow tendency due to autoconversion from cloud ice.
This is a simplified version for use in simulations without supersaturation
(e.g., with saturation adjustment).

# Arguments
- `acnv`: autoconversion parameters (contains `τ`, `q_threshold`, `k`)
- `q_icl`: cloud ice specific content
- `smooth_transition`: flag to switch on smoothing
"""
conv_q_icl_to_q_sno_no_supersat(
    (; τ, q_threshold, k)::CMP.Acnv1M{FT},
    q_icl::FT,
    smooth_transition::Bool = false,
) where {FT} =
    smooth_transition ?
    CO.logistic_function_integral(q_icl, q_threshold, k) / τ :
    max(0, q_icl - q_threshold) / τ

"""
    conv_q_icl_to_q_sno(ice::CloudIce, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)

Returns the snow tendency due to autoconversion from ice.
Parameterized following:
  - Harrington et al. (1995), https://doi.org/10.1175/1520-0469(1995)052<4344:POICCP>2.0.CO;2  
  - Kaul et al. (2015), https://doi.org/10.1175/MWR-D-14-00319.1

# Arguments
- `ice`: ice parameters (contains `r_ice_snow`, `pdf`, `mass`)
- `aps`: air properties struct  
- `tps`: thermodynamics parameters struct
- `q_tot`: total water specific content
- `q_lcl`: cloud liquid water specific content
- `q_icl`: cloud ice specific content
- `q_rai`: rain specific content
- `q_sno`: snow specific content
- `ρ`: air density
- `T`: air temperature
"""
function conv_q_icl_to_q_sno(
    (; r_ice_snow, pdf, mass)::CMP.CloudIce{FT},
    aps::CMP.AirProperties{FT},
    tps::TDI.PS,
    q_tot::FT,
    q_lcl::FT,
    q_icl::FT,
    q_rai::FT,
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT}
    acnv_rate = FT(0)
    S = TDI.supersaturation_over_ice(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)

    if (q_icl > CO.ϵ_numerics(FT) && S > FT(0))
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
    accretion(cloud::CloudCondensateType, precip::PrecipitationType, vel, ce, q_clo, q_pre, ρ)

Returns the source of precipitating water (rain or snow) due to collisions
with cloud water (liquid or ice).

# Arguments
- `cloud`: type for cloud water or cloud ice
- `precip`: type for rain or snow
- `vel`: terminal velocity parameters (Blk1MVelTypeRain or Blk1MVelTypeSnow, contains `gamma_accr`)
- `ce`: collision efficiency parameters
- `q_clo`: cloud liquid water or cloud ice specific content
- `q_pre`: rain or snow specific content
- `ρ`: air density
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
    if (q_clo > CO.ϵ_numerics(FT) && q_pre > CO.ϵ_numerics(FT))

        n0::FT = get_n0(precip.pdf, q_pre, ρ)
        v0::FT = get_v0(vel, ρ)

        (; r0) = precip.mass
        (; χv, ve, Δv, gamma_accr) = vel
        (; a0, ae, χa, Δa) = precip.area

        λ_inv = lambda_inverse(precip.pdf, precip.mass, q_pre, ρ)
        E = Ec(cloud, precip, ce)

        # gamma_accr = SF.gamma(ae + ve + Δa + Δv + 1) (pre-computed in vel)
        accr_rate =
            q_clo * E * n0 * a0 * v0 * χa * χv * λ_inv *
            gamma_accr / (r0 / λ_inv)^(ae + ve + Δa + Δv)
    end
    return accr_rate
end

"""
    accretion_rain_sink(rain::Rain, ice::CloudIce, vel::Blk1MVelTypeRain, ce, q_icl, q_rai, ρ)

Returns the sink of rain water (partial source of snow) due to collisions
with cloud ice.

# Arguments
- `rain`: rain type parameters
- `ice`: ice type parameters
- `vel`: terminal velocity parameters for rain
- `ce`: collision efficiency parameters
- `q_icl`: cloud ice specific content
- `q_rai`: rain water specific content
- `ρ`: air density
"""
function accretion_rain_sink(
    rain::CMP.Rain{FT},
    ice::CMP.CloudIce{FT},
    vel::CMP.Blk1MVelTypeRain{FT},
    ce::CMP.CollisionEff,
    q_icl::FT,
    q_rai::FT,
    ρ::FT,
) where {FT}
    accr_rate = FT(0)
    if (q_icl > CO.ϵ_numerics(FT) && q_rai > CO.ϵ_numerics(FT))

        n0_ice = get_n0(ice.pdf)
        λ_ice_inv = lambda_inverse(ice.pdf, ice.mass, q_icl, ρ)

        n0 = get_n0(rain.pdf, q_rai, ρ)
        v0 = get_v0(vel, ρ)
        (; r0, m0, me, Δm, χm) = rain.mass
        (; χv, ve, Δv) = vel
        (; a0, ae, χa, Δa) = rain.area

        E = Ec(ice, rain, ce)

        λ_inv = lambda_inverse(rain.pdf, rain.mass, q_rai, ρ)

        accr_rate =
            E / ρ * n0 * n0_ice * m0 * a0 * v0 * χm * χa * χv * λ_ice_inv * λ_inv *
            SF.gamma(me + ae + ve + Δm + Δa + Δv + 1) /
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

# Arguments
- `type_i`: snow (T < T_freeze) or rain (T > T_freeze)
- `type_j`: rain (T < T_freeze) or snow (T > T_freeze)  
- `blk1mveltype_ti`, `blk1mveltype_tj`: 1M terminal velocity parameters
- `ce`: collision efficiency parameters (contains `e_rai_sno`, `coeff_disp`)
- `q_i`, `q_j`: specific contents of snow or rain [kg/kg]
- `ρ`: air density [kg/m³]

# Returns
- Accretion rate [kg/kg/s]
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
    if (q_i > CO.ϵ_numerics(FT) && q_j > CO.ϵ_numerics(FT))

        n0_i = get_n0(type_i.pdf, q_i, ρ)
        n0_j = get_n0(type_j.pdf, q_j, ρ)

        (; r0, m0, me, Δm, χm, gamma_coeff) = type_j.mass
        δ = me + Δm

        E_ij = Ec(type_i, type_j, ce)

        λ_i_inv = lambda_inverse(type_i.pdf, type_i.mass, q_i, ρ)
        λ_j_inv = lambda_inverse(type_j.pdf, type_j.mass, q_j, ρ)

        v_ti = terminal_velocity(type_i, blk1mveltype_ti, ρ, q_i)
        v_tj = terminal_velocity(type_j, blk1mveltype_tj, ρ, q_j)

        # Add simple parameterization for velocity dispersion, assuming that fall velocity 
        # standard deviations are proportional to the mean fall velocities, with coefficient 
        # ce.coeff_disp
        Δv_eff = sqrt((v_ti - v_tj)^2 + ce.coeff_disp * (v_ti^2 + v_tj^2))

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
    evaporation_sublimation(rain::Rain, vel::Blk1MVelTypeRain, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)
    evaporation_sublimation(snow::Snow, vel::Blk1MVelTypeSnow, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, T)

Returns the tendency due to rain evaporation or snow sublimation/deposition.
Ventilation factor parameterization follows Seifert and Beheng (2006),
https://doi.org/10.1007/s00703-005-0112-4.

**Rain**: Only evaporation is considered (S < 0), result is clamped to be ≤ 0.
**Snow**: Both sublimation (S < 0) and deposition (S > 0) are considered.

# Arguments
- `rain` or `snow`: particle parameters (contains `pdf`, `mass`, `vent`)
- `vel`: terminal velocity parameters (Blk1MVelTypeRain or Blk1MVelTypeSnow)
- `aps`: air properties struct
- `tps`: thermodynamics parameters struct
- `q_tot`: total water specific content [kg/kg]
- `q_lcl`: cloud liquid water specific content [kg/kg]
- `q_icl`: cloud ice specific content [kg/kg]
- `q_rai`: rain specific content [kg/kg]
- `q_sno`: snow specific content [kg/kg]
- `ρ`: air density [kg/m³]
- `T`: air temperature [K]

# Returns
- Evaporation/sublimation/deposition rate [kg/kg/s]
"""
function evaporation_sublimation(
    (; pdf, mass, vent)::CMP.Rain{FT},
    vel::CMP.Blk1MVelTypeRain{FT},
    aps::CMP.AirProperties{FT},
    tps::TDI.PS,
    q_tot::FT,
    q_lcl::FT,
    q_icl::FT,
    q_rai::FT,
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT}
    evap_subl_rate = FT(0)
    S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)

    if (q_rai > CO.ϵ_numerics(FT) && S < FT(0))

        (; ν_air, D_vapor) = aps
        G = CO.G_func_liquid(aps, tps, T)
        n0 = get_n0(pdf, q_rai, ρ)
        v0 = get_v0(vel, ρ)
        (; χv, ve, Δv, gamma_vent) = vel
        (; r0) = mass
        a_vent = vent.a
        b_vent = vent.b

        λ_inv = lambda_inverse(pdf, mass, q_rai, ρ)

        # Schmidt number
        Sc = ν_air / D_vapor

        evap_subl_rate =
            4 * FT(π) * n0 / ρ * S * G * λ_inv^2 *
            (
                a_vent +
                b_vent * cbrt(Sc) /
                (r0 / λ_inv)^((ve + Δv) / 2) *
                sqrt(2 * v0 * χv / ν_air * λ_inv) *
                gamma_vent
            )
    end
    # only evaporation is considered for rain
    return min(0, evap_subl_rate)
end

function evaporation_sublimation(
    (; pdf, mass, vent)::CMP.Snow{FT},
    vel::CMP.Blk1MVelTypeSnow{FT},
    aps::CMP.AirProperties{FT},
    tps::TDI.PS,
    q_tot::FT,
    q_lcl::FT,
    q_icl::FT,
    q_rai::FT,
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT}
    evap_subl_rate = FT(0)
    if q_sno > CO.ϵ_numerics(FT)
        (; ν_air, D_vapor) = aps

        S = TDI.supersaturation_over_ice(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)
        G = CO.G_func_ice(aps, tps, T)

        n0 = get_n0(pdf, q_sno, ρ)
        v0 = get_v0(vel, ρ)
        (; r0) = mass
        (; χv, ve, Δv, gamma_vent) = vel

        a_vent = vent.a
        b_vent = vent.b

        λ_inv = lambda_inverse(pdf, mass, q_sno, ρ)

        # Schmidt number
        Sc = ν_air / D_vapor

        evap_subl_rate =
            4 * FT(π) * n0 / ρ * S * G * λ_inv^2 *
            (
                a_vent +
                b_vent * cbrt(Sc) /
                (r0 / λ_inv)^((ve + Δv) / 2) *
                sqrt(2 * v0 * χv / ν_air * λ_inv) *
                gamma_vent
            )
    end
    # both sublimation (S < 0) and deposition (S > 0) are considered for snow
    return evap_subl_rate
end

"""
    snow_melt(snow::Snow, vel::Blk1MVelTypeSnow, aps, tps, q_sno, ρ, T)

Returns the tendency due to snow melt.

# Arguments
- `snow`: snow parameters (contains `T_freeze`, `pdf`, `mass`, `vent`)
- `vel`: terminal velocity parameters
- `aps`: air properties struct
- `tps`: thermodynamics parameters struct
- `q_sno`: snow water specific content
- `ρ`: air density
- `T`: air temperature
"""
function snow_melt(
    (; T_freeze, pdf, mass, vent)::CMP.Snow{FT},
    vel::CMP.Blk1MVelTypeSnow{FT},
    aps::CMP.AirProperties{FT},
    tps::TDI.PS,
    q_sno::FT,
    ρ::FT,
    T::FT,
) where {FT}
    snow_melt_rate = FT(0)

    if (q_sno > CO.ϵ_numerics(FT) && T > T_freeze)
        (; ν_air, D_vapor, K_therm) = aps

        L = TDI.Lf(tps, T)

        n0 = get_n0(pdf, q_sno, ρ)
        v0 = get_v0(vel, ρ)
        (; r0) = mass
        (; χv, ve, Δv, gamma_vent) = vel

        a_vent = vent.a
        b_vent = vent.b

        λ_inv = lambda_inverse(pdf, mass, q_sno, ρ)

        # Schmidt number
        Sc = ν_air / D_vapor

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
