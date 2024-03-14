"""
    Module for functions shared by different parameterizations.
"""
module Common

import SpecialFunctions as SF

import Thermodynamics as TD
const TPS = TD.Parameters.ThermodynamicsParameters

import ..Parameters as CMP
const HPS = CMP.H2SO4SolutionParameters

export G_func
export H2SO4_soln_saturation_vapor_pressure
export a_w_xT
export a_w_eT
export a_w_ice
export Chen2022_vel_add
export Chen2022_vel_coeffs_small

"""
    G_func(air_props, tps, T, Liquid())
    G_func(air_props, tps, T, Ice())

 - `air_props` - struct with air parameters
 - `tps` - struct with thermodynamics parameters
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Utility function combining thermal conductivity and vapor diffusivity effects.
"""
function G_func(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TPS,
    T::FT,
    ::TD.Liquid,
) where {FT}
    R_v = TD.Parameters.R_v(tps)
    L = TD.latent_heat_vapor(tps, T)
    p_vs = TD.saturation_vapor_pressure(tps, T, TD.Liquid())

    return FT(1) /
           (L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs)
end
function G_func(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TPS,
    T::FT,
    ::TD.Ice,
) where {FT}
    R_v = TD.Parameters.R_v(tps)
    L = TD.latent_heat_sublim(tps, T)
    p_vs = TD.saturation_vapor_pressure(tps, T, TD.Ice())

    return FT(1) /
           (L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs)
end

"""
    A Heaviside step function
"""
heaviside(x::FT) where {FT} = FT(x > 0)

"""
    logistic_function(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the curve, characterizing steepness of the transition

Returns the value of the logistic function for smooth transitioning at thresholds. This is
a normalized curve changing from 0 to 1 while x varies from 0 to Inf (for positive k). For
x < 0 the value at x = 0 (zero) is returned. For x_0 = 0 H(x) is returned.
"""
function logistic_function(x::FT, x_0::FT, k::FT) where {FT}

    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return FT(1)
    end

    return FT(1) / (FT(1) + exp(-k * (x / x_0 - x_0 / x)))
end

"""
    logistic_function_integral(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the logistic function, characterizing steepness of the transition

Returns the value of the indefinite integral of the logistic function, for smooth transitioning
of piecewise linear profiles at thresholds. This curve smoothly transition from y = 0
for 0 < x < x_0 to y = x - x_0 for x_0 < x.
"""
function logistic_function_integral(x::FT, x_0::FT, k::FT) where {FT}
    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return x
    end

    # translation of the curve in x and y to enforce zero at x = 0
    trnslt = -log(FT(1) - exp(-k)) / k

    kt = k * (x / x_0 - FT(1) + trnslt)
    return (kt > FT(40)) ? x - x_0 : (log(FT(1) + exp(kt)) / k - trnslt) * x_0
end

"""
    H2SO4_soln_saturation_vapor_pressure(prs, x, T)

 - `prs` - a struct with H2SO4 solution free parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns the saturation vapor pressure above a sulphuric acid solution droplet in Pa.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight
"""
function H2SO4_soln_saturation_vapor_pressure(
    (;
        T_max,
        T_min,
        w_2,
        c1,
        c2,
        c3,
        c4,
        c5,
        c6,
        c7,
    )::CMP.H2SO4SolutionParameters{FT},
    x::FT,
    T::FT,
) where {FT}

    @assert T < T_max
    @assert T > T_min

    w_h = w_2 * x
    p_sol =
        exp(
            c1 - c2 * x + c3 * x * w_h - c4 * x * w_h^2 +
            (c5 + c6 * x - c7 * x * w_h) / T,
        ) * 100 # * 100 converts mbar --> Pa
    return p_sol
end

"""
    a_w_xT(H2SO4_prs, tps, x, T)

 - `H2SO4_prs` - a struct with H2SO4 solution free parameters
 - `tps` - a struct with thermodynamics parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns water activity of H2SO4 containing droplet.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight.
"""
function a_w_xT(H2SO4_prs::HPS, tps::TPS, x::FT, T::FT) where {FT}

    p_sol = H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x, T)
    p_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())

    return p_sol / p_sat
end

"""
    a_w_eT(tps, e, T)

 - `tps` - struct with thermodynamics parameters
 - `e` - partial pressure of water [Pa]
 - `T` - air temperature [K].

Returns water activity of pure water droplet.
Valid when droplet is in equilibrium with surroundings.
"""
function a_w_eT(tps::TPS, e::FT, T::FT) where {FT}
    # RH
    return e / TD.saturation_vapor_pressure(tps, T, TD.Liquid())
end

"""
    a_w_ice(tps, T)

- `tps` - struct with thermodynamics parameters
 - `T` - air temperature [K].

Returns water activity of ice.
"""
function a_w_ice(tps::TPS, T::FT) where {FT}

    return TD.saturation_vapor_pressure(tps, T, TD.Ice()) /
           TD.saturation_vapor_pressure(tps, T, TD.Liquid())
end

"""
    Chen2022_vel_coeffs_small(precip_type, velo_scheme, ρ)

 - velo_scheme - type for terminal velocity scheme (contains free parameters)
 - ρ - air density

Returns the coefficients from Appendix B in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function Chen2022_vel_coeffs_small(
    velo_scheme::CMP.Chen2022VelTypeRain,
    ρ::FT,
) where {FT}

    (; ρ0, a, a3_pow, b, b_ρ, c) = velo_scheme

    q = exp(ρ0 * ρ)
    ai = (a[1] * q, a[2] * q, a[3] * q * ρ^a3_pow)
    bi = (b[1] - b_ρ * ρ, b[2] - b_ρ * ρ, b[3] - b_ρ * ρ)
    ci = (c[1], c[2], c[3])

    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end
function Chen2022_vel_coeffs_small(
    velo_scheme::CMP.Chen2022VelTypeSnowIce{FT},
    ρ::FT,
) where {FT}
    (; As, Bs, Cs, Es, Fs, Gs) = velo_scheme

    ai = (Es * ρ^As, Fs * ρ^As)
    bi = (Bs + ρ * Cs, Bs + ρ * Cs)
    ci = (FT(0), Gs)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end

function Chen2022_vel_coeffs_large(
    velo_scheme::CMP.Chen2022VelTypeSnowIce{FT},
    ρ::FT,
) where {FT}
    (; Al, Bl, Cl, El, Fl, Gl, Hl) = velo_scheme

    ai = (Bl * ρ^Al, El * ρ^Al * exp(Hl * ρ))
    bi = (Cl, Fl)
    ci = (FT(0), Gl)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end

"""
    Chen2022_vel_coeffs_large(velo_scheme, ρ)

 - velo_scheme - type for terminal velocity scheme (contains free parameters)
 - ρ - air density

Returns the coefficients from Appendix B (table B4 for large particles) in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function Chen2022_vel_coeffs_large(
    velo_scheme::CMP.Chen2022VelTypeSnowIce{FT},
    ρ::FT,
) where {FT}
    (; Al, Bl, Cl, El, Fl, Gl, Hl) = velo_scheme

    ai = (Bl * ρ^Al, El * ρ^Al * exp(Hl * ρ))
    bi = (Cl, Fl)
    ci = (FT(0), Gl)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end

"""
    Chen2022_vel_add(a, b, c, λ, k)

 - a, b, c, - free parameters defined in Chen etl al 2022
 - λ - size distribution parameter
 - k - size distribution moment for which we compute the bulk fall speed

Returns the addends of the bulk fall speed of rain or ice particles
following Chen et al 2022 DOI: 10.1016/j.atmosres.2022.106171 in [m/s].
We are assuming exponential size distribution and hence μ=0.
"""
function Chen2022_vel_add(a::FT, b::FT, c::FT, λ::FT, k::Int) where {FT}
    μ = 0 # Exponential instead of gamma distribution
    δ = FT(μ + k + 1)
    return a * λ^δ * SF.gamma(b + δ) / (λ + c)^(b + δ) / SF.gamma(δ)
end
end
