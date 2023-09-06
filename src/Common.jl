"""
    Module for functions shared by different parameterizations.
"""
module Common

import SpecialFunctions as SF

import Thermodynamics as TD

import ..Parameters as CMP
const APS = CMP.AbstractCloudMicrophysicsParameters

import ..CommonTypes as CT

export G_func
export H2SO4_soln_saturation_vapor_pressure
export Chen2022_vel_add
export Chen2022_vel_coeffs

"""
    G_func(param_set, T, Liquid())
    G_func(param_set, T, Ice())

 - `param_set` - abstract set with earth parameters
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Utility function combining thermal conductivity and vapor diffusivity effects.
"""
function G_func(param_set::APS, T::FT, ::TD.Liquid) where {FT <: Real}

    thermo_params = CMP.thermodynamics_params(param_set)
    _K_therm::FT = CMP.K_therm(param_set)
    _R_v::FT = CMP.R_v(param_set)
    _D_vapor::FT = CMP.D_vapor(param_set)

    L = TD.latent_heat_vapor(thermo_params, T)
    p_vs = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end
function G_func(param_set::APS, T::FT, ::TD.Ice) where {FT <: Real}

    thermo_params = CMP.thermodynamics_params(param_set)
    _K_therm::FT = CMP.K_therm(param_set)
    _R_v::FT = CMP.R_v(param_set)
    _D_vapor::FT = CMP.D_vapor(param_set)

    L = TD.latent_heat_sublim(thermo_params, T)
    p_vs = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end

"""
    A Heaviside step function
"""
function heaviside(x::FT) where {FT <: Real}
    return FT(x > 0)
end

"""
    logistic_function(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the curve, characterizing steepness of the transition

Returns the value of the logistic function for smooth transitioning at thresholds. This is
a normalized curve changing from 0 to 1 while x varies from 0 to Inf (for positive k). For
x < 0 the value at x = 0 (zero) is returned. For x_0 = 0 H(x) is returned.
"""
function logistic_function(x::FT, x_0::FT, k::FT) where {FT <: Real}

    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return FT(1)
    end

    return 1 / (1 + exp(-k * (x / x_0 - x_0 / x)))
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
function logistic_function_integral(x::FT, x_0::FT, k::FT) where {FT <: Real}

    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return x
    end

    # translation of the curve in x and y to enforce zero at x = 0
    _trnslt::FT = -log(1 - exp(-k)) / k

    _kt::FT = k * (x / x_0 - 1 + _trnslt)
    _result::FT =
        (_kt > 40.0) ? x - x_0 : (log(1 + exp(_kt)) / k - _trnslt) * x_0
    return _result
end

"""
    H2SO4_soln_saturation_vapor_pressure(prs, x, T)

 - `prs` - a set with free parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns the saturation vapor pressure above a sulphuric acid solution droplet in Pa.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight
"""
function H2SO4_soln_saturation_vapor_pressure(
    prs::APS,
    x::FT,
    T::FT,
) where {FT <: Real}

    T_max::FT = CMP.H2SO4_sol_T_max(prs)
    T_min::FT = CMP.H2SO4_sol_T_min(prs)
    w_2::FT = CMP.H2SO4_sol_w_2(prs)

    c1::FT = CMP.H2SO4_sol_c1(prs)
    c2::FT = CMP.H2SO4_sol_c2(prs)
    c3::FT = CMP.H2SO4_sol_c3(prs)
    c4::FT = CMP.H2SO4_sol_c4(prs)
    c5::FT = CMP.H2SO4_sol_c5(prs)
    c6::FT = CMP.H2SO4_sol_c6(prs)
    c7::FT = CMP.H2SO4_sol_c7(prs)

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
    Delta_a_w(prs, x, T)

 - `prs` - set with model parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns the change in water activity when droplet undergoes immersion freezing.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight.
"""
function Delta_a_w(prs::APS, x::FT, T::FT) where {FT <: Real}

    thermo_params = CMP.thermodynamics_params(prs)

    p_sol = H2SO4_soln_saturation_vapor_pressure(prs, x, T)
    p_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    p_ice = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())

    a_w = p_sol / p_sat
    a_w_ice = p_ice / p_sat
    Δa_w = a_w - a_w_ice

    return min(Δa_w, FT(1))
end

"""
    Chen2022_snow_ice_coeffs(prs, ρ_i)

 - prs - set with model parameters
 - ρ_i - cloud ice density

Returns the coefficients from Appendix B, Table B3 in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
needed for snow and ice terminal velocity
"""
function Chen2022_snow_ice_coeffs(prs::APS, ρ_i::FT) where {FT <: Real}

    As_1::FT = CMP.As_coeff_1_Ch2022(prs)
    As_2::FT = CMP.As_coeff_2_Ch2022(prs)
    As_3::FT = CMP.As_coeff_3_Ch2022(prs)
    Bs_1::FT = CMP.Bs_coeff_1_Ch2022(prs)
    Bs_2::FT = CMP.Bs_coeff_2_Ch2022(prs)
    Bs_3::FT = CMP.Bs_coeff_3_Ch2022(prs)
    Cs_1::FT = CMP.Cs_coeff_1_Ch2022(prs)
    Cs_2::FT = CMP.Cs_coeff_2_Ch2022(prs)
    Cs_3::FT = CMP.Cs_coeff_3_Ch2022(prs)
    Cs_4::FT = CMP.Cs_coeff_4_Ch2022(prs)
    Es_1::FT = CMP.Es_coeff_1_Ch2022(prs)
    Es_2::FT = CMP.Es_coeff_2_Ch2022(prs)
    Es_3::FT = CMP.Es_coeff_3_Ch2022(prs)
    Fs_1::FT = CMP.Fs_coeff_1_Ch2022(prs)
    Fs_2::FT = CMP.Fs_coeff_2_Ch2022(prs)
    Fs_3::FT = CMP.Fs_coeff_3_Ch2022(prs)
    Gs_1::FT = CMP.Gs_coeff_1_Ch2022(prs)
    Gs_2::FT = CMP.Gs_coeff_2_Ch2022(prs)
    Gs_3::FT = CMP.Gs_coeff_3_Ch2022(prs)

    As = As_1 * (log(ρ_i))^2 − As_2 * log(ρ_i) - As_3
    Bs = FT(1) / (Bs_1 + Bs_2 * log(ρ_i) + Bs_3 / sqrt(ρ_i))
    Cs = Cs_1 + Cs_2 * exp(Cs_3 * ρ_i) + Cs_4 * sqrt(ρ_i)
    Es = Es_1 - Es_2 * (log(ρ_i))^2 + Es_3 * sqrt(ρ_i)
    Fs = -exp(Fs_1 - Fs_2 * (log(ρ_i))^2 + Fs_3 * log(ρ_i))
    Gs = FT(1) / (Gs_1 + Gs_2 / (log(ρ_i)) - Gs_3 * log(ρ_i) / ρ_i)

    return (As, Bs, Cs, Es, Fs, Gs)
end

"""
    Chen2022_vel_coeffs(prs, precip_type, ρ)

 - prs - set with free parameters
 - precip_type - type for ice, rain or snow
 - ρ - air density

Returns the coefficients from Appendix B in Chen et al 2022
DOI: 10.1016/j.atmosres.2022.106171
"""
function Chen2022_vel_coeffs(prs::APS, ::CT.RainType, ρ::FT) where {FT <: Real}

    ρ0::FT = CMP.q_coeff_rain_Ch2022(prs)
    a1::FT = CMP.a1_coeff_rain_Ch2022(prs)
    a2::FT = CMP.a2_coeff_rain_Ch2022(prs)
    a3::FT = CMP.a3_coeff_rain_Ch2022(prs)
    a3_pow::FT = CMP.a3_pow_coeff_rain_Ch2022(prs)
    b1::FT = CMP.b1_coeff_rain_Ch2022(prs)
    b2::FT = CMP.b2_coeff_rain_Ch2022(prs)
    b3::FT = CMP.b3_coeff_rain_Ch2022(prs)
    b_ρ::FT = CMP.b_rho_coeff_rain_Ch2022(prs)
    c1::FT = CMP.c1_coeff_rain_Ch2022(prs)
    c2::FT = CMP.c2_coeff_rain_Ch2022(prs)
    c3::FT = CMP.c3_coeff_rain_Ch2022(prs)

    q = exp(ρ0 * ρ)
    ai = (a1 * q, a2 * q, a3 * q * ρ^a3_pow)
    bi = (b1 - b_ρ * ρ, b2 - b_ρ * ρ, b3 - b_ρ * ρ)
    ci = (c1, c2, c3)

    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000

    return (aiu, bi, ciu)
end
function Chen2022_vel_coeffs(
    prs::APS,
    ::Union{CT.IceType, CT.SnowType},
    ρ::FT,
) where {FT <: Real}

    ρ_i::FT = CMP.ρ_cloud_ice(prs)

    _As, _Bs, _Cs, _Es, _Fs, _Gs = Chen2022_snow_ice_coeffs(prs, ρ_i)

    ai = (_Es * ρ^_As, _Fs * ρ^_As)
    bi = (_Bs + ρ * _Cs, _Bs + ρ * _Cs)
    ci = (FT(0), _Gs)
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
function Chen2022_vel_add(a::FT, b::FT, c::FT, λ::FT, k::Int) where {FT <: Real}
    μ = 0 # Exponantial instaed of gamma distribution
    δ = μ + k + 1
    return a * λ^δ * SF.gamma(b + δ) / (λ + c)^(b + δ) / SF.gamma(δ)
end
end
