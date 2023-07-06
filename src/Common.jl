"""
    Module for functions shared by different parameterizations.
"""
module Common

import Thermodynamics
const TD = Thermodynamics

import ..Parameters
const CMP = Parameters
const APS = Parameters.AbstractCloudMicrophysicsParameters

export G_func
export H2SO4_soln_saturation_vapor_pressure
export ABIFM_Delta_a_w

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
    H2SO4_soln_saturation_vapor_pressure(x, T)

 - `x` - wt percent sulphuric acid [unitless] 
 - `T` - air temperature [K].

Returns the saturation vapor pressure above a sulphuric acid solution droplet in Pa.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight
"""
function H2SO4_soln_saturation_vapor_pressure(x::FT, T::FT) where {FT <: Real}

    @assert T < FT(235)
    @assert T > FT(185)

    w_h = 1.4408 * x
    p_sol =
        exp(
            23.306 - 5.3465 * x + 12 * x * w_h - 8.19 * x * w_h^2 +
            (-5814 + 928.9 * x - 1876.7 * x * w_h) / T,
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

    p_sol = H2SO4_soln_saturation_vapor_pressure(x, T)
    p_sat = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    p_ice = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())

    a_w = p_sol / p_sat
    a_w_ice = p_ice / p_sat
    Δa_w = a_w - a_w_ice

    return min(Δa_w, FT(1))
end

end
