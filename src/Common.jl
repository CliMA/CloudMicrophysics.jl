"""
    Module for functions shared by different parameterizations.
"""
module Common

using Thermodynamics

using CLIMAParameters
using CLIMAParameters.Planet: R_v
using CLIMAParameters.Atmos.Microphysics

const APS = AbstractParameterSet

export G_func

"""
    G_func(param_set, T, Liquid())
    G_func(param_set, T, Ice())

 - `param_set` - abstract set with earth parameters
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Utility function combining thermal conductivity and vapor diffusivity effects.
"""
function G_func(param_set::APS, T::FT, ::Liquid) where {FT <: Real}

    _K_therm::FT = K_therm(param_set)
    _R_v::FT = R_v(param_set)
    _D_vapor::FT = D_vapor(param_set)

    L = latent_heat_vapor(param_set, T)
    p_vs = saturation_vapor_pressure(param_set, T, Liquid())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end
function G_func(param_set::APS, T::FT, ::Ice) where {FT <: Real}

    _K_therm::FT = K_therm(param_set)
    _R_v::FT = R_v(param_set)
    _D_vapor::FT = D_vapor(param_set)

    L = latent_heat_sublim(param_set, T)
    p_vs = saturation_vapor_pressure(param_set, T, Ice())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end

end
