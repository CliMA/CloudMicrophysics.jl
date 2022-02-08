"""
    Module for functions shared by different parameterizations.
"""
module Common

import Thermodynamics
const TD = Thermodynamics

import CloudMicrophysics.CloudMicrophysicsParameters

export G_func

"""
    G_func(param_set, T, Liquid())
    G_func(param_set, T, Ice())

 - `param_set` - abstract set with earth parameters
 - `T` - air temperature
 - `Liquid()`, `Ice()` - liquid or ice phase to dispatch over.

Utility function combining thermal conductivity and vapor diffusivity effects.
"""
function G_func(param_set::CloudMicrophysicsParameters, T::FT, ::TD.Liquid) where {FT <: Real}

    K_therm = param_set.K_therm
    R_v = param_set.R_v
    D_vapor = param_set.D_vapor

    L = TD.latent_heat_vapor(param_set.ThermodynamicsParameters, T)
    p_vs = TD.saturation_vapor_pressure(param_set.ThermodynamicsParameters, T, TD.Liquid())
    
    return FT(1) / (
        L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs
    )
end
function G_func(param_set::CloudMicrophysicsParameters, T::FT, ::TD.Ice) where {FT <: Real}

    K_therm = param_set.K_therm
    R_v = param_set.R_v
    D_vapor = param_set.D_vapor

    L = TD.latent_heat_sublim(param_set.ThermodynamicsParameters, T)
    p_vs = TD.saturation_vapor_pressure(param_set.ThermodynamicsParameters, T, TD.Ice())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end








end
