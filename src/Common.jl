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

end
