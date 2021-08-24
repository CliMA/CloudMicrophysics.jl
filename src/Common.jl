"""
    Module for functions shared by different parameterizations.
"""
module Common

import Thermodynamics
import CLIMAParameters
import CLIMAParameters.Planet
import CLIMAParameters.Atmos.Microphysics

const TD = Thermodynamics
const CP = CLIMAParameters
const CP_planet = CLIMAParameters.Planet
const CP_micro = CLIMAParameters.Atmos.Microphysics
const APS = CP.AbstractParameterSet

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

    _K_therm::FT = CP_micro.K_therm(param_set)
    _R_v::FT = CP_planet.R_v(param_set)
    _D_vapor::FT = CP_micro.D_vapor(param_set)

    L = TD.latent_heat_vapor(param_set, T)
    p_vs = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end
function G_func(param_set::APS, T::FT, ::TD.Ice) where {FT <: Real}

    _K_therm::FT = CP_micro.K_therm(param_set)
    _R_v::FT = CP_planet.R_v(param_set)
    _D_vapor::FT = CP_micro.D_vapor(param_set)

    L = TD.latent_heat_sublim(param_set, T)
    p_vs = TD.saturation_vapor_pressure(param_set, T, TD.Ice())

    return FT(1) / (
        L / _K_therm / T * (L / _R_v / T - FT(1)) + _R_v * T / _D_vapor / p_vs
    )
end

end
