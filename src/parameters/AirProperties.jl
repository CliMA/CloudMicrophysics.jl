export AirProperties

"""
    AirProperties{FT}

Parameters with air properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AirProperties{FT} <: ParametersType
    "thermal conductivity of air [w/m/K]"
    K_therm::FT
    "diffusivity of water vapor [m²/s]"
    D_vapor::FT
    "kinematic viscosity of air [m²/s]"
    ν_air::FT
end

ShowMethods.field_units(::AirProperties) =
    (; K_therm = "W/m/K", D_vapor = "m²/s", ν_air = "m²/s")

function AirProperties(td::CP.ParamDict)
    name_map = (;
        :thermal_conductivity_of_air => :K_therm,
        :diffusivity_of_water_vapor => :D_vapor,
        :kinematic_viscosity_of_air => :ν_air,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return AirProperties(; parameters...)
end
