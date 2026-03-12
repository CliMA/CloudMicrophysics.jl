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
    "diffusivity of water vapor [m2/s]"
    D_vapor::FT
    "kinematic viscosity of air [m2/s]"
    ν_air::FT
end


function AirProperties(td::CP.ParamDict)
    name_map = (;
        :thermal_conductivity_of_air => :K_therm,
        :diffusivity_of_water_vapor => :D_vapor,
        :kinematic_viscosity_of_air => :ν_air,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return AirProperties(; parameters...)
end
