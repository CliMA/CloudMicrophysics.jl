export AirProperties

"""
    AirProperties{FT}

Parameters with air properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AirProperties{FT} <: ParametersType{FT}
    "thermal conductivity of air [w/m/K]"
    K_therm::FT
    "diffusivity of water vapor [m2/s]"
    D_vapor::FT
    "kinematic viscosity of air [m2/s]"
    ν_air::FT
end

AirProperties(::Type{FT}) where {FT <: AbstractFloat} =
    AirProperties(CP.create_toml_dict(FT))

function AirProperties(td::CP.AbstractTOMLDict)
    name_map = (;
        :thermal_conductivity_of_air => :K_therm,
        :diffusivity_of_water_vapor => :D_vapor,
        :kinematic_viscosity_of_air => :ν_air,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AirProperties{FT}(; parameters...)
end
