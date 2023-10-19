export AirProperties

"""
    AirProperties{FT}

Parameters with air properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AirProperties{FT} <: ParametersType{FT}
    "thermal conductivity of air [w/m/K]"
    K_therm::FT
    "diffusivity of water vapor [m2/s]"
    D_vapor::FT
    "kinematic viscosity of air [m2/s]"
    ν_air::FT
end

function AirProperties(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return AirProperties(
        FT(data["thermal_conductivity_of_air"]["value"]),
        FT(data["diffusivity_of_water_vapor"]["value"]),
        FT(data["kinematic_viscosity_of_air"]["value"]),
    )
end
