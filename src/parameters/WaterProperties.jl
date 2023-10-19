export WaterProperties

"""
    WaterProperties{FT}

Parameters with water properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct WaterProperties{FT} <: ParametersType{FT}
    "density of liquid water [kg/m3]"
    ρw::FT
    "density of ice [kg/m3]"
    ρi::FT
end

function WaterProperties(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return WaterProperties(
        FT(data["density_liquid_water"]["value"]),
        FT(data["density_ice_water"]["value"]),
    )
end
