export WaterProperties

"""
    WaterProperties{FT}

Parameters with water properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
@params struct WaterProperties{FT} <: ParametersType
    "density of liquid water [kg/m3]"
    ρw::FT
    "density of ice [kg/m3]"
    ρi::FT
end
