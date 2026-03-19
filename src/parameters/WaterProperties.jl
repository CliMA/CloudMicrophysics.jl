export WaterProperties

"""
    WaterProperties{FT}

Parameters with water properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct WaterProperties{FT} <: ParametersType
    "density of liquid water [kg/m3]"
    ρw::FT
    "density of ice [kg/m3]"
    ρi::FT
end

function WaterProperties(td::CP.ParamDict)
    name_map = (; :density_liquid_water => :ρw, :density_ice_water => :ρi)
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return WaterProperties(; parameters...)
end
