export WaterProperties

"""
    WaterProperties{FT}

Parameters with water properties.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct WaterProperties{FT} <: ParametersType{FT}
    "density of liquid water [kg/m3]"
    ρw::FT
    "density of ice [kg/m3]"
    ρi::FT
end

WaterProperties(::Type{FT}) where {FT <: AbstractFloat} =
    WaterProperties(CP.create_toml_dict(FT))

function WaterProperties(td::CP.AbstractTOMLDict)
    name_map = (; :density_liquid_water => :ρw, :density_ice_water => :ρi)
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return WaterProperties{FT}(; parameters...)
end
