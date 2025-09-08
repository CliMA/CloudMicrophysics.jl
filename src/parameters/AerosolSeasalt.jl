export Seasalt

"""
    Seasalt{FT}

Parameters for seasalt

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Seasalt{FT} <: AerosolType{FT}
    "molar mass [kg/mol]"
    M::FT
    "density [kg/m3]"
    ρ::FT
    "osmotic coefficient [-]"
    ϕ::FT
    "ion number [-]"
    ν::FT
    "water soluble mass fraction [-]"
    ϵ::FT
    "hygroscopicity parameter [-]"
    κ::FT
end

Seasalt(::Type{FT}) where {FT <: AbstractFloat} =
    Seasalt(CP.create_toml_dict(FT))

function Seasalt(td::CP.ParamDict)
    name_map = (;
        :seasalt_aerosol_molar_mass => :M,
        :seasalt_aerosol_density => :ρ,
        :seasalt_aerosol_osmotic_coefficient => :ϕ,
        :seasalt_aerosol_ion_number => :ν,
        :seasalt_aerosol_water_soluble_mass_fraction => :ϵ,
        :seasalt_aerosol_kappa => :κ,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Seasalt{FT}(; parameters...)
end
