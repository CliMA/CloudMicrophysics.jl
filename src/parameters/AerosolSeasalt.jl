export Seasalt

"""
    Seasalt{FT}

Parameters for seasalt

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Seasalt{FT} <: AerosolType
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
    return Seasalt(; parameters...)
end
