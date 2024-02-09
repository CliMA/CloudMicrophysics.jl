export Sulfate

"""
    Sulfate{FT}

Parameters for sulfate aerosol

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Sulfate{FT} <: AerosolType{FT}
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

Sulfate(::Type{FT}) where {FT <: AbstractFloat} =
    Sulfate(CP.create_toml_dict(FT))

function Sulfate(td::CP.AbstractTOMLDict)
    name_map = (;
        :sulfate_aerosol_molar_mass => :M,
        :sulfate_aerosol_density => :ρ,
        :sulfate_aerosol_osmotic_coefficient => :ϕ,
        :sulfate_aerosol_ion_number => :ν,
        :sulfate_aerosol_water_soluble_mass_fraction => :ϵ,
        :sulfate_aerosol_kappa => :κ,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Sulfate{FT}(; parameters...)
end
