export Seasalt

"""
    Seasalt{FT}

Parameters for seasalt

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Seasalt{FT} <: AerosolType{FT}
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

function Seasalt(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Seasalt(
        FT(data["seasalt_aerosol_molar_mass"]["value"]),
        FT(data["seasalt_aerosol_density"]["value"]),
        FT(data["seasalt_aerosol_osmotic_coefficient"]["value"]),
        FT(data["seasalt_aerosol_ion_number"]["value"]),
        FT(data["seasalt_aerosol_water_soluble_mass_fraction"]["value"]),
        FT(data["seasalt_aerosol_kappa"]["value"]),
    )
end
