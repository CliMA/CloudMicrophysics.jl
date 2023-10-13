"""
    SeasaltParameters{FT}

Parameters for seasalt

# Fields
$(DocStringExtensions.FIELDS)
"""
struct SeasaltParameters{FT} <: AbstractAerosolProperties
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
Base.broadcastable(x::SeasaltParameters) = tuple(x)

function SeasaltParameters(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return SeasaltParameters(
        FT(data["seasalt_aerosol_molar_mass"]["value"]),
        FT(data["seasalt_aerosol_density"]["value"]),
        FT(data["seasalt_aerosol_osmotic_coefficient"]["value"]),
        FT(data["seasalt_aerosol_ion_number"]["value"]),
        FT(data["seasalt_aerosol_water_soluble_mass_fraction"]["value"]),
        FT(data["seasalt_aerosol_kappa"]["value"]),
    )
end
