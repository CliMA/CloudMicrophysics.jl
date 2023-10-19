export AerosolActivationParameters

"""
    AerosolActivationParameters{FT}

Parameters for Abdul-Razzak and Ghan 2000 aerosol activation scheme
DOI: 10.1029/1999JD901161

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AerosolActivationParameters{FT} <: ParametersType{FT}
    "molar mass of water [kg/mol]"
    M_w::FT
    "gas constant [J/mol/K]"
    R::FT
    "cloud water density [kg/m3]"
    ρ_w::FT
    "surface tension of water [N/m]"
    σ::FT
    "gravitational_acceleration [m/s2]"
    g::FT
end

function AerosolActivationParameters(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return AerosolActivationParameters(
        FT(data["molar_mass_water"]["value"]),
        FT(data["gas_constant"]["value"]),
        FT(data["density_liquid_water"]["value"]),
        FT(data["surface_tension_water"]["value"]),
        FT(data["gravitational_acceleration"]["value"]),
    )
end
