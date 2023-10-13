"""
    AerosolActivationParameters{FT}

Parameters for aerosol activation

# Fields
$(DocStringExtensions.FIELDS)
"""
struct AerosolActivationParameters{FT} <:
       AbstractAerosolActivationParameters{FT}
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

function AerosolActivationParameters(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return AerosolActivationParameters(
        FT(data["molar_mass_water"]["value"]),
        FT(data["gas_constant"]["value"]),
        FT(data["density_liquid_water"]["value"]),
        FT(data["surface_tension_water"]["value"]),
        FT(data["gravitational_acceleration"]["value"]),
    )
end
