export Ferrihydrite

"""
    Ferrihydrite{FT}

Parameters for Ferrihydrite from Alpert et al 2022
DOI: 10.1039/D1EA00077B

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Ferrihydrite{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end

function Ferrihydrite(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Ferrihydrite(
        FT(data["Alpert2022_J_deposition_m_Ferrihydrite"]["value"]),
        FT(data["Alpert2022_J_deposition_c_Ferrihydrite"]["value"]),
    )
end
