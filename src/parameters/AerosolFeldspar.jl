export Feldspar

"""
    Feldspar{FT}

Parameters for Feldspar from Alpert et al 2022
DOI: 10.1039/D1EA00077B

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Feldspar{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end

function Feldspar(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Feldspar(
        FT(data["Alpert2022_J_deposition_m_Feldspar"]["value"]),
        FT(data["Alpert2022_J_deposition_c_Feldspar"]["value"]),
    )
end
