export Illite

"""
    Illite{FT}

Parameters for illite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Illite{FT} <: AerosolType{FT}
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function Illite(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Illite(
        FT(data["KnopfAlpert2013_J_ABIFM_m_Illite"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_c_Illite"]["value"]),
    )
end
