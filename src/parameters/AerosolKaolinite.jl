export Kaolinite

"""
    Kaolinite{FT}

Parameters for kaolinite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Kaolinite{FT} <: AerosolType{FT}
    "m coefficient for immersion freezing J [-]"
    m::FT
    "c coefficient for immersion freezing J [-]"
    c::FT
end

function Kaolinite(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Kaolinite(
        FT(data["KnopfAlpert2013_J_ABIFM_m_Kaolinite"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_c_Kaolinite"]["value"]),
    )
end
