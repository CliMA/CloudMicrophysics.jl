export Kaolinite

"""
    Kaolinite{FT}

Parameters for kaolinite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D and China et al 2017
DOI: 10.1002/2016JD025817

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Kaolinite{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function Kaolinite(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Kaolinite(
        FT(data["China2017_J_deposition_m_Kaolinite"]["value"]),
        FT(data["China2017_J_deposition_c_Kaolinite"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_m_Kaolinite"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_c_Kaolinite"]["value"]),
    )
end
