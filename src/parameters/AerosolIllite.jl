export Illite

"""
    Illite{FT}

Parameters for illite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Illite{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function Illite(td::CP.ParamDict)
    name_map = (;
        :J_ABDINM_m_Illite => :deposition_m,
        :J_ABDINM_c_Illite => :deposition_c,
        :KnopfAlpert2013_J_ABIFM_m_Illite => :ABIFM_m,
        :KnopfAlpert2013_J_ABIFM_c_Illite => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Illite(; parameters...)
end
