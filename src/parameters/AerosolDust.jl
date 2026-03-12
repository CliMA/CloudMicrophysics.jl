export Dust

"""
    Dust{FT}

Parameters for generic dust.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Dust{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function Dust(td::CP.ParamDict)
    name_map = (;
        :J_ABDINM_m_Dust => :deposition_m,
        :J_ABDINM_c_Dust => :deposition_c,
        :J_ABIFM_m_Dust => :ABIFM_m,
        :J_ABIFM_c_Dust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Dust(; parameters...)
end
