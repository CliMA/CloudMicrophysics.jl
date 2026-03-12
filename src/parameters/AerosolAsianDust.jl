export AsianDust

"""
    AsianDust{FT}

Parameters for Asian Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AsianDust{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function AsianDust(td::CP.ParamDict)
    name_map = (;
        :J_ABDINM_m_AsianDust => :deposition_m,
        :J_ABDINM_c_AsianDust => :deposition_c,
        :J_ABIFM_m_AsianDust => :ABIFM_m,
        :J_ABIFM_c_AsianDust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return AsianDust(; parameters...)
end
