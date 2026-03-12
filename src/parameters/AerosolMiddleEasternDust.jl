export MiddleEasternDust

"""
    MiddleEasternDust{FT}

Parameters for Middle Eastern Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct MiddleEasternDust{FT} <: AerosolType
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function MiddleEasternDust(td::CP.ParamDict)
    name_map = (;
        :J_ABIFM_m_MiddleEasternDust => :ABIFM_m,
        :J_ABIFM_c_MiddleEasternDust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return MiddleEasternDust(; parameters...)
end
