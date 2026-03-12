export SaharanDust

"""
    SaharanDust{FT}

Parameters for Saharan Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SaharanDust{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end

function SaharanDust(td::CP.ParamDict)
    name_map = (;
        :J_ABDINM_m_SaharanDust => :deposition_m,
        :J_ABDINM_c_SaharanDust => :deposition_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return SaharanDust(; parameters...)
end
