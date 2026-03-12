export ArizonaTestDust

"""
    ArizonaTestDust{FT}

Parameters for Arizona Test Dust from
Mohler et al, 2006. DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ArizonaTestDust{FT} <: AerosolType
    "S₀ for T > T_thr [-]"
    S₀_warm::FT
    "S₀ for T < T_thr [-]"
    S₀_cold::FT
    "a for T > T_thr [-]"
    a_warm::FT
    "a for T < T_thr [-]"
    a_cold::FT
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function ArizonaTestDust(td::CP.ParamDict)
    name_map = (;
        :Mohler2006_S0_warm_ArizonaTestDust => :S₀_warm,
        :Mohler2006_S0_cold_ArizonaTestDust => :S₀_cold,
        :Mohler2006_a_warm_ArizonaTestDust => :a_warm,
        :Mohler2006_a_cold_ArizonaTestDust => :a_cold,
        :J_ABDINM_m_ArizonaTestDust => :deposition_m,
        :J_ABDINM_c_ArizonaTestDust => :deposition_c,
        :J_ABIFM_m_ArizonaTestDust => :ABIFM_m,
        :J_ABIFM_c_ArizonaTestDust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return ArizonaTestDust(; parameters...)
end
