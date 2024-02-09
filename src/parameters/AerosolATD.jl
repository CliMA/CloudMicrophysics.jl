export ArizonaTestDust

"""
    ArizonaTestDust{FT}

Parameters for Arizona Test Dust from
Mohler et al, 2006. DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ArizonaTestDust{FT} <: AerosolType{FT}
    "S₀ for T > T_thr [-]"
    S₀_warm::FT
    "S₀ for T < T_thr [-]"
    S₀_cold::FT
    "a for T > T_thr [-]"
    a_warm::FT
    "a for T < T_thr [-]"
    a_cold::FT
end

ArizonaTestDust(::Type{FT}) where {FT <: AbstractFloat} =
    ArizonaTestDust(CP.create_toml_dict(FT))

function ArizonaTestDust(td::CP.AbstractTOMLDict)
    name_map = (;
        :Mohler2006_S0_warm_ArizonaTestDust => :S₀_warm,
        :Mohler2006_S0_cold_ArizonaTestDust => :S₀_cold,
        :Mohler2006_a_warm_ArizonaTestDust => :a_warm,
        :Mohler2006_a_cold_ArizonaTestDust => :a_cold,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return ArizonaTestDust{FT}(; parameters...)
end
