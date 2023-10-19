export ArizonaTestDust

"""
    ArizonaTestDust{FT}

Parameters for Arizona Test Dust from
Mohler et al, 2006. DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ArizonaTestDust{FT} <: AerosolType{FT}
    "S₀ for T > T_thr [-]"
    S₀_warm::FT
    "S₀ for T < T_thr [-]"
    S₀_cold::FT
    "a for T > T_thr [-]"
    a_warm::FT
    "a for T < T_thr [-]"
    a_cold::FT
end

function ArizonaTestDust(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return ArizonaTestDust(
        FT(data["Mohler2006_S0_warm_ArizonaTestDust"]["value"]),
        FT(data["Mohler2006_S0_cold_ArizonaTestDust"]["value"]),
        FT(data["Mohler2006_a_warm_ArizonaTestDust"]["value"]),
        FT(data["Mohler2006_a_cold_ArizonaTestDust"]["value"]),
    )
end
