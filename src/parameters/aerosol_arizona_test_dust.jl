"""
    ArizonaTestDust{FT}

Parameters for Arizona Test Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ArizonaTestDust{FT} <: AbstractAerosolProperties
    "S₀ for T > T_thr [-]"
    S₀_warm::FT
    "S₀ for T < T_thr [-]"
    S₀_cold::FT
    "a for T > T_thr [-]"
    a_warm::FT
    "a for T < T_thr [-]"
    a_cold::FT
end
Base.broadcastable(x::ArizonaTestDust) = tuple(x)

function ArizonaTestDust(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return ArizonaTestDust(
        FT(data["Mohler2006_S0_warm_ArizonaTestDust"]["value"]),
        FT(data["Mohler2006_S0_cold_ArizonaTestDust"]["value"]),
        FT(data["Mohler2006_a_warm_ArizonaTestDust"]["value"]),
        FT(data["Mohler2006_a_cold_ArizonaTestDust"]["value"]),
    )
end
