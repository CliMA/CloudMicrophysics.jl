"""
    DesertDust{FT}

Parameters for desert dust

# Fields
$(DocStringExtensions.FIELDS)
"""
struct DesertDust{FT} <: AbstractAerosolProperties
    "S₀ for T > T_thr [-]"
    S₀_warm::FT
    "S₀ for T < T_thr [-]"
    S₀_cold::FT
    "a for T > T_thr [-]"
    a_warm::FT
    "a for T < T_thr [-]"
    a_cold::FT
    "m coefficient for immersion freezing J [-]"
    m::FT
    "c coefficient for immersion freezing J [-]"
    c::FT
end
Base.broadcastable(x::DesertDust) = tuple(x)

function DesertDust(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return DesertDust(
        FT(data["Mohler2006_S0_warm_DesertDust"]["value"]),
        FT(data["Mohler2006_S0_cold_DesertDust"]["value"]),
        FT(data["Mohler2006_a_warm_DesertDust"]["value"]),
        FT(data["Mohler2006_a_cold_DesertDust"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_m_DesertDust"]["value"]),
        FT(data["KnopfAlpert2013_J_ABIFM_c_DesertDust"]["value"]),
    )
end
