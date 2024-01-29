export DesertDust

"""
    DesertDust{FT}

Parameters for desert dust
from Knopf and Alpert 2013 DOI: 10.1039/C3FD00035D
and from Mohler et al, 2006 DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
struct DesertDust{FT} <: AerosolType{FT}
    "S₀ for T > T_thr [-]"
    S₀_warm::FT
    "S₀ for T < T_thr [-]"
    S₀_cold::FT
    "a for T > T_thr [-]"
    a_warm::FT
    "a for T < T_thr [-]"
    a_cold::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

function DesertDust(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
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
