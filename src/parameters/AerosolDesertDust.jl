export DesertDust

"""
    DesertDust{FT}

Parameters for desert dust
from Knopf and Alpert 2013 DOI: 10.1039/C3FD00035D
and from Mohler et al, 2006 DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct DesertDust{FT} <: AerosolType{FT}
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

DesertDust(::Type{FT}) where {FT <: AbstractFloat} =
    DesertDust(CP.create_toml_dict(FT))

function DesertDust(td::CP.AbstractTOMLDict)
    name_map = (;
        :Mohler2006_S0_warm_DesertDust => :S₀_warm,
        :Mohler2006_S0_cold_DesertDust => :S₀_cold,
        :Mohler2006_a_warm_DesertDust => :a_warm,
        :Mohler2006_a_cold_DesertDust => :a_cold,
        :AlpertKnopf2016_J_ABIFM_m_DesertDust => :ABIFM_m,
        :AlpertKnopf2016_J_ABIFM_c_DesertDust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return DesertDust{FT}(; parameters...)
end
