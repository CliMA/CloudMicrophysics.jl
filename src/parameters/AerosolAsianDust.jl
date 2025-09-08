export AsianDust

"""
    AsianDust{FT}

Parameters for Asian Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AsianDust{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

AsianDust(::Type{FT}) where {FT <: AbstractFloat} =
    AsianDust(CP.create_toml_dict(FT))

function AsianDust(td::CP.ParamDict)
    name_map = (;
        :J_ABDINM_m_AsianDust => :deposition_m,
        :J_ABDINM_c_AsianDust => :deposition_c,
        :J_ABIFM_m_AsianDust => :ABIFM_m,
        :J_ABIFM_c_AsianDust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AsianDust{FT}(; parameters...)
end
