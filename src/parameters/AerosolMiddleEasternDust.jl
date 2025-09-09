export MiddleEasternDust

"""
    MiddleEasternDust{FT}

Parameters for Middle Eastern Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MiddleEasternDust{FT} <: AerosolType{FT}
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

MiddleEasternDust(::Type{FT}) where {FT <: AbstractFloat} =
    MiddleEasternDust(CP.create_toml_dict(FT))

function MiddleEasternDust(td::CP.ParamDict)
    name_map = (;
        :J_ABIFM_m_MiddleEasternDust => :ABIFM_m,
        :J_ABIFM_c_MiddleEasternDust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return MiddleEasternDust{FT}(; parameters...)
end
