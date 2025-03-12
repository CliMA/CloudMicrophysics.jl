export SaharanDust

"""
    SaharanDust{FT}

Parameters for Saharan Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SaharanDust{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end

SaharanDust(::Type{FT}) where {FT <: AbstractFloat} =
    SaharanDust(CP.create_toml_dict(FT))

function SaharanDust(td::CP.AbstractTOMLDict)
    name_map = (;
        :J_ABDINM_m_SaharanDust => :deposition_m,
        :J_ABDINM_c_SaharanDust => :deposition_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return SaharanDust{FT}(; parameters...)
end
