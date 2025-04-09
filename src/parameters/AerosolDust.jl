export Dust

"""
    Dust{FT}

Parameters for generic dust.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Dust{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

Dust(::Type{FT}) where {FT <: AbstractFloat} = Dust(CP.create_toml_dict(FT))

function Dust(td::CP.AbstractTOMLDict)
    name_map = (;
        :J_ABDINM_m_Dust => :deposition_m,
        :J_ABDINM_c_Dust => :deposition_c,
        :J_ABIFM_m_Dust => :ABIFM_m,
        :J_ABIFM_c_Dust => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Dust{FT}(; parameters...)
end
