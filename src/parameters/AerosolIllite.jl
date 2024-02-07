export Illite

"""
    Illite{FT}

Parameters for illite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Illite{FT} <: AerosolType{FT}
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

Illite(::Type{FT}) where {FT <: AbstractFloat} = Illite(CP.create_toml_dict(FT))

function Illite(td::CP.AbstractTOMLDict)
    name_map = (;
        :KnopfAlpert2013_J_ABIFM_m_Illite => :ABIFM_m,
        :KnopfAlpert2013_J_ABIFM_c_Illite => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Illite{FT}(; parameters...)
end
