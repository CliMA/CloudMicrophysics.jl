export Kaolinite

"""
    Kaolinite{FT}

Parameters for kaolinite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D and China et al 2017
DOI: 10.1002/2016JD025817

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct Kaolinite{FT} <: AerosolType{FT}
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end

Kaolinite(::Type{FT}) where {FT <: AbstractFloat} =
    Kaolinite(CP.create_toml_dict(FT))

function Kaolinite(td::CP.ParamDict)
    name_map = (;
        :China2017_J_deposition_m_Kaolinite => :deposition_m,
        :China2017_J_deposition_c_Kaolinite => :deposition_c,
        :KnopfAlpert2013_J_ABIFM_m_Kaolinite => :ABIFM_m,
        :KnopfAlpert2013_J_ABIFM_c_Kaolinite => :ABIFM_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return Kaolinite{FT}(; parameters...)
end
