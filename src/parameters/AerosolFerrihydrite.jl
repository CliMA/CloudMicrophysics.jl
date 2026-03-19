export Ferrihydrite

"""
    Ferrihydrite{FT}

Parameters for Ferrihydrite from Alpert et al 2022
DOI: 10.1039/D1EA00077B

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Ferrihydrite{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end

function Ferrihydrite(td::CP.ParamDict)
    name_map = (;
        :Alpert2022_J_deposition_m_Ferrihydrite => :deposition_m,
        :Alpert2022_J_deposition_c_Ferrihydrite => :deposition_c,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Ferrihydrite(; parameters...)
end
