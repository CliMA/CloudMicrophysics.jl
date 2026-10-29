export Kaolinite

"""
    Kaolinite{FT}

Parameters for kaolinite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D and China et al 2017
DOI: 10.1002/2016JD025817

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Kaolinite{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end
