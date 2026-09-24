export Illite

"""
    Illite{FT}

Parameters for illite from Knopf and Alpert 2013
DOI: 10.1039/C3FD00035D

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Illite{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end
