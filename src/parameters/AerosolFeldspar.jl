export Feldspar

"""
    Feldspar{FT}

Parameters for Feldspar from Alpert et al 2022
DOI: 10.1039/D1EA00077B

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Feldspar{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end
