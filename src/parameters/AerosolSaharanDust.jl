export SaharanDust

"""
    SaharanDust{FT}

Parameters for Saharan Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SaharanDust{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
end
