export Dust

"""
    Dust{FT}

Parameters for generic dust.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Dust{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end
