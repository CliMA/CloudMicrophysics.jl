export AsianDust

"""
    AsianDust{FT}

Parameters for Asian Dust

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AsianDust{FT} <: AerosolType
    "m coefficient for deposition nucleation J [-]"
    deposition_m::FT
    "c coefficient for deposition nucleation J [-]"
    deposition_c::FT
    "m coefficient for immersion freezing J [-]"
    ABIFM_m::FT
    "c coefficient for immersion freezing J [-]"
    ABIFM_c::FT
end
