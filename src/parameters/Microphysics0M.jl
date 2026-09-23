export Parameters0M

"""
    Parameters0M{FT}

Parameters for zero-moment bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Parameters0M{FT} <: ParametersType
    "precipitation timescale [s]"
    τ_precip::FT
    "condensate specific content precipitation threshold [-]"
    qc_0::FT
    "supersaturation precipitation threshold [-]"
    S_0::FT
end
