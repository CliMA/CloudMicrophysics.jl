"""
    CloudMicrophysicsParameters0M{FT}

Parameters for zero-moment bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudMicrophysicsParameters0M{FT} <: AbstractCloudMicrophysicsParameters
    "precipitation timescale [s]"
    τ_precip::FT
    "specific humidity precipitation threshold [-]"
    qc_0::FT
    "supersaturation precipitation threshold [-]"
    S_0::FT
end

CloudMicrophysicsParameters0M(param_set) = CloudMicrophysicsParameters0M(
    param_set.τ_precip,
    param_set.qc_0,
    param_set.S_0,
)
