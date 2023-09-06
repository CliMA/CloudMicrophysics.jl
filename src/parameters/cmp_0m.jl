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

function CloudMicrophysicsParameters0M(::Type{FT}) where {FT}
    toml_dict = CP.create_toml_dict(FT)
    (; data) = toml_dict
    return CloudMicrophysicsParameters0M(
        FT(data["precipitation_timescale"]["value"]),
        FT(data["specific_humidity_precipitation_threshold"]["value"]),
        FT(data["supersaturation_precipitation_threshold"]["value"]),
    )
end
