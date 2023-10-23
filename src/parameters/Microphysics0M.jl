export Parameters0M

"""
    Parameters0M{FT}

Parameters for zero-moment bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Parameters0M{FT} <: ParametersType{FT}
    "precipitation timescale [s]"
    τ_precip::FT
    "specific humidity precipitation threshold [-]"
    qc_0::FT
    "supersaturation precipitation threshold [-]"
    S_0::FT
end

function Parameters0M(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return Parameters0M(
        FT(data["precipitation_timescale"]["value"]),
        FT(data["specific_humidity_precipitation_threshold"]["value"]),
        FT(data["supersaturation_precipitation_threshold"]["value"]),
    )
end
