export Microphysics0MParams

"""
    Microphysics0MParams{P, PU, E}

Unified parameter container for 0-moment microphysics.

# Fields
$(DocStringExtensions.FIELDS)

# Constructors

    Microphysics0MParams(FT)
    Microphysics0MParams(toml_dict::ClimaParams.ParamDict)
    Microphysics0MParams(; precip, precip_updraft = precip, evaporation = nothing)

where
- `FT::Type`: A floating point type, e.g. `Float32`
- `toml_dict::ParamDict`: A ClimaParams parameter TOML dictionary (`ParamDict`)

The TOML constructor builds `precip_updraft` with the same thresholds as `precip`
and the timescale `precipitation_timescale_updraft_ratio * precipitation_timescale`.
"""
@kwdef struct Microphysics0MParams{P, PU, E} <: ParametersType
    "Precipitation removal parameters"
    precip::P
    "Precipitation removal parameters for convective updrafts"
    precip_updraft::PU = precip
    "Precipitation evaporation parameters (`nothing` if not provided)"
    evaporation::E = nothing
end
function Microphysics0MParams(toml_dict::CP.ParamDict)
    precip = Parameters0M(toml_dict)
    (; precipitation_timescale_updraft_ratio) = CP.get_parameter_values(
        toml_dict,
        "precipitation_timescale_updraft_ratio",
        "CloudMicrophysics",
    )
    precip_updraft = Parameters0M(;
        precip.qc_0,
        precip.S_0,
        τ_precip = precipitation_timescale_updraft_ratio * precip.τ_precip,
    )
    return Microphysics0MParams(;
        precip,
        precip_updraft,
        evaporation = Evaporation0M(toml_dict),
    )
end
