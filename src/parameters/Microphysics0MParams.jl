export Microphysics0MParams

"""
    Microphysics0MParams{P, PU, E}

Unified parameter container for 0-moment microphysics.

# Fields
$(DocStringExtensions.FIELDS)

# Constructors

    Microphysics0MParams(FT; with_evaporation = false)
    Microphysics0MParams(toml_dict::ClimaParams.ParamDict; with_evaporation = false)
    Microphysics0MParams(; precip, precip_updraft = precip, evaporation = nothing)

where
- `FT::Type`: A floating point type, e.g. `Float32`
- `toml_dict::ParamDict`: A ClimaParams parameter TOML dictionary (`ParamDict`)
- `with_evaporation`: Include the precipitation evaporation parameters
  ([`Evaporation0M`](@ref)); otherwise `evaporation` is `nothing` (default: false)

The TOML constructor builds `precip_updraft` with the same thresholds as `precip`
and the timescale `precipitation_timescale_updraft_ratio * precipitation_timescale`.

# Example
```julia
import CloudMicrophysics.Parameters as CMP

# Precipitation removal only
mp = CMP.Microphysics0MParams(Float64)

# Precipitation removal + evaporation of the precipitation flux
mp_evap = CMP.Microphysics0MParams(Float64; with_evaporation = true)
```
"""
@kwdef struct Microphysics0MParams{P, PU, E} <: ParametersType
    "Precipitation removal parameters"
    precip::P
    "Precipitation removal parameters for convective updrafts"
    precip_updraft::PU = precip
    "Precipitation evaporation parameters (`nothing` without evaporation)"
    evaporation::E = nothing
end
function Microphysics0MParams(toml_dict::CP.ParamDict; with_evaporation = false)
    precip = Parameters0M(toml_dict)
    (; τ_precip_updraft_ratio) = make_params(toml_dict, name_map(Microphysics0MParams))
    precip_updraft = Parameters0M(;
        precip.qc_0,
        precip.S_0,
        τ_precip = τ_precip_updraft_ratio * precip.τ_precip,
    )
    return Microphysics0MParams(;
        precip,
        precip_updraft,
        evaporation = with_evaporation ? Evaporation0M(toml_dict) : nothing,
    )
end
