export Microphysics0MParams

"""
    Microphysics0MParams{P}

Unified parameter container for 0-moment microphysics.

# Fields
$(DocStringExtensions.FIELDS)

# Constructors

    Microphysics0MParams(FT)
    Microphysics0MParams(toml_dict::ClimaParams.ParamDict)

where
- `FT::Type`: A floating point type, e.g. `Float32`
- `toml_dict::ParamDict`: A ClimaParams parameter TOML dictionary (`ParamDict`)
"""
@kwdef struct Microphysics0MParams{P} <: ParametersType
    "Precipitation removal parameters"
    precip::P
end
Microphysics0MParams(toml_dict::CP.ParamDict) =
    Microphysics0MParams(; precip = Parameters0M(toml_dict))
