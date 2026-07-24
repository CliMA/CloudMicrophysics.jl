export Microphysics1MParams, CloudPhaseParams1M, PrecipPhaseParams1M

"""
    CloudPhaseParams1M{LCL, ICL}

Parameters for cloud-phase (non-precipitating) hydrometeors in 1-moment scheme.

# Fields
- `liquid::LCL`: CloudLiquid — cloud liquid water parameters
- `ice::ICL`: CloudIce — cloud ice parameters
"""
@kwdef struct CloudPhaseParams1M{LCL, ICL} <: ParametersType
    liquid::LCL
    ice::ICL
end

Base.show(io::IO, mime::MIME"text/plain", x::CloudPhaseParams1M) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    PrecipPhaseParams1M{RAI, SNO}

Parameters for precipitating hydrometeors in 1-moment scheme.

# Fields
- `rain::RAI`: Rain — rain parameters
- `snow::SNO`: Snow — snow parameters
"""
@kwdef struct PrecipPhaseParams1M{RAI, SNO} <: ParametersType
    rain::RAI
    snow::SNO
end

Base.show(io::IO, mime::MIME"text/plain", x::PrecipPhaseParams1M) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    Microphysics1MParams{OPT, PPR, CP, PP, AP, VL}

Unified parameter container for 1-moment bulk microphysics.

`options` selects which variant of each process runs. The parameter values each
selected variant needs live in `process_params`, whose fields mirror `options`
one-to-one. For example, `options.rain_autoconversion = Kessler1M()` picks the
Kessler scheme, and its `τ`, threshold, and `k` live in
`process_params.rain_autoconversion`. Turning a process off with `nothing` (e.g.
`options.cloud_ice_melt = nothing`) gives it a `nothing` slot in
`process_params`, as does any option that needs no parameters. Shared parameters
(particle size distributions, air properties, terminal velocities) are stored
directly.

# Fields
- `options::OPT`: Microphysics1MOptions — process selection
- `process_params::PPR`: parameter data for each selected process, mirroring
  `options` (`nothing` when a process is disabled with `nothing` or needs no
  parameters)
- `cloud::CP`: CloudPhaseParams1M — cloud liquid and ice parameters
- `precip::PP`: PrecipPhaseParams1M — rain and snow parameters
- `air_properties::AP`: AirProperties — air properties (diffusivities, thermal conductivity)
- `terminal_velocity::VL`: Blk1MVelType — terminal velocity parameters for rain and snow

# Constructors

    Microphysics1MParams(::Type{FT}; options_kwargs...)
    Microphysics1MParams(toml_dict::CP.ParamDict; options_kwargs...)

Create a `Microphysics1MParams` from a float type or a ClimaParams TOML dictionary.
Pass keyword arguments to override individual options.

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Create 1M parameters with default (baseline) settings
mp = CMP.Microphysics1MParams(Float64)

# Create with temperature-dependent ice formation and cloud ice melt disabled
mp = CMP.Microphysics1MParams(Float64;
    cloud_ice_formation = CMP.TemperatureDependent(),
    cloud_ice_melt = nothing,
)
```
"""
@kwdef struct Microphysics1MParams{OPT, PPR, CP, PP, AP, VL} <: ParametersType
    options::OPT
    process_params::PPR
    cloud::CP
    precip::PP
    air_properties::AP
    terminal_velocity::VL
end
Base.show(io::IO, mime::MIME"text/plain", x::Microphysics1MParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    Microphysics1MParams(toml_dict::CP.ParamDict; options_kwargs...)

Create a `Microphysics1MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `options_kwargs...`: Keyword arguments forwarded to `Microphysics1MOptions`
"""
function Microphysics1MParams(toml_dict::CP.ParamDict; options_kwargs...)
    options = Microphysics1MOptions(; options_kwargs...)
    return Microphysics1MParams(;
        options,
        process_params = microphysics_1m_process_params(toml_dict, options),
        cloud = CloudPhaseParams1M(;
            liquid = CloudLiquid(toml_dict),
            ice = CloudIce(toml_dict),
        ),
        precip = PrecipPhaseParams1M(;
            rain = Rain(toml_dict),
            snow = Snow(toml_dict),
        ),
        air_properties = AirProperties(toml_dict),
        terminal_velocity = Blk1MVelType(toml_dict),
    )
end
