# Parameters Interface

All free parameters of `CloudMicrophysics.jl` live outside the source code,
  in the `parameters.toml` file of
  [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl).
The library never hard-codes a coefficient: every process function receives
  the parameters it needs through a typed parameter struct, and every such
  struct is filled from a ClimaParams dictionary.
This page explains how those structs are organized, how to build and override
  them, and how to add a new parameter.
The [Parameter reference](generated/ParametersReference.md) page, generated from
  the source code at every documentation build, lists all ClimaParams parameters
  used by the library and the struct fields they map to.

## From a TOML entry to a struct field

A ClimaParams entry stores a value, a type and a description:

```toml
[rain_autoconversion_timescale]
value = 1000.0
type = "float"
description = "Rain formation timescale for the 1-moment microphysics scheme (s)."
```

On the `CloudMicrophysics.jl` side, the parameters of one parameterization are
  grouped in a struct that subtypes `ParametersType`, with a docstring per
  field, and a constructor that reads a `ClimaParams.ParamDict`.
The mapping from the long, unique ClimaParams name to the short field name is
  a *name map*, and all name maps live in one file,
  `src/parameters/name_maps.jl`, as methods of
  [`name_map`](@ref CloudMicrophysics.Parameters.name_map):

```julia
@kwdef struct CloudLiquid{FT} <: ParametersType
    "liquid water density [kg/m3]"
    ρw::FT
    "prescribed liquid cloud droplet effective radius [m]"
    r_eff::FT
    "prescribed cloud droplet number concentration [1/m3]"
    N_0::FT
end

# src/parameters/name_maps.jl
name_map(::Type{CloudLiquid}) = (;
    :density_liquid_water => :ρw,
    :liquid_cloud_effective_radius => :r_eff,
    :cloud_liquid_sedimentation_number_concentration => :N_0,
)
```

For a struct like this one, whose constructor only reads the map and fills the
  fields, the `ParamDict` constructor is generated from the map
  (`CloudLiquid(td) = CloudLiquid(; make_params(td, name_map(CloudLiquid))...)`;
  the struct is listed in `PLAIN_PARAMETER_TYPES`).
Constructors that validate values, pre-compute derived fields (gamma-function
  factors, unit rescaling) or compose other structs are written by hand next to
  their struct, and read their parameters with
  [`make_params`](@ref CloudMicrophysics.Parameters.make_params)`(td, name_map(T))`,
  so that the name map stays the single list of what the struct reads.
`make_params` (a thin wrapper of `CP.get_parameter_values`) returns a
  `NamedTuple` with the values converted to the float type of the dictionary,
  and tags each parameter as used by the `"CloudMicrophysics"` component.
That tag is what lets ClimaParams report unused overrides and write a log of
  the parameters a simulation actually consumed
  (see [Logging and reproducibility](@ref)), and what the test
  `test/parameter_name_maps_tests.jl` uses to check that every name map matches
  what its constructor reads.
The [Parameter reference](generated/ParametersReference.md) is rendered from
  the same name maps.

## Building parameter structs

Every struct with a name map has two constructors, one taking a float type
  and one taking a `ClimaParams.ParamDict` (small building blocks such as
  `Ventilation` or `ParticlePDFSnow` are only constructed by their containers):

```julia
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP

rain = CMP.Rain(Float32)                     # default ClimaParams values, Float32

toml_dict = CP.create_toml_dict(Float64)     # explicit dictionary ...
rain = CMP.Rain(toml_dict)                   # ... shared between structs
```

The `T(FT; kwargs...)` form is defined once for all subtypes of
  `ParametersType` (see the end of `src/parameters/Parameters.jl`) and simply
  calls `T(CP.create_toml_dict(FT); kwargs...)`.
Use the `ParamDict` form when you build several structs, so that all of them
  see the same overrides and the same usage log.

Small structs are composed into the containers that the process functions
  take as their `mp` argument:

- [`Microphysics0MParams`](@ref CloudMicrophysics.Parameters.Microphysics0MParams)
  wraps the 0-moment removal parameters.
- [`Microphysics1MParams`](@ref CloudMicrophysics.Parameters.Microphysics1MParams)
  stores the cloud and precipitation particle properties, air properties and
  terminal velocity parameters, plus the process options and their parameters
  (see below).
- [`Microphysics2MParams`](@ref CloudMicrophysics.Parameters.Microphysics2MParams)
  stores the Seifert & Beheng warm-rain parameters and, optionally, the
  [`P3IceParams`](@ref CloudMicrophysics.Parameters.P3IceParams) for the ice phase.
- [`TerminalVelocityParams`](@ref CloudMicrophysics.Parameters.TerminalVelocityParams)
  gathers the different terminal velocity parameterizations.
- Aerosol activation, ice nucleation and aerosol nucleation schemes take their
  parameter structs directly (e.g.
  [`AerosolActivationParameters`](@ref CloudMicrophysics.Parameters.AerosolActivationParameters),
  [`Illite`](@ref CloudMicrophysics.Parameters.Illite)).

The [Parameter reference](generated/ParametersReference.md) lists, for each
  container, the structs it is built from.

### Process options in the 1-moment scheme

`Microphysics1MParams` separates *which* variant of a process runs from the
  parameters that variant needs.
`mp.processes` is a
  [`Microphysics1MOptions`](@ref CloudMicrophysics.Parameters.Microphysics1MOptions)
  whose fields hold empty singleton option types (`Kessler1M()`,
  `PrescribedNd()`, ...) or `nothing` to switch a process off.
`mp.process_params` mirrors it field by field and holds the parameters the
  selected option needs, built by
  [`process_params_for`](@ref CloudMicrophysics.Parameters.process_params_for):

```julia
mp = CMP.Microphysics1MParams(Float64;
    rain_autoconversion = CMP.PrescribedNd(),   # Kessler1M() by default
    cloud_ice_melt = nothing,                   # disable a process
)
mp.processes.rain_autoconversion       # PrescribedNd()
mp.process_params.rain_autoconversion  # VarTimescaleAcnv{Float64}(τ = ..., α = ..., Nc = ...)
```

Process functions dispatch on the option and read their parameters from the
  matching `process_params` slot, e.g.
  `conv_q_lcl_to_q_rai(mp.processes.rain_autoconversion, mp, tps, micro, thermo)`.
Because the option is passed separately from `mp`, the two can disagree if a
  caller passes an option that `mp` was not built with.
The process functions guard against this with `Utilities.consistent_params`
  and raise an `ArgumentError` that names the option and the constructor call
  that fixes it, instead of failing with an obscure field-access error.

## Overriding parameter values

Overrides are applied when the `ParamDict` is created, so a single override
  file (or dictionary) affects every struct built from that dictionary.
For a parameter that exists in ClimaParams only the `value` needs to be given;
  `type` and `description` are taken from the defaults
  (a parameter that is not in ClimaParams yet also needs its `type`):

```julia
toml_dict = CP.create_toml_dict(Float64; override_file = "my_calibration.toml")

# or, without a file
toml_dict = CP.create_toml_dict(Float64;
    override_file = Dict("rain_autoconversion_timescale" => Dict("value" => 1200.0)),
)
mp = CMP.Microphysics1MParams(toml_dict)
```

Several override files can be layered, later files taking precedence:
  `CP.create_toml_dict(Float64, "base.toml", "experiment.toml")`.

The package ships three calibrated or literature parameter sets in
  `src/parameters/toml/` that can be used as override files:

- `CliMA_1M.toml`: calibrated 1-moment rain formation, evaporation and
  terminal velocity coefficients,
- `ARG2000.toml`: calibrated Abdul-Razzak & Ghan (2000) activation coefficients,
- `SB2006_limiters.toml`: the raindrop size distribution limiters and minimal
  raindrop mass from the original Seifert & Beheng (2006) paper, used in the
  tests that compare against the paper.

```julia
import CloudMicrophysics as CM
override_file = joinpath(pkgdir(CM), "src", "parameters", "toml", "CliMA_1M.toml")
mp = CMP.Microphysics1MParams(CP.create_toml_dict(Float64; override_file))
```

A worked example, including how the overridden value propagates to a process
  rate, is in the [Handling parameters](guides/literated/Parameters.md) guide.

### Logging and reproducibility

Because every read is tagged, ClimaParams can tell which parameters a run
  used and whether an override was never picked up (typically a typo in the
  parameter name):

```julia
CP.log_parameter_information(toml_dict, "parameters_used.toml"; strict = true)
```

This writes all used parameters (with their final values) to a TOML file that
  can be passed back as an override file to reproduce the run, and errors if
  an override was not used.

## Adding a new parameter

1. Add the entry to `parameters.toml` in ClimaParams, with `value`, `type` and
   a `description` that states the units and the source.
   Release ClimaParams and bump the `ClimaParams` compat entry in
   `Project.toml`.
   While prototyping it is fine to keep the value in a local override file
   and move it to ClimaParams as the last step.
2. Add a field with a docstring to the parameter struct, and the
   `:clima_name => :field` pair to `name_map(::Type{T})` in
   `src/parameters/name_maps.jl`.
   Prefer a long, descriptive ClimaParams name (they must be unique across all
   CliMA packages) and a short field name that matches the notation of the
   documentation.
   A new struct gets a `name_map` method and an entry in `parameter_groups()`;
   add it to `PLAIN_PARAMETER_TYPES` if its constructor needs nothing beyond
   the map, otherwise write the constructor by hand using
   `make_params(td, name_map(T))`.
3. If the parameter belongs to a new variant of a 1-moment process, add an
   option type in `src/parameters/Microphysics1MOptions.jl`, its `name_map`
   (and an entry in `NAMEDTUPLE_OPTION_TYPES`, or a hand-written
   `process_params_for` method), and a method of the process function
   dispatching on the option.

The [Parameter reference](generated/ParametersReference.md) picks up the new
  entry automatically at the next documentation build, and
  `test/parameter_name_maps_tests.jl` fails if a name map lists a parameter the
  constructor does not read, or the constructor reads one that no name map lists.
