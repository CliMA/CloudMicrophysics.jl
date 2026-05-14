export Microphysics1MParams, CloudPhaseParams1M, PrecipPhaseParams1M

"""
    CloudPhaseParams1M{FT, LCL, ICL}

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
    PrecipPhaseParams1M{FT, RAI, SNO}

Parameters for precipitating hydrometeors in 1-moment scheme.

# Fields
- `rain::RAI`: Rain — rain parameters (includes autoconversion)
- `snow::SNO`: Snow — snow parameters (includes autoconversion)
"""
@kwdef struct PrecipPhaseParams1M{RAI, SNO} <: ParametersType
    rain::RAI
    snow::SNO
end

Base.show(io::IO, mime::MIME"text/plain", x::PrecipPhaseParams1M) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    Microphysics1MParams{FT, OPT, CP, PP, CE, AP, VL, VA, FR}

Unified parameter container for 1-moment bulk microphysics.

# Fields
- `options::OPT`: Microphysics1MOptions — process configuration (selects parameterization variants)
- `cloud::CP`: CloudPhaseParams1M — cloud liquid and ice parameters
- `precip::PP`: PrecipPhaseParams1M — rain and snow parameters
- `collision::CE`: CollisionEff — collision efficiencies between species
- `air_properties::AP`: AirProperties — air properties (diffusivities, thermal conductivity)
- `terminal_velocity::VL`: Blk1MVelType — terminal velocity parameters for rain and snow
- `autoconv_2M::VA`: VarTimescaleAcnv or Nothing — 2M autoconversion parameters (built when `LiquidAutoconv2M` is selected)
- `prescribed_Nc::FT`: Prescribed cloud droplet number concentration [1/m³]
- `frostenberg2023::FR`: Frostenberg 2023 INP parameters or Nothing (built when `INPDependentIceFormation` is selected)

# Constructors

    Microphysics1MParams(::Type{FT}; options = Microphysics1MOptions())
    Microphysics1MParams(toml_dict::CP.ParamDict; options = Microphysics1MOptions())

Create a `Microphysics1MParams` from a float type or a ClimaParams TOML dictionary.
The `options` keyword argument selects parameterization variants; sub-parameters
(`autoconv_2M`, `frostenberg2023`) are conditionally built based on the chosen options.

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Create 1M parameters with default (baseline) settings
mp = CMP.Microphysics1MParams(Float64)

# Create with INP-dependent ice formation and cloud ice melt
mp = CMP.Microphysics1MParams(Float64;
    options = CMP.Microphysics1MOptions(
        cloud_ice_formation  = CMP.TemperatureDependentCloudIceFormation(),
        cloud_ice_melt = CMP.CloudIceMeltToLiquid(),
    ),
)

# Create with 2M autoconversion
mp = CMP.Microphysics1MParams(Float64;
    options = CMP.Microphysics1MOptions(
        cloud_liquid_autoconversion = CMP.LiquidAutoconv2M(),
    ),
)
```
"""
@kwdef struct Microphysics1MParams{FT, OPT, CP, PP, CE, AP, VL, VA, FR} <: ParametersType
    options::OPT
    cloud::CP
    precip::PP
    collision::CE
    air_properties::AP
    terminal_velocity::VL
    autoconv_2M::VA
    prescribed_Nc::FT
    frostenberg2023::FR
end
Base.show(io::IO, mime::MIME"text/plain", x::Microphysics1MParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    Microphysics1MParams(toml_dict::CP.ParamDict; options = Microphysics1MOptions())

Create a `Microphysics1MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `options`: Process configuration (default: baseline behavior)
"""
function Microphysics1MParams(toml_dict::CP.ParamDict; options = Microphysics1MOptions())
    (; prescribed_cloud_droplet_number_concentration) = CP.get_parameter_values(
        toml_dict, "prescribed_cloud_droplet_number_concentration", "CloudMicrophysics",
    )
    # Conditional construction driven by options
    autoconv_2M = options.cloud_liquid_autoconversion isa LiquidAutoconv2M ?
                  VarTimescaleAcnv(toml_dict) : nothing
    frostenberg2023 =
        options.cloud_ice_formation isa TemperatureDependentCloudIceFormation ?
        Frostenberg2023(toml_dict) : nothing

    return Microphysics1MParams(;
        # Process options
        options,
        # Cloud phase parameters
        cloud = CloudPhaseParams1M(;
            liquid = CloudLiquid(toml_dict),
            ice = CloudIce(toml_dict),
        ),
        # Precipitation phase parameters
        precip = PrecipPhaseParams1M(;
            rain = Rain(toml_dict),
            snow = Snow(toml_dict),
        ),
        # Shared physics parameters
        collision = CollisionEff(toml_dict),
        air_properties = AirProperties(toml_dict),
        terminal_velocity = Blk1MVelType(toml_dict),
        # Conditionally built parameters
        autoconv_2M,
        # Prescribed cloud droplet number concentration
        prescribed_Nc = prescribed_cloud_droplet_number_concentration,
        frostenberg2023,
    )
end
