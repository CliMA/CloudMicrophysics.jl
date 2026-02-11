export Microphysics1MParams, CloudPhaseParams1M, PrecipPhaseParams1M

"""
    CloudPhaseParams1M{FT, LCL, ICL}

Parameters for cloud-phase (non-precipitating) hydrometeors in 1-moment scheme.

# Fields
- `liquid::LCL`: CloudLiquid — cloud liquid water parameters
- `ice::ICL`: CloudIce — cloud ice parameters
"""
struct CloudPhaseParams1M{FT, LCL, ICL} <: ParametersType{FT}
    liquid::LCL
    ice::ICL
end

"""
    PrecipPhaseParams1M{FT, RAI, SNO}

Parameters for precipitating hydrometeors in 1-moment scheme.

# Fields
- `rain::RAI`: Rain — rain parameters (includes autoconversion)
- `snow::SNO`: Snow — snow parameters (includes autoconversion)
"""
struct PrecipPhaseParams1M{FT, RAI, SNO} <: ParametersType{FT}
    rain::RAI
    snow::SNO
end

"""
    Microphysics1MParams{FT, CP, PP, CE, AP, VL, VA}

Unified parameter container for 1-moment bulk microphysics.

# Fields
- `cloud::CP`: CloudPhaseParams1M — cloud liquid and ice parameters
- `precip::PP`: PrecipPhaseParams1M — rain and snow parameters
- `collision::CE`: CollisionEff — collision efficiencies between species
- `air_properties::AP`: AirProperties — air properties (diffusivities, thermal conductivity)
- `terminal_velocity::VL`: Blk1MVelType — terminal velocity parameters for rain and snow
- `autoconv_2M::VA`: VarTimescaleAcnv or Nothing — optional 2M autoconversion fallback
- `prescribed_Nc::FT`: Prescribed cloud droplet number concentration [1/m³]

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Create 1M parameters with default settings
mp = CMP.Microphysics1MParams(Float64)

# Create with optional 2M autoconversion fallback
mp_with_2M = CMP.Microphysics1MParams(Float64; with_2M_autoconv = true)
```
"""
struct Microphysics1MParams{FT, CP, PP, CE, AP, VL, VA} <: ParametersType{FT}
    cloud::CP
    precip::PP
    collision::CE
    air_properties::AP
    terminal_velocity::VL
    autoconv_2M::VA
    prescribed_Nc::FT
end

"""
    Microphysics1MParams(::Type{FT}; with_2M_autoconv = false) where {FT <: AbstractFloat}

Create a `Microphysics1MParams` object from a floating point type.

# Arguments
- `FT`: Floating point type (e.g., Float64, Float32)
- `with_2M_autoconv`: Include 2-moment autoconversion parameters (default: false)
"""
Microphysics1MParams(::Type{FT}; with_2M_autoconv = false) where {FT <: AbstractFloat} =
    Microphysics1MParams(CP.create_toml_dict(FT); with_2M_autoconv)

"""
    Microphysics1MParams(toml_dict::CP.ParamDict; with_2M_autoconv = false)

Create a `Microphysics1MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `with_2M_autoconv`: Include 2-moment autoconversion parameters (default: false)
"""
function Microphysics1MParams(toml_dict::CP.ParamDict; with_2M_autoconv = false)
    FT = CP.float_type(toml_dict)

    # Cloud phase parameters
    liquid = CloudLiquid(toml_dict)
    ice = CloudIce(toml_dict)
    cloud = CloudPhaseParams1M{FT, typeof(liquid), typeof(ice)}(liquid, ice)

    # Precipitation phase parameters
    rain = Rain(toml_dict)
    snow = Snow(toml_dict)
    precip = PrecipPhaseParams1M{FT, typeof(rain), typeof(snow)}(rain, snow)

    # Shared physics parameters
    collision = CollisionEff(toml_dict)
    air_properties = AirProperties(toml_dict)
    terminal_velocity = Blk1MVelType(toml_dict)

    # Optional 2M autoconversion
    autoconv_2M = with_2M_autoconv ? VarTimescaleAcnv(toml_dict) : nothing


    # Prescribed cloud droplet number concentration
    (; prescribed_cloud_droplet_number_concentration) =
        CP.get_parameter_values(
            toml_dict,
            "prescribed_cloud_droplet_number_concentration",
            "CloudMicrophysics",
        )
    prescribed_Nc = FT(prescribed_cloud_droplet_number_concentration)
    return Microphysics1MParams{
        FT,
        typeof(cloud),
        typeof(precip),
        typeof(collision),
        typeof(air_properties),
        typeof(terminal_velocity),
        typeof(autoconv_2M),
    }(
        cloud,
        precip,
        collision,
        air_properties,
        terminal_velocity,
        autoconv_2M,
        prescribed_Nc,
    )
end
