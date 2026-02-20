export Microphysics2MParams, WarmRainParams2M, P3IceParams

"""
    WarmRainParams2M{FT, SB, AP}

Parameters for 2-moment warm rain processes (Seifert-Beheng 2006).

# Fields
- `seifert_beheng::SB`: SB2006 — all warm rain parameters (autoconversion, accretion, etc.)
- `air_properties::AP`: AirProperties — air properties for evaporation
"""
@kwdef struct WarmRainParams2M{FT, SB, AP} <: ParametersType{FT}
    seifert_beheng::SB
    air_properties::AP
end

function WarmRainParams2M(toml_dict::CP.ParamDict; is_limited = true)
    FT = CP.float_type(toml_dict)
    seifert_beheng = SB2006(toml_dict, is_limited)
    air_properties = AirProperties(toml_dict)
    return WarmRainParams2M{FT}(; seifert_beheng, air_properties)
end

"""
    P3IceParams{FT, P3, VL, PDc, PDr}

Parameters for P3 ice-phase processes (optional).

# Fields
- `scheme::P3`: ParametersP3 — P3 scheme parameters
- `terminal_velocity::VL`: Chen2022VelType — terminal velocity for ice
- `cloud_pdf::PDc`: CloudParticlePDF_SB2006 — cloud droplet size distribution
- `rain_pdf::PDr`: RainParticlePDF_SB2006 — rain drop size distribution
- `heterogeneous::HET`: Frostenberg2023 — heterogeneous ice nucleation parameters
- `deposition_condfreeze::DEP`: MorrisonMilbrandt2014 — deposition ice nucleation and
    condensation freezing parameters
"""
@kwdef struct P3IceParams{FT, P3, VL, PDc, PDr, HET, DEP} <: ParametersType{FT}
    scheme::P3
    terminal_velocity::VL
    cloud_pdf::PDc
    rain_pdf::PDr
    heterogeneous::HET
    deposition_condfreeze::DEP  # TODO: Combine or split these structs???
end

P3IceParams(toml_dict::CP.ParamDict; is_limited = true) = 
    P3IceParams{CP.float_type(toml_dict)}(; 
        scheme = ParametersP3(toml_dict), 
        terminal_velocity = Chen2022VelType(toml_dict), 
        cloud_pdf = CloudParticlePDF_SB2006(toml_dict), 
        rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited), 
        heterogeneous = Frostenberg2023(toml_dict),
        deposition_condfreeze = MorrisonMilbrandt2014(toml_dict),
    )

"""
    Microphysics2MParams{FT, WR, ICE}

Unified parameter container for 2-moment microphysics.

Supports:
- **Warm rain only** (SB2006): when `ice` is `nothing`
- **Warm rain + P3 ice**: when `ice` is `P3IceParams`

# Fields
- `warm_rain::WR`: WarmRainParams2M — SB2006 parameters + air properties
- `ice::ICE`: P3IceParams or Nothing — optional P3 ice parameters

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Warm rain only
mp_warm = CMP.Microphysics2MParams(Float64; with_ice = false)

# Warm rain + P3 ice
mp_p3 = CMP.Microphysics2MParams(Float64; with_ice = true)
```
"""
@kwdef struct Microphysics2MParams{FT, WR, ICE} <: ParametersType{FT}
    warm_rain::WR
    ice::ICE
end

"""
    Microphysics2MParams(::Type{FT}; with_ice = false, is_limited = true)

Create a `Microphysics2MParams` object from a floating point type.

# Arguments
- `FT`: Floating point type (e.g., `Float64`, `Float32`)
- `with_ice`: Include P3 ice-phase parameters (default: `false`)
- `is_limited`: Use limited rain size distribution parameters (default: `true`)
"""
Microphysics2MParams(::Type{FT}; with_ice = false, is_limited = true) where {FT} =
    Microphysics2MParams(CP.create_toml_dict(FT); with_ice, is_limited)

"""
    Microphysics2MParams(toml_dict::CP.ParamDict; with_ice = false, is_limited = true)

Create a `Microphysics2MParams` object from a `ClimaParams` TOML dictionary.

# Arguments
- `toml_dict`: `ClimaParams` parameter dictionary
- `with_ice`: Include P3 ice-phase parameters (default: `false`)
- `is_limited`: Use limited rain size distribution parameters (default: `true`)
"""
function Microphysics2MParams(toml_dict::CP.ParamDict; with_ice = false, is_limited = true)
    FT = CP.float_type(toml_dict)

    # Warm rain parameters (always present)
    warm_rain = WarmRainParams2M(toml_dict; is_limited)

    # Optional ice phase parameters
    ice = with_ice ? P3IceParams(toml_dict; is_limited) : nothing

    return Microphysics2MParams{FT}(; warm_rain, ice)
end
