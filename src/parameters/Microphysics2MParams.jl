export Microphysics2MParams, WarmRainParams2M, P3IceParams

"""
    WarmRainParams2M

Parameters for 2-moment warm rain processes (Seifert-Beheng 2006).

# Fields
- `seifert_beheng::SB`: SB2006 — all warm rain parameters (autoconversion, accretion, etc.)
- `air_properties::AP`: AirProperties — air properties for evaporation
- `condevap::CE`: MM2015 cond-evap relaxation timescale
- `subdep::SD`: MM2015 sub-dep relaxation timescale
- `cloud_velocity::CV`: [`StokesRegimeVelType`](@ref), the cloud droplet fall speed
- `rain_velocity::RV`: [`Chen2022VelTypeRain`](@ref), the raindrop fall speed

# Why the fall speeds live here

A sedimentation velocity is a microphysics parameter of the species it falls, so the warm species
carry theirs beside the processes that act on them. Before this they were reachable only from a
bundle the HOST assembled separately, which left the kernel able to compute every warm rate and
unable to compute the speed at which the result falls, and left the raindrop speed hanging off the
ice parameters, where it is neither an ice quantity nor available at all without ice.
"""
@kwdef struct WarmRainParams2M{SB, AP, CE, SD, CV, RV} <: ParametersType
    seifert_beheng::SB
    air_properties::AP
    condevap::CE
    subdep::SD
    cloud_velocity::CV
    rain_velocity::RV
end
# Construct WarmRainParams2M from a ClimaParams TOML dictionary
WarmRainParams2M(toml_dict::CP.ParamDict; is_limited = true) =
    WarmRainParams2M(;
        seifert_beheng = SB2006(toml_dict; is_limited),
        air_properties = AirProperties(toml_dict),
        condevap = CondEvap2M(toml_dict),
        subdep = SubDep2M(toml_dict),
        cloud_velocity = StokesRegimeVelType(toml_dict),
        rain_velocity = Chen2022VelType(toml_dict).rain,
    )

Base.show(io::IO, mime::MIME"text/plain", x::WarmRainParams2M) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    P3IceParams

Parameters for P3 ice-phase processes.

# Fields
$(DocStringExtensions.FIELDS)

# Constructor

The main constructor is
```
P3IceParams(toml_dict::CP.ParamDict; is_limited = true, slope_law = DEFAULT_SLOPE_LAW,
    aspect_ratio = DEFAULT_ASPECT_RATIO, quadrature_order = 6,
    quad = Quadrature.GaussLegendre(FT, quadrature_order))
```
which constructs the parameterization with components:
- `scheme` = [`ParametersP3`](@ref), built with `slope_law` and `aspect_ratio`
- `terminal_velocity` = [`Chen2022VelType`](@ref)
- `cloud_pdf` = [`CloudParticlePDF_SB2006`](@ref)
- `rain_pdf` = [`RainParticlePDF_SB2006`](@ref)
- `ice_nucleation` = [`Frostenberg2023`](@ref)
- `rain_freezing` = [`RainFreezing`](@ref)

# Keyword arguments
- `is_limited`: use limited rain size-distribution parameters (default: true)
- `slope_law`, `aspect_ratio`: forwarded to [`ParametersP3`](@ref). Defaults are
  [`DEFAULT_SLOPE_LAW`](@ref) and [`DEFAULT_ASPECT_RATIO`](@ref), named once there.
- `quadrature_order`: order of the default `Quadrature.GaussLegendre` rule (default: 6). The
  order was set by a three-level error study against a 128-node reference over seventeen regime
  states: the transport components converge to at most 0.45 percent at order 6, and every
  larger error, at most 1.4 percent, sits in a component proportional to the collision
  efficiency, whose `E = 1` parameterization uncertainty exceeds it by an order of magnitude.
  Order 5 exceeds 2 percent in the transport family. The two dominant ice kernels are nested
  quadratures and cost the square of the order.
- `quad`: the size-distribution `Quadrature.QuadratureRule` (default:
  `Quadrature.GaussLegendre(FT, quadrature_order)`). Pass this to use a rule other than
  Gauss-Legendre.

"""
@kwdef struct P3IceParams{P3, VL, PDc, PDr, HET, RF, INPDM, Q} <: ParametersType
    "The core P3 scheme parameters"
    scheme::P3
    "The terminal velocity parameterization"
    terminal_velocity::VL
    "The cloud droplet size distribution"
    cloud_pdf::PDc
    "The rain drop size distribution"
    rain_pdf::PDr
    "The ice nucleation parameters (empirical INP closure)"
    ice_nucleation::HET
    "The rain freezing parameters (Bigg-type immersion freezing)"
    rain_freezing::RF
    "Model for F23 INP-activation depletion. Currently only
    [`NIceProxyDepletion`](@ref) (n_ice-as-proxy form) is provided;
    it sets the value subtracted from `INPC(T)/ρ` in the F23 deposition +
    immersion-cap rates."
    inp_depletion_model::INPDM = NIceProxyDepletion()
    "Quadrature rule for the size-distribution integrals
    (deposition / sublimation, melting, riming, ice-rain collection,
    sedimentation). See also [`Quadrature.GaussLegendre`](@ref)."
    quad::Q = QUAD.GaussLegendre(Float64, 6)
end
Base.show(io::IO, mime::MIME"text/plain", x::P3IceParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

P3IceParams(toml_dict::CP.ParamDict;
    is_limited = true,
    quadrature_order = 6,
    quad = QUAD.GaussLegendre(CP.float_type(toml_dict), quadrature_order),
    inp_depletion_model = NIceProxyDepletion(τ_act = 300),
    slope_law = DEFAULT_SLOPE_LAW,
    aspect_ratio = DEFAULT_ASPECT_RATIO,
) = P3IceParams(;
    # Forwarded rather than left at the `ParametersP3` default, because a
    # `ParametersP3` built here is the only one a host model ever sees. Without
    # these an alternative slope law could not be selected from ClimaAtmos at
    # all, while `P3_mu_smoothing_sharpness` ships as a tunable TOML key, which
    # is a configurable law with no way to configure it. Defaults named once, at
    # [`DEFAULT_SLOPE_LAW`](@ref) and [`DEFAULT_ASPECT_RATIO`](@ref).
    scheme = ParametersP3(toml_dict; slope_law, aspect_ratio),
    terminal_velocity = Chen2022VelType(toml_dict),
    cloud_pdf = CloudParticlePDF_SB2006(toml_dict),
    rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited),
    ice_nucleation = Frostenberg2023(toml_dict),
    rain_freezing = RainFreezing(toml_dict),
    inp_depletion_model,
    quad,
)

"""
    Microphysics2MParams{WR, ICE}

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
@kwdef struct Microphysics2MParams{WR, ICE} <: ParametersType
    warm_rain::WR
    ice::ICE
end
Base.show(io::IO, mime::MIME"text/plain", x::Microphysics2MParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

"""
    Microphysics2MParams(toml_dict::CP.ParamDict; with_ice = false, is_limited = true,
        slope_law = DEFAULT_SLOPE_LAW, aspect_ratio = DEFAULT_ASPECT_RATIO,
        quadrature_order = 6, quad = Quadrature.GaussLegendre(FT, quadrature_order))

Create a `Microphysics2MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `with_ice`: Include P3 ice-phase parameters (default: false)
- `is_limited`: Use limited rain size distribution parameters (default: true)
- `slope_law`, `aspect_ratio`: passed to [`P3IceParams`](@ref) when `with_ice`, which
  forwards them to [`ParametersP3`](@ref). Defaults are [`DEFAULT_SLOPE_LAW`](@ref)
  and [`DEFAULT_ASPECT_RATIO`](@ref), named once there.
- `quadrature_order`: order of the default `Quadrature.GaussLegendre` rule passed to
  [`P3IceParams`](@ref) when `with_ice` (default: 6)
- `quad`: the size-distribution `Quadrature.QuadratureRule` passed to
  [`P3IceParams`](@ref) when `with_ice` (default: `Quadrature.GaussLegendre(FT, quadrature_order)`)
"""
Microphysics2MParams(toml_dict::CP.ParamDict;
    with_ice = false, is_limited = true,
    quadrature_order = 6,
    quad = QUAD.GaussLegendre(CP.float_type(toml_dict), quadrature_order),
    inp_depletion_model = NIceProxyDepletion(τ_act = 300),
    slope_law = DEFAULT_SLOPE_LAW,
    aspect_ratio = DEFAULT_ASPECT_RATIO,
) = Microphysics2MParams(;
    # Warm rain parameters (always present)
    warm_rain = WarmRainParams2M(toml_dict; is_limited),
    # Optional ice phase parameters
    ice = with_ice ?
          P3IceParams(toml_dict; is_limited, quad, inp_depletion_model,
        slope_law, aspect_ratio) :
          nothing,
)
