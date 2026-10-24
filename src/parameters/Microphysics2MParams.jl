export Microphysics2MParams, WarmRainParams2M, P3IceParams

"""
    WarmRainParams2M

Parameters for 2-moment warm rain processes (Seifert-Beheng 2006).

# Fields
- `seifert_beheng::SB`: SB2006 — all warm rain parameters (autoconversion, accretion, etc.)
- `air_properties::AP`: AirProperties — air properties for evaporation
- `condevap::CE`: MM2015 cond-evap relaxation timescale
- `subdep::SD`: MM2015 sub-dep relaxation timescale
"""
@kwdef struct WarmRainParams2M{SB, AP, CE, SD} <: ParametersType
    seifert_beheng::SB
    air_properties::AP
    condevap::CE
    subdep::SD
end
# Construct WarmRainParams2M from a ClimaParams TOML dictionary
WarmRainParams2M(toml_dict::CP.ParamDict; is_limited = true) =
    WarmRainParams2M(;
        seifert_beheng = SB2006(toml_dict; is_limited),
        air_properties = AirProperties(toml_dict),
        condevap = CondEvap2M(toml_dict),
        subdep = SubDep2M(toml_dict),
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
P3IceParams(toml_dict::CP.ParamDict; is_limited = true, quadrature_order = 16,
    integration = Quadrature.GaussLegendre(FT, quadrature_order),
    ice_nucleation = Frostenberg2023(toml_dict),
    inp_depletion_model = NIceProxyDepletion(), freezing_quadrature_order = 8)
```
which constructs the parameterization with components:
- `scheme` = [`ParametersP3`](@ref)
- `terminal_velocity` = [`Chen2022VelType`](@ref)
- `cloud_pdf` = [`CloudParticlePDF_SB2006`](@ref)
- `rain_pdf` = [`RainParticlePDF_SB2006`](@ref)
- `rain_freezing` = [`RainFreezing`](@ref)
- `homogeneous` = [`Koop2000`](@ref)
- `dendritic_freezing` = [`DendriticFreezing`](@ref)

# Keyword arguments
- `is_limited`: use limited rain size-distribution parameters (default: true)
- `quadrature_order`: order of the default `Quadrature.GaussLegendre` rule (default: 16)
- `integration`: how the size-distribution integrals are evaluated, a `Quadrature.QuadratureRule`
  (default: `Quadrature.GaussLegendre(FT, quadrature_order)`). Pass a quadrature rule to use a rule
  other than Gauss-Legendre; for lookup tables, build the parameters with `P3Scheme.with_tables`.
  Tables built from other parameters throw an `ArgumentError`, see [`check_integration`](@ref).
- `ice_nucleation`: the deposition nucleation target spectrum, an [`INPTargetSpectrum`](@ref):
  [`Frostenberg2023`](@ref) (the default), the Cooper (1986) spectrum [`ExponentialSupercoolingINP`](@ref), or a
  [`StochasticINPSpectrum`](@ref)
- `inp_depletion_model`: the INP depletion model (default: `NIceProxyDepletion()`)
- `freezing_quadrature_order`: number of nodes of the [`Quadrature.GaussLaguerre`](@ref) rules of the rain and
  cloud droplet freezing integrals (default: 8)

"""
@kwdef struct P3IceParams{P3, VL, PDc, PDr, HET, RF, HOM, DF, INPDM, Q, QR, QC} <: ParametersType
    "The core P3 scheme parameters"
    scheme::P3
    "The terminal velocity parameterization"
    terminal_velocity::VL
    "The cloud droplet size distribution"
    cloud_pdf::PDc
    "The rain drop size distribution"
    rain_pdf::PDr
    "The deposition nucleation target spectrum, an [`INPTargetSpectrum`](@ref)"
    ice_nucleation::HET
    "The rain freezing parameters (Bigg-type immersion freezing)"
    rain_freezing::RF
    "The homogeneous freezing parameters, [`Koop2000`](@ref)"
    homogeneous::HOM
    "The growth of the freezing front across a nucleated drop, [`DendriticFreezing`](@ref)"
    dendritic_freezing::DF
    "Model for F23 INP-activation depletion. Currently only
    [`NIceProxyDepletion`](@ref) (n_ice-as-proxy form) is provided;
    it sets the value subtracted from `INPC(T)/ρ` in the F23 deposition +
    immersion-cap rates."
    inp_depletion_model::INPDM = NIceProxyDepletion()
    "How the size-distribution integrals are evaluated, a quadrature rule such as [`Quadrature.GaussLegendre`](@ref), or lookup tables, [`P3Scheme.P3Tables`](@ref CloudMicrophysics.P3Scheme.P3Tables)"
    integration::Q = QUAD.GaussLegendre(Float64, 16)
    "Quadrature rule of the raindrop freezing integrals, a [`Quadrature.GaussLaguerre`](@ref) rule in the diameter
    over the mean diameter"
    rain_freezing_quad::QR
    "Quadrature rule of the cloud droplet freezing integrals, a [`Quadrature.GaussLaguerre`](@ref) rule in
    `λc D^μcD` with `α = (νc + 1)/μc - 1`"
    cloud_freezing_quad::QC
end
Base.show(io::IO, mime::MIME"text/plain", x::P3IceParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
Adapt.@adapt_structure P3IceParams

function P3IceParams(toml_dict::CP.ParamDict;
    is_limited = true,
    quadrature_order = 16,
    integration = QUAD.GaussLegendre(CP.float_type(toml_dict), quadrature_order),
    ice_nucleation = Frostenberg2023(toml_dict),
    inp_depletion_model = NIceProxyDepletion(),
    freezing_quadrature_order = 8,
)
    FT = CP.float_type(toml_dict)
    cloud_pdf = CloudParticlePDF_SB2006(toml_dict)
    (; νc, μc) = cloud_pdf
    ice_params = P3IceParams(;
        scheme = ParametersP3(toml_dict),
        terminal_velocity = Chen2022VelType(toml_dict),
        cloud_pdf,
        rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited),
        ice_nucleation,
        rain_freezing = RainFreezing(toml_dict),
        homogeneous = Koop2000(toml_dict),
        dendritic_freezing = DendriticFreezing(toml_dict),
        inp_depletion_model,
        integration,
        rain_freezing_quad = QUAD.GaussLaguerre(FT, freezing_quadrature_order),
        cloud_freezing_quad = QUAD.GaussLaguerre(FT, freezing_quadrature_order; α = (νc + 1) / μc - 1),
    )
    check_integration(integration, ice_params)
    return ice_params
end

"""
    check_integration(integration, ice_params)

Throw an `ArgumentError` if `integration` cannot evaluate the size-distribution integrals of the [`P3IceParams`](@ref) `ice_params`.

A quadrature rule evaluates them for any parameters. Lookup tables evaluate them only for the
parameters they were built from, see `P3Scheme.with_tables`.
"""
check_integration(integration, ice_params) = nothing

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
Adapt.@adapt_structure Microphysics2MParams

"""
    Microphysics2MParams(toml_dict::CP.ParamDict; with_ice = false, is_limited = true,
        quadrature_order = 16, integration = Quadrature.GaussLegendre(FT, quadrature_order),
        ice_nucleation = Frostenberg2023(toml_dict),
        inp_depletion_model = NIceProxyDepletion())

Create a `Microphysics2MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `with_ice`: Include P3 ice-phase parameters (default: false)
- `is_limited`: Use limited rain size distribution parameters (default: true)
- `quadrature_order`: order of the default `Quadrature.GaussLegendre` rule passed to
  [`P3IceParams`](@ref) when `with_ice` (default: 16)
- `integration`: the evaluation of the size-distribution integrals passed to
  [`P3IceParams`](@ref) when `with_ice` (default: `Quadrature.GaussLegendre(FT, quadrature_order)`)
- `ice_nucleation`: the deposition nucleation target spectrum passed to [`P3IceParams`](@ref)
  when `with_ice` (default: `Frostenberg2023(toml_dict)`)
- `inp_depletion_model`: the INP depletion model passed to [`P3IceParams`](@ref) when
  `with_ice` (default: `NIceProxyDepletion()`)
"""
Microphysics2MParams(toml_dict::CP.ParamDict;
    with_ice = false, is_limited = true,
    quadrature_order = 16,
    integration = QUAD.GaussLegendre(CP.float_type(toml_dict), quadrature_order),
    ice_nucleation = Frostenberg2023(toml_dict),
    inp_depletion_model = NIceProxyDepletion(),
) = Microphysics2MParams(;
    # Warm rain parameters (always present)
    warm_rain = WarmRainParams2M(toml_dict; is_limited),
    # Optional ice phase parameters
    ice = with_ice ?
          P3IceParams(toml_dict; is_limited, integration, ice_nucleation, inp_depletion_model) :
          nothing,
)
