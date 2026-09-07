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
- `activation::AC`: [`AerosolActivationParameters`](@ref), the ARG2000 fit
- `aerosol::AE`: the aerosol population droplet activation draws on, a
  [`PrescribedAerosol`](@ref) or `nothing`

# Why the fall speeds live here

A sedimentation velocity is a microphysics parameter of the species it falls, so the warm species
carry theirs beside the processes that act on them. Before this they were reachable only from a
bundle the HOST assembled separately, which left the kernel able to compute every warm rate and
unable to compute the speed at which the result falls, and left the raindrop speed hanging off the
ice parameters, where it is neither an ice quantity nor available at all without ice.

# Constructor keyword arguments
- `aerosol`: the prescribed aerosol population. By default `nothing`, which computes zero
  activation: the supply of cloud condensation nuclei is a per-configuration statement, so a run
  that specifies it passes a [`PrescribedAerosol`](@ref) here rather than inheriting one.
"""
@kwdef struct WarmRainParams2M{SB, AP, CE, SD, CV, RV, AC, AE} <: ParametersType
    seifert_beheng::SB
    air_properties::AP
    condevap::CE
    subdep::SD
    cloud_velocity::CV
    rain_velocity::RV
    activation::AC
    aerosol::AE = nothing
end
# Construct WarmRainParams2M from a ClimaParams TOML dictionary
WarmRainParams2M(toml_dict::CP.ParamDict; is_limited = true, aerosol = nothing,
    rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited)) =
    WarmRainParams2M(;
        seifert_beheng = SB2006(toml_dict; is_limited, rain_pdf),
        air_properties = AirProperties(toml_dict),
        condevap = CondEvap2M(toml_dict),
        subdep = SubDep2M(toml_dict),
        cloud_velocity = StokesRegimeVelType(toml_dict),
        rain_velocity = Chen2022VelType(toml_dict).rain,
        activation = AerosolActivationParameters(toml_dict),
        aerosol,
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
    quad = Quadrature.GaussLegendre(FT, quadrature_order), liqice_partition = nothing,
    inp_depletion_model)
```
which constructs the parameterization with components:
- `scheme` = [`ParametersP3`](@ref), built with `slope_law` and `aspect_ratio`
- `terminal_velocity` = [`Chen2022VelType`](@ref)
- `cloud_pdf` = [`CloudParticlePDF_SB2006`](@ref)
- `rain_pdf` = [`RainParticlePDF_SB2006`](@ref)
- `ice_nucleation` = [`ExponentialSupercoolingINP`](@ref), the Cooper (1986) deposition
  target spectrum; [`Frostenberg2023`](@ref) is selectable as a non-default option
- `rain_freezing` = [`RainFreezing`](@ref)
- `homogeneous` = [`Koop2000`](@ref)

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
- `liqice_partition`: the liquid-ice collision closure (default: `nothing`, the per-particle
  partition). Pass `P3Scheme.BulkPartition()` for the reference P3 code's bulk partition.

"""
@kwdef struct P3IceParams{P3, VL, PDc, PDr, HET, RF, HOM, INPDM, Q, LIA} <: ParametersType
    "The core P3 scheme parameters"
    scheme::P3
    "The terminal velocity parameterization"
    terminal_velocity::VL
    "The cloud droplet size distribution"
    cloud_pdf::PDc
    "The rain drop size distribution"
    rain_pdf::PDr
    "The deposition nucleation target spectrum, an [`AbstractINPTargetSpectrum`](@ref)"
    ice_nucleation::HET
    "The rain freezing parameters (Bigg-type immersion freezing), also used for the
    cloud-droplet immersion line"
    rain_freezing::RF
    "The homogeneous freezing parameters, [`Koop2000`](@ref), composed alongside the
    heterogeneous coefficient in [`HetIceNucleation.rain_freezing_rate`](@ref) and
    [`HetIceNucleation.cloud_freezing_rate`](@ref)"
    homogeneous::HOM
    "Depletion proxy model for the deposition nucleation target. Currently only
    [`NIceProxyDepletion`](@ref) (legacy n_ice-as-proxy form) is provided; it sets the
    value subtracted from the target concentration in
    [`HetIceNucleation.deposition_rate`](@ref). (A prognostic activation-memory model is
    deferred to a follow-up PR.)"
    inp_depletion_model::INPDM = NIceProxyDepletion()
    "Quadrature rule for the size-distribution integrals
    (deposition / sublimation, melting, riming, ice-rain collection,
    sedimentation). See also [`Quadrature.GaussLegendre`](@ref).
    A [`P3Scheme.P3TabulatedQuadrature`](@ref) may be passed instead, to replace
    selected size-distribution integrals with a lookup table."
    quad::Q = QUAD.GaussLegendre(Float64, 6)
    "How the liquid-ice collision entry partitions collected mass between freezing and
    shedding, or `nothing` to let the quadrature decide. `nothing` selects a per-particle
    partition, [`P3Scheme.PartitionedOuter`](@ref) with a plain rule and
    [`P3Scheme.SplitCorrection`](@ref) with a carrier holding the liquid-ice tables, which are
    two forms of one closure. [`P3Scheme.BulkPartition`](@ref) selects the reference P3 code's
    partition instead, keeping this scheme's own freezing capacity: it compares the collected mass
    with that capacity once for the
    whole population rather than at every ice diameter. That is a choice of physics rather than of
    numerics, and is made
    here rather than through `quad`, which carries the numerics."
    liqice_partition::LIA = nothing
end
Base.show(io::IO, mime::MIME"text/plain", x::P3IceParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

# `quad` may be a `P3Scheme.P3TabulatedQuadrature`, which holds device arrays. It reaches a
# kernel inside this struct rather than as a broadcast argument of its own, captured by
# closures at many call sites; `Adapt` returns a struct with no rule unchanged and does not
# visit its interior, so every layer between `ClimaAtmosParameters` and the carrier needs one.
# Adapting a rule with no rule of its own is the identity, so this is free for every other
# quadrature.
Adapt.@adapt_structure P3IceParams

P3IceParams(toml_dict::CP.ParamDict;
    is_limited = true,
    quadrature_order = 6,
    quad = QUAD.GaussLegendre(CP.float_type(toml_dict), quadrature_order),
    liqice_partition = nothing,
    inp_depletion_model = NIceProxyDepletion(),
    slope_law = DEFAULT_SLOPE_LAW,
    aspect_ratio = DEFAULT_ASPECT_RATIO,
    rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited),
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
    rain_pdf,
    ice_nucleation = ExponentialSupercoolingINP(toml_dict),
    rain_freezing = RainFreezing(toml_dict),
    homogeneous = Koop2000(toml_dict),
    inp_depletion_model,
    quad,
    liqice_partition,
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

# The middle link of the adaptation chain described at `P3IceParams` above.
Adapt.@adapt_structure Microphysics2MParams

"""
    Microphysics2MParams(toml_dict::CP.ParamDict; with_ice = false, is_limited = true,
        slope_law = DEFAULT_SLOPE_LAW, aspect_ratio = DEFAULT_ASPECT_RATIO,
        quadrature_order = 6, quad = Quadrature.GaussLegendre(FT, quadrature_order),
        liqice_partition = nothing, inp_depletion_model, aerosol = nothing)

Create a `Microphysics2MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `with_ice`: Include P3 ice-phase parameters (default: false)
- `is_limited`: Use limited rain size distribution parameters (default: true)
- `inp_depletion_model`: the F23 INP-activation depletion model passed to
  [`P3IceParams`](@ref) when `with_ice`. By default, [`NIceProxyDepletion`](@ref).
- `slope_law`, `aspect_ratio`: passed to [`P3IceParams`](@ref) when `with_ice`, which
  forwards them to [`ParametersP3`](@ref). Defaults are [`DEFAULT_SLOPE_LAW`](@ref)
  and [`DEFAULT_ASPECT_RATIO`](@ref), named once there.
- `quadrature_order`: order of the default `Quadrature.GaussLegendre` rule passed to
  [`P3IceParams`](@ref) when `with_ice` (default: 6)
- `quad`: the size-distribution `Quadrature.QuadratureRule` passed to
  [`P3IceParams`](@ref) when `with_ice` (default: `Quadrature.GaussLegendre(FT, quadrature_order)`)
- `liqice_partition`: the liquid-ice collision closure passed to [`P3IceParams`](@ref) when
  `with_ice` (default: `nothing`, the per-particle partition)
- `aerosol`: the aerosol population droplet activation draws on, passed to
  [`WarmRainParams2M`](@ref). By default `nothing`, which computes zero activation; a
  configuration that specifies cloud condensation nuclei passes a [`PrescribedAerosol`](@ref).
"""
Microphysics2MParams(toml_dict::CP.ParamDict;
    with_ice = false, is_limited = true, aerosol = nothing,
    quadrature_order = 6,
    quad = QUAD.GaussLegendre(CP.float_type(toml_dict), quadrature_order),
    liqice_partition = nothing,
    inp_depletion_model = NIceProxyDepletion(),
    slope_law = DEFAULT_SLOPE_LAW,
    aspect_ratio = DEFAULT_ASPECT_RATIO,
    rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited),
) = Microphysics2MParams(;
    # One `rain_pdf` object reaches both halves rather than each building its own from the
    # same TOML. They consume the same size-distribution inversion, so they must not be able
    # to disagree, and a caller selecting `RainParticlePDF_SB2006_windowed` must reach both.
    warm_rain = WarmRainParams2M(toml_dict; is_limited, aerosol, rain_pdf),
    # Optional ice phase parameters
    ice = with_ice ?
          P3IceParams(toml_dict; is_limited, quad, liqice_partition, inp_depletion_model,
        slope_law, aspect_ratio, rain_pdf) :
          nothing,
)
