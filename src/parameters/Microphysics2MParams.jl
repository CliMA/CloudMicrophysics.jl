export Microphysics2MParams, WarmRainParams2M, P3IceParams

"""
    WarmRainParams2M

Parameters for 2-moment warm rain processes (Seifert-Beheng 2006).

# Fields
- `seifert_beheng::SB`: SB2006 — all warm rain parameters (autoconversion, accretion, etc.)
- `air_properties::AP`: AirProperties — air properties for evaporation
- `condevap::CE`: MM2015 cond-evap relaxation timescale
- `subdep::SD`: MM2015 sub-dep relaxation timescale
- `activation_scheme::AS`: aerosol → cloud-droplet activation scheme
  (defaults to [`NoActivation`](@ref); set to e.g. `DiagnosticNc`,
  `TwomeyActivation`, or `FixedARGActivation` to enable activation in
  the fused 2M tendency).
"""
@kwdef struct WarmRainParams2M{SB, AP, CE, SD, AS <: AbstractActivationScheme} <: ParametersType
    seifert_beheng::SB
    air_properties::AP
    condevap::CE
    subdep::SD
    activation_scheme::AS = NoActivation()
end
# Construct WarmRainParams2M from a ClimaParams TOML dictionary
WarmRainParams2M(toml_dict::CP.ParamDict; is_limited = true,
                 activation_scheme::AbstractActivationScheme = NoActivation()) =
    WarmRainParams2M(;
        seifert_beheng = SB2006(toml_dict; is_limited),
        air_properties = AirProperties(toml_dict),
        condevap = CondEvap2M(toml_dict),
        subdep = SubDep2M(toml_dict),
        activation_scheme,
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
P3IceParams(toml_dict::CP.ParamDict; is_limited = true)
```
which constructs the parameterization with components:
- `scheme` = [`ParametersP3`](@ref)
- `terminal_velocity` = [`Chen2022VelType`](@ref)
- `cloud_pdf` = [`CloudParticlePDF_SB2006`](@ref)
- `rain_pdf` = [`RainParticlePDF_SB2006`](@ref)
- `ice_nucleation` = [`Frostenberg2023`](@ref)
- `rain_freezing` = [`RainFreezing`](@ref)

"""
@kwdef struct P3IceParams{P3, VL, PDc, PDr, HET, RF, INPDM} <: ParametersType
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
    "Model for F23 INP-activation memory: choose
    [`NIceProxyDepletion`](@ref) (default; legacy n_ice-as-proxy form)
    or [`PrognosticINPDecay`](@ref) (Phillips-style activation memory
    via a host-supplied `n_INP_used` prognostic). Selects what value
    is subtracted from `INPC(T)/ρ` in the F23 deposition + immersion-
    cap rates, and whether the host needs to maintain a `n_INP_used`
    tracer."
    inp_depletion_model::INPDM = NIceProxyDepletion()
    "Number of Chebyshev–Gauss nodes used for size-distribution integrals
    (deposition / sublimation, melting, riming, ice-rain collection,
    sedimentation). Lower → faster, slightly less accurate. Default 100
    matches the original P3 paper sensitivity studies; n_elem=128 KiD runs
    show ~5× speed-up at qorder=40 with negligible bulk error."
    quadrature_order::Int = 100
end
Base.show(io::IO, mime::MIME"text/plain", x::P3IceParams) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)

P3IceParams(toml_dict::CP.ParamDict;
    is_limited = true, quadrature_order::Int = 100,
    inp_depletion_model = NIceProxyDepletion(τ_act = 300),
) = P3IceParams(;
        scheme = ParametersP3(toml_dict),
        terminal_velocity = Chen2022VelType(toml_dict),
        cloud_pdf = CloudParticlePDF_SB2006(toml_dict),
        rain_pdf = RainParticlePDF_SB2006(toml_dict; is_limited),
        ice_nucleation = Frostenberg2023(toml_dict),
        rain_freezing = RainFreezing(toml_dict),
        inp_depletion_model,
        quadrature_order,
    )

"""

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
    Microphysics2MParams(toml_dict::CP.ParamDict;
                         with_ice = false, is_limited = true,
                         activation_scheme = NoActivation())

Create a `Microphysics2MParams` object from a ClimaParams TOML dictionary.

# Arguments
- `toml_dict`: ClimaParams parameter dictionary
- `with_ice`: Include P3 ice-phase parameters (default: false)
- `is_limited`: Use limited rain size distribution parameters (default: true)
- `activation_scheme`: aerosol → cloud-droplet activation scheme to bake
  into `warm_rain.activation_scheme`. Default `NoActivation()`. Pass e.g.
  `DiagnosticNc{FT}(; N_c, q_thresh, τ_relax)` to enable activation in
  the fused 2M tendency without mutating the host's BMT call sites.
"""
Microphysics2MParams(toml_dict::CP.ParamDict;
    with_ice = false, is_limited = true,
    activation_scheme::AbstractActivationScheme = NoActivation(),
    quadrature_order::Int = 100,
    inp_depletion_model = NIceProxyDepletion(τ_act = 300),
) = Microphysics2MParams(;
        # Warm rain parameters (always present)
        warm_rain = WarmRainParams2M(toml_dict; is_limited, activation_scheme),
        # Optional ice phase parameters
        ice = with_ice ?
              P3IceParams(toml_dict; is_limited, quadrature_order, inp_depletion_model) :
              nothing,
    )
