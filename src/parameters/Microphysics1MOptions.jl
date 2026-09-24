export Microphysics1MOptions,
    MicrophysicsOption,
    CloudLiquidFormation,
    CloudIceFormation,
    ConstantTimescale,
    PrescribedIceNumber,
    TemperatureDependent,
    TemperatureDependentIceNumber,
    CloudIceMelt,
    RainAutoconversion,
    Kessler1M,
    PrescribedNd,
    SnowAutoconversion,
    NoSupersaturation,
    WithSupersaturation,
    RainEvaporation,
    SnowDepositionSublimation,
    SublimationOnly,
    DepositionAndSublimation,
    SnowMelt,
    CloudLiquidRainAccretion,
    CloudLiquidSnowAccretion,
    CloudIceRainAccretion,
    CloudIceSnowAccretion,
    RainSnowAccretion,
    Homogeneous,
    Heterogeneous,
    HomogeneousAndHeterogeneous

"""
    MicrophysicsOption

Abstract type for all microphysics process options.

Option types are empty singletons that select which variant of a process runs.
The parameter values a variant needs live in the `process_params` field of
[`Microphysics1MParams`](@ref), built by [`process_params_for`](@ref).
"""
abstract type MicrophysicsOption end

# ═══════════════════════════════════════════════════════════════════
# Multi-variant processes: abstract type + concrete subtypes
# ═══════════════════════════════════════════════════════════════════

"""
    CloudIceFormation <: MicrophysicsOption

Abstract type for cloud ice formation (deposition/sublimation) methods.
See subtypes: [`ConstantTimescale`](@ref), [`PrescribedIceNumber`](@ref),
[`TemperatureDependent`](@ref), [`TemperatureDependentIceNumber`](@ref).
"""
abstract type CloudIceFormation <: MicrophysicsOption end

"""
    RainAutoconversion <: MicrophysicsOption

Abstract type for rain autoconversion methods.
See subtypes: [`Kessler1M`](@ref), [`PrescribedNd`](@ref).
"""
abstract type RainAutoconversion <: MicrophysicsOption end

"""
    SnowAutoconversion <: MicrophysicsOption

Abstract type for snow autoconversion methods.
See subtypes: [`NoSupersaturation`](@ref), [`WithSupersaturation`](@ref).
"""
abstract type SnowAutoconversion <: MicrophysicsOption end

"""
    SnowDepositionSublimation <: MicrophysicsOption

Abstract type for snow deposition/sublimation methods.
See subtypes: [`SublimationOnly`](@ref), [`DepositionAndSublimation`](@ref).
"""
abstract type SnowDepositionSublimation <: MicrophysicsOption end

# ═══════════════════════════════════════════════════════════════════
# Cloud liquid formation (single variant → concrete type = process name)
# ═══════════════════════════════════════════════════════════════════

"""
    CloudLiquidFormation <: MicrophysicsOption

Constant relaxation timescale for liquid condensation and evaporation.
Parameters (`τ_relax`) are stored in `process_params.cloud_liquid_formation` in
[`Microphysics1MParams`](@ref).
"""
struct CloudLiquidFormation <: MicrophysicsOption end

# ═══════════════════════════════════════════════════════════════════
# Cloud ice formation variants
# ═══════════════════════════════════════════════════════════════════

"""
    ConstantTimescale <: CloudIceFormation

Constant relaxation timescale for ice deposition and sublimation.
Parameters (`τ_relax`) are stored in `process_params.cloud_ice_formation` in
[`Microphysics1MParams`](@ref).
"""
struct ConstantTimescale <: CloudIceFormation end

"""
    PrescribedIceNumber <: CloudIceFormation

Ice deposition/sublimation timescale derived from the prescribed cloud-ice
number concentration `N_0` in `CloudIce` (the same `N_0` used for
sedimentation).  Both deposition and sublimation use the dynamically
computed timescale.  No additional process parameters are required.
"""
struct PrescribedIceNumber <: CloudIceFormation end

"""
    TemperatureDependent <: CloudIceFormation

INP-dependent Frostenberg (2023) timescale for deposition,
with constant timescale for sublimation. Parameters (`τ_relax`, `frostenberg`)
are stored in `process_params.cloud_ice_formation` in
[`Microphysics1MParams`](@ref).
"""
struct TemperatureDependent <: CloudIceFormation end

"""
    TemperatureDependentIceNumber <: CloudIceFormation

Ice deposition/sublimation timescale derived, as for [`PrescribedIceNumber`](@ref),
from a cloud-ice number concentration and the resulting mean crystal radius, but
with the number concentration a prescribed function of temperature,
`N_ice(T) = min(N_ref exp(a + b max(T_freeze - T, 0)), N_max)`, instead of the constant
sedimentation number `N_0`. Both deposition and sublimation use the same
timescale. Parameters ([`IceNumberTemperatureFit`](@ref)) are stored in
`process_params.cloud_ice_formation` in [`Microphysics1MParams`](@ref).
"""
struct TemperatureDependentIceNumber <: CloudIceFormation end

# ═══════════════════════════════════════════════════════════════════
# Rain autoconversion variants
# ═══════════════════════════════════════════════════════════════════

"""
    Kessler1M <: RainAutoconversion

1-moment Kessler autoconversion of cloud liquid to rain: the smooth logistic transition of
`q_lcl` across a threshold, divided by a timescale. The threshold and the timescale may each
depend on the air vertical velocity `w`.
Parameters (a [`KesslerAcnv`](@ref) with `τ_slow`, `τ_fast`, `q_threshold_slow`,
`q_threshold_fast`, `w_0`, `k`) are stored in `process_params.rain_autoconversion` in
[`Microphysics1MParams`](@ref).
"""
struct Kessler1M <: RainAutoconversion end

"""
    PrescribedNd <: RainAutoconversion

Variable-timescale autoconversion using prescribed cloud droplet number Nc.
Parameters (a `VarTimescaleAcnv` with `τ`, `α`, `Nc`) are stored in
`process_params.rain_autoconversion` in [`Microphysics1MParams`](@ref).
"""
struct PrescribedNd <: RainAutoconversion end

# ═══════════════════════════════════════════════════════════════════
# Snow autoconversion variants
# ═══════════════════════════════════════════════════════════════════

"""
    NoSupersaturation <: SnowAutoconversion

Simplified autoconversion of cloud ice to snow without supersaturation dependence.
Parameters (an `Acnv1M` with `τ`, `q_threshold`, `k`) are stored in
`process_params.snow_autoconversion` in [`Microphysics1MParams`](@ref).
"""
struct NoSupersaturation <: SnowAutoconversion end

"""
    WithSupersaturation <: SnowAutoconversion

Harrington/Kaul autoconversion of cloud ice to snow with supersaturation dependence.
Parameters (`r_ice_snow`) are stored in `process_params.snow_autoconversion` in
[`Microphysics1MParams`](@ref).
"""
struct WithSupersaturation <: SnowAutoconversion end

# ═══════════════════════════════════════════════════════════════════
# Accretion (single variant each → concrete type = process name)
# ═══════════════════════════════════════════════════════════════════

"""
    CloudLiquidRainAccretion <: MicrophysicsOption

Cloud liquid + rain → rain (Marshall-Palmer kernel).
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_liquid_rain_accretion` in [`Microphysics1MParams`](@ref).
"""
struct CloudLiquidRainAccretion <: MicrophysicsOption end

"""
    CloudLiquidSnowAccretion <: MicrophysicsOption

Cloud liquid + snow → snow/rain depending on temperature
(includes warm-rain melt contribution).
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_liquid_snow_accretion` in [`Microphysics1MParams`](@ref).
"""
struct CloudLiquidSnowAccretion <: MicrophysicsOption end

"""
    CloudIceRainAccretion <: MicrophysicsOption

Cloud ice + rain → snow (Marshall-Palmer kernel).
The coupled rain-sink arm (rain + cloud ice → snow) is toggled automatically.
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_ice_rain_accretion` in [`Microphysics1MParams`](@ref).
"""
struct CloudIceRainAccretion <: MicrophysicsOption end

"""
    CloudIceSnowAccretion <: MicrophysicsOption

Cloud ice + snow → snow (Marshall-Palmer kernel).
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_ice_snow_accretion` in [`Microphysics1MParams`](@ref).
"""
struct CloudIceSnowAccretion <: MicrophysicsOption end

"""
    RainSnowAccretion <: MicrophysicsOption

Snow-rain collisions: both temperature pathways
(cold: rain→snow, warm: snow→rain) plus thermal melt.
Parameters (collision efficiency `e`, velocity dispersion `coeff_disp`) are
stored in `process_params.rain_snow_accretion` in [`Microphysics1MParams`](@ref).
"""
struct RainSnowAccretion <: MicrophysicsOption end

# ═══════════════════════════════════════════════════════════════════
# Snow deposition/sublimation variants
# ═══════════════════════════════════════════════════════════════════

"""Only sublimation (S < 0 over ice); deposition handled separately by non-equilibrium."""
struct SublimationOnly <: SnowDepositionSublimation end

"""Both sublimation (S < 0) and deposition (S > 0) in the Marshall-Palmer integral."""
struct DepositionAndSublimation <: SnowDepositionSublimation end

# ═══════════════════════════════════════════════════════════════════
# Single-variant on/off processes (concrete type = process name)
# ═══════════════════════════════════════════════════════════════════

"""Rain evaporation (sub-saturated conditions over liquid)."""
struct RainEvaporation <: MicrophysicsOption end

"""Cloud ice melts to cloud liquid above freezing."""
struct CloudIceMelt <: MicrophysicsOption end

"""Snow melts to rain above freezing."""
struct SnowMelt <: MicrophysicsOption end

"""
    Homogeneous <: MicrophysicsOption

All cloud liquid freezes to ice below the homogeneous nucleation temperature
(T < T_hom ≈ 233 K) on a short relaxation timescale.
"""
struct Homogeneous <: MicrophysicsOption end

"""
    Heterogeneous <: MicrophysicsOption

Cloud liquid freezing to cloud ice via Bigg (1953) immersion freezing
for T < T_freeze, using the Reisner et al. (1998) parameterization.
Droplet volume is based on prescribed cloud droplet number concentration.
"""
struct Heterogeneous <: MicrophysicsOption end

"""
    HomogeneousAndHeterogeneous <: MicrophysicsOption

Both homogeneous and heterogeneous cloud liquid freezing are active.
"""
struct HomogeneousAndHeterogeneous <: MicrophysicsOption end

# ═══════════════════════════════════════════════════════════════════
# Options struct
# ═══════════════════════════════════════════════════════════════════

"""
    Microphysics1MOptions{CLF, CIF, CIM, CLFr, RA, SA, RCE, SDS, SM, CLRA, CLSA, CIRA, CISA, RSA}

Process configuration for 1-moment microphysics.

Each field selects a process variant (a concrete `MicrophysicsOption` subtype).
Set any field to `nothing` to disable that process entirely.

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Default options
opts = CMP.Microphysics1MOptions()

# Disable cloud ice melt and snow melt:
opts = CMP.Microphysics1MOptions(;
    cloud_ice_melt = nothing,
    snow_melt = nothing,
)

# Switch rain autoconversion to the prescribed-Nd variant:
opts = CMP.Microphysics1MOptions(; rain_autoconversion = CMP.PrescribedNd())
```
"""
@kwdef struct Microphysics1MOptions{
    CLF, CIF, CIM, CLFr, RA, SA, RCE, SDS, SM, CLRA, CLSA, CIRA, CISA, RSA,
}
    "cloud liquid formation option"
    cloud_liquid_formation::CLF = CloudLiquidFormation()
    "cloud ice formation option"
    cloud_ice_formation::CIF = ConstantTimescale()
    "cloud ice melting option"
    cloud_ice_melt::CIM = CloudIceMelt()
    "cloud liquid freezing option"
    cloud_liquid_freezing::CLFr = HomogeneousAndHeterogeneous()
    "rain autoconversion option"
    rain_autoconversion::RA = Kessler1M()
    "cloud ice to snow autoconversion option"
    snow_autoconversion::SA = NoSupersaturation()
    "rain condensation/evaporation option"
    rain_condensation_evaporation::RCE = RainEvaporation()
    "snow sublimation/deposition option"
    snow_deposition_sublimation::SDS = DepositionAndSublimation()
    "snow melting option"
    snow_melt::SM = SnowMelt()
    "cloud liquid + rain accretion option"
    cloud_liquid_rain_accretion::CLRA = CloudLiquidRainAccretion()
    "cloud liquid + snow accretion option"
    cloud_liquid_snow_accretion::CLSA = CloudLiquidSnowAccretion()
    "cloud ice + rain accretion option (also toggles rain sink arm)"
    cloud_ice_rain_accretion::CIRA = CloudIceRainAccretion()
    "cloud ice + snow accretion option"
    cloud_ice_snow_accretion::CISA = CloudIceSnowAccretion()
    "rain-snow collisions option"
    rain_snow_accretion::RSA = RainSnowAccretion()
end

# ═══════════════════════════════════════════════════════════════════
# Process parameters: option → parameter data
# ═══════════════════════════════════════════════════════════════════

"""
    process_params_for(option, toml_dict)

Return the parameter data that the selected process `option` needs, read from
`toml_dict`. Returns `nothing` for disabled processes (`option === nothing`)
and for options that carry no parameters. The result is stored in the matching
field of `Microphysics1MParams.process_params` and read back at call time by
the process's tendency function.
"""
process_params_for(::Nothing, ::CP.ParamDict) = nothing
process_params_for(::MicrophysicsOption, ::CP.ParamDict) = nothing

function process_params_for(::TemperatureDependent, td::CP.ParamDict)
    p = make_params(td, name_map(TemperatureDependent))
    return (; τ_relax = p.τ_relax, frostenberg = Frostenberg2023(td))
end

process_params_for(::TemperatureDependentIceNumber, td::CP.ParamDict) =
    IceNumberTemperatureFit(td)

function process_params_for(::HomogeneousAndHeterogeneous, td::CP.ParamDict)
    return (;
        process_params_for(Homogeneous(), td)...,
        process_params_for(Heterogeneous(), td)...,
    )
end

process_params_for(::Kessler1M, td::CP.ParamDict) = KesslerAcnv(td)

process_params_for(::PrescribedNd, td::CP.ParamDict) = VarTimescaleAcnv(td)

function process_params_for(::NoSupersaturation, td::CP.ParamDict)
    p = make_params(td, name_map(NoSupersaturation))
    return Acnv1M(p.τ, p.q_threshold, p.k)
end

"""
    microphysics_1m_process_params(toml_dict, options)

Assemble the `process_params` container for [`Microphysics1MParams`](@ref) by
mapping each field of `options` through [`process_params_for`](@ref). The result
mirrors the `options` fields one-to-one.
"""
microphysics_1m_process_params(td::CP.ParamDict, o::Microphysics1MOptions) = (;
    cloud_liquid_formation = process_params_for(o.cloud_liquid_formation, td),
    cloud_ice_formation = process_params_for(o.cloud_ice_formation, td),
    cloud_ice_melt = process_params_for(o.cloud_ice_melt, td),
    cloud_liquid_freezing = process_params_for(o.cloud_liquid_freezing, td),
    rain_autoconversion = process_params_for(o.rain_autoconversion, td),
    snow_autoconversion = process_params_for(o.snow_autoconversion, td),
    rain_condensation_evaporation = process_params_for(o.rain_condensation_evaporation, td),
    snow_deposition_sublimation = process_params_for(o.snow_deposition_sublimation, td),
    snow_melt = process_params_for(o.snow_melt, td),
    cloud_liquid_rain_accretion = process_params_for(o.cloud_liquid_rain_accretion, td),
    cloud_liquid_snow_accretion = process_params_for(o.cloud_liquid_snow_accretion, td),
    cloud_ice_rain_accretion = process_params_for(o.cloud_ice_rain_accretion, td),
    cloud_ice_snow_accretion = process_params_for(o.cloud_ice_snow_accretion, td),
    rain_snow_accretion = process_params_for(o.rain_snow_accretion, td),
)

# ═══════════════════════════════════════════════════════════════════
# Deprecated TOML-dict constructors
#
# Option types carry no parameters, so they take no TOML dict.
# These shims keep downstream callers that still pass `toml_dict`
# working (the argument is ignored). Drop the argument to migrate.
# ═══════════════════════════════════════════════════════════════════
@deprecate CloudLiquidFormation(::CP.ParamDict) CloudLiquidFormation() false
@deprecate ConstantTimescale(::CP.ParamDict) ConstantTimescale() false
@deprecate PrescribedIceNumber(::CP.ParamDict) PrescribedIceNumber() false
@deprecate TemperatureDependent(::CP.ParamDict) TemperatureDependent() false
@deprecate Kessler1M(::CP.ParamDict) Kessler1M() false
@deprecate PrescribedNd(::CP.ParamDict) PrescribedNd() false
@deprecate NoSupersaturation(::CP.ParamDict) NoSupersaturation() false
@deprecate WithSupersaturation(::CP.ParamDict) WithSupersaturation() false
@deprecate CloudLiquidRainAccretion(::CP.ParamDict) CloudLiquidRainAccretion() false
@deprecate CloudLiquidSnowAccretion(::CP.ParamDict) CloudLiquidSnowAccretion() false
@deprecate CloudIceRainAccretion(::CP.ParamDict) CloudIceRainAccretion() false
@deprecate CloudIceSnowAccretion(::CP.ParamDict) CloudIceSnowAccretion() false
@deprecate RainSnowAccretion(::CP.ParamDict) RainSnowAccretion() false
