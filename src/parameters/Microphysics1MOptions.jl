export Microphysics1MOptions,
    MicrophysicsOption,
    CloudLiquidFormation,
    CloudIceFormation,
    ConstantTimescale,
    TemperatureDependent,
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
    RainSnowAccretion

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
See subtypes: [`ConstantTimescale`](@ref), [`TemperatureDependent`](@ref).
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
Parameters (`τ_relax`) are stored in `process_params.cloud_liquid_formation`.
"""
struct CloudLiquidFormation <: MicrophysicsOption end

# ═══════════════════════════════════════════════════════════════════
# Cloud ice formation variants
# ═══════════════════════════════════════════════════════════════════

"""
    ConstantTimescale <: CloudIceFormation

Constant relaxation timescale for ice deposition and sublimation.
Parameters (`τ_relax`) are stored in `process_params.cloud_ice_formation`.
"""
struct ConstantTimescale <: CloudIceFormation end

"""
    TemperatureDependent <: CloudIceFormation

INP-dependent Frostenberg (2023) timescale for deposition,
with constant timescale for sublimation. Parameters (`τ_relax`, `frostenberg`)
are stored in `process_params.cloud_ice_formation`.
"""
struct TemperatureDependent <: CloudIceFormation end

# ═══════════════════════════════════════════════════════════════════
# Rain autoconversion variants
# ═══════════════════════════════════════════════════════════════════

"""
    Kessler1M <: RainAutoconversion

1-moment Kessler autoconversion of cloud liquid to rain.
Parameters (an `Acnv1M` with `τ`, `q_threshold`, `k`) are stored in
`process_params.rain_autoconversion`.
"""
struct Kessler1M <: RainAutoconversion end

"""
    PrescribedNd <: RainAutoconversion

Variable-timescale autoconversion using prescribed cloud droplet number Nc.
Parameters (a `VarTimescaleAcnv` with `τ`, `α`, `Nc`) are stored in
`process_params.rain_autoconversion`.
"""
struct PrescribedNd <: RainAutoconversion end

# ═══════════════════════════════════════════════════════════════════
# Snow autoconversion variants
# ═══════════════════════════════════════════════════════════════════

"""
    NoSupersaturation <: SnowAutoconversion

Simplified autoconversion of cloud ice to snow without supersaturation dependence.
Parameters (an `Acnv1M` with `τ`, `q_threshold`, `k`) are stored in
`process_params.snow_autoconversion`.
"""
struct NoSupersaturation <: SnowAutoconversion end

"""
    WithSupersaturation <: SnowAutoconversion

Harrington/Kaul autoconversion of cloud ice to snow with supersaturation dependence.
Parameters (`r_ice_snow`) are stored in `process_params.snow_autoconversion`.
"""
struct WithSupersaturation <: SnowAutoconversion end

# ═══════════════════════════════════════════════════════════════════
# Accretion (single variant each → concrete type = process name)
# ═══════════════════════════════════════════════════════════════════

"""
    CloudLiquidRainAccretion <: MicrophysicsOption

Cloud liquid + rain → rain (Marshall-Palmer kernel).
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_liquid_rain_accretion`.
"""
struct CloudLiquidRainAccretion <: MicrophysicsOption end

"""
    CloudLiquidSnowAccretion <: MicrophysicsOption

Cloud liquid + snow → snow/rain depending on temperature
(includes warm-rain melt contribution).
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_liquid_snow_accretion`.
"""
struct CloudLiquidSnowAccretion <: MicrophysicsOption end

"""
    CloudIceRainAccretion <: MicrophysicsOption

Cloud ice + rain → snow (Marshall-Palmer kernel).
The coupled rain-sink arm (rain + cloud ice → snow) is toggled automatically.
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_ice_rain_accretion`.
"""
struct CloudIceRainAccretion <: MicrophysicsOption end

"""
    CloudIceSnowAccretion <: MicrophysicsOption

Cloud ice + snow → snow (Marshall-Palmer kernel).
Parameters (collision efficiency `e`) are stored in
`process_params.cloud_ice_snow_accretion`.
"""
struct CloudIceSnowAccretion <: MicrophysicsOption end

"""
    RainSnowAccretion <: MicrophysicsOption

Snow-rain collisions: both temperature pathways
(cold: rain→snow, warm: snow→rain) plus thermal melt.
Parameters (collision efficiency `e`, velocity dispersion `coeff_disp`) are
stored in `process_params.rain_snow_accretion`.
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

# ═══════════════════════════════════════════════════════════════════
# Options struct
# ═══════════════════════════════════════════════════════════════════

"""
    Microphysics1MOptions{CLF, CIF, CIM, RA, SA, RCE, SDS, SM, CLRA, CLSA, CIRA, CISA, RSA}

Process configuration for 1-moment microphysics.
Each field selects a variant for one microphysical process.
Set a field to `nothing` to disable the process entirely.

# Fields
$(DocStringExtensions.FIELDS)

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
    CLF, CIF, CIM, RA, SA, RCE, SDS, SM, CLRA, CLSA, CIRA, CISA, RSA,
}
    "cloud liquid formation option"
    cloud_liquid_formation::CLF = CloudLiquidFormation()
    "cloud ice formation option"
    cloud_ice_formation::CIF = ConstantTimescale()
    "cloud ice melting option"
    cloud_ice_melt::CIM = CloudIceMelt()
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

function process_params_for(::CloudLiquidFormation, td::CP.ParamDict)
    name_map = (; :condensation_evaporation_timescale => :τ_relax)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::ConstantTimescale, td::CP.ParamDict)
    name_map = (; :sublimation_deposition_timescale => :τ_relax)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::TemperatureDependent, td::CP.ParamDict)
    name_map = (; :sublimation_deposition_timescale => :τ_relax)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return (; τ_relax = p.τ_relax, frostenberg = Frostenberg2023(td))
end

function process_params_for(::Kessler1M, td::CP.ParamDict)
    name_map = (;
        :rain_autoconversion_timescale => :τ,
        :cloud_liquid_water_specific_humidity_autoconversion_threshold => :q_threshold,
        :threshold_smooth_transition_steepness => :k,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Acnv1M(p.τ, p.q_threshold, p.k)
end

process_params_for(::PrescribedNd, td::CP.ParamDict) = VarTimescaleAcnv(td)

function process_params_for(::NoSupersaturation, td::CP.ParamDict)
    name_map = (;
        :snow_autoconversion_timescale => :τ,
        :cloud_ice_specific_humidity_autoconversion_threshold => :q_threshold,
        :threshold_smooth_transition_steepness => :k,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Acnv1M(p.τ, p.q_threshold, p.k)
end

function process_params_for(::WithSupersaturation, td::CP.ParamDict)
    name_map = (; :ice_snow_threshold_radius => :r_ice_snow)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::CloudLiquidRainAccretion, td::CP.ParamDict)
    name_map = (; :cloud_liquid_rain_collision_efficiency => :e)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::CloudLiquidSnowAccretion, td::CP.ParamDict)
    name_map = (; :cloud_liquid_snow_collision_efficiency => :e)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::CloudIceRainAccretion, td::CP.ParamDict)
    name_map = (; :cloud_ice_rain_collision_efficiency => :e)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::CloudIceSnowAccretion, td::CP.ParamDict)
    name_map = (; :cloud_ice_snow_collision_efficiency => :e)
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
end

function process_params_for(::RainSnowAccretion, td::CP.ParamDict)
    name_map = (;
        :rain_snow_collision_efficiency => :e,
        :rain_snow_velocity_dispersion_coefficient => :coeff_disp,
    )
    return CP.get_parameter_values(td, name_map, "CloudMicrophysics")
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
