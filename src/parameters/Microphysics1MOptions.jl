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
    CloudLiquidFormation{FT} <: MicrophysicsOption

Constant relaxation timescale for liquid condensation and evaporation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudLiquidFormation{FT} <: MicrophysicsOption
    "condensation/evaporation non-equilibrium relaxation timescale [s]"
    τ_relax::FT
end

function CloudLiquidFormation(td::CP.ParamDict)
    name_map = (; :condensation_evaporation_timescale => :τ_relax)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return CloudLiquidFormation(p.τ_relax)
end

# ═══════════════════════════════════════════════════════════════════
# Cloud ice formation variants
# ═══════════════════════════════════════════════════════════════════

"""
    ConstantTimescale{FT} <: CloudIceFormation

Constant relaxation timescale for ice deposition and sublimation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ConstantTimescale{FT} <: CloudIceFormation
    "deposition/sublimation non-equilibrium relaxation timescale [s]"
    τ_relax::FT
end

function ConstantTimescale(td::CP.ParamDict)
    name_map = (; :sublimation_deposition_timescale => :τ_relax)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return ConstantTimescale(p.τ_relax)
end

"""
    TemperatureDependent{FT, FR} <: CloudIceFormation

INP-dependent Frostenberg (2023) timescale for deposition,
with constant timescale for sublimation.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct TemperatureDependent{FT, FR} <: CloudIceFormation
    "sublimation arm: constant relaxation timescale [s]"
    τ_relax::FT
    "Frostenberg (2023) INP parameters for deposition"
    frostenberg::FR
end

function TemperatureDependent(td::CP.ParamDict)
    name_map = (; :sublimation_deposition_timescale => :τ_relax)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return TemperatureDependent(p.τ_relax, Frostenberg2023(td))
end

# ═══════════════════════════════════════════════════════════════════
# Rain autoconversion variants
# ═══════════════════════════════════════════════════════════════════

"""
    Kessler1M{AC} <: RainAutoconversion

1-moment Kessler autoconversion of cloud liquid to rain.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct Kessler1M{AC} <: RainAutoconversion
    "autoconversion parameters (τ, q_threshold, k)"
    acnv1M::AC
end

function Kessler1M(td::CP.ParamDict)
    name_map = (;
        :rain_autoconversion_timescale => :τ,
        :cloud_liquid_water_specific_humidity_autoconversion_threshold => :q_threshold,
        :threshold_smooth_transition_steepness => :k,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Kessler1M(Acnv1M(p.τ, p.q_threshold, p.k))
end

"""
    PrescribedNd{VA} <: RainAutoconversion

Variable-timescale autoconversion using prescribed cloud droplet number Nc.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct PrescribedNd{VA} <: RainAutoconversion
    "variable-timescale autoconversion parameters (τ, α, Nc)"
    autoconv::VA
end

function PrescribedNd(td::CP.ParamDict)
    return PrescribedNd(VarTimescaleAcnv(td))
end

# ═══════════════════════════════════════════════════════════════════
# Snow autoconversion variants
# ═══════════════════════════════════════════════════════════════════

"""
    NoSupersaturation{AC} <: SnowAutoconversion

Simplified autoconversion of cloud ice to snow without supersaturation dependence.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct NoSupersaturation{AC} <: SnowAutoconversion
    "autoconversion parameters (τ, q_threshold, k)"
    acnv1M::AC
end

function NoSupersaturation(td::CP.ParamDict)
    name_map = (;
        :snow_autoconversion_timescale => :τ,
        :cloud_ice_specific_humidity_autoconversion_threshold => :q_threshold,
        :threshold_smooth_transition_steepness => :k,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return NoSupersaturation(Acnv1M(p.τ, p.q_threshold, p.k))
end

"""
    WithSupersaturation{FT} <: SnowAutoconversion

Harrington/Kaul autoconversion of cloud ice to snow with supersaturation dependence.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct WithSupersaturation{FT} <: SnowAutoconversion
    "ice-snow threshold radius [m]"
    r_ice_snow::FT
end

function WithSupersaturation(td::CP.ParamDict)
    name_map = (; :ice_snow_threshold_radius => :r_ice_snow)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return WithSupersaturation(p.r_ice_snow)
end

# ═══════════════════════════════════════════════════════════════════
# Accretion (single variant each → concrete type = process name)
# ═══════════════════════════════════════════════════════════════════

"""
    CloudLiquidRainAccretion{FT} <: MicrophysicsOption

Cloud liquid + rain → rain (Marshall-Palmer kernel).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudLiquidRainAccretion{FT} <: MicrophysicsOption
    "collision efficiency [-]"
    e::FT
end

function CloudLiquidRainAccretion(td::CP.ParamDict)
    name_map = (; :cloud_liquid_rain_collision_efficiency => :e)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return CloudLiquidRainAccretion(p.e)
end

"""
    CloudLiquidSnowAccretion{FT} <: MicrophysicsOption

Cloud liquid + snow → snow/rain depending on temperature
(includes warm-rain melt contribution).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudLiquidSnowAccretion{FT} <: MicrophysicsOption
    "collision efficiency [-]"
    e::FT
end

function CloudLiquidSnowAccretion(td::CP.ParamDict)
    name_map = (; :cloud_liquid_snow_collision_efficiency => :e)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return CloudLiquidSnowAccretion(p.e)
end

"""
    CloudIceRainAccretion{FT} <: MicrophysicsOption

Cloud ice + rain → snow (Marshall-Palmer kernel).
The coupled rain-sink arm (rain + cloud ice → snow) is toggled automatically.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudIceRainAccretion{FT} <: MicrophysicsOption
    "collision efficiency [-]"
    e::FT
end

function CloudIceRainAccretion(td::CP.ParamDict)
    name_map = (; :cloud_ice_rain_collision_efficiency => :e)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return CloudIceRainAccretion(p.e)
end

"""
    CloudIceSnowAccretion{FT} <: MicrophysicsOption

Cloud ice + snow → snow (Marshall-Palmer kernel).

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudIceSnowAccretion{FT} <: MicrophysicsOption
    "collision efficiency [-]"
    e::FT
end

function CloudIceSnowAccretion(td::CP.ParamDict)
    name_map = (; :cloud_ice_snow_collision_efficiency => :e)
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return CloudIceSnowAccretion(p.e)
end

"""
    RainSnowAccretion{FT} <: MicrophysicsOption

Snow-rain collisions: both temperature pathways
(cold: rain→snow, warm: snow→rain) plus thermal melt.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct RainSnowAccretion{FT} <: MicrophysicsOption
    "collision efficiency [-]"
    e::FT
    "velocity dispersion coefficient [-]"
    coeff_disp::FT
end

function RainSnowAccretion(td::CP.ParamDict)
    name_map = (;
        :rain_snow_collision_efficiency => :e,
        :rain_snow_velocity_dispersion_coefficient => :coeff_disp,
    )
    p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return RainSnowAccretion(p.e, p.coeff_disp)
end

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

# Default options require TOML dict for parametric types
(i.e. for populating the free parameter values)
opts = CMP.Microphysics1MOptions(toml_dict)

# Disable cloud ice melt and snow melt:
opts = CMP.Microphysics1MOptions(toml_dict;
    cloud_ice_melt = nothing,
    snow_melt = nothing,
)
```
"""
@kwdef struct Microphysics1MOptions{CLF, CIF, CIM, RA, SA, RCE, SDS, SM, CLRA, CLSA, CIRA, CISA, RSA}
    "cloud liquid formation option"
    cloud_liquid_formation::CLF
    "cloud ice formation option"
    cloud_ice_formation::CIF
    "cloud ice melting option"
    cloud_ice_melt::CIM
    "rain autoconversion option"
    rain_autoconversion::RA
    "cloud ice to snow autoconversion option"
    snow_autoconversion::SA
    "rain condensation/evaporation option"
    rain_condensation_evaporation::RCE
    "snow sublimation/deposition option"
    snow_deposition_sublimation::SDS
    "snow melting option"
    snow_melt::SM
    "cloud liquid + rain accretion option"
    cloud_liquid_rain_accretion::CLRA
    "cloud liquid + snow accretion option"
    cloud_liquid_snow_accretion::CLSA
    "cloud ice + rain accretion option (also toggles rain sink arm)"
    cloud_ice_rain_accretion::CIRA
    "cloud ice + snow accretion option"
    cloud_ice_snow_accretion::CISA
    "rain-snow collisions option"
    rain_snow_accretion::RSA
end

"""
    Microphysics1MOptions(toml_dict::CP.ParamDict; kwargs...)

Create a `Microphysics1MOptions` with default options and their parameters populated
from the TOML dictionary. Override individual options via keyword arguments.
"""
function Microphysics1MOptions(toml_dict::CP.ParamDict; kwargs...)
    defaults = (;
        cloud_liquid_formation = CloudLiquidFormation(toml_dict),
        cloud_ice_formation = ConstantTimescale(toml_dict),
        cloud_ice_melt = CloudIceMelt(),
        rain_autoconversion = Kessler1M(toml_dict),
        snow_autoconversion = NoSupersaturation(toml_dict),
        rain_condensation_evaporation = RainEvaporation(),
        snow_deposition_sublimation = DepositionAndSublimation(),
        snow_melt = SnowMelt(),
        cloud_liquid_rain_accretion = CloudLiquidRainAccretion(toml_dict),
        cloud_liquid_snow_accretion = CloudLiquidSnowAccretion(toml_dict),
        cloud_ice_rain_accretion = CloudIceRainAccretion(toml_dict),
        cloud_ice_snow_accretion = CloudIceSnowAccretion(toml_dict),
        rain_snow_accretion = RainSnowAccretion(toml_dict),
    )
    merged = merge(defaults, NamedTuple(kwargs))
    return Microphysics1MOptions(; merged...)
end
