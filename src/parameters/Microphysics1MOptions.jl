export Microphysics1MOptions,
    MicrophysicsOption,
    NoCloudLiquidFormation,
    ConstantTimescaleCloudLiquidFormation,
    NoCloudIceFormation,
    ConstantTimescaleCloudIceFormation,
    TemperatureDependentCloudIceFormation,
    NoCloudIceMelt,
    CloudIceMelt,
    NoRainAutoconversion,
    RainAutoconversion1M,
    RainAutoconversionPrescribedNd,
    NoSnowAutoconversion,
    SnowAutoconversionNoSupersaturation,
    SnowAutoconversionWithSupersaturation,
    NoRainCondensationEvaporation,
    RainEvaporation,
    NoSnowDepositionSublimation,
    SnowSublimation,
    SnowDepositionSublimation,
    NoSnowMelt,
    SnowMelt,
    NoCloudLiquidRainAccretion,
    CloudLiquidRainAccretion,
    NoCloudLiquidSnowAccretion,
    CloudLiquidSnowAccretion,
    NoCloudIceRainAccretion,
    CloudIceRainAccretion,
    NoCloudIceSnowAccretion,
    CloudIceSnowAccretion,
    NoRainSnowAccretion,
    RainSnowAccretion
"""
    MicrophysicsOption

Abstract type for all microphysics process options.
"""
abstract type MicrophysicsOption end

"""No cloud liquid condensation/evaporation."""
struct NoCloudLiquidFormation <: MicrophysicsOption end

"""Use a constant relaxation timescale for liquid condensation and evaporation."""
struct ConstantTimescaleCloudLiquidFormation <: MicrophysicsOption end

"""No cloud ice deposition/sublimation."""
struct NoCloudIceFormation <: MicrophysicsOption end

"""Use a constant relaxation timescale for ice deposition and sublimation."""
struct ConstantTimescaleCloudIceFormation <: MicrophysicsOption end

"""Use the INP-dependent Frostenberg (2023) timescale for deposition,
   with constant timescale for sublimation."""
struct TemperatureDependentCloudIceFormation <: MicrophysicsOption end

"""No cloud ice melting process."""
struct NoCloudIceMelt <: MicrophysicsOption end

"""Cloud ice melts to cloud liquid above freezing."""
struct CloudIceMelt <: MicrophysicsOption end

"""No rain autoconversion."""
struct NoRainAutoconversion <: MicrophysicsOption end

"""1-moment Kessler autoconversion of cloud liquid to rain."""
struct RainAutoconversion1M <: MicrophysicsOption end

"""Variable-timescale 2-moment autoconversion (uses prescribed Nc)."""
struct RainAutoconversionPrescribedNd <: MicrophysicsOption end

"""No snow autoconversion."""
struct NoSnowAutoconversion <: MicrophysicsOption end

"""Simplified autoconversion without supersaturation dependence."""
struct SnowAutoconversionNoSupersaturation <: MicrophysicsOption end

"""Harrington/Kaul autoconversion with supersaturation dependence."""
struct SnowAutoconversionWithSupersaturation <: MicrophysicsOption end

"""No rain condensation/evaporation."""
struct NoRainCondensationEvaporation <: MicrophysicsOption end

"""Rain evaporation (sub-saturated conditions over liquid)."""
struct RainEvaporation <: MicrophysicsOption end

"""No snow deposition or sublimation."""
struct NoSnowDepositionSublimation <: MicrophysicsOption end

"""Only sublimation (S < 0 over ice); deposition handled separately by non-equilibrium."""
struct SnowSublimation <: MicrophysicsOption end

"""Both sublimation (S < 0) and deposition (S > 0) in the Marshall-Palmer integral."""
struct SnowDepositionSublimation <: MicrophysicsOption end

"""No snow melting."""
struct NoSnowMelt <: MicrophysicsOption end

"""Snow melts to rain above freezing."""
struct SnowMelt <: MicrophysicsOption end

"""No cloud liquid + rain accretion."""
struct NoCloudLiquidRainAccretion <: MicrophysicsOption end

"""Cloud liquid + rain → rain (Marshall-Palmer kernel).
Also triggers the coupled rain-sink arm of cloud ice + rain collisions."""
struct CloudLiquidRainAccretion <: MicrophysicsOption end

"""No cloud liquid + snow accretion."""
struct NoCloudLiquidSnowAccretion <: MicrophysicsOption end

"""Cloud liquid + snow → snow/rain depending on temperature (includes warm-rain melt contribution)."""
struct CloudLiquidSnowAccretion <: MicrophysicsOption end

"""No cloud ice + rain accretion."""
struct NoCloudIceRainAccretion <: MicrophysicsOption end

"""Cloud ice + rain → snow (Marshall-Palmer kernel).
The coupled rain-sink arm (rain + cloud ice → snow) is toggled automatically."""
struct CloudIceRainAccretion <: MicrophysicsOption end

"""No cloud ice + snow accretion."""
struct NoCloudIceSnowAccretion <: MicrophysicsOption end

"""Cloud ice + snow → snow (Marshall-Palmer kernel)."""
struct CloudIceSnowAccretion <: MicrophysicsOption end

"""No rain-snow collisions."""
struct NoRainSnowAccretion <: MicrophysicsOption end

"""Snow-rain collisions: both temperature pathways (cold: rai→sno, warm: sno→rai) plus thermal melt."""
struct RainSnowAccretion <: MicrophysicsOption end

"""
    Microphysics1MOptions{CLF, CIF, CIM, RA, SA, RCE, SDS, SM, CLRA, CLSA, CIRA, CISA, RSA}

Process configuration for 1-moment microphysics.
Each field selects a variant for one microphysical process.
The `No*` variants disable a process entirely (return zero tendency).

Defaults match the baseline (old main branch) behavior.

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Default (baseline) options:
opts = CMP.Microphysics1MOptions()

# Disable all accretion, enable INP-dependent ice and cloud ice melt:
opts = CMP.Microphysics1MOptions(
    cloud_ice_formation           = CMP.TemperatureDependentCloudIceFormation(),
    cloud_ice_melt                = CMP.CloudIceMelt(),
    cloud_liquid_rain_accretion   = CMP.NoCloudLiquidRainAccretion(),
    cloud_liquid_snow_accretion   = CMP.NoCloudLiquidSnowAccretion(),
    cloud_ice_rain_accretion      = CMP.NoCloudIceRainAccretion(),
    cloud_ice_snow_accretion      = CMP.NoCloudIceSnowAccretion(),
    rain_snow_accretion           = CMP.NoRainSnowAccretion(),
    snow_deposition_sublimation   = CMP.SnowDepositionSublimation(),
)
```
"""
@kwdef struct Microphysics1MOptions{
    CLF <: MicrophysicsOption,
    CIF <: MicrophysicsOption,
    CIM <: MicrophysicsOption,
    RA <: MicrophysicsOption,
    SA <: MicrophysicsOption,
    RCE <: MicrophysicsOption,
    SDS <: MicrophysicsOption,
    SM <: MicrophysicsOption,
    CLRA <: MicrophysicsOption,
    CLSA <: MicrophysicsOption,
    CIRA <: MicrophysicsOption,
    CISA <: MicrophysicsOption,
    RSA <: MicrophysicsOption,
} <: ParametersType
    "cloud liquid formation timescale option"
    cloud_liquid_formation::CLF = ConstantTimescaleCloudLiquidFormation()
    "cloud ice formation timescale option"
    cloud_ice_formation::CIF = ConstantTimescaleCloudIceFormation()
    "cloud ice melting option"
    cloud_ice_melt::CIM = CloudIceMelt()
    "rain autoconversion option"
    rain_autoconversion::RA = RainAutoconversion1M()
    "cloud ice to snow autoconversion option"
    snow_autoconversion::SA = SnowAutoconversionNoSupersaturation()
    "rain condensation/evaporation option"
    rain_condensation_evaporation::RCE = RainEvaporation()
    "snow sublimation/deposition option"
    snow_deposition_sublimation::SDS = SnowDepositionSublimation()
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
