export Microphysics1MOptions,
    MicrophysicsOption,
    ConstantTimescaleCloudLiquidFormation,
    ConstantTimescaleCloudIceFormation,
    TemperatureDependentCloudIceFormation,
    NoCloudIceMelt,
    CloudIceMeltToLiquid,
    LiquidAutoconv1M,
    LiquidAutoconv2M,
    SnowAutoconvNoSupersat,
    SnowAutoconvWithSupersat,
    EvaporationOnly,
    SublimationOnly,
    DepositionSublimation,
    SnowMelt
"""
    MicrophysicsOption

Abstract type for all microphysics process options.
"""
abstract type MicrophysicsOption end

"""Use a constant relaxation timescale for liquid condensation and evaporation."""
struct ConstantTimescaleCloudLiquidFormation <: MicrophysicsOption end

"""Use a constant relaxation timescale for ice deposition and sublimation."""
struct ConstantTimescaleCloudIceFormation <: MicrophysicsOption end

"""Use the INP-dependent Frostenberg (2023) timescale for deposition,
   with constant timescale for sublimation."""
struct TemperatureDependentCloudIceFormation <: MicrophysicsOption end

"""No cloud ice melting process."""
struct NoCloudIceMelt <: MicrophysicsOption end

"""Cloud ice melts to cloud liquid above freezing."""
struct CloudIceMeltToLiquid <: MicrophysicsOption end

"""1-moment Kessler autoconversion only."""
struct LiquidAutoconv1M <: MicrophysicsOption end

"""Variable-timescale 2-moment autoconversion (uses prescribed Nc)."""
struct LiquidAutoconv2M <: MicrophysicsOption end

"""Simplified autoconversion without supersaturation dependence."""
struct SnowAutoconvNoSupersat <: MicrophysicsOption end

"""Harrington/Kaul autoconversion with supersaturation dependence."""
struct SnowAutoconvWithSupersat <: MicrophysicsOption end

"""Evaporation only (sub-saturated conditions over liquid)."""
struct EvaporationOnly <: MicrophysicsOption end

"""Only sublimation (S < 0); deposition handled separately by non-equilibrium."""
struct SublimationOnly <: MicrophysicsOption end

"""Both sublimation (S < 0) and deposition (S > 0) in the Marshall-Palmer integral."""
struct DepositionSublimation <: MicrophysicsOption end

"""Snow melts to rain above freezing."""
struct SnowMelt <: MicrophysicsOption end

"""
    Microphysics1MOptions{CLF, CIF, CIM, CLA, SA, RE, SS, SM}

Process configuration for 1-moment microphysics.
Each field selects a variant for one microphysical process.

Defaults match the baseline (old main branch) behavior.

# Fields
$(DocStringExtensions.FIELDS)

# Example
```julia
using CloudMicrophysics.Parameters as CMP

# Default (baseline) options:
opts = CMP.Microphysics1MOptions()

# New branch options (INP-dependent ice, cloud ice melt, deposition+sublimation):
opts = CMP.Microphysics1MOptions(
    cloud_liquid_formation  = CMP.ConstantTimescaleCloudLiquidFormation(),
    cloud_ice_formation     = CMP.TemperatureDependentCloudIceFormation(),
    cloud_ice_melt          = CMP.CloudIceMeltToLiquid(),
    snow_sublimation    = CMP.DepositionSublimation(),
)
```
"""
@kwdef struct Microphysics1MOptions{
    CLF <: MicrophysicsOption,
    CIF <: MicrophysicsOption,
    CIM <: MicrophysicsOption,
    CLA <: MicrophysicsOption,
    SA <: MicrophysicsOption,
    RE <: MicrophysicsOption,
    SS <: MicrophysicsOption,
    SM <: MicrophysicsOption,
} <: ParametersType
    "cloud liquid formation timescale option"
    cloud_liquid_formation::CLF = ConstantTimescaleCloudLiquidFormation()
    "cloud ice formation timescale option"
    cloud_ice_formation::CIF = ConstantTimescaleCloudIceFormation()
    "cloud ice melting option"
    cloud_ice_melt::CIM = CloudIceMeltToLiquid()
    "cloud liquid autoconversion option"
    cloud_liquid_autoconversion::CLA = LiquidAutoconv1M()
    "cloud ice to snow autoconversion option"
    snow_autoconversion::SA = SnowAutoconvNoSupersat()
    "rain evaporation option"
    rain_evaporation::RE = EvaporationOnly()
    "snow sublimation/deposition option"
    snow_sublimation::SS = DepositionSublimation()
    "snow melting option"
    snow_melt::SM = SnowMelt()
end
