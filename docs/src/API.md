```@meta
CurrentModule = CloudMicrophysics
```
# .

# Non-equilibrium cloud formation

```@docs
MicrophysicsNonEq
MicrophysicsNonEq.τ_relax
MicrophysicsNonEq.conv_q_vap_to_q_liq_ice
MicrophysicsNonEq.conv_q_vap_to_q_liq_ice_MM2015
MicrophysicsNonEq.terminal_velocity
```

# 0-moment precipitation microphysics

```@docs
Microphysics0M
Microphysics0M.remove_precipitation
```

# 1-moment precipitation microphysics

```@docs
Microphysics1M
Microphysics1M.get_v0
Microphysics1M.get_n0
Microphysics1M.lambda_inverse
Microphysics1M.terminal_velocity
Microphysics1M.conv_q_liq_to_q_rai
Microphysics1M.conv_q_ice_to_q_sno_no_supersat
Microphysics1M.conv_q_ice_to_q_sno
Microphysics1M.accretion
Microphysics1M.accretion_rain_sink
Microphysics1M.accretion_snow_rain
Microphysics1M.evaporation_sublimation
Microphysics1M.snow_melt
```

# 2-moment precipitation microphysics
```@docs
Microphysics2M
```
## Size distributions

### Parameters
```@docs
Microphysics2M.pdf_cloud_parameters
Microphysics2M.pdf_rain_parameters
Microphysics2M.pdf_cloud_parameters_mass
Microphysics2M.pdf_rain_parameters_mass
Microphysics2M.log_pdf_cloud_parameters_mass
```

### Size distributions
```@docs
Microphysics2M.size_distribution(::CMP.RainParticlePDF_SB2006, _, _, _)
Microphysics2M.size_distribution(::CMP.CloudParticlePDF_SB2006, _, _, _)
Microphysics2M.size_distribution_value
Microphysics2M.get_size_distribution_bounds
```

## Rates
```@docs
Microphysics2M.LiqRaiRates
Microphysics2M.autoconversion
Microphysics2M.accretion
Microphysics2M.liquid_self_collection
Microphysics2M.autoconversion_and_liquid_self_collection
Microphysics2M.rain_self_collection
Microphysics2M.rain_breakup
Microphysics2M.rain_self_collection_and_breakup
Microphysics2M.rain_terminal_velocity
Microphysics2M.rain_evaporation
Microphysics2M.conv_q_liq_to_q_rai
Microphysics2M.number_increase_for_mass_limit
Microphysics2M.number_decrease_for_mass_limit
```

## Distribution tools for 2-moment microphysics

```@docs
DistributionTools
DistributionTools.generalized_gamma_quantile
DistributionTools.generalized_gamma_cdf
DistributionTools.exponential_cdf
DistributionTools.exponential_quantile
```

# P3 scheme

```@docs
P3Scheme
```

## Construct parameterization set
```@docs
CMP.ParametersP3{Float64}
CMP.ParametersP3(::Float64)
CMP.ParametersP3(::Parameters.CP.AbstractTOMLDict)
```

### Sub-parameterizations
```@docs
CMP.MassPowerLaw
CMP.AreaPowerLaw
CMP.SlopeLaw
CMP.SlopePowerLaw{Float64}
CMP.SlopeConstant{Float64}
CMP.VentilationFactor
CMP.LocalRimeDensity
```

## Obtain particle state

```@docs
P3Scheme.P3State
P3Scheme.get_state
```

### State relationships

```@docs
P3Scheme.get_thresholds_ρ_g
P3Scheme.get_ρ_d
P3Scheme.get_ρ_g
P3Scheme.get_D_th
P3Scheme.get_D_gr
P3Scheme.get_D_cr
```

### Derived quantities

#### Main particle methods
These methods are a function of the particle diameter, `D`.
```@docs
P3Scheme.ice_mass
P3Scheme.ice_density
P3Scheme.∂ice_mass_∂D
P3Scheme.ice_area
P3Scheme.ϕᵢ
```

## Distribution parameters

```@docs
P3Scheme.get_μ
P3Scheme.get_logN₀
P3Scheme.get_distribution_logλ
```

### Distribution relationships

```@docs
P3Scheme.logN′ice
P3Scheme.size_distribution(::P3Scheme.P3State, _)
P3Scheme.loggamma_inc_moment
P3Scheme.loggamma_moment
P3Scheme.logmass_gamma_moment
P3Scheme.logLdivN
```

### Derived integral quantities
These methods integrate over the particle size distribution.
```@docs
P3Scheme.D_m
P3Scheme.ice_particle_terminal_velocity
P3Scheme.ice_terminal_velocity_number_weighted
P3Scheme.ice_terminal_velocity_mass_weighted
```

### Processes

#### Heterogeneous ice nucleation
```@docs
P3Scheme.het_ice_nucleation
```

#### Ice melting
```@docs
P3Scheme.ice_melt
```

#### Collisions with liquid droplets
```@docs
P3Scheme.bulk_liquid_ice_collision_sources
```
Supporting methods:
```@docs
P3Scheme.volumetric_collision_rate_integrand
P3Scheme.compute_max_freeze_rate
P3Scheme.compute_local_rime_density
P3Scheme.get_liquid_integrals
P3Scheme.∫liquid_ice_collisions
```

### Supporting integral methods

```@docs
P3Scheme.∫fdD
P3Scheme.∫fdD_error
P3Scheme.integral_bounds
```

# Aerosol model

```@docs
AerosolModel
AerosolModel.Mode_B
AerosolModel.Mode_κ
AerosolModel.AerosolDistribution
```

# Aerosol activation

```@docs
AerosolActivation
AerosolActivation.mean_hygroscopicity_parameter
AerosolActivation.max_supersaturation
AerosolActivation.N_activated_per_mode
AerosolActivation.M_activated_per_mode
AerosolActivation.total_N_activated
AerosolActivation.total_M_activated
```

# Artifact calling

```@docs
ArtifactCalling
ArtifactCalling.AIDA_ice_nucleation
```

# Heterogeneous ice nucleation
```@docs
HetIceNucleation
HetIceNucleation.dust_activated_number_fraction
HetIceNucleation.MohlerDepositionRate
HetIceNucleation.deposition_J
HetIceNucleation.ABIFM_J
HetIceNucleation.P3_deposition_N_i
HetIceNucleation.P3_het_N_i
HetIceNucleation.INP_concentration_frequency
HetIceNucleation.INP_concentration_mean
```

# Homogeneous ice nucleation
```@docs
HomIceNucleation
HomIceNucleation.homogeneous_J_cubic
HomIceNucleation.homogeneous_J_linear
```

# Cloud diagnostics
```@docs
CloudDiagnostics
CloudDiagnostics.radar_reflectivity_1M
CloudDiagnostics.radar_reflectivity_2M
CloudDiagnostics.effective_radius_const
CloudDiagnostics.effective_radius_Liu_Hallet_97
CloudDiagnostics.effective_radius_2M
```

# Common utility functions

```@docs
Common
Common.G_func_liquid
Common.G_func_ice
Common.logistic_function
Common.logistic_function_integral
Common.H2SO4_soln_saturation_vapor_pressure
Common.a_w_xT
Common.a_w_eT
Common.a_w_ice
Common.Chen2022_vel_coeffs
Common.Chen2022_monodisperse_pdf
Common.Chen2022_exponential_pdf
Common.volume_sphere_D
Common.volume_sphere_R
Common.ventilation_factor
```

# Parameters

```@docs
Parameters
Parameters.ParametersType
Parameters.AerosolType
Parameters.AerosolDistributionType
Parameters.CloudCondensateType
Parameters.PrecipitationType
Parameters.TerminalVelocityType
Parameters.Precipitation2MType
Parameters.AirProperties
Parameters.WaterProperties
Parameters.ArizonaTestDust
Parameters.DesertDust
Parameters.AsianDust
Parameters.MiddleEasternDust
Parameters.SaharanDust
Parameters.Illite
Parameters.Kaolinite
Parameters.Feldspar
Parameters.Ferrihydrite
Parameters.Dust
Parameters.Seasalt
Parameters.Sulfate
Parameters.AerosolActivationParameters
Parameters.IceNucleationParameters
Parameters.Frostenberg2023
Parameters.H2SO4SolutionParameters
Parameters.Mohler2006
Parameters.Koop2000
Parameters.H2S04NucleationParameters
Parameters.OrganicNucleationParameters
Parameters.MixedNucleationParameters
Parameters.Parameters0M
Parameters.ParticlePDFSnow
Parameters.ParticlePDFIceRain
Parameters.ParticleMass
Parameters.ParticleArea
Parameters.SnowAspectRatio
Parameters.Ventilation
Parameters.Acnv1M
Parameters.CloudLiquid
Parameters.CloudIce
Parameters.Rain
Parameters.Snow
Parameters.CollisionEff
Parameters.KK2000
Parameters.AcnvKK2000
Parameters.AccrKK2000
Parameters.B1994
Parameters.AcnvB1994
Parameters.AccrB1994
Parameters.TC1980
Parameters.AcnvTC1980
Parameters.AccrTC1980
Parameters.LD2004
Parameters.VarTimescaleAcnv
Parameters.SB2006
Parameters.RainParticlePDF_SB2006
Parameters.RainParticlePDF_SB2006_limited
Parameters.RainParticlePDF_SB2006_notlimited
Parameters.CloudParticlePDF_SB2006
Parameters.AcnvSB2006
Parameters.AccrSB2006
Parameters.SelfColSB2006
Parameters.BreakupSB2006
Parameters.EvaporationSB2006
Parameters.NumberAdjustmentHorn2012
Parameters.Blk1MVelType
Parameters.Blk1MVelTypeRain
Parameters.Blk1MVelTypeSnow
Parameters.SB2006VelType
Parameters.Chen2022VelType
Parameters.Chen2022VelTypeSmallIce
Parameters.Chen2022VelTypeLargeIce
Parameters.Chen2022VelTypeRain
```

# Precipitation susceptibility

```@docs
PrecipitationSusceptibility.precipitation_susceptibility_autoconversion
PrecipitationSusceptibility.precipitation_susceptibility_accretion
```
