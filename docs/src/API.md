```@meta
CurrentModule = CloudMicrophysics
```
# .

# Non-equilibrium cloud formation

```@docs
MicrophysicsNonEq
MicrophysicsNonEq.τ_relax
MicrophysicsNonEq.conv_q_vap_to_q_liq_ice
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
Microphysics1M.lambda
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
Microphysics2M.LiqRaiRates
Microphysics2M.raindrops_limited_vars
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
```

# P3 scheme

```@docs
P3Scheme
P3Scheme.thresholds
P3Scheme.p3_mass
P3Scheme.distribution_parameter_solver
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

# Heterogeneous ice nucleation
```@docs
HetIceNucleation
HetIceNucleation.dust_activated_number_fraction
HetIceNucleation.deposition_J
HetIceNucleation.ABIFM_J
```

# Homogeneous ice nucleation
```@docs
HomIceNucleation
HomIceNucleation.homogeneous_J
```

# Common utility functions

```@docs
Common
Common.G_func
Common.logistic_function
Common.logistic_function_integral
Common.H2SO4_soln_saturation_vapor_pressure
Common.a_w_xT
Common.a_w_eT
Common.a_w_ice
Common.Chen2022_vel_add
Common.Chen2022_vel_coeffs
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
Parameters.Illite
Parameters.Kaolinite
Parameters.Feldspar
Parameters.Ferrihydrite
Parameters.Seasalt
Parameters.Sulfate
Parameters.AerosolActivationParameters
Parameters.IceNucleationParameters
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
Parameters.ParticlePDF_SB2006
Parameters.AcnvSB2006
Parameters.AccrSB2006
Parameters.SelfColSB2006
Parameters.BreakupSB2006
Parameters.EvaporationSB2006
Parameters.ParametersP3
Parameters.Blk1MVelType
Parameters.Blk1MVelTypeRain
Parameters.Blk1MVelTypeSnow
Parameters.SB2006VelType
Parameters.Chen2022VelType
Parameters.Chen2022VelTypeSnowIce
Parameters.Chen2022VelTypeRain
```

# Precipitation susceptibility

```@docs
PrecipitationSusceptibility.precipitation_susceptibility_autoconversion
PrecipitationSusceptibility.precipitation_susceptibility_accretion
```
