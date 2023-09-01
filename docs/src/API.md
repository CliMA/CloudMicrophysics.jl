```@meta
CurrentModule = CloudMicrophysics
```
# .

# Parameters

```@docs
Parameters.CloudMicrophysicsParameters0M
```

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
Microphysics1M.v0
Microphysics1M.n0
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
P3Scheme.α_va
P3Scheme.D_th
P3Scheme.thresholds
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
Common.Delta_a_w
Common.Chen2022_snow_ice_coeffs
Common.Chen2022_vel_add
Common.Chen2022_vel_coeffs
```

# Common utility types

```@docs
CommonTypes
CommonTypes.AbstractAerosolDistribution
CommonTypes.AbstractCloudType
CommonTypes.AbstractPrecipType
CommonTypes.LiquidType
CommonTypes.IceType
CommonTypes.RainType
CommonTypes.SnowType
CommonTypes.Abstract2MPrecipType
CommonTypes.KK2000Type
CommonTypes.B1994Type
CommonTypes.TC1980Type
CommonTypes.LD2004Type
CommonTypes.SB2006Type
CommonTypes.AbstractTerminalVelocityType
CommonTypes.Blk1MVelType
CommonTypes.SB2006VelType
CommonTypes.Chen2022Type
CommonTypes.AbstractAerosolType
CommonTypes.ArizonaTestDustType
CommonTypes.DesertDustType
CommonTypes.KaoliniteType
CommonTypes.IlliteType
```

# Precipitation susceptibility

```@docs
PrecipitationSusceptibility.precipitation_susceptibility_autoconversion
PrecipitationSusceptibility.precipitation_susceptibility_accretion
```
