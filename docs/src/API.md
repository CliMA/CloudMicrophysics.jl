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
Microphysics1M.v0_rai
Microphysics1M.n0_sno
Microphysics1M.lambda
Microphysics1M.unpack_params
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
Microphysics2M.conv_q_liq_to_q_rai_KK2000
Microphysics2M.conv_q_liq_to_q_rai_B1994
Microphysics2M.conv_q_liq_to_q_rai_TC1980
Microphysics2M.conv_q_liq_to_q_rai_LD2004
Microphysics2M.accretion_KK2000
Microphysics2M.accretion_B1994
Microphysics2M.accretion_TC1980
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

# Common utility functions

```@docs
Common
Common.G_func
Common.logistic_function
Common.logistic_function_integral
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
```
