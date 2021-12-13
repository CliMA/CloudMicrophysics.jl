```@meta
CurrentModule = CloudMicrophysics
```
# .

# 0-moment microphysics

```@docs
Microphysics_0M
Microphysics_0M.remove_precipitation
```

# 1-moment microphysics

```@docs
Microphysics_1M
Microphysics_1M.AbstractCloudType
Microphysics_1M.AbstractPrecipType
Microphysics_1M.LiquidType
Microphysics_1M.IceType
Microphysics_1M.RainType
Microphysics_1M.SnowType
Microphysics_1M.v0_rai
Microphysics_1M.n0_sno
Microphysics_1M.τ_relax
Microphysics_1M.lambda
Microphysics_1M.unpack_params
Microphysics_1M.terminal_velocity
Microphysics_1M.conv_q_vap_to_q_liq_ice
Microphysics_1M.conv_q_liq_to_q_rai
Microphysics_1M.conv_q_ice_to_q_sno_no_supersat
Microphysics_1M.conv_q_ice_to_q_sno
Microphysics_1M.accretion
Microphysics_1M.accretion_rain_sink
Microphysics_1M.accretion_snow_rain
Microphysics_1M.evaporation_sublimation
Microphysics_1M.snow_melt
```

# Aerosol model

```@docs
AerosolModel
AerosolModel.Mode_B
AerosolModel.Mode_κ
AerosolModel.AbstractAerosolDistribution
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
```
