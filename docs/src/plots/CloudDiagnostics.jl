import CairoMakie as PL
PL.activate!(type = "svg")

import ClimaParams as CP

import CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.CloudDiagnostics as CMD

FT = Float64

# parameters
wtr = CMP.WaterProperties(FT)
rain = CMP.Rain(FT)
cloud_liquid = CMP.CloudLiquid(FT)
cloud_ice = CMP.CloudIce(FT)

override_file = joinpath(
    pkgdir(CloudMicrophysics),
    "src",
    "parameters",
    "toml",
    "SB2006_limiters.toml",
)
toml_dict = CP.create_toml_dict(FT; override_file)
SB = CMP.SB2006(toml_dict)
SB_no_limiters = CMP.SB2006(toml_dict, false)

ρ_air = FT(1)

# tested rain range for reflectivity plots
q_rain_range = range(0, stop = 5e-3, length = 1000)

#! format: off

Z_1M = [CMD.radar_reflectivity_1M(rain, q_rai, ρ_air) for q_rai in q_rain_range]

Z_2M_100 = [CMD.radar_reflectivity_2M(SB, FT(0), q_rai, FT(0), FT(100), ρ_air) for q_rai in q_rain_range]
Z_2M_10  = [CMD.radar_reflectivity_2M(SB, FT(0), q_rai, FT(0), FT(10), ρ_air) for q_rai in q_rain_range]

Z_2M_100_nolim = [CMD.radar_reflectivity_2M(SB_no_limiters, FT(0), q_rai, FT(0), FT(100), ρ_air) for q_rai in q_rain_range]
Z_2M_10_nolim  = [CMD.radar_reflectivity_2M(SB_no_limiters, FT(0), q_rai, FT(0), FT(10), ρ_air) for q_rai in q_rain_range]

# tested cloud range for effective radiius plots
q_liq_range = range(0, stop = 5e-3, length = 1000)

reff_2M_100  = [CMD.effective_radius_2M(SB, q_liq, FT(0), FT(100), FT(0), ρ_air) for q_liq in q_liq_range]
reff_2M_1000 = [CMD.effective_radius_2M(SB, q_liq, FT(0), FT(1000), FT(0), ρ_air) for q_liq in q_liq_range]

reff_2M_100_nolim  = [CMD.effective_radius_2M(SB_no_limiters, q_liq, FT(0), FT(100), FT(0), ρ_air) for q_liq in q_liq_range]
reff_2M_1000_nolim = [CMD.effective_radius_2M(SB_no_limiters, q_liq, FT(0), FT(1000), FT(0), ρ_air) for q_liq in q_liq_range]

reff_2M_100_LH  = [CMD.effective_radius_Liu_Hallet_97(cloud_liquid, ρ_air, q_liq, FT(100), FT(0), FT(0)) for q_liq in q_liq_range]
reff_2M_1000_LH = [CMD.effective_radius_Liu_Hallet_97(cloud_liquid, ρ_air, q_liq, FT(1000), FT(0), FT(0)) for q_liq in q_liq_range]

# plotting
fig = PL.Figure(size = (1100, 1000), fontsize=22, linewidth=3)

ax1 = PL.Axis(fig[1, 1])
ax2 = PL.Axis(fig[2, 1])

PL.ylims!(ax2, [-10, 70])

ax1.xlabel = "q_liq [g/kg]"
ax1.ylabel = "effective radius [um]"
ax2.xlabel = "q_rai [g/kg]"
ax2.ylabel = "radar reflectivity [dBZ]"

#p_reff_2M_100        = PL.lines!(ax1, q_liq_range * 1e3,  reff_2M_100 * 1e6,  color = :blue, linestyle = :dot)
#p_reff_2M_1000       = PL.lines!(ax1, q_liq_range * 1e3,  reff_2M_1000 * 1e6,  color = :skyblue1, linestyle = :dot)
p_reff_2M_100_nolim  = PL.lines!(ax1, q_liq_range * 1e3,  reff_2M_100_nolim * 1e6,  color = :blue)
p_reff_2M_100_LH     = PL.lines!(ax1, q_liq_range * 1e3,  reff_2M_100_LH * 1e6,     color = :crimson)

p_reff_2M_1000_nolim = PL.lines!(ax1, q_liq_range * 1e3,  reff_2M_1000_nolim * 1e6,  color = :skyblue1)
p_reff_2M_1000_LH    = PL.lines!(ax1, q_liq_range * 1e3,  reff_2M_1000_LH * 1e6,     color = :orange)

p_Z_1M           = PL.lines!(ax2, q_rain_range * 1e3,  Z_1M,  color = :green)

p_Z_2M_100       = PL.lines!(ax2, q_rain_range * 1e3,  Z_2M_100,        color = :skyblue1, linestyle = :dot)
p_Z_2M_100_nolim = PL.lines!(ax2, q_rain_range * 1e3,  Z_2M_100_nolim,  color = :skyblue1)

p_Z_2M_10        = PL.lines!(ax2, q_rain_range * 1e3,  Z_2M_10,        color = :blue, linestyle = :dot)
p_Z_2M_10_nolim  = PL.lines!(ax2, q_rain_range * 1e3,  Z_2M_10_nolim,  color = :blue)

PL.Legend(
    fig[1, 2],
    [p_reff_2M_100_nolim, p_reff_2M_100_LH, p_reff_2M_1000_nolim, p_reff_2M_1000_LH],
    [
       "SB2006 N = 100  1/cm3",
       "LH1997 N = 100  1/cm3",
       "SB2006 N = 1000 1/cm3",
       "LH1997 N = 1000 1/cm3",
    ],
    framevisible = false,
)
PL.Legend(
    fig[2, 2],
    [p_Z_1M, p_Z_2M_100, p_Z_2M_100_nolim, p_Z_2M_10, p_Z_2M_10_nolim],
    [
       "1M",
       "SB2006 N = 100 1/cm3",
       "SB2006 no limiters N = 100 1/cm3",
       "SB2006 N = 10 1/cm3",
       "SB2006 no limiters N = 10 1/cm3",
    ],
    framevisible = false,
)

#! format: on

PL.save("CloudDiagnostics.svg", fig)
