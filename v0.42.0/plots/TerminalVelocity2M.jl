import CairoMakie as MK
import Random

import ClimaParams as CP
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP

FT = Float64

const rain = CMP.Rain(FT)
const SB2006 = CMP.SB2006(FT)
const SB2006_no_lim = CMP.SB2006(FT; is_limited = false)
const SB2006Vel = CMP.SB2006VelType(FT)

ρ_air = 1.2
# random q and N values for plotting vt vs. mean radius
Random.seed!(1234)
q_rain_random = rand(10000) .* 5e-3
N_rain_random = 10.0 .^ (1.0 .+ 6.0 .* rand(10000))
r_mean =
    ((ρ_air .* q_rain_random ./ N_rain_random) / 1000.0 * 3 / 4 / pi) .^
    (1.0 / 3) .* 1e6

#! format: off

# SB2006 terminal velocities for random combinations of q and N
SB_rain_bN_rm = [CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rain_random[i], ρ_air, N_rain_random[i])[1] for i in 1:length(q_rain_random)]
SB_rain_bM_rm = [CM2.rain_terminal_velocity(SB2006, SB2006Vel, q_rain_random[i], ρ_air, N_rain_random[i])[2] for i in 1:length(q_rain_random)]
SB_rain_bN_rm_nolim = [CM2.rain_terminal_velocity(SB2006_no_lim, SB2006Vel, q_rain_random[i], ρ_air, N_rain_random[i])[1] for i in 1:length(q_rain_random)]
SB_rain_bM_rm_nolim = [CM2.rain_terminal_velocity(SB2006_no_lim, SB2006Vel, q_rain_random[i], ρ_air, N_rain_random[i])[2] for i in 1:length(q_rain_random)]

#! format: on

# SB2006 group velocities vs. mean radius
fig = MK.Figure(size = (900, 450))
ax1 = MK.Axis(fig[1, 1]; xlabel = "r_mean [mm]", ylabel = "terminal velocity [m/s]", title = "Number-weighted")
ax2 = MK.Axis(fig[1, 2]; xlabel = "r_mean [mm]", ylabel = "terminal velocity [m/s]", title = "Mass-weighted")
# rasterize the dense scatters to keep the SVG small
scatter_style = (; markersize = 3, strokewidth = 0, rasterize = 2)
legend_marker = (; markersize = 10)
MK.scatter!(ax1, r_mean ./ 1000, SB_rain_bN_rm; label = "Rain-SB2006 [ND]" => legend_marker, scatter_style...)
MK.scatter!(
    ax1,
    r_mean ./ 1000,
    SB_rain_bN_rm_nolim;
    label = "Rain-SB2006 [ND], no lim" => legend_marker,
    scatter_style...,
)
MK.scatter!(ax2, r_mean ./ 1000, SB_rain_bM_rm; label = "Rain-SB2006 [M]" => legend_marker, scatter_style...)
MK.scatter!(
    ax2,
    r_mean ./ 1000,
    SB_rain_bM_rm_nolim;
    label = "Rain-SB2006 [M], no lim" => legend_marker,
    scatter_style...,
)
MK.limits!(ax1, 0, 2, 0, 5.5)
MK.limits!(ax2, 0, 2, 0, 10)
MK.axislegend(ax1; position = :rb)
MK.axislegend(ax2; position = :rb)
MK.save("2M_terminal_velocity_comparisons.svg", fig)
nothing
