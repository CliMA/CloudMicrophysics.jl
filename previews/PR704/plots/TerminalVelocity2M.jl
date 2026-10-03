import Plots as PL
using Measures

import ClimaParams as CP
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP

FT = Float64

const rain = CMP.Rain(FT)
const SB2006 = CMP.SB2006(FT)
const SB2006_no_lim = CMP.SB2006(FT, false)
const SB2006Vel = CMP.SB2006VelType(FT)

ρ_air = 1.2
# random q and N values for plotting vt vs. mean radius
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

# SB2006 group velocities vs. mean radius
p3 = PL.scatter( r_mean./1000, SB_rain_bN_rm,       ms = 1.0, markerstrokewidth=0, xlabel = "r_mean [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-SB2006 [ND]", margin=10mm)
p3 = PL.scatter!(r_mean./1000, SB_rain_bN_rm_nolim, ms = 1.0, markerstrokewidth=0, xlabel = "r_mean [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-SB2006 [ND], no lim")
p3 = PL.plot!(xlim = [0, 2], ylim = [0, 5.5])

p4 = PL.scatter( r_mean./1000, SB_rain_bM_rm,       ms = 1.0, markerstrokewidth=0, xlabel = "r_mean [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-SB2006 [M]")
p4 = PL.scatter!(r_mean./1000, SB_rain_bM_rm_nolim, ms = 1.0, markerstrokewidth=0, xlabel = "r_mean [mm]", ylabel = "terminal velocity [m/s]", label = "Rain-SB2006 [M], no lim")
p4 = PL.plot!(xlim = [0, 2], ylim = [0, 10])
# save plot
PL.plot(p3, p4, layout = (1, 2), size = (900, 450), dpi = 300)
PL.savefig("2M_terminal_velocity_comparisons.svg")

#! format: on
