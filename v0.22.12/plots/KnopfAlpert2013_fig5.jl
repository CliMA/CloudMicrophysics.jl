import CairoMakie as MK

import ClimaParams
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI
import CloudMicrophysics.Parameters as CMP

FT = Float64
tps = TD.Parameters.ThermodynamicsParameters(FT)
H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
illite = CMP.Illite(FT)    # dust type

# Knopf and Alpert 2013 Figure 5A: Cirrus
# https://doi.org/10.1039/C3FD00035D
temp = [
    228.20357,
    228.33571,
    228.50357,
    228.75000,
    228.92143,
    229.16786,
    229.39643,
    229.52143,
]
KA13_fig5A_J = [
    70170382.86704,
    12101528.74384,
    1277935.32665,
    52710.05764,
    6040.39151,
    293.39697,
    18.97171,
    4.18121,
]

# Our parameterization
T_range = range(228.2, stop = 229.6, length = 100)  # air temperature
T_dew = FT(228.0)               # dew point temperature
x_sulph = FT(0)                 # sulphuric acid concentration in droplets
a_sol = [                       # water activity of solution droplet at equilibrium
    CMO.H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x_sulph, T_dew) /
    TD.saturation_vapor_pressure(tps, T, TD.Liquid()) for T in T_range
]
# water activity of ice
a_ice = [CMO.a_w_ice(tps, T) for T in T_range]

Δa_w = @. max(abs(a_sol - a_ice), FT(0.0))
J_ABIFM = @. CMI.ABIFM_J(illite, Δa_w) * 1e-4 # converted from SI units to cm^-2 s^-1

# Plot results
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(
    fig[1, 1],
    ylabel = "J_het [cm^-2 s^-1]",
    xlabel = "Temp [K]",
    yscale = log10,
)

MK.ylims!(ax1, FT(1), FT(1e8))
MK.lines!(
    ax1,
    temp,
    KA13_fig5A_J,
    label = "KA13 Figure 5A",
    linestyle = :dash,
    color = :green,
)
MK.lines!(ax1, T_range, J_ABIFM, label = "CliMA", color = :green)
fig[1, 2] = MK.Legend(fig, ax1)

MK.save("KnopfAlpert2013_fig5.svg", fig)
