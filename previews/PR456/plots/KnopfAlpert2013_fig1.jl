import Plots as PL

import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Common as CO
import CloudMicrophysics.HetIceNucleation as IN
import CloudMicrophysics.Parameters as CMP

FT = Float64
const tps = TD.Parameters.ThermodynamicsParameters(FT)
const H2SO4_prs = CMP.H2SO4SolutionParameters(FT)
const kaolinite = CMP.Kaolinite(FT) # dust type

# Initializing
T_range = range(210, stop = 232, length = 100)  # air temperature
x = 0.1                                         # wt% sulphuric acid in droplets
Δa_w = [CO.a_w_xT(H2SO4_prs, tps, x, T) - CO.a_w_ice(tps, T) for T in T_range]    # difference in solution and ice water activity
J = @. IN.ABIFM_J(kaolinite, Δa_w)      # J in SI units
log10J_converted = @. log10(J * 1e-4)           # converts J into cm^-2 s^-1 and takes log

# Knopf and Alpert 2013 Figure 4A
# https://doi.org/10.1039/C3FD00035D

#! format: off
# data read from Fig 4 in Knopf & Alpert 2013
# using https://automeris.io/WebPlotDigitizer/
KA13_Delta_a_obs = [
    0.13641, 0.16205, 0.21538, 0.23897, 0.24513, 0.24718, 0.25026,
    0.25128, 0.25231, 0.25333, 0.25538, 0.25744, 0.25846, 0.25949,
    0.26051, 0.26051, 0.26462, 0.26462, 0.26872, 0.26974, 0.27077,
    0.27077, 0.27179, 0.27385, 0.27692, 0.27795, 0.27795, 0.27795,
    0.28308, 0.28410, 0.28410, 0.28615, 0.28718, 0.28718, 0.29128,
    0.29128, 0.29231, 0.29333, 0.29744, 0.29744, 0.29744, 0.29949,
    0.30359, 0.30462, 0.30564, 0.30667, 0.31077, 0.31077, 0.31077,
]
KA13_log10J_obs = [
    -3.51880, -3.20301, 2.21053, 2.57143, 2.25564, 3.56391, 3.20301, 2.25564,
    3.78947, 4.42105, 3.51880, 2.84211, 4.15038, 3.24812, 3.78947, 4.37594,
    3.38346, 4.46617, 4.06015, 4.73684, 4.06015, 3.60902, 6.13534, 4.51128,
    4.37594, 4.82707, 4.96241, 5.23308, 3.92481, 5.36842, 5.63910, 5.81955,
    4.60150, 4.96241, 5.50376, 6.00000, 5.14286, 5.77444, 5.41353, 6.09023,
    5.77444, 5.14286, 6.18045, 5.86466, 5.54887, 5.27820, 6.09023, 5.77444, 5.54887,
]
#! format: on

KA13_Delta_a_param = [0.10256, 0.35692, 0.21949]
KA13_log10J_param = [-4.91729, 8.97744, 1.44361]

PL.plot(
    Δa_w,
    log10J_converted,
    label = "CliMA",
    xlabel = "Delta a_w [unitless]",
    ylabel = "log10(J) [cm^-2 s^-1]",
)
PL.scatter!(
    KA13_Delta_a_obs,
    KA13_log10J_obs,
    markercolor = :black,
    label = "paper observations",
)
PL.plot!(
    KA13_Delta_a_param,
    KA13_log10J_param,
    linecolor = :red,
    label = "paper parameterization",
)

PL.savefig("Knopf_Alpert_fig_1.svg")
