import RootSolvers as RS
import CairoMakie as MK

import CloudMicrophysics as CM
import Thermodynamics as TD
import CLIMAParameters as CP

const CMO = CM.Common
const CMP = CM.Parameters

FT = Float64
tps = TD.Parameters.ThermodynamicsParameters(FT)
H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

# Baumgartner at al 2022 Figure 5
# https://acp.copernicus.org/articles/22/65/2022/acp-22-65-2022.pdf
#! format: off
BG2022_fig5_T = [182.903, 189.941, 199.9706, 208.4164, 216.3929, 223.2551, 227.5366]
BG2022_fig5_aw = [0.7930, 0.8129, 0.8416, 0.8679, 0.89781, 0.928486, 0.95039]
#! format: on

T_range = range(190, stop = 234, length = 100)
# sat vap pressure over pure liq water using TD package
p_sat_liq = [TD.saturation_vapor_pressure(tps, T, TD.Liquid()) for T in T_range]
# sat vap pressure over ice using TD package
p_sat_ice = [TD.saturation_vapor_pressure(tps, T, TD.Ice()) for T in T_range]
a_w_ice = p_sat_ice ./ p_sat_liq

radius = 1.0e-4                             # cm
J_crit = 1 / (4 / 3 * pi * radius^3) / 60.0 # cm^-3 s^-1
# homogeneous J from Koop 2000
fun(Δa) = -906.7 + 8502 * Δa - 26924 * Δa^2 + 29180 * Δa^3 - log10(J_crit)
sol = RS.find_zero(fun, RS.SecantMethod(FT(0), FT(0.5)), RS.CompactSolution())
Δa_crit = sol.root
a_w = Δa_crit .+ a_w_ice

# various ways Baumgartner calculates vapor pressures
Baumgartner_p_ice(T) =
    exp(9.550426 - 5723.265 / T + 3.53068 * log(T) - 0.00728332 * T)
MK_p_liq(T) = exp(
    54.842763 - 6763.22 / T - 4.21 * log(T) +
    0.000367 * T +
    tanh(0.0415 * (T - 218.8)) *
    (53.878 - 1331.22 / T - 9.44523 * log(T) + 0.014025 * T),
)
Tab_p_liq(T) = exp(
    100 * (18.4524 - 3505.15788 / T - 330918.55 / (T^2) + 12725068.26 / (T^3)),
)
Nach_p_liq(T) = exp(74.8727 - 7167.405 / T - 7.77107 * log(T) + 0.00505 * T)

a_w_MK = [Δa_crit + (Baumgartner_p_ice(T)) / (MK_p_liq(T)) for T in T_range]
a_w_Tab = [Δa_crit + (Baumgartner_p_ice(T)) / (Tab_p_liq(T)) for T in T_range]
a_w_Nach = [Δa_crit + (Baumgartner_p_ice(T)) / (Nach_p_liq(T)) for T in T_range]
a_w_Luo = [
    Δa_crit +
    (Baumgartner_p_ice(T)) /
    (CMO.H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, 0.0, T)) for
    T in T_range
]

#! format: off
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], title = "Water Activity vs Temperature", ylabel = "a_w [-]", xlabel = "T [K]")
MK.lines!(ax1, BG2022_fig5_T, BG2022_fig5_aw, label = "Baumgartner2022", linestyle = :dash, color = :blue)
MK.lines!(ax1, T_range, a_w, label = "CloudMicrophysics", color = :blue)
MK.lines!(ax1, T_range, a_w_MK, label = "Baum_MK", color = :green)
MK.lines!(ax1, T_range, a_w_Nach, label = "Baum_Nach", color = :lightgreen)
MK.lines!(ax1, T_range, a_w_Luo, label = "Baum_Luo", color = :darkgreen)
MK.lines!(ax1, T_range, a_w_ice, label = "a_w_ice", color = :red)
MK.axislegend(position = :lt)
MK.save("Baumgartner2022_fig5.svg", fig)

fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], title = "Saturated Vapor Pressure vs Temperature", ylabel = "Pressure [Pa]", xlabel = "T [K]")
MK.lines!(ax1, T_range, [Baumgartner_p_ice(T) for T in T_range], label = "Baumgartner's p_ice", color = :blue, linestyle = :dash)
MK.lines!(ax1, T_range, p_sat_ice, label = "CM's p_ice", color = :blue)
MK.lines!(ax1, T_range, [MK_p_liq(T) for T in T_range], label = "MK_p_liq", color = :green, linestyle = :dash)
MK.lines!(ax1, T_range, [Nach_p_liq(T) for T in T_range], label = "Nach_p_liq", color = :lightgreen, linestyle = :dash)
MK.lines!(ax1, T_range, p_sat_liq, label = "CM's p_liq", color = :green)
MK.axislegend(position = :lt)
MK.save("vap_pressure_vs_T.svg", fig)
#! format: on
