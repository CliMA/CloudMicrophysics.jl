import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

const CMO = CM.Common
const CMP = CM.Parameters

FT = Float64
tps = TD.Parameters.ThermodynamicsParameters(FT)
H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

T_range = range(190, stop = 234, length = 100)
x = FT(0.1)
#! format: off
p_sol_1 = [CMO.H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x, T) for T in T_range]     # p_sol for concentration x
p_sol_0 = [CMO.H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, 0.0, T) for T in T_range]   # sat vap pressure over pure liq water using p_sol eqn
p_sat_liq = [TD.saturation_vapor_pressure(tps, T, TD.Liquid()) for T in T_range]  # sat vap pressure over pure liq water using TD package
p_sat_ice = [TD.saturation_vapor_pressure(tps, T, TD.Ice()) for T in T_range]     # sat vap pressure over ice using TD package

a_w = p_sol_1 ./ p_sat_liq                 # a_w current parameterization
a_w_alternate = p_sol_1 ./ p_sol_0         # a_w if sat vapor pressure over pure liq water was calculated with p_sol eqn
a_w_ice = p_sat_ice ./ p_sat_liq           # a_ice current parameterization
a_w_ice_alternate = p_sat_ice ./ p_sol_0   # a_ice if sat vapor pressure over pure liq water was calculated with p_sol eqn
a_w_ice_μ = [exp((210368 + 131.438*T - (3.32373e6 /T) - 41729.1*log(T))/(8.31441*T)) for T in T_range]  # a_ice using chemical potential parameterization

# Plotting. NOTE: all ".* -1" is only to flip x axis
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], title = "Temperature vs Water Activity", ylabel = "T [K]", xlabel = "-a_w", limits = ((-1, -0.4), nothing))
MK.lines!(ax1, a_w .* -1, T_range, label = "CM default a_w", color = :blue)
MK.lines!(ax1, a_w_alternate .* -1, T_range, label = "a_w using p(0,T)", linestyle = :dash, color = :blue)
MK.lines!(ax1, a_w_ice .* -1, T_range, label = "CM default a_w_ice", color = :green)
MK.lines!(ax1, a_w_ice_alternate .* -1, T_range, label = "a_w_ice using p(0,T)", color = :green, linestyle = :dash)
MK.lines!(ax1, a_w_ice_μ .* -1, T_range, label = "a_w_ice using μ", color = :lightgreen)
MK.axislegend()

MK.save("T_vs_wateractivity.svg", fig)
#! format: on
