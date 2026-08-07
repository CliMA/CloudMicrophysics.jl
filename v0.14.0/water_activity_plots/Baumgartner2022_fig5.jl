using Roots

import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const CMO = CM.Common
const CMI = CM.HetIceNucleation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
# Boiler plate code to have access to model parameters and constants
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
prs = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(prs)

# Baumgartner at al 2022 Figure 5
# https://acp.copernicus.org/articles/22/65/2022/acp-22-65-2022.pdf
#! format: off
Baumgartner2022_fig5_temp = [182.903, 189.941, 199.9706, 208.4164, 216.3929, 223.2551, 227.5366]
Baumgartner2022_fig5_water_activity = [0.7930, 0.8129, 0.8416, 0.8679, 0.89781, 0.928486, 0.95039]

T_range = range(190, stop = 234, length = 100)

p_sat_liq = [TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) for T in T_range]  # sat vap pressure over pure liq water using TD package
p_sat_ice = [TD.saturation_vapor_pressure(thermo_params, T, TD.Ice()) for T in T_range]     # sat vap pressure over ice using TD package
a_w_ice = p_sat_ice ./ p_sat_liq

radius = 1.0e-4                             # cm
J_crit = 1 / (4 / 3 * pi * radius^3) / 60.0 # cm^-3 s^-1
fun(Delta_a)= - 906.7 + 8502 * Delta_a - 26924 * Delta_a^2 + 29180 * Delta_a^3 - log10(J_crit) # homogeneous J from Koop 2000
initial_guess_Delta_a = 0.34
Delta_a_crit = find_zero(fun, initial_guess_Delta_a)

a_w = [Delta_a_crit + (TD.saturation_vapor_pressure(thermo_params, T, TD.Ice()))/(TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())) for T in T_range]

# various ways Baumgartner calculates vapor pressures
Baumgartner_p_ice(T) = exp(9.550426 - 5723.265 / T + 3.53068*log(T) - 0.00728332*T)
MK_p_liq(T) = exp(54.842763 - 6763.22 / T - 4.21*log(T) + 0.000367*T + tanh(0.0415*(T - 218.8)) * (53.878 - 1331.22 / T - 9.44523*log(T) + 0.014025*T))
Tab_p_liq(T) = exp(100 * (18.4524 - 3505.15788/T - 330918.55/(T^2) + 12725068.26/(T^3)))
Nach_p_liq(T) = exp(74.8727 - 7167.405/T - 7.77107*log(T) + 0.00505*T)

a_w_MK = [Delta_a_crit + (Baumgartner_p_ice(T))/(MK_p_liq(T)) for T in T_range]
a_w_Tab = [Delta_a_crit + (Baumgartner_p_ice(T))/(Tab_p_liq(T)) for T in T_range]
a_w_Nach = [Delta_a_crit + (Baumgartner_p_ice(T))/(Nach_p_liq(T)) for T in T_range]
a_w_Luo = [Delta_a_crit + (Baumgartner_p_ice(T))/(CMO.H2SO4_soln_saturation_vapor_pressure(prs, 0.0, T)) for T in T_range]

fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], title = "Water Activity vs Temperature", ylabel = "a_w [-]", xlabel = "T [K]")
MK.lines!(ax1, Baumgartner2022_fig5_temp, Baumgartner2022_fig5_water_activity, label = "Baumgartner2022", linestyle = :dash, color = :blue)
MK.lines!(ax1, T_range, a_w, label = "CloudMicrophysics", color = :blue)
MK.lines!(ax1, T_range, a_w_MK, label = "Baum_MK", color = :green)
MK.lines!(ax1, T_range, a_w_Nach, label = "Baum_Nach", color = :lightgreen)
MK.lines!(ax1, T_range, a_w_Luo, label = "Baum_Luo", color = :darkgreen)
MK.lines!(ax1, T_range, a_w_ice, label = "a_w_ice", color = :red)
MK.axislegend(position = :lt)
MK.save("Baumgartner2022_fig5.svg", fig)
#! format: on
