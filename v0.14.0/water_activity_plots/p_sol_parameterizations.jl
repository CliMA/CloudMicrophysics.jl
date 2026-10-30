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

# our vapor pressure parameterizations
p_sat_liq = [TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) for T in T_range]  # sat vap pressure over pure liq water using TD package
p_sat_ice = [TD.saturation_vapor_pressure(thermo_params, T, TD.Ice()) for T in T_range]     # sat vap pressure over ice using TD package

# vapor pressure parameterizations mentioned in Baumgartner 2022
Baumgartner_p_ice(T) = exp(9.550426 - 5723.265 / T + 3.53068*log(T) - 0.00728332*T)
MK_p_liq(T) = exp(54.842763 - 6763.22 / T - 4.21*log(T) + 0.000367*T + tanh(0.0415*(T - 218.8)) * (53.878 - 1331.22 / T - 9.44523*log(T) + 0.014025*T))
Tab_p_liq(T) = exp(100 * (18.4524 - 3505.15788/T - 330918.55/(T^2) + 12725068.26/(T^3)))
Nach_p_liq(T) = exp(74.8727 - 7167.405/T - 7.77107*log(T) + 0.00505*T)

fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], title = "Saturated Vapor Pressure vs Temperature", ylabel = "Pressure [Pa]", xlabel = "T [K]")
MK.lines!(ax1,T_range, [Baumgartner_p_ice(T) for T in T_range], label = "Baumgartner's p_ice", color = :blue, linestyle = :dash)
MK.lines!(ax1,T_range, p_sat_ice, label = "CM's p_ice", color = :blue)
MK.lines!(ax1,T_range, [MK_p_liq(T) for T in T_range], label = "MK_p_liq", color = :green, linestyle = :dash)
MK.lines!(ax1,T_range, [Nach_p_liq(T) for T in T_range], label = "Nach_p_liq", color = :lightgreen, linestyle = :dash)
MK.lines!(ax1, T_range, [TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid()) for T in T_range], label = "CM's p_liq", color = :green)
MK.axislegend(position = :lt)
MK.save("vap_pressure_vs_T.svg", fig)
#! format: on
