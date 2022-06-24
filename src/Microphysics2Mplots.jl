import Plots
using GLMakie

import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

const PL = Plots
const CMT = CloudMicrophysics.CommonTypes
const CM1 = CloudMicrophysics.Microphysics1M
const CM2 = CloudMicrophysics.Microphysics2M
const TD = Thermodynamics
const CP = CLIMAParameters
const CP_planet = CLIMAParameters.Planet

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

const liquid = CMT.LiquidType()
const ice = CMT.IceType()
const rain = CMT.RainType()
const snow = CMT.SnowType()

include(joinpath(pkgdir(CloudMicrophysics), "docs", "src", "Wooddata.jl"))

# Example values
q_liq_range  = range(1e-8, stop=1e-3, length=100)

yvalsKK2000 = [CM2.conv_q_liq_to_q_rai_KK2000(param_set, q_liq, 1.0) for q_liq in q_liq_range]
yvalsB1994 = [CM2.conv_q_liq_to_q_rai_B1994(param_set, q_liq) for q_liq in q_liq_range]
yvalsTC1980 = [CM2.conv_q_liq_to_q_rai_TC1980(param_set, q_liq) for q_liq in q_liq_range]
yvalsLD2004 = [CM2.conv_q_liq_to_q_rai_LD2004(param_set, q_liq) for q_liq in q_liq_range]
yvalsK1969 = [CM1.conv_q_liq_to_q_rai(param_set, q_liq) for q_liq in q_liq_range]

fig = Figure()
ax1 = Axis(fig[1, 1]; yscale = log10)
ylims!(ax1, [1e-13, 1e-5])
l1 = lines!(ax1, q_liq_range * 1e3, yvalsKK2000, color = :red)
l2 = lines!(ax1, q_liq_range * 1e3, yvalsB1994, color = :green)
l3 = lines!(ax1, q_liq_range * 1e3, yvalsTC1980, color = :blue)
l4 = lines!(ax1, q_liq_range * 1e3, yvalsLD2004, color = :purple)
l5 = lines!(ax1, q_liq_range * 1e3, yvalsK1969, color = :black)
l6 = lines!(ax1, KK2000_x, KK2000_y, color = :red, linestyle = :dash)
l7 = lines!(ax1, B1994_x, B1994_y, color = :green, linestyle = :dash)
l8 = lines!(ax1, TC1980_x, TC1980_y, color = :blue, linestyle = :dash)
l9 = lines!(ax1, LD2004_x, LD2004_y, color = :purple, linestyle = :dash)
ax1.xlabel = "q_liq [g/kg]"
ax1.ylabel = "autoconversion rate [kg m^-3 s^-1]"
Legend(fig[1, 2], 
[l1, l2, l3, l4, l5, l6, l7, l8, l9], 
["KK2000", "B1994", "TC1980", "LD2004 (2.25 not 7.5)", "K1969", "Wood_KK2000", "Wood_B1994", "Wood_TC1980", "Wood_LD2004"])
display(fig)
