import CairoMakie as PL
PL.activate!(type = "svg")

import Thermodynamics as TD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

# parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
pps = CMP.ParametersP3(FT)
vel = CMP.Chen2022VelType(FT)
aps = CMP.AirProperties(FT)
tps = TD.Parameters.ThermodynamicsParameters(FT)

# model time step (for limiting)
dt = FT(1)

# initial ice content and number concentration
Lᵢ = FT(1e-4)
Nᵢ = FT(2e5)
F_liq = FT(0)

# tested temperature range
ΔT_range = range(1e-4, stop = 0.025, length = 1000)

# limiters to not melt more mass and number than we have
max_dLdt = [Lᵢ / dt for ΔT in ΔT_range]
max_dNdt = [Nᵢ / dt for ΔT in ΔT_range]

#! format: off

ρₐ1 = FT(1.2)
Fᵣ1 = FT(0.8)
ρᵣ1 = FT(800)
dLdt1_rai = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ1, ρᵣ1, F_liq, dt).dLdt_rai for ΔT in ΔT_range]
dLdt1_liq = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ1, ρᵣ1, F_liq, dt).dLdt_liq for ΔT in ΔT_range]
dNdt1 = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ1, ρᵣ1, F_liq, dt).dNdt_ice for ΔT in ΔT_range]

Fᵣ2 = FT(0.2)
dLdt2_rai = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ2, ρᵣ1, F_liq, dt).dLdt_rai for ΔT in ΔT_range]
dLdt2_liq = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ2, ρᵣ1, F_liq, dt).dLdt_liq for ΔT in ΔT_range]
dNdt2 = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ2, ρᵣ1, F_liq, dt).dNdt_ice for ΔT in ΔT_range]

ρᵣ2 = FT(200)
dLdt3_rai = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ2, ρᵣ2, F_liq, dt).dLdt_rai for ΔT in ΔT_range]
dLdt3_liq = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ2, ρᵣ2, F_liq, dt).dLdt_liq for ΔT in ΔT_range]
dNdt3 = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ1, Fᵣ2, ρᵣ2, F_liq, dt).dNdt_ice for ΔT in ΔT_range]

ρₐ2 = FT(0.5)
dLdt4_rai = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ2, Fᵣ2, ρᵣ2, F_liq, dt).dLdt_rai for ΔT in ΔT_range]
dLdt4_liq = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ2, Fᵣ2, ρᵣ2, F_liq, dt).dLdt_liq for ΔT in ΔT_range]
dNdt4 = [P3.ice_melt(pps, vel, aps, tps, Lᵢ, Nᵢ, pps.T_freeze .+ ΔT, ρₐ2, Fᵣ2, ρᵣ2, F_liq, dt).dNdt_ice for ΔT in ΔT_range]

# plotting
fig = PL.Figure(size = (2000, 1500), fontsize=24, linewidth=3)

ax1 = PL.Axis(fig[1, 1]; yscale = log10)
ax2 = PL.Axis(fig[1, 2]; yscale = log10)
ax3 = PL.Axis(fig[2, 1]; yscale = log10)
ax4 = PL.Axis(fig[2, 2]; yscale = log10)

ax1.xlabel = "T [C]"
ax1.ylabel = "total ice mass melting rate [g/m3/s]"
ax2.xlabel = "T [C]"
ax2.ylabel = "ice number melting rate [1/cm3/s]"
ax3.xlabel = "T [C]"
ax3.ylabel = "rain mass source [g/m3/s]"
ax4.xlabel = "T [C]"
ax4.ylabel = "liquid mass content source [g/m3/s]"

for ax in [ax1, ax3, ax4]
    l_max_dLdt = PL.lines!(ax, ΔT_range,  max_dLdt * 1e3,  color = :thistle)
end
l_max_dNdt = PL.lines!(ax2, ΔT_range,  max_dNdt * 1e-6, color = :thistle)

l_dLdt1 = PL.lines!(ax1, ΔT_range,  (dLdt1_rai + dLdt1_liq) * 1e3,  color = :skyblue)
l_dNdt1 = PL.lines!(ax2, ΔT_range,  dNdt1 * 1e-6, color = :skyblue)
l_dLdt1_rai = PL.lines!(ax3, ΔT_range,  (dLdt1_rai) * 1e3,  color = :skyblue)
l_dLdt1_liq = PL.lines!(ax4, ΔT_range,  (dLdt1_liq) * 1e3,  color = :skyblue)


l_dLdt2 = PL.lines!(ax1, ΔT_range,  (dLdt2_rai + dLdt2_liq) * 1e3,  color = :blue3)
l_dNdt2 = PL.lines!(ax2, ΔT_range,  dNdt2 * 1e-6, color = :blue3)
l_dLdt2_rai = PL.lines!(ax3, ΔT_range,  (dLdt2_rai) * 1e3,  color = :blue3)
l_dLdt2_liq = PL.lines!(ax4, ΔT_range,  (dLdt2_liq) * 1e3,  color = :blue3)

l_dLdt3 = PL.lines!(ax1, ΔT_range,  (dLdt3_rai + dLdt3_liq) * 1e3,  color = :orchid)
l_dNdt3 = PL.lines!(ax2, ΔT_range,  dNdt3 * 1e-6, color = :orchid)
l_dLdt3_rai = PL.lines!(ax3, ΔT_range,  (dLdt3_rai) * 1e3,  color = :orchid)
l_dLdt3_liq = PL.lines!(ax4, ΔT_range,  (dLdt3_liq) * 1e3,  color = :orchid)

l_dLdt4 = PL.lines!(ax1, ΔT_range,  (dLdt4_rai + dLdt4_liq) * 1e3,  color = :purple)
l_dNdt4 = PL.lines!(ax2, ΔT_range,  dNdt4 * 1e-6, color = :purple)
l_dLdt4_rai = PL.lines!(ax3, ΔT_range,  (dLdt4_rai) * 1e3,  color = :purple)
l_dLdt4_liq = PL.lines!(ax4, ΔT_range,  (dLdt4_liq) * 1e3,  color = :purple)


PL.Legend(
    fig[1:2, 3],
    [l_max_dNdt, l_dNdt1, l_dNdt2, l_dNdt3, l_dNdt4],
    [
       "limit",
       "ρₐ=1.2 kg/m3, Fᵣ=0.8, ρᵣ=800kg/m3",
       "ρₐ=1.2 kg/m3, Fᵣ=0.2, ρᵣ=800kg/m3",
       "ρₐ=1.2 kg/m3, Fᵣ=0.2, ρᵣ=200kg/m3",
       "ρₐ=0.5 kg/m3, Fᵣ=0.2, ρᵣ=200kg/m3",
    ],
    framevisible = false,
)
PL.resize_to_layout!(fig)
PL.save("P3_ice_melt.svg", fig)
