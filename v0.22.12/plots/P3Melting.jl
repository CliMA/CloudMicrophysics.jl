import CairoMakie: CairoMakie, Makie
CairoMakie.activate!(type = "svg")

import Thermodynamics as TD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
FT = Float64

# parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
params = CMP.ParametersP3(FT; slope_law = :constant)
vel = CMP.Chen2022VelType(FT)
aps = CMP.AirProperties(FT)
tps = TD.Parameters.ThermodynamicsParameters(FT)

# model time step (for limiting)
dt = FT(1)

# initial ice content and number concentration
Lᵢ = FT(1e-4)
Nᵢ = FT(2e5)

# tested temperature range
ΔT_range = range(1e-4, stop = 0.025, length = 1000)

# limiters to not melt more mass and number than we have
max_dLdt = [Lᵢ / dt for ΔT in ΔT_range]
max_dNdt = [Nᵢ / dt for ΔT in ΔT_range]

#! format: off

ρₐ1 = FT(1.2)
Fᵣ1 = FT(0.8)
ρᵣ1 = FT(800)
state = P3.get_state(params; F_rim = Fᵣ1, ρ_r = ρᵣ1)
dist = P3.get_distribution_parameters(state; L = Lᵢ, N = Nᵢ)
dLdt1 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ1, dt).dLdt for ΔT in ΔT_range]
dNdt1 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ1, dt).dNdt for ΔT in ΔT_range]

Fᵣ2 = FT(0.2)
state = P3.get_state(params; F_rim = Fᵣ2, ρ_r = ρᵣ1)
dist = P3.get_distribution_parameters(state; L = Lᵢ, N = Nᵢ)
dLdt2 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ1, dt).dLdt for ΔT in ΔT_range]
dNdt2 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ1, dt).dNdt for ΔT in ΔT_range]

ρᵣ2 = FT(200)
state = P3.get_state(params; F_rim = Fᵣ2, ρ_r = ρᵣ2)
dist = P3.get_distribution_parameters(state; L = Lᵢ, N = Nᵢ)
dLdt3 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ1, dt).dLdt for ΔT in ΔT_range]
dNdt3 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ1, dt).dNdt for ΔT in ΔT_range]

ρₐ2 = FT(0.5)
state = P3.get_state(params; F_rim = Fᵣ2, ρ_r = ρᵣ2)
dist = P3.get_distribution_parameters(state; L = Lᵢ, N = Nᵢ)
dLdt4 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ2, dt).dLdt for ΔT in ΔT_range]
dNdt4 = [P3.ice_melt(dist, vel, aps, tps, params.T_freeze .+ ΔT, ρₐ2, dt).dNdt for ΔT in ΔT_range]

# plotting
fig = Makie.Figure(size = (1500, 500), fontsize=22, linewidth=3)

ax1 = Makie.Axis(fig[1, 1]; yscale = log10)
ax2 = Makie.Axis(fig[1, 2]; yscale = log10)

ax1.xlabel = "T [C]"
ax1.ylabel = "ice mass melting rate [g/m3/s]"
ax2.xlabel = "T [C]"
ax2.ylabel = "ice number melting rate [1/cm3/s]"

l_max_dLdt = Makie.lines!(ax1, ΔT_range,  max_dLdt * 1e3,  color = :thistle)
l_max_dNdt = Makie.lines!(ax2, ΔT_range,  max_dNdt * 1e-6, color = :thistle)

l_dLdt1 = Makie.lines!(ax1, ΔT_range,  dLdt1 * 1e3,  color = :skyblue)
l_dNdt1 = Makie.lines!(ax2, ΔT_range,  dNdt1 * 1e-6, color = :skyblue)

l_dLdt2 = Makie.lines!(ax1, ΔT_range,  dLdt2 * 1e3,  color = :blue3)
l_dNdt2 = Makie.lines!(ax2, ΔT_range,  dNdt2 * 1e-6, color = :blue3)

l_dLdt3 = Makie.lines!(ax1, ΔT_range,  dLdt3 * 1e3,  color = :orchid)
l_dNdt3 = Makie.lines!(ax2, ΔT_range,  dNdt3 * 1e-6, color = :orchid)

l_dLdt4 = Makie.lines!(ax1, ΔT_range,  dLdt4 * 1e3,  color = :purple)
l_dNdt4 = Makie.lines!(ax2, ΔT_range,  dNdt4 * 1e-6, color = :purple)

Makie.Legend(
    fig[1, 3],
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
Makie.save("P3_ice_melt.svg", fig)
