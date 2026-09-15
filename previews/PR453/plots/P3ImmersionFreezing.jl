import CairoMakie as PL
PL.activate!(type = "svg")

import Thermodynamics as TD
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

# thermodynamics parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)

# helper functions
function RH2qₜ(T, RH)
    eᵥ_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    eᵥ = RH * eᵥ_sat
    qᵥ = 1 / (1 - tps.molmass_dryair / tps.molmass_water * (eᵥ - p) / eᵥ)
    qₜ = qᵥ + qₗ + qᵢ
    return qₜ
end
function p2ρ(T, RH)
    return TD.air_density(tps, T, p, TD.PhasePartition(RH2qₜ(T, RH), qₗ, qᵢ))
end

# ambient conditions
Nₗ = FT(500 * 1e6)
qᵢ = FT(1 * 1e-3)
qₗ = FT(1 * 1e-3)
p = FT(800 * 1e2)

# supported aerosol types
dd = CMP.DesertDust(FT)
il = CMP.Illite(FT)

# model time step (for limiting)
dt = FT(1)

# plot data
RH_range = range(0.8, stop = 1.2, length = 1000)
T1 = FT(273.15 - 15)
T2 = FT(273.15 - 35)

#! format: off

# limiters to not nucleate more mass and number than we have in liquid phase
max_dLdt_T1 = [qₗ* p2ρ(T1, RH) / dt for RH in RH_range]
max_dLdt_T2 = [qₗ* p2ρ(T2, RH) / dt for RH in RH_range]
max_dNdt =    [Nₗ              / dt for RH in RH_range]

dLdt_dd_T1 = [P3.het_ice_nucleation(dd, tps, TD.PhasePartition(RH2qₜ(T1, RH), qₗ, qᵢ), Nₗ, RH, T1, p2ρ(T1, RH), dt).dLdt for RH in RH_range]
dNdt_dd_T1 = [P3.het_ice_nucleation(dd, tps, TD.PhasePartition(RH2qₜ(T1, RH), qₗ, qᵢ), Nₗ, RH, T1, p2ρ(T1, RH), dt).dNdt for RH in RH_range]

dLdt_il_T1 = [P3.het_ice_nucleation(il, tps, TD.PhasePartition(RH2qₜ(T1, RH), qₗ, qᵢ), Nₗ, RH, T1, p2ρ(T1, RH), dt).dLdt for RH in RH_range]
dNdt_il_T1 = [P3.het_ice_nucleation(il, tps, TD.PhasePartition(RH2qₜ(T1, RH), qₗ, qᵢ), Nₗ, RH, T1, p2ρ(T1, RH), dt).dNdt for RH in RH_range]


dLdt_dd_T2 = [P3.het_ice_nucleation(dd, tps, TD.PhasePartition(RH2qₜ(T2, RH), qₗ, qᵢ), Nₗ, RH, T2, p2ρ(T2, RH), dt).dLdt for RH in RH_range]
dNdt_dd_T2 = [P3.het_ice_nucleation(dd, tps, TD.PhasePartition(RH2qₜ(T2, RH), qₗ, qᵢ), Nₗ, RH, T2, p2ρ(T2, RH), dt).dNdt for RH in RH_range]

dLdt_il_T2 = [P3.het_ice_nucleation(il, tps, TD.PhasePartition(RH2qₜ(T2, RH), qₗ, qᵢ), Nₗ, RH, T2, p2ρ(T2, RH), dt).dLdt for RH in RH_range]
dNdt_il_T2 = [P3.het_ice_nucleation(il, tps, TD.PhasePartition(RH2qₜ(T2, RH), qₗ, qᵢ), Nₗ, RH, T2, p2ρ(T2, RH), dt).dNdt for RH in RH_range]

# plotting
fig = PL.Figure(size = (1500, 500), fontsize=22, linewidth=3)

ax1 = PL.Axis(fig[1, 1]; yscale = log10)
ax2 = PL.Axis(fig[1, 2]; yscale = log10)

ax1.xlabel = "RH [%]"
ax1.ylabel = "ice mass nucleation rate [g/m3/s]"
ax2.xlabel = "RH [%]"
ax2.ylabel = "ice number nucleation rate [1/cm3/s]"

l_max_dLdt_T1 = PL.lines!(ax1, RH_range * 1e2,  max_dLdt_T1 * 1e3,  color = :thistle)
l_max_dLdt_T2 = PL.lines!(ax1, RH_range * 1e2,  max_dLdt_T2 * 1e3,  color = :thistle)
l_max_dNdt    = PL.lines!(ax2, RH_range * 1e2,  max_dNdt    * 1e-6, color = :thistle)


l_dLdt_dd_T1 = PL.lines!(ax1, RH_range * 1e2,  dLdt_dd_T1 * 1e3,  color = :skyblue)
l_dNdt_dd_T1 = PL.lines!(ax2, RH_range * 1e2,  dNdt_dd_T1 * 1e-6, color = :skyblue)

l_dLdt_dd_T2 = PL.lines!(ax1, RH_range * 1e2,  dLdt_dd_T2 * 1e3,  color = :blue3)
l_dNdt_dd_T2 = PL.lines!(ax2, RH_range * 1e2,  dNdt_dd_T2 * 1e-6, color = :blue3)


l_dLdt_il_T1 = PL.lines!(ax1, RH_range * 1e2,  dLdt_il_T1 * 1e3,  color = :orchid)
l_dNdt_il_T1 = PL.lines!(ax2, RH_range * 1e2,  dNdt_il_T1 * 1e-6, color = :orchid)

l_dLdt_il_T2 = PL.lines!(ax1, RH_range * 1e2,  dLdt_il_T2 * 1e3,  color = :purple)
l_dNdt_il_T2 = PL.lines!(ax2, RH_range * 1e2,  dNdt_il_T2 * 1e-6, color = :purple)

PL.Legend(
    fig[1, 3],
    [l_max_dNdt,
     l_dNdt_dd_T1, l_dNdt_dd_T2,
     l_dNdt_il_T1, l_dNdt_il_T2],
    [
       "limit",
       "T=-15C, desert dust",
       "T=-35C, desert dust",
       "T=-15C, illite",
       "T=-35C, illite",
    ],
    framevisible = false,
)
PL.save("P3_het_ice_nucleation.svg", fig)
