import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Distributions as DS

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32

# Get free parameters
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)

# Initial conditions
#Nₐ = FT(300 * 1e6)
Nₐ = FT(0)
Nₗ = FT(300 * 1e6)
Nᵢ = FT(0)
T₀ = FT(190)
cᵥ₀ = FT(5 * 1e-6)
ln_INPC = FT(0)

# Constants
ρₗ = wps.ρw
ϵₘ = TDI.Rd_over_Rv(tps)
eₛ = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
qᵥ = ϵₘ / (ϵₘ - 1 + 1 / cᵥ₀)
# Compute qₗ assuming that initially droplets are lognormally distributed
# with N(r₀, σ). We are not keeping that size distribution assumption
# in the simulation
r₀ = FT(25 * 1e-9)
σ = FT(2)
qₗ = Nₗ * FT(4 / 3 * π) * exp((6 * log(r₀) + 9 * σ^2) / (2))
qᵢ = FT(0)
Sᵢ = FT(1.55)
Sₗ = Sᵢ / ξ(tps, T₀)
e = Sₗ * eₛ
p₀ = e / cᵥ₀
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

# Simulation parameters passed into ODE solver
w = FT(1)           # updraft speed
const_dt = FT(0.01) # model timestep
t_max = FT(120)     # total time
homogeneous = "ABHOM"     # homogeneous freezing only
deposition_growth = "Deposition"

# Data from Jensen(2022) Figure 1
# https://agupubs.onlinelibrary.wiley.com/doi/10.1029/2022JD036535
#! format: off
Jensen_t_sat = [0, 62.71, 70.52, 76.87, 82.4, 84.84, 88.1, 92, 96.07, 100.63, 105.35, 112.51, 119.83]
Jensen_sat = [1.55, 1.694, 1.7107, 1.7208, 1.725, 1.726, 1.7259, 1.722, 1.715, 1.702, 1.686, 1.653, 1.6126]
Jensen_t_T = [0, 120]
Jensen_T = [190, 189]
Jensen_t_ICNC = [0.217, 42.69, 50.02, 54.41, 58.97, 65.316, 72.477, 82.08, 92.658, 94.123, 95.5877, 119.84]
Jensen_ICNC = [0, 0, 0.282, 0.789, 1.804, 4.1165, 7.218, 12.12, 16.35, 16.8, 16.97, 17.086]
#! format: on

fig = MK.Figure(size = (800, 300), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "Saturation", xlabel = "Time [s]")
ax2 = MK.Axis(fig[1, 2], ylabel = "ICNC [cm^-3]", xlabel = "Time [s]")

MK.lines!(
    ax1,
    Jensen_t_sat,
    Jensen_sat,
    label = "Jensen et al 2022",
    color = :green,
    linewidth = 2,
)
MK.lines!(
    ax2,
    Jensen_t_ICNC,
    Jensen_ICNC,
    color = :green,
    label = "Jensen et al 2022",
    linewidth = 2,
)

params = parcel_params{FT}(
    w = w,
    const_dt = const_dt,
    homogeneous = homogeneous,
    deposition_growth = deposition_growth,
)

# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)

# Plot results
MK.lines!(
    ax1,
    sol.t,
    S_i.(tps, sol[3, :], (sol[1, :])),
    label = "ice",
    color = :blue,
    linewidth = 2,
)
MK.lines!(
    ax1,
    sol.t,
    (sol[1, :]),
    label = "liquid",
    color = :red,
    linewidth = 2,
) # liq saturation
MK.lines!(ax2, sol.t, sol[9, :] * 1e-6, label = "CM.jl", linewidth = 2) # ICNC

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 16,
    orientation = :vertical,
    position = :lc,
)
MK.axislegend(
    ax2,
    framevisible = false,
    labelsize = 16,
    orientation = :vertical,
    position = :lc,
)

MK.save("Jensen_et_al_2022.svg", fig)
nothing
