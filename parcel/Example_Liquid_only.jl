import OrdinaryDiffEq as ODE
import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams as CP

include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

FT = Float32

# Get free parameters
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
R_v = TDI.Rᵥ(tps)
R_d = TDI.Rd(tps)

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(200 * 1e6)
Nᵢ = FT(0)
r₀ = FT(8e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(273.15 + 7.0)
ln_INPC = FT(0)
e = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
Sₗ = FT(1)
md_v = (p₀ - e) / R_d / T₀
mv_v = e / R_v / T₀
ml_v = Nₗ * 4 / 3 * FT(π) * ρₗ * r₀^3
qᵥ = mv_v / (md_v + mv_v + ml_v)
qₗ = ml_v / (md_v + mv_v + ml_v)
qᵢ = FT(0)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

# Simulation parameters passed into ODE solver
w = FT(10)                                 # updraft speed
const_dt = FT(0.5)                         # model timestep
t_max = FT(20)
liq_size_distribution_list = ["Monodisperse", "Gamma"]
condensation_growth = "Condensation"

# Data from Rogers(1975) Figure 1
# https://www.tandfonline.com/doi/abs/10.1080/00046973.1975.9648397
#! format: off
Rogers_time_supersat = [0.0645, 0.511, 0.883, 1.4, 2.07, 2.72, 3.24, 3.89, 4.53, 5.87, 7.16, 9.79, 16.0, 19.8]
Rogers_supersat = [0.0268, 0.255, 0.393, 0.546, 0.707, 0.805, 0.863, 0.905, 0.938, 0.971, 0.978, 0.963, 0.910, 0.885]
Rogers_time_radius = [0.561, 2, 3.99, 10.7, 14.9, 19.9]
Rogers_radius = [8.0, 8.08, 8.26, 8.91, 9.26, 9.68]

# Setup the plots
fig = MK.Figure(size = (800, 300))
ax1 = MK.Axis(fig[1, 1], xlabel = "Time [s]", ylabel = "Supersaturation [%]")
ax2 = MK.Axis(fig[1, 2], xlabel = "Time [s]", ylabel = "radius [μm]")
ax3 = MK.Axis(fig[1, 3], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
MK.lines!(ax1, Rogers_time_supersat, Rogers_supersat, label = "Rogers_1975", color = :red)
MK.lines!(ax2, Rogers_time_radius, Rogers_radius, color = :red)
#! format: on

for DSD in liq_size_distribution_list
    local params = parcel_params{FT}(
        liq_size_distribution = DSD,
        condensation_growth = condensation_growth,
        const_dt = const_dt,
        w = w,
    )
    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)

    # Plot results
    MK.lines!(ax1, sol.t, (sol[1, :] .- 1) * 100.0, label = DSD)
    MK.lines!(ax3, sol.t, sol[5, :] * 1e3)

    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    # Compute the current air density
    sol_T = sol[3, :]
    sol_p = sol[2, :]
    sol_qᵥ = sol[4, :]
    sol_qₗ = sol[5, :]
    sol_qᵢ = sol[6, :]
    local ρₐ = TDI.air_density.(tps, sol_T, sol_p, sol_qᵥ + sol_qₗ + sol_qᵢ, sol_qₗ, sol_qᵢ)
    # Compute the mean particle size based on the distribution
    distr = sol.prob.p.liq_distr
    moms = distribution_moments.(distr, sol_qₗ, sol_Nₗ, ρₗ, ρₐ)
    local r = similar(sol_T)
    for it in range(1, length(sol_T))
        r[it] = moms[it].r
    end
    MK.lines!(ax2, sol.t, r * 1e6)
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :vertical,
    position = :rb,
)

MK.save("liquid_only_parcel.svg", fig)
nothing
