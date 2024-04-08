
import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)

# Initial conditions
ρₗ = wps.ρw
Nₐ = FT(0)
Nₗ = FT(2000)
Nᵢ = FT(0)
r₀ = FT(1e-6)       # radius of droplets
p₀ = FT(800 * 1e2)
T₀ = FT(251)
qᵥ = FT(8.1e-4)
qₗ = Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2) # 1.2 should be ρₐ
qᵢ = FT(0)
ln_INPC = FT(0)

# Moisture dependent initial conditions
q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
R_v = TD.Parameters.R_v(tps)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

# Simulation parameters passed into ODE solver
w = FT(0.4)                                # updraft speed
const_dt = FT(1)                           # model timestep
t_max = FT(600)
aerosol = CMP.Illite(FT)
heterogeneous = "ABIFM"
condensation_growth = "Condensation"
deposition_growth = "Deposition"
size_distribution_list = ["Monodisperse", "Gamma"]

# Plotting
fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Supersaturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_ice [g/kg]")
ax4 = MK.Axis(fig[2, 2], ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[3, 1], xlabel = "Time [min]", ylabel = "N_liq")
ax6 = MK.Axis(fig[3, 2], xlabel = "Time [min]", ylabel = "N_ice / N_tot")
MK.ylims!(ax6, 1e-6, 1)
# ax7 = MK.Axis(fig[4, 1], ylabel = "r_liq [m]", yscale = log10)
# MK.ylims!(ax7, 1e-6, 0.01)

for DSD in size_distribution_list
    local params = parcel_params{FT}(
        const_dt = const_dt,
        w = w,
        aerosol = aerosol,
        heterogeneous = heterogeneous,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        size_distribution = DSD,
    )
    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)

    # Plot results
    MK.lines!(
        ax1,
        sol.t / 60,
        S_i.(tps, sol[3, :], sol[1, :]) .- 1,
        label = DSD,
    )
    MK.lines!(ax2, sol.t / 60, sol[3, :])
    MK.lines!(ax3, sol.t / 60, sol[6, :] * 1e3)
    MK.lines!(ax4, sol.t / 60, sol[5, :] * 1e3)
    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    MK.lines!(ax5, sol.t / 60, sol_Nₗ)
    MK.lines!(ax6, sol.t / 60, sol_Nᵢ ./ (sol_Nₗ .+ sol_Nᵢ))
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rt,
)

MK.save("immersion_freezing.svg", fig)
