import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)

# Initial conditions
Nₐ = FT(1 * 1e6)
Nₗ = FT(0)
Nᵢ = FT(0)
T₀ = FT(230)
cᵥ₀ = FT(5 * 1e-5)
ln_INPC = FT(0)

# Constants
ρₗ = wps.ρw
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
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
const_dt = FT(1)    # model timestep
t_max = FT(30)      # total time
aerosol = CMP.Sulfate(FT)

aerosol_act = "AeroAct"     # turn on aerosol activation
aero_σ_g = FT(2.3)

params = parcel_params{FT}(
    w = w,
    const_dt = const_dt,
    aerosol_act = aerosol_act,
    aerosol = aerosol,
    aero_σ_g = aero_σ_g,
)

# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)

# Plot results
fig = MK.Figure(size = (1000, 1000), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "Saturation")
ax2 = MK.Axis(fig[1, 2], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_vap [g/kg]")
ax4 = MK.Axis(fig[2, 2], ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[3, 1], ylabel = "N_aero [m^-3]", xlabel = "Time [s]")
ax6 = MK.Axis(fig[3, 2], ylabel = "N_liq [m^-3]", xlabel = "Time [s]")

MK.lines!(
    ax1,
    sol.t,
    S_i.(tps, sol[3, :], (sol[1, :])),
    label = "ice",
    color = :blue,
)
MK.lines!(ax1, sol.t, (sol[1, :]), label = "liquid", color = :red) # liq saturation
MK.lines!(ax2, sol.t, sol[3, :])                                   # temperature
MK.lines!(ax3, sol.t, sol[4, :] * 1e3)                             # q_vap
MK.lines!(ax4, sol.t, sol[5, :] * 1e3)                             # q_liq
MK.lines!(ax5, sol.t, sol[7, :])                                   # N_aero
MK.lines!(ax6, sol.t, sol[8, :] * 1e3)                             # q_ice

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 20,
    orientation = :horizontal,
    nbanks = 2,
    position = :rc,
)

MK.save("Parcel_Aerosol_Activation.svg", fig)
