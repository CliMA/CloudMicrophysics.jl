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
Nₐ = FT(5e8)
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
r₀ = FT(3e-7)
σ = FT(2)
qₗ = Nₗ * FT(4 / 3 * π) * exp((6 * log(r₀) + 9 * σ^2) / (2))
qᵢ = FT(0)
Sₗ = FT(0.99)
e = Sₗ * eₛ
p₀ = e / cᵥ₀
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]

# Simulation parameters passed into ODE solver
w = FT(1.2)         # updraft speed
const_dt = FT(1)    # model timestep
t_max = FT(100)     # total time
aerosol = CMP.Sulfate(FT)

condensation_growth = "Condensation"
aerosol_act = "AeroAct"     # turn on aerosol activation
aero_σ_g = FT(2.3)
r_nuc = r₀

params = parcel_params{FT}(
    w = w,
    const_dt = const_dt,
    aerosol_act = aerosol_act,
    aerosol = aerosol,
    aero_σ_g = aero_σ_g,
    r_nuc = r_nuc,
    condensation_growth = condensation_growth,
)

# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)

# Plot results
fig = MK.Figure(size = (1000, 800), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "Liquid Saturation [-]")
ax2 = MK.Axis(fig[1, 2], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "N_aero [m^-3]", xlabel = "Time [s]")
ax4 = MK.Axis(fig[2, 2], ylabel = "N_liq [m^-3]", xlabel = "Time [s]")
ax5 = MK.Axis(fig[3, 1], ylabel = "q_vap [g/kg]", xlabel = "Time [s]")
ax6 = MK.Axis(fig[3, 2], ylabel = "q_liq [g/kg]", xlabel = "Time [s]")

MK.lines!(ax1, sol.t, (sol[1, :]))            # liq saturation
MK.lines!(ax2, sol.t, sol[3, :])              # temperature
MK.lines!(ax3, sol.t, sol[7, :])              # N_aero
MK.lines!(ax4, sol.t, sol[8, :])              # N_liq
MK.lines!(ax5, sol.t, sol[4, :] .* 1e3)       # q_vap
MK.lines!(ax6, sol.t, sol[5, :] .* 1e3)       # q_liq


MK.save("Parcel_Aerosol_Activation.svg", fig)
