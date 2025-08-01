import OrdinaryDiffEq as ODE
import CairoMakie as MK

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
Nₐ = FT(5e8)
Nₗ = FT(0)
Nᵢ = FT(0)
T₀ = FT(230)
cᵥ₀ = FT(5 * 1e-5)
ln_INPC = FT(0)

# Constants
ρₗ = wps.ρw
ϵₘ = TDI.Rd_over_Rv(tps)
eₛ = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
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
t_max = FT(35)     # total time
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
    Nₐ = Nₐ,
)

# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)

# Plot results
fig = MK.Figure(size = (800, 300), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "Liquid Saturation [-]", xlabel = "Time [s]")
ax2 = MK.Axis(fig[1, 2], ylabel = "N [m^-3]", xlabel = "Time [s]")

MK.lines!(ax1, sol.t, (sol[1, :]), linewidth = 2)
MK.lines!(ax2, sol.t, sol[7, :], label = "N_aero", linewidth = 2, color = :red)
MK.lines!(ax2, sol.t, sol[8, :], label = "N_lcl", linewidth = 2, color = :blue)

MK.axislegend(ax2, framevisible = false, labelsize = 16, position = :rc)

MK.save("Parcel_Aerosol_Activation.svg", fig)
nothing
