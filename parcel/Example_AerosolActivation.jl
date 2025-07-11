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
Nₗ = FT(1e8)
Nᵢ = FT(1e8)
T₀ = FT(230)
cᵥ₀ = FT(5 * 1e-5)
ln_INPC = FT(0)

# Constants
ρₗ = wps.ρw
ρᵢ = wps.ρi
ϵₘ = TDI.Rd_over_Rv(tps)
eₛ = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
qᵥ = ϵₘ / (ϵₘ - 1 + 1 / cᵥ₀)
Sₗ = FT(1.0)
e = Sₗ * eₛ
p₀ = e / cᵥ₀
ρ_air = TDI.air_density(tps, T₀, p₀, qᵥ, FT(0), FT(0))
rₗ = FT(6.333333e-6)
qₗ = Nₗ * FT(4 / 3 * π) * (rₗ)^3 * ρₗ / ρ_air
rᵢ = FT(10e-6)
qᵢ = Nᵢ * FT(4 / 3 * π) * (rᵢ)^3 * ρᵢ / ρ_air
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC, FT(0)]

# Simulation parameters passed into ODE solver
w = FT(1.2)         # updraft speed
const_dt = FT(1e-2)    # model timestep
t_max = FT(60)     # total time
aerosol = CMP.Sulfate(FT)

condensation_growth = "Condensation"
deposition_growth = "Deposition"
aerosol_act = "AeroAct"     # turn on aerosol activation
aero_σ_g = FT(2.3)
r_nuc = FT(4e-8)

params = parcel_params{FT}(
    w = w,
    const_dt = const_dt,
    aerosol_act = aerosol_act,
    aerosol = aerosol,
    aero_σ_g = aero_σ_g,
    r_nuc = r_nuc,
    condensation_growth = condensation_growth,
    deposition_growth = deposition_growth,
    Nₐ = Nₐ,
    qₗ₀ = qₗ,
    Nₗ₀ = Nₗ,
    liq_size_distribution = "MonodisperseMix",
)

# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)

# Plot results
fig = MK.Figure(size = (800, 600), fontsize = 20)
ax1 = MK.Axis(fig[1, 1], ylabel = "Liquid Saturation [-]", xlabel = "Time [s]")
ax2 = MK.Axis(fig[1, 2], ylabel = "N [m^-3]", xlabel = "Time [s]")
ax3 = MK.Axis(fig[2, 1], ylabel = "T [k]", xlabel = "Time [s]")
ax4 = MK.Axis(fig[2, 2], ylabel = "q_liq [g kg^-1]", xlabel = "Time [s]")

MK.lines!(ax1, sol.t, sol[1, :], linewidth = 2)
MK.lines!(ax2, sol.t, sol[7, :], label = "N_aero", linewidth = 2, color = :red)
MK.lines!(ax2, sol.t, sol[8, :], label = "N_liq", linewidth = 2, color = :blue)
MK.lines!(ax3, sol.t, sol[3, :], label = "T", linewidth = 2, color = :blue)
MK.lines!(ax4, sol.t, sol[5, :] * 1e3, label = "q_liq", linewidth = 2, color = :blue)

MK.axislegend(ax2, framevisible = false, labelsize = 16, position = :rc)

MK.save("Parcel_Aerosol_Activation.png", fig)

# ARG evaluation
S_max = maximum(sol[1, :])-FT(1)
ad = AM.Mode_κ(
        params.r_nuc,
        params.aero_σ_g,
        params.Nₐ,
        (FT(1.0),),
        (FT(1.0),),
        (params.aerosol.M,),
        (params.aerosol.κ,),
    )
all_ad = AM.AerosolDistribution((ad,))
S_max_ARG = AA.max_supersaturation(params.aap, all_ad, params.aps, params.tps, T₀, p₀, w, qᵥ+qₗ+qᵢ, qₗ, qᵢ)
error_ARG = abs(S_max_ARG - S_max) / S_max_ARG * 100
S_max_mod = AA.max_supersaturation(params.aap, all_ad, params.aps, params.tps, T₀, p₀, w, qᵥ+qₗ+qᵢ, qₗ, qᵢ, Nₗ, Nᵢ)
error_mod = abs(S_max_mod - S_max) / S_max_ARG * 100

@show S_max
@show S_max_ARG, error_ARG
@show S_max_mod, error_mod

nothing
