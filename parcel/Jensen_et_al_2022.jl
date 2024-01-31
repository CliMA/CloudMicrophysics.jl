import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP
import Distributions as DS
import CloudMicrophysics.Parameters as CMP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))

FT = Float64

# Get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)
aps = CMP.AirProperties(FT)
ip = CMP.IceNucleationParameters(FT)
H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

# Constants
ρₗ = wps.ρw
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v

# Saturation conversion equations
ξ(T) =
    TD.saturation_vapor_pressure(tps, T, TD.Liquid()) /
    TD.saturation_vapor_pressure(tps, T, TD.Ice())
S_i(T, S_liq) = @. S_liq * ξ(T) # saturation ratio

# Initial conditions
Nₐ = FT(300 * 1e6)
Nₗ = FT(0)
Nᵢ = FT(0)
r₀ = FT(25 * 1e-9)
T₀ = FT(190)
cᵥ₀ = FT(5 * 1e-6)
x_sulph = FT(0)

# saturation and partial pressure
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
qᵥ = ϵₘ / (ϵₘ - 1 + 1 / cᵥ₀)
qₗ = Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / 1.2
qᵢ = FT(0)
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
R_a = TD.gas_constant_air(tps, q)
Sᵢ = FT(1.55)
Sₗ = Sᵢ / ξ(T₀)
e = Sₗ * eₛ
p₀ = e / cᵥ₀

## Moisture dependent initial conditions
ts = TD.PhaseNonEquil_pTq(tps, p₀, T₀, q)

IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
r_nuc = r₀          # assumed size of nucleated particles
w = FT(1)           # updraft speed
α_m = FT(1)         # accomodation coefficient
const_dt = FT(0.01) # model timestep
t_max = FT(120)     # total time
aerosol = []        # aerosol type. fill in if immersion or deposition freezing
R_mode_liq = r_nuc  # lognormal distribution mode radius
σ = FT(2)           # lognormal distribution std deviation
ice_nucleation_modes = ["HomogeneousFreezing"]     # homogeneous freezing only
growth_modes = ["Deposition"]
droplet_size_distribution_list = [["Jensen_Example"]]

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

fig = MK.Figure(resolution = (1000, 1000))
ax1 = MK.Axis(fig[1, 1], ylabel = "Saturation")
ax2 = MK.Axis(fig[3, 1], xlabel = "Time [s]", ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_vap [g/kg]")
ax4 = MK.Axis(fig[2, 2], xlabel = "Time [s]", ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[1, 2], ylabel = "ICNC [cm^-3]")
ax6 = MK.Axis(fig[3, 2], ylabel = "q_ice [g/kg]")

MK.ylims!(ax2, 188.5, 190)
MK.xlims!(ax2, -5, 150)
MK.xlims!(ax3, -5, 150)
MK.xlims!(ax4, -5, 150)

#! format: off
MK.lines!(ax1, Jensen_t_sat, Jensen_sat, label = "Jensen et al 2022", color = :green)
MK.lines!(ax2, Jensen_t_T, Jensen_T, color = :green, label = "Jensen et al 2022")
MK.lines!(ax5, Jensen_t_ICNC, Jensen_ICNC, color = :green, label = "Jensen et al 2022")

for droplet_size_distribution in droplet_size_distribution_list
    p = (;
        wps,
        aps,
        tps,
        ip,
        const_dt,
        r_nuc,
        w,
        α_m,
        aerosol,
        ice_nucleation_modes,
        growth_modes,
        droplet_size_distribution,
        R_mode_liq,
        σ,
    )
    # solve ODE
    sol = run_parcel(IC, FT(0), t_max, p)

    DSD = droplet_size_distribution[1]

    # Plot results
    MK.lines!(ax1, sol.t, S_i(sol[3, :], (sol[1, :])), label = "ice", color = :blue)
    MK.lines!(ax1, sol.t, (sol[1, :]), label = "liquid", color = :red) # liq saturation
    MK.lines!(ax2, sol.t, sol[3, :], label = "CM.jl")                  # temperature
    MK.lines!(ax3, sol.t, sol[4, :] * 1e3) # q_vap
    MK.lines!(ax4, sol.t, sol[5, :] * 1e3) # q_liq
    MK.lines!(ax5, sol.t, sol[9, :] * 1e-6, label = "CM.jl")           # ICNC
    MK.lines!(ax6, sol.t, sol[6, :] * 1e3) # q_ice

    MK.axislegend(
        ax1,
        framevisible = false,
        labelsize = 12,
        orientation = :horizontal,
        nbanks = 2,
        position = :lc,
    )
    MK.axislegend(
        ax5,
        framevisible = false,
        labelsize = 12,
        orientation = :horizontal,
        nbanks = 2,
        position = :lc,
    )
    MK.axislegend(
        ax2,
        framevisible = false,
        labelsize = 12,
        orientation = :horizontal,
        nbanks = 2,
        position = :lc,
    )
#! format: on
end

MK.save("Jensen_et_al_2022.svg", fig)
