import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))
FT = Float64
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
aps = CMP.AirProperties(FT)
wps = CMP.WaterProperties(FT)
ip = CMP.IceNucleationParameters(FT)

# Constants
ρₗ = wps.ρw
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)

# Initial conditions
Nₐ = FT(0)
Nₗ = FT(500 * 1e3)
Nᵢ = FT(0)
r₀ = FT(1e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(251)
qᵥ = FT(8.1e-4)
qₗ = Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / 1.2 # 1.2 should be ρₐ
qᵢ = FT(0)
x_sulph = FT(0.01)

# Moisture dependent initial conditions
q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
ts = TD.PhaseNonEquil_pTq(tps, p₀, T₀, q)
ρₐ = TD.air_density(tps, ts)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = qᵥ * p₀ * R_v / Rₐ

Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
r_nuc = FT(0.5 * 1.e-4 * 1e-6)             # assumed size of nucleated particles
w = FT(0.7)                                # updraft speed
α_m = FT(0.5)                              # accomodation coefficient
const_dt = FT(1)                           # model timestep
t_max = FT(1200)
aerosol = CMP.Illite(FT)
ice_nucleation_modes = ["ImmersionFreezing"]
growth_modes = ["Condensation", "Deposition"]
droplet_size_distribution_list = [["Monodisperse"], ["Gamma"]]

# Plotting
fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Supersaturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_ice [g/kg]")
ax4 = MK.Axis(fig[2, 2], ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[3, 1], xlabel = "Time [min]", ylabel = "N_liq")
ax6 = MK.Axis(
    fig[3, 2],
    xlabel = "Time [min]",
    ylabel = "N_ice / N_tot",
    yscale = log10,
)
MK.ylims!(ax6, 1e-6, 1)

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
    )
    # solve ODE
    sol = run_parcel(IC, FT(0), t_max, p)

    DSD = droplet_size_distribution[1]

    # Plot results
    ξ(T) =
        TD.saturation_vapor_pressure(tps, T, TD.Liquid()) /
        TD.saturation_vapor_pressure(tps, T, TD.Ice())
    S_i(T, S_liq) = ξ(T) * S_liq - 1

    MK.lines!(ax1, sol.t / 60, S_i.(sol[3, :], sol[1, :]), label = DSD)
    MK.lines!(ax2, sol.t / 60, sol[3, :])
    MK.lines!(ax3, sol.t / 60, sol[6, :] * 1e3)
    MK.lines!(ax4, sol.t / 60, sol[5, :] * 1e3)

    sol_Nₗ = sol[8, :]
    sol_Nᵢ = sol[9, :]
    sol_qₗ = sol[5, :]
    if DSD == "Monodisperse"
        rₗ = cbrt.(sol_qₗ ./ sol_Nₗ ./ (4 / 3 * FT(π)) / ρₗ * ρₐ)
    end
    if DSD == "Gamma"
        λ = cbrt.(32 .* FT(π) .* sol_Nₗ ./ sol_qₗ * ρₗ / ρₐ)
        rₗ = 2 ./ λ
    end
    MK.lines!(ax5, sol.t / 60, sol_Nₗ)
    MK.lines!(ax6, sol.t / 60, sol_Nᵢ ./ (sol_Nₗ .+ sol_Nᵢ))
end

MK.axislegend(
    ax1,
    framevisible = false,
    labelsize = 12,
    orientation = :horizontal,
    nbanks = 2,
    position = :rb,
)

MK.save("immersion_freezing.svg", fig)
