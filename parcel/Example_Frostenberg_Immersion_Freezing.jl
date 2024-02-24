import OrdinaryDiffEq as ODE
import CairoMakie as MK
import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP
import Random as RAND
random_seeds = [0, 1234, 5678, 3443]

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TD.Parameters.ThermodynamicsParameters(FT)
wps = CMP.WaterProperties(FT)

# Initial conditions
ρₗ = wps.ρw
Nₐ = FT(0)
Nₗ = FT(500 * 1e3)
Nᵢ = FT(0)
r₀ = FT(1e-6)
p₀ = FT(800 * 1e2)
T₀ = FT(251)
qᵥ = FT(8.1e-4)
qₗ = Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2) # 1.2 should be ρₐ
qᵢ = FT(0)
x_sulph = FT(0.01)

# Moisture dependent initial conditions
q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
R_v = TD.Parameters.R_v(tps)
Rₐ = TD.gas_constant_air(tps, q)
eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
e = eᵥ(qᵥ, p₀, Rₐ, R_v)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, x_sulph]

# Simulation parameters passed into ODE solver
w = FT(0.7)                                # updraft speed
const_dt = FT(1)                           # model timestep
t_max = FT(1200)
aerosol = CMP.Illite(FT)
condensation_growth = "Condensation"
deposition_growth = "Deposition"
DSD = "Monodisperse"

# Plotting
fig = MK.Figure(resolution = (900, 700))
ax1 = MK.Axis(fig[1, 1], ylabel = "Ice Supersaturation [-]")
ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
ax3 = MK.Axis(fig[2, 1], ylabel = "q_ice [g/kg]")
ax4 = MK.Axis(fig[2, 2], ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[3, 1], xlabel = "Time [min]", ylabel = "N_liq")
ax6 = MK.Axis(fig[3, 2], xlabel = "Time [min]", ylabel = "N_ice")

colors = [:blue, :green, :darkorange]

function plot_results(
    sol_t,
    sol,
    variable,
    label = "",
    unit = "",
    linestyle = :solid,
    color = :black,
)

    MK.lines!(
        ax1,
        sol_t / 60,
        S_i.(tps, sol[3, :], sol[1, :]) .- 1,
        label = label * string(variable) * unit,
        linestyle = linestyle,
        color = color,
    )
    MK.lines!(
        ax2,
        sol_t / 60,
        sol[3, :],
        label = label * string(variable) * unit,
        linestyle = linestyle,
        color = color,
    )
    MK.lines!(
        ax3,
        sol_t / 60,
        sol[6, :] * 1e3,
        label = label * string(variable) * unit,
        linestyle = linestyle,
        color = color,
    )
    MK.lines!(
        ax4,
        sol_t / 60,
        sol[5, :] * 1e3,
        label = label * string(variable) * unit,
        linestyle = linestyle,
        color = color,
    )
    MK.lines!(
        ax5,
        sol_t / 60,
        sol[8, :],
        label = label * string(variable) * unit,
        linestyle = linestyle,
        color = color,
    )
    MK.lines!(
        ax6,
        sol_t / 60,
        sol[9, :],
        label = label * string(variable) * unit,
        linestyle = linestyle,
        color = color,
    )

    MK.axislegend(ax2, position = :lt)
end

#Frostenberg_mean
params = parcel_params{FT}(
    const_dt = const_dt,
    w = w,
    aerosol = aerosol,
    heterogeneous = "Frostenberg_mean",
    condensation_growth = condensation_growth,
    deposition_growth = deposition_growth,
    size_distribution = DSD,
)
# solve ODE
sol = run_parcel(IC, FT(0), t_max, params)
plot_results(sol.t, sol, "mean")


# Frostenberg_random with different drawing frequencies
drawing_interval_range = range(FT(1), stop = FT(5), length = 3)

for (drawing_interval, color) in zip(drawing_interval_range, colors)

    # creating an ensamble of solutions
    solutions = []
    for random_seed in random_seeds

        RAND.seed!(trunc(Int, random_seed)) #set the random seed

        params = parcel_params{FT}(
            const_dt = const_dt,
            w = w,
            aerosol = aerosol,
            heterogeneous = "Frostenberg_random",
            condensation_growth = condensation_growth,
            deposition_growth = deposition_growth,
            size_distribution = DSD,
            drawing_interval = drawing_interval,
        )
        # solve ODE
        sol = run_parcel(IC, FT(0), t_max, params)
        push!(solutions, sol)
    end
    mean_sol = sum(solutions) / length(solutions)
    plot_results(sol.t, mean_sol, drawing_interval, "t_d=", " s", :dot, color)
end


# Frostenberg_stochastic with different timescales γ
γ_range = range(FT(1), stop = FT(5), length = 3)

for (γ, color) in zip(γ_range, colors)

    # creating an ensamble of solutions
    solutions = []
    for random_seed in random_seeds

        RAND.seed!(trunc(Int, random_seed)) #set the random seed

        params = parcel_params{FT}(
            const_dt = const_dt,
            w = w,
            aerosol = aerosol,
            heterogeneous = "Frostenberg_stochastic",
            condensation_growth = condensation_growth,
            deposition_growth = deposition_growth,
            size_distribution = DSD,
            γ = γ,
        )
        # solve ODE
        sol = run_parcel(IC, FT(0), t_max, params)
        push!(solutions, sol)
    end
    mean_sol = sum(solutions) / length(solutions)
    plot_results(sol.t, mean_sol, γ, "γ=", " s", :solid, color)
end

MK.save("frostenberg_immersion_freezing.svg", fig)
