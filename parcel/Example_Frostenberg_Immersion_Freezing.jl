import OrdinaryDiffEq as ODE
import CairoMakie as MK

import CloudMicrophysics as CM
import CloudMicrophysics.HetIceNucleation as CMI
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams as CP

import Random as RAND
RAND.seed!(44)
N_ensemble = 32

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
FT = Float32
# get free parameters
tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
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
qₗ = Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2) # ρₐ ~ 1.2
qᵢ = FT(0)
qₜ = qᵥ + qₗ + qᵢ
ln_INPC₀ = FT(CMI_het.INP_concentration_mean(T₀))

# Moisture dependent initial conditions
Rᵥ = TDI.Rᵥ(tps)
Rₐ = TDI.Rₘ(tps, qₜ, qₗ, qᵢ)
eₛ = TDI.saturation_vapor_pressure_over_liquid(tps, T₀)
e = eᵥ(qᵥ, p₀, Rₐ, Rᵥ)
Sₗ = FT(e / eₛ)
IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC₀]

# Simulation parameters passed into ODE solver
w = FT(0.7)
t_max = FT(1200)
condensation_growth = "Condensation"
deposition_growth = "Deposition"
DSD = "Monodisperse"

# Model time step
const_dt_range = [FT(1), FT(0.25)]
# Every how many time steps will I sample for the "random" option
n_dt_range_12 = [[1, 60], [4, 240]]
# Assumed process timescale for the stochastic option
τ_range = [1, 100]

results_mean = []
labels_mean = []
time_mean = []
for const_dt in const_dt_range
    # Frostenberg_mean
    local params = parcel_params{FT}(
        const_dt = const_dt,
        w = w,
        heterogeneous = "Frostenberg_mean",
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        liq_size_distribution = DSD,
    )
    # solve ODE
    sol = run_parcel(IC, FT(0), t_max, params)
    push!(time_mean, sol.t)
    push!(results_mean, sol)
    push!(labels_mean, "dt = " * string(const_dt))
end

results_random = []
labels_random = []
time_random = []
for (const_dt, n_dt_range) in zip(const_dt_range, n_dt_range_12)
    # Frostenberg_random with different sampling frequencies
    for n_dt in n_dt_range
        # creating an ensamble of solutions
        sampling_interval = FT(n_dt * const_dt)
        solutions_rd = []
        for ensemble_member in range(1, length = N_ensemble)
            local params = parcel_params{FT}(
                const_dt = const_dt,
                w = w,
                heterogeneous = "Frostenberg_random",
                condensation_growth = condensation_growth,
                deposition_growth = deposition_growth,
                liq_size_distribution = DSD,
                sampling_interval = sampling_interval,
            )
            # solve ODE
            sol_rd = run_parcel(IC, FT(0), t_max, params)
            push!(solutions_rd, sol_rd)
        end
        mean_sol = sum(solutions_rd) / length(solutions_rd)
        push!(time_random, solutions_rd[1].t)
        push!(results_random, mean_sol)
        push!(
            labels_random,
            "dt = " * string(const_dt) * " Δt = " * string(sampling_interval),
        )
    end
end

results_stoch = []
labels_stoch = []
time_stoch = []
for const_dt in const_dt_range
    # Frostenberg_stochastic with different timescales γ
    for τ in τ_range
        # creating an ensamble of solutions
        γ = FT(1) / τ
        solutions_st = []
        for ensemble_member in range(1, length = N_ensemble)
            local params = parcel_params{FT}(
                const_dt = const_dt,
                w = w,
                heterogeneous = "Frostenberg_stochastic",
                condensation_growth = condensation_growth,
                deposition_growth = deposition_growth,
                liq_size_distribution = DSD,
                γ = γ,
            )
            # solve ODE
            sol_st = run_parcel(IC, FT(0), t_max, params)
            push!(solutions_st, sol_st)
        end
        mean_sol = sum(solutions_st) / length(solutions_st)
        push!(time_stoch, solutions_st[1].t)
        push!(results_stoch, mean_sol)
        push!(labels_stoch, "dt = " * string(const_dt) * " τ = " * string(τ))
    end
end

# Plotting
fig = MK.Figure(size = (1200, 900))
plot_theme = MK.Theme(Axis = (; xgridvisible = false, ygridvisible = false))
MK.set_theme!(plot_theme)
ax1 = MK.Axis(fig[1, 1]; ylabel = "Ice Supersaturation [-]")
ax2 = MK.Axis(fig[1, 2]; ylabel = "T [K]")
ax3 = MK.Axis(fig[2, 1]; ylabel = "q_ice [g/kg]")
ax4 = MK.Axis(fig[2, 2]; ylabel = "q_liq [g/kg]")
ax5 = MK.Axis(fig[3, 1]; xlabel = "Time [min]", ylabel = "N_ice [1/cm3]")
ax6 = MK.Axis(fig[3, 2]; xlabel = "Time [min]", ylabel = "N_liq [1/cm3]")
MK.xlims!.(fig.content, 0, t_max / 60)
ax7 = MK.Axis(
    fig[4, 1:2];
    xlabel = "T [K]",
    ylabel = "INPC [1/m3]",
    yscale = log10,
)
# top axis with time
ax7_top = MK.Axis(fig[4, 1:2]; xaxisposition = :top, xlabel = "Time [min]")
MK.hidespines!(ax7_top)
MK.hideydecorations!(ax7_top)
MK.xlims!(ax7_top, 0, t_max / 60)

colors = [:red, :green, :orange, :limegreen]
colors_mean = [:black, :gray]

function plot_results(sol, t, ll, cl, ls)
    MK.lines!(
        ax1,
        t / 60,
        S_i.(tps, sol[3, :], sol[1, :]) .- 1,
        linestyle = ls,
        color = cl,
    )
    MK.lines!(ax2, t / 60, sol[3, :], linestyle = ls, color = cl)
    MK.lines!(ax3, t / 60, sol[6, :] * 1e3, linestyle = ls, color = cl)
    MK.lines!(ax4, t / 60, sol[5, :] * 1e3, linestyle = ls, color = cl)
    MK.lines!(ax5, t / 60, sol[9, :] * 1e-6, linestyle = ls, color = cl)
    MK.lines!(
        ax6,
        t / 60,
        sol[8, :] * 1e-6,
        linestyle = ls,
        color = cl,
        label = ll,
    )
end
for (sol, t, ll, cl) in zip(results_stoch, time_stoch, labels_stoch, colors)
    ls = :solid
    plot_results(sol, t, ll, cl, ls)
    MK.lines!(ax7, sol[3, :], exp.(sol[10, :]), linestyle = ls, color = cl)
end
for (sol, t, ll, cl) in zip(results_random, time_random, labels_random, colors)
    ls = :dot
    plot_results(sol, t, ll, cl, ls)
end
for (sol, t, ll, cl) in zip(results_mean, time_mean, labels_mean, colors_mean)
    ls = :solid
    plot_results(sol, t, ll, cl, ls)
    MK.lines!(
        ax7,
        sol[3, :],
        exp.(CMI_het.INP_concentration_mean.(sol[3, :])),
        linestyle = ls,
        color = cl,
    )
end
MK.xlims!(ax7, extrema(results_mean[2][3, :]) |> reverse)
fig[1:2, 3] = MK.Legend(fig, ax6, framevisible = false)
MK.save("frostenberg_immersion_freezing.svg", fig)
nothing
