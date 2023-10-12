import CloudMicrophysics as CM
import CairoMakie as MK
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het

FT = Float64

# AK 2016 data and the frozen fraction fit to data
include(joinpath(pkgdir(CM), "box", "Alpert_Knopf_2016_data.jl"))

# N_frz, N_ufrz - numbers of frozen and unforzen droplets
function frozen_frac(N_frz, N_ufrz)
    return N_frz / (N_frz + N_ufrz)
end

# definition of the ODE problem for box model
include(joinpath(pkgdir(CM), "box", "box.jl"))

# initial conditions following KA16 figure 4 (Cr1 case)
A_droplet = FT(1e-5) * 1e-4 # INP surface area, m^2
σg = 10
N₀ = 1000
N_ice = 0
T_initial = FT(256)         # initial temperature, K
cooling_rate = FT(0.5 / 60) # prescribed cooling rate K s^-1
aerosol = CMP.Illite(FT)    # aerosol free parameters
tps = CMP.ThermodynamicsParameters(FT) # thermodynamics free parameters
t_0 = 0
t_end = 3310
dt = 10
A_sum = N₀ * A_droplet # Total surface area m^2 when all droplets are equal

# A vector with surface areas sampled from a lognormal distribution and sorted
A_distr = RD.rand(DS.LogNormal(log(A_droplet), log(σg)), N₀)
Aj_unsorted = zeros(N₀)
for it in range(start = 1, stop = N₀, step = 1)
    Aj_unsorted[it] = A_distr[it]
end
Aj_sorted = zeros(N₀)
Aj_sorted = sort(Aj_unsorted, rev = true)

# initial condition for the ODE problem
IC = [T_initial, A_sum, FT(N₀), FT(N_ice)]
# additional model parameters
p_def = (; tps, A_droplet, aerosol, cooling_rate, N₀)

# run the box model without and with distribution sampling
sol_cst = run_box(IC, t_0, t_end, (p_def..., const_dt = dt, flag = false))
sol_var = run_box(
    IC,
    t_0,
    t_end,
    (p_def..., const_dt = dt, flag = true, Aj_sorted = Aj_sorted),
)

# Compute the immersion freezing rate for the no-sampling case
Δa = @. FT(1) - CMO.a_w_ice(tps, sol_cst[1, :])
J_immer = @. CMI_het.ABIFM_J(aerosol, Δa) # m^-2 s^-1

# Compute the frozen fraction for all four runs
ff_cst = @. frozen_frac(sol_cst[4, :], sol_cst[3, :])
ff_var = @. frozen_frac(sol_var[4, :], sol_var[3, :])

#! format: off
# Plot results
fig = MK.Figure(resolution = (1200, 400))
ax1 = MK.Axis(fig[1, 1], ylabel = "J_immer [cm^-2 s^-1]", xlabel = "Temperature [K]", yscale = log10)
ax2 = MK.Axis(fig[1, 2], ylabel = "frozen fraction",      xlabel = "Temperature [K]")
ax3 = MK.Axis(fig[1, 3], ylabel = "total A [cm2]",        xlabel = "Temperature [K]", yscale = log10)

MK.lines!(ax1, sol_cst[1, :], J_immer .* 1e-4, color = :rosybrown2, label = "CliMA ABIFM", linewidth = 3)
MK.lines!(ax1, AK16_T_J,      AK16_J,          color = :red,        label = "ABIFM",       linewidth = 3)
MK.ylims!(ax1, 1e-4, 1e10)

MK.lines!(ax2, AK16_T_frozen_fraction, AK16_frozen_fraction, color = :black,  linewidth = 3, label = "prescribed")
MK.lines!(ax2, sol_cst[1, :],          ff_cst,               color = :orange, linewidth = 3, label = "const A")
MK.lines!(ax2, sol_var[1, :],          ff_var,               color = :green3, linewidth = 3, label = "variable A")

idx = (sol_cst[3, :] .> 0)
MK.lines!(ax3, sol_cst[1, idx], sol_cst[3, idx] .* A_droplet .* 1e4, color = :orange, linewidth = 3, label = "const A")
idx = (sol_var[2, :] .> 0)
MK.lines!(ax3, sol_var[1, idx], sol_var[2, idx] .* 1e4,              color = :green3, linewidth = 3, label = "variable A")
MK.ylims!(ax3, 0.75e-6, 2e-1)

MK.axislegend(ax1, framevisible = false, labelsize = 14, orientation = :horizontal, nbanks = 2, position = :rt)
MK.axislegend(ax2, framevisible = false, labelsize = 14, orientation = :horizontal, nbanks = 3, position = :lb)
MK.axislegend(ax3, framevisible = false, labelsize = 14, orientation = :horizontal, nbanks = 2, position = :lt)

MK.save("Alpert_Knopf_2016_forward.svg", fig)
#! format: on
