import CairoMakie as MK
import CloudMicrophysics as CM

# Float32 does not work because inconsistent types using EKP
FT = Float64
include(
    joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "perfect_model_J.jl"),
)

ABDINM_output = calibrate_J_parameters(FT, "ABDINM")
ABIFM_output = calibrate_J_parameters(FT, "ABIFM")
ABHOM_output = calibrate_J_parameters(FT, "ABHOM")
iterations = collect(1:size(ABHOM_output[5])[1])

ABDINM_calibrated_parameters = [ABDINM_output[1], ABDINM_output[2]]
ABDINM_coeff_true = [ABDINM_output[3], ABDINM_output[4]]
ABIFM_calibrated_parameters = [ABIFM_output[1], ABIFM_output[2]]
ABIFM_coeff_true = [ABIFM_output[3], ABIFM_output[4]]
ABHOM_calibrated_parameters = [ABHOM_output[1], ABHOM_output[2]]
ABHOM_coeff_true = [ABHOM_output[3], ABHOM_output[4]]

# ensemble_means(ϕ_n_values, N_iterations, ensembles)
ABDINM_ensemble_means = ensemble_means(
    ABDINM_output[5],
    size(ABDINM_output[5])[1],
    size(ABDINM_output[5][1])[2],
)
ABIFM_ensemble_means = ensemble_means(
    ABIFM_output[5],
    size(ABIFM_output[5])[1],
    size(ABIFM_output[5][1])[2],
)
ABHOM_ensemble_means = ensemble_means(
    ABHOM_output[5],
    size(ABHOM_output[5])[1],
    size(ABHOM_output[5][1])[2],
)

# Plotting calibrated parameters per iteration
#! format: off
ABDINM_fig = MK.Figure(size = (800, 600))
ax1 = MK.Axis(ABDINM_fig[1, 1], ylabel = "m coefficient [-]", xlabel = "iteration number", title = "ABDINM")
ax2 = MK.Axis(ABDINM_fig[2, 1], ylabel = "c coefficient [-]", xlabel = "iteration number")

MK.lines!(ax1, iterations, ABDINM_ensemble_means[1],                                    label = "ensemble mean", color = :orange)
MK.lines!(ax1, iterations, zeros(length(ABDINM_ensemble_means[1])) .+ ABDINM_output[3], label = "default value")
MK.lines!(ax2, iterations, ABDINM_ensemble_means[2],                                    label = "ensemble mean", color = :orange)
MK.lines!(ax2, iterations, zeros(length(ABDINM_ensemble_means[1])) .+ ABDINM_output[4], label = "default value")

ABIFM_fig = MK.Figure(size = (800, 600))
ax3 = MK.Axis(ABIFM_fig[1, 1], ylabel = "m coefficient [-]", xlabel = "iteration number", title = "ABIFM")
ax4 = MK.Axis(ABIFM_fig[2, 1], ylabel = "c coefficient [-]", xlabel = "iteration number")

MK.lines!(ax3, iterations, ABIFM_ensemble_means[1],                                   label = "ensemble mean", color = :orange)
MK.lines!(ax3, iterations, zeros(length(ABIFM_ensemble_means[1])) .+ ABIFM_output[3], label = "default value")
MK.lines!(ax4, iterations, ABIFM_ensemble_means[2],                                   label = "ensemble mean", color = :orange)
MK.lines!(ax4, iterations, zeros(length(ABIFM_ensemble_means[1])) .+ ABIFM_output[4], label = "default value")

ABHOM_fig = MK.Figure(size = (800, 600))
ax5 = MK.Axis(ABHOM_fig[1, 1], ylabel = "m coefficient [-]", xlabel = "iteration number", title = "ABHOM")
ax6 = MK.Axis(ABHOM_fig[2, 1], ylabel = "c coefficient [-]", xlabel = "iteration number")

MK.lines!(ax5, iterations, ABHOM_ensemble_means[1],                                   label = "ensemble mean", color = :orange)
MK.lines!(ax5, iterations, zeros(length(ABHOM_ensemble_means[1])) .+ ABHOM_output[3], label = "default value")
MK.lines!(ax6, iterations, ABHOM_ensemble_means[2],                                   label = "ensemble mean", color = :orange)
MK.lines!(ax6, iterations, zeros(length(ABHOM_ensemble_means[1])) .+ ABHOM_output[4], label = "default value")
#! format: on

MK.axislegend(ax1, framevisible = true, labelsize = 12, position = :rc)
MK.axislegend(ax2, framevisible = true, labelsize = 12)
MK.axislegend(ax3, framevisible = true, labelsize = 12, position = :rc)
MK.axislegend(ax4, framevisible = true, labelsize = 12)
MK.axislegend(ax5, framevisible = true, labelsize = 12, position = :rc)
MK.axislegend(ax6, framevisible = true, labelsize = 12)

MK.save("perfect_calibration_dep.svg", ABDINM_fig)
MK.save("perfect_calibration_imm.svg", ABIFM_fig)
MK.save("perfect_calibration_hom.svg", ABHOM_fig)
