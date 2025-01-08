import CairoMakie as MK
import CloudMicrophysics as CM

# Float32 does not work because inconsistent types using EKP
FT = Float64
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"))

IN_mode_list = ["ABDINM", "ABIFM", "ABHOM"]

#! format: off
for IN_mode in IN_mode_list
    params = perf_model_params(FT, IN_mode)
    IC = perf_model_IC(FT, IN_mode)
    end_sim = 25

    pseudo_data = perf_model_pseudo_data(FT, IN_mode, params, IC, end_sim)
    Γ = pseudo_data[2]
    y_truth = pseudo_data[1]
    coeff_true = pseudo_data[3]

    output = calibrate_J_parameters_EKI(
        FT,
        IN_mode,
        params,
        IC,
        y_truth,
        end_sim,
        Γ,
        perfect_model = true,
    )

    iterations = collect(1:size(output[2])[1])
    calibrated_parameters = output[1]
    ϕ_n_values = output[2]
    ABDINM_m_mean = []
    ABDINM_c_mean = []
    ABIFM_m_mean = []
    ABIFM_c_mean = []
    ABHOM_m_mean = []
    ABHOM_c_mean = []

    ABDINM_m_mean = ensemble_means(
        output[2],
        size(ϕ_n_values)[1],
        size(ϕ_n_values[1])[2],
    )[1]
    ABDINM_c_mean = ensemble_means(
        output[2],
        size(ϕ_n_values)[1],
        size(ϕ_n_values[1])[2],
    )[2]
    ABIFM_m_mean = ensemble_means(
        output[2],
        size(ϕ_n_values)[1],
        size(ϕ_n_values[1])[2],
    )[3]
    ABIFM_c_mean = ensemble_means(
        output[2],
        size(ϕ_n_values)[1],
        size(ϕ_n_values[1])[2],
    )[4]
    ABHOM_m_mean = ensemble_means(
        output[2],
        size(ϕ_n_values)[1],
        size(ϕ_n_values[1])[2],
    )[5]
    ABHOM_c_mean = ensemble_means(
        output[2],
        size(ϕ_n_values)[1],
        size(ϕ_n_values[1])[2],
    )[6]

    # Plotting calibrated parameters per iteration
    fig = MK.Figure(size = (1100, 900))
    ax1 = MK.Axis(fig[1, 1], ylabel = "ABDINM m coefficient [-]", title = IN_mode)
    ax2 = MK.Axis(fig[1, 2], ylabel = "ABDINM c coefficient [-]", xlabel = "iteration #")
    ax3 = MK.Axis(fig[2, 1], ylabel = "ABIFM m coefficient [-]", title = IN_mode)
    ax4 = MK.Axis(fig[2, 2], ylabel = "ABIFM c coefficient [-]", xlabel = "iteration #")
    ax5 = MK.Axis(fig[3, 1], ylabel = "ABHOM m coefficient [-]", title = IN_mode)
    ax6 = MK.Axis(fig[3, 2], ylabel = "ABHOM c coefficient [-]", xlabel = "iteration #")
    

    MK.lines!(ax1, iterations, ABDINM_m_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax1, iterations, zeros(length(ABDINM_m_mean)) .+ coeff_true[1], label = "default value")
    MK.lines!(ax2, iterations, ABDINM_c_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax2, iterations, zeros(length(ABDINM_c_mean)) .+ coeff_true[2], label = "default value")
    MK.lines!(ax3, iterations, ABIFM_m_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax3, iterations, zeros(length(ABIFM_m_mean)) .+ coeff_true[3], label = "default value")
    MK.lines!(ax4, iterations, ABIFM_c_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax4, iterations, zeros(length(ABIFM_c_mean)) .+ coeff_true[4], label = "default value")
    MK.lines!(ax5, iterations, ABHOM_m_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax5, iterations, zeros(length(ABHOM_m_mean)) .+ coeff_true[5], label = "default value")
    MK.lines!(ax6, iterations, ABHOM_c_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax6, iterations, zeros(length(ABHOM_c_mean)) .+ coeff_true[6], label = "default value")

    MK.axislegend(ax1, framevisible = true, labelsize = 12, position = :rc)
    MK.axislegend(ax2, framevisible = true, labelsize = 12)

    if IN_mode == "ABDINM"
        mode_label = "dep"
    elseif IN_mode == "ABIFM"
        mode_label = "imm"
    elseif IN_mode == "ABHOM"
        mode_label = "hom"
    end

    plot_name = "perfect_calibration_$mode_label.svg"
    MK.save(plot_name, fig)

    # Plotting calibrated parcel's ICNC
    fig2 = MK.Figure(size = (800, 600))
    ax3 = MK.Axis(fig2[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = IN_mode)
    
    soln = run_model(params, calibrated_parameters, FT, IC, end_sim)
    soln_dflt = run_model(params, coeff_true, FT, IC, end_sim)

    MK.lines!(ax3, soln.t, soln[9,:]./ (IC[7] + IC[8] + IC[9]), label = "calibrated")
    MK.lines!(ax3, soln_dflt.t, soln_dflt[9,:]./ (IC[7] + IC[8] + IC[9]), label = "default", color = :orange, linestyle = :dot)
    plot_name = "perfect_calibration_ICNC_$mode_label.svg"
    MK.save(plot_name, fig2)
end
#! format: on
