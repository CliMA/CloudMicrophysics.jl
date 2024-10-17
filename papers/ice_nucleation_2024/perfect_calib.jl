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

    pseudo_data = perf_model_pseudo_data(FT, IN_mode, params, IC)
    Γ = pseudo_data[2]
    y_truth = pseudo_data[1]
    coeff_true = pseudo_data[3]
    end_sim = 25

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
    m_mean = []
    m_mean = ensemble_means(
        output[2],
        size(output[2])[1],
        size(output[2][1])[2],
    )[1]
    c_mean = []
    c_mean = ensemble_means(
        output[2],
        size(output[2])[1],
        size(output[2][1])[2],
    )[2]

    # Plotting calibrated parameters per iteration
    fig = MK.Figure(size = (800, 600))
    ax1 = MK.Axis(fig[1, 1], ylabel = "m coefficient [-]", xlabel = "iteration number", title = IN_mode)
    ax2 = MK.Axis(fig[2, 1], ylabel = "c coefficient [-]", xlabel = "iteration number")
    
    MK.lines!(ax1, iterations, m_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax1, iterations, zeros(length(m_mean)) .+ coeff_true[1], label = "default value")
    MK.lines!(ax2, iterations, c_mean,                                 label = "ensemble mean", color = :orange)
    MK.lines!(ax2, iterations, zeros(length(c_mean)) .+ coeff_true[2], label = "default value")
    
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
    
    soln = run_calibrated_model(FT, IN_mode, calibrated_parameters, params, IC)
    soln_dflt = run_calibrated_model(FT, IN_mode, coeff_true, params, IC)

    MK.lines!(ax3, soln.t, soln[9,:]./ (IC[7] + IC[8] + IC[9]), label = "calibrated")
    MK.lines!(ax3, soln_dflt.t, soln_dflt[9,:]./ (IC[7] + IC[8] + IC[9]), label = "default", color = :orange)
    plot_name = "perfect_calibration_ICNC_$mode_label.svg"
    MK.save(plot_name, fig2)
end
#! format: on
