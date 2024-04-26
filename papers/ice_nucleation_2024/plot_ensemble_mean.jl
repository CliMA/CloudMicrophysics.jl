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

    output = calibrate_J_parameters(
        FT,
        IN_mode,
        params,
        IC,
        y_truth,
        Γ,
        perfect_model = true,
    )

    iterations = collect(1:size(output[3])[1])
    calibrated_parameters = [output[1], output[2]]
    m_mean = []
    m_mean = ensemble_means(
        output[3],
        size(output[3])[1],
        size(output[3][1])[2],
    )[1]
    c_mean = []
    c_mean = ensemble_means(
        output[3],
        size(output[3])[1],
        size(output[3][1])[2],
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

end
#! format: on
