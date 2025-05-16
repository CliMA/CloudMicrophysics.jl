import CairoMakie as MK

include(joinpath(pkgdir(CM), "parcel", "ParcelCommon.jl"))

#! format: off

"""
    plot_AIDA_ICNC_data(
        exp_name,
        AIDA_t_profile, AIDA_ICNC_profile,
        t_profile, ICNC_moving_avg, frozen_frac_moving_mean,
    )

Plots raw AIDA ICNC data. Plot is saved and returned.
"""
function plot_AIDA_ICNC_data(
    exp_name, start_time,
    AIDA_t_profile, AIDA_ICNC_profile,
    t_profile, ICNC_moving_avg, frozen_frac_moving_mean,
)

    AIDA_data_fig = MK.Figure(size = (1000, 600), fontsize = 24)
    data_ax1 = MK.Axis(AIDA_data_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "AIDA data $exp_name")
    data_ax2 = MK.Axis(
        AIDA_data_fig[1, 2],
        ylabel = "Frozen Frac Moving Mean [-]",
        xlabel = "time [s]",
        title = "AIDA data $exp_name",
    )
    MK.lines!(
        data_ax1,
        AIDA_t_profile,
        AIDA_ICNC_profile,
        label = "Raw AIDA",
        color = :blue,
        linestyle = :dash,
        linewidth = 2,
    )
    MK.lines!(data_ax1, t_profile .+ start_time, ICNC_moving_avg, label = "AIDA ICNC moving mean", linewidth = 2.5, color = :blue)
    MK.lines!(data_ax2, t_profile .+ start_time, frozen_frac_moving_mean, linewidth = 2.5, color = :blue)
    MK.axislegend(data_ax1, framevisible = true, labelsize = 12, position = :rc)
    MK.save("$exp_name" * "_ICNC.svg", AIDA_data_fig)

    return AIDA_data_fig

end

"""
    plot_calibrated_coeffs(
        batch_name,
        UKI_n_iterations, UKI_n_ensembles, UKI_all_params,
        UKI_mean_each_iter,
    )

Plots evolution of calibrated coefficients over the UKI iterations. Plot is saved and returned.
"""
function plot_calibrated_coeffs(
    batch_name,
    UKI_n_iterations, UKI_n_ensembles, UKI_all_params,
    UKI_mean_each_iter,
)

    # Every ensemble member at every iteration
    UKI_m = []
    UKI_c = []
    UKI_x_axis_ensembles = []
    for iter in 1:UKI_n_iterations
        for ensemble_n in 1:UKI_n_ensembles
            append!(UKI_x_axis_ensembles, iter)
            append!(UKI_m, UKI_all_params[iter][1, ensemble_n])
            append!(UKI_c, UKI_all_params[iter][2, ensemble_n])
        end
    end

    # Mean of each iteration
    UKI_m_mean = []
    UKI_c_mean = []
    UKI_x_axis_mean = []
    for iter in 1:UKI_n_iterations
        append!(UKI_x_axis_mean, iter)
        append!(UKI_m_mean, UKI_mean_each_iter[iter][1])
        append!(UKI_c_mean, UKI_mean_each_iter[iter][2])
    end

    calibrated_coeffs_fig = MK.Figure(size = (900, 450), fontsize = 24)
    m_coeff_ax = MK.Axis(
        calibrated_coeffs_fig[1, 1],
        ylabel = "m coefficient [-]",
        xlabel = "Iteration #",
        yticklabelcolor = :blue,
        )
    c_coeff_ax = MK.Axis(
        calibrated_coeffs_fig[1, 1],
        ylabel = "c coefficient [-]",
        xlabel = "Iteration #",
        yticklabelcolor = :red,
        yaxisposition = :right,
        )

    MK.Label(
        calibrated_coeffs_fig[0, :], 
        "$batch_name: Calibrating Coefficients",
        font = :bold,
        fontsize = 30,
        halign = :center,   # Horizontal alignment
        valign = :top       # Vertical alignment (relative to its grid cell)
    )
    MK.lines!(
        m_coeff_ax,
        UKI_x_axis_mean,
        UKI_m_mean,
        label = "UKI Mean",
        color = (:blue, 0.6),
        linewidth = 3,
    )
    # MK.scatter!(
    #     m_coeff_ax,
    #     UKI_x_axis_ensembles,
    #     UKI_m,
    #     label = "UKI Ensemble",
    #     color = (:deepskyblue, 0.5),
    # )
    MK.lines!(
        c_coeff_ax,
        UKI_x_axis_mean,
        UKI_c_mean,
        label = "UKI Mean",
        color = (:red, 0.6),
        linewidth = 3,
    )
    # MK.scatter!(
    #     c_coeff_ax,
    #     UKI_x_axis_ensembles,
    #     UKI_c,
    #     label = "UKI Ensemble",
    #     color = (:deepskyblue, 0.5),
    # )

    # MK.axislegend(m_coeff_ax, framevisible = false, labelsize = 18, position = :rt)
    MK.save("$batch_name" * "_calibrated_coeffs_fig.svg", calibrated_coeffs_fig)

    return calibrated_coeffs_fig

end

"""
    plot_calibrated_parcels(
        exp_name, Nₜ,
        UKI_parcel,
        t_profile, T_profile, P_profile, ICNC_profile,
    )

Plots parcel simulations ran with calibrated parameters. Plot is saved and returned.
"""
function plot_calibrated_parcels(
    exp_name, Nₜ,
    UKI_parcel,
    t_profile, T_profile, P_profile, ICNC_profile,
    Sₗ_profile,
)

    calibrated_parcel_fig = MK.Figure(size = (1500, 800), fontsize = 20)

    ax_parcel_1 =
        MK.Axis(calibrated_parcel_fig[1, 1], ylabel = "S [-]")
    ax_parcel_2 = MK.Axis(calibrated_parcel_fig[1, 3], ylabel = "qₗ [kg/kg]", )
    ax_parcel_3 = MK.Axis(calibrated_parcel_fig[2, 2], ylabel = "T [K]")
    ax_parcel_4 = MK.Axis(calibrated_parcel_fig[3, 1], ylabel = "qᵢ [×10⁻⁶ kg/kg]", xlabel = "time [s]")
    ax_parcel_5 = MK.Axis(calibrated_parcel_fig[2, 3], ylabel = "Nₗ [×10⁶ m⁻³]")
    ax_parcel_6 = MK.Axis(calibrated_parcel_fig[3, 2], ylabel = "Nᵢ [×10⁶ m⁻³]", xlabel = "time [s]")
    ax_parcel_7 = MK.Axis(calibrated_parcel_fig[2, 1], ylabel = "p [hPa]",)
    ax_parcel_8 = MK.Axis(calibrated_parcel_fig[1, 2], ylabel = "Nₐ [×10⁶ m⁻³]")

    MK.Label(
        calibrated_parcel_fig[0, 1:3], 
        "$exp_name: Calibrated Parcel",
        font = :bold,
        fontsize = 30,
        halign = :center,   # Horizontal alignment
        valign = :top       # Vertical alignment (relative to its grid cell)
    )

    MK.lines!(ax_parcel_1, UKI_parcel.t, UKI_parcel[1, :], label = "Liquid", color = :dodgerblue, linewidth = 2.5) # label = "liquid"
    MK.lines!(
        ax_parcel_1,
        UKI_parcel.t,
        S_i.(tps, UKI_parcel[3, :], UKI_parcel[1, :]),
        label = "Ice",
        color = :dodgerblue,
        linestyle = :dash,
        linewidth = 2.5,
    )
    MK.lines!(ax_parcel_1, t_profile, Sₗ_profile, label = "AIDA", color = :black, linewidth = 2.5)

    MK.lines!(ax_parcel_2, UKI_parcel.t, UKI_parcel[5, :], color = :dodgerblue, linewidth = 2.5)

    MK.lines!(ax_parcel_3, UKI_parcel.t, UKI_parcel[3, :], color = :dodgerblue, linewidth = 2.5)
    # MK.lines!(ax_parcel_3, t_profile, T_profile, color = :black, linestyle = :dash, linewidth = 2.5)

    MK.lines!(ax_parcel_4, UKI_parcel.t, UKI_parcel[6, :] .* 10^6, color = :dodgerblue, linewidth = 2.5)

    MK.lines!(ax_parcel_5, UKI_parcel.t, UKI_parcel[8, :] ./ 10^6, color = :dodgerblue, linewidth = 2.5)

    MK.lines!(ax_parcel_6, UKI_parcel.t, UKI_parcel[9, :] ./ 10^6, color = :dodgerblue, label = "UKI", linewidth = 2.5)
    MK.lines!(ax_parcel_6, t_profile, ICNC_profile ./ 10^6 , color = :black, label = "AIDA", linewidth = 2.5)

    error = (ICNC_profile ./ 10^6) .* 0.1
    MK.errorbars!(ax_parcel_6, t_profile, ICNC_profile ./ 10^6, error, color = (:gray, 0.3))

    MK.lines!(ax_parcel_7, UKI_parcel.t, UKI_parcel[2, :] ./ 10^2, color = :dodgerblue, linewidth = 2.5)
    # MK.lines!(ax_parcel_7, t_profile, P_profile, color = :black, linestyle = :dash, linewidth = 2.5)

    MK.lines!(ax_parcel_8, UKI_parcel.t, UKI_parcel[7, :] ./ 10^6, color = :dodgerblue, linewidth = 2.5)

    MK.axislegend(ax_parcel_1, framevisible = false, labelsize = 16, position = :lt)
    MK.axislegend(ax_parcel_6, framevisible = false, labelsize = 16, position = :rb)
    MK.save("$exp_name" * "_calibrated_parcel_fig.svg", calibrated_parcel_fig)

    return calibrated_parcel_fig

end

"""
    plot_compare_ICNC(
        exp_name,
        Nₜ, UKI_parcel,
        t_profile, frozen_frac_moving_mean, frozen_frac,
    )

Plots frozen fraction evolution for UKI calibrated parcel simulation and AIDA data.
Plot is saved and returned.
"""
function plot_compare_ICNC(
    exp_name,
    Nₜ, UKI_parcel,
    t_profile, frozen_frac_moving_mean, frozen_frac,
)

    ICNC_comparison_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_compare =
        MK.Axis(ICNC_comparison_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")
    MK.lines!(
        ax_compare,
        UKI_parcel.t,
        UKI_parcel[9, :] ./ Nₜ,
        label = "UKI Calibrated",
        linewidth = 2.5,
        color = :blue,
    )
    error = frozen_frac_moving_mean .* 0.1
    MK.errorbars!(ax_compare, t_profile, frozen_frac_moving_mean, error, color = (:gray, 0.3))
    MK.lines!(
        ax_compare,
        t_profile,
        frozen_frac,
        label = "Raw AIDA",
        linewidth = 2,
        color = :black,
        linestyle = :dash,
    )
    MK.lines!(
        ax_compare,
        t_profile,
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color = :black,
    )

    MK.axislegend(ax_compare, framevisible = false, labelsize = 20, position = :rb)
    MK.save("$exp_name" * "_ICNC_comparison_fig.svg", ICNC_comparison_fig)

    return ICNC_comparison_fig
end

"""
    plot_UKI_spread(
        exp_name,
        t_profile, frozen_frac_moving_mean, Nₜ,
        UKI_parcel_1, UKI_parcel_2, UKI_parcel_3, UKI_parcel_4, UKI_parcel_5, UKI_parcel,
    )

Plots UKI spread of the 5 ensembles.
"""
function plot_UKI_spread(
    exp_name,
    t_profile, frozen_frac_moving_mean, Nₜ,
    UKI_parcel_1, UKI_parcel_2, UKI_parcel_3, UKI_parcel_4, UKI_parcel_5, UKI_parcel,
)

    UKI_spread_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_spread = MK.Axis(UKI_spread_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")
    MK.lines!(
        ax_spread,
        t_profile,
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color = :blue,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_1.t,
        UKI_parcel_1[9, :] ./ Nₜ,
        linewidth = 2.5,
        color = :grey80,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_2.t,
        UKI_parcel_2[9, :] ./ Nₜ,
        linewidth = 2.5,
        color = :grey60,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_3.t,
        UKI_parcel_3[9, :] ./ Nₜ,
        linewidth = 2.5,
        color = :grey40,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_4.t,
        UKI_parcel_4[9, :] ./ Nₜ,
        linewidth = 2.5,
        color = :grey25,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_5.t,
        UKI_parcel_5[9, :] ./ Nₜ,
        linewidth = 2.5,
        color = :grey12,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel.t,
        UKI_parcel[9, :] ./ Nₜ,
        label = "CM.jl Parcel (UKI Calibrated)",
        linewidth = 2.5,
        color = :fuchsia,
        linestyle = :dash,
    )
    error = frozen_frac_moving_mean .* 0.1
    MK.errorbars!(ax_spread, t_profile, frozen_frac_moving_mean, error, color = (:blue, 0.3))
    MK.save("$exp_name" * "_UKI_spread_fig.svg", UKI_spread_fig)

    return UKI_spread_fig

end

"""
    plot_ICNC_overview(overview_data)

Plots frozen fraction evolution for UKI calibrated parcel simulation
and AIDA data for all batch experiments. Plot is saved and returned.
"""
function plot_ICNC_overview(overview_data)

    UKI_calibrated_parcel = overview_data.UKI_calibrated_parcel
    P3_parcel_list = overview_data.P3_parcel
    Nₜ_list = overview_data.Nₜ_list
    t_profile_list = overview_data.t_profile_list
    frozen_frac_moving_mean_list = overview_data.frozen_frac_moving_mean_list
    frozen_frac_list = overview_data.frozen_frac_list
    uncertainty = 0.1

    overview_fig = MK.Figure(size = (1000, 800), fontsize = 24)
    MK.Label(
        overview_fig[0, 1], 
        "         DEP         ",
        font = :bold,
        fontsize = 30,
        halign = :center,   # Horizontal alignment
        valign = :top       # Vertical alignment (relative to its grid cell)
    )
    MK.Label(
        overview_fig[0, 2], 
        "         IMM         ",
        font = :bold,
        fontsize = 30,
        halign = :center,   # Horizontal alignment
        valign = :top       # Vertical alignment (relative to its grid cell)
    )
    MK.Label(
        overview_fig[0, 3], 
        "         HOM         ",
        font = :bold,
        fontsize = 30,
        halign = :center,   # Horizontal alignment
        valign = :top       # Vertical alignment (relative to its grid cell)
    )

    exp_name_list = [
        "IN05_17", "IN05_18",
        "IN07_01", "IN07_19", "EXP45",
        "ACI04_22", "EXP19",
    ]
    position = [
        overview_fig[1, 3], overview_fig[2, 3],
        overview_fig[1, 1], overview_fig[2, 1], overview_fig[3, 1],
        overview_fig[1, 2], overview_fig[2, 2],
    ]

    for (i, exp_name) in enumerate(exp_name_list)

        UKI_parcel = UKI_calibrated_parcel[i]
        P3_parcel = P3_parcel_list[i]
        Nₜ = Nₜ_list[i]
        frozen_frac_moving_mean = frozen_frac_moving_mean_list[i]
        t_profile = t_profile_list[i]
        frozen_frac = frozen_frac_list[i]

        ax = MK.Axis(position[i], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")

        up_bound_uncertain = frozen_frac_moving_mean .* (1 + uncertainty)
        low_bound_uncertain = frozen_frac_moving_mean .* (1 - uncertainty)
        MK.lines!(
            ax,
            t_profile,
            up_bound_uncertain,
            linewidth = 1,
            color = (:gray, 0.5),
        )
        MK.lines!(
            ax,
            t_profile,
            low_bound_uncertain,
            linewidth = 1,
            color = (:gray, 0.5),
        )
        MK.band!(t_profile, low_bound_uncertain, up_bound_uncertain, color = (:gray, 0.5))

        MK.lines!(
            ax,
            UKI_parcel.t,
            UKI_parcel[9, :] ./ Nₜ,
            label = "UKI Calibrated",
            linewidth = 3,
            color = :dodgerblue,
        )
        # MK.lines!(
        #     ax,
        #     P3_parcel.t,
        #     P3_parcel[9, :] ./ Nₜ,
        #     label = "CM.jl Parcel (P3)",
        #     linewidth = 3,
        #     color = :red,
        # )
        MK.lines!(
            ax,
            t_profile,
            frozen_frac,
            label = "Raw AIDA",
            linewidth = 2,
            color = :black,
            linestyle = :dash,
        )
        MK.lines!(
            ax,
            t_profile,
            frozen_frac_moving_mean,
            label = "AIDA Moving Avg",
            linewidth = 3,
            color = :black,
        )
    end

    legend_axis = MK.Axis(overview_fig[3, 2:3])
    MK.lines!(legend_axis, 1, 1, label = "AIDA Raw Data", color = :black, linestyle = :dash, linewidth = 2)
    MK.lines!(legend_axis, 1, 1, label = "AIDA Moving Average", color = :black, linewidth = 3)
    MK.lines!(legend_axis, 1, 1, label = "UKI Calibrated Simulation", color = :dodgerblue, linewidth = 3)
    # MK.lines!(legend_axis, 1, 1, label = "CM.jl Parcel (P3)", color = :red)
    MK.axislegend(legend_axis, "Legend", framevisible = false, labelsize = 18, position = :cc)

    MK.hidespines!(legend_axis)
    MK.hidedecorations!(legend_axis)

    MK.save("Overview_ICNC_fig.svg", overview_fig)
    return overview_fig
end

"""
    plot_loss_func(batch_name, UKI_n_iterations, UKI_error)

Plots the difference in loss functions over each iteration of
    calibration for a batch. Plot is saved and returned.
"""
function plot_loss_func(batch_name, UKI_n_iterations, UKI_error)

    loss_fig = MK.Figure(size = (600, 530), fontsize = 24)
    ax_loss = MK.Axis(loss_fig[1, 1], ylabel = "L(θ) [-]", xlabel = "Iteration # [-]", title = "$batch_name: Minimizing Loss")
    MK.lines!(
        ax_loss,
        collect(1:(UKI_n_iterations - 1)),
        UKI_error,
        linewidth = 2.5,
        color = :red,
    )
    # MK.axislegend(ax_loss, framevisible = false, labelsize = 18, position = :rc)
    MK.save("$batch_name" * "_loss_fig.svg", loss_fig)

    return loss_fig

end
#! format: on
