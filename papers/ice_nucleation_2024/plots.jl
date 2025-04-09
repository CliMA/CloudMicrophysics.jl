import CairoMakie as MK

include(joinpath(pkgdir(CM), "parcel", "ParcelCommon.jl"))

"""
    plot_AIDA_ICNC_data(
        exp_name,
        AIDA_t_profile, AIDA_ICNC_profile,
        t_profile, ICNC_moving_avg, frozen_frac_moving_mean,
    )

Plots raw AIDA ICNC data. Plot is saved and returned.
"""
function plot_AIDA_ICNC_data(
    exp_name,
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
    MK.lines!(data_ax1, t_profile, ICNC_moving_avg, label = "AIDA ICNC moving mean", linewidth = 2.5, color = :blue)
    MK.lines!(data_ax2, t_profile, frozen_frac_moving_mean, linewidth = 2.5, color = :blue)
    MK.axislegend(data_ax1, framevisible = true, labelsize = 12, position = :rc)
    MK.save("$exp_name" * "_ICNC.svg", AIDA_data_fig)

    return AIDA_data_fig

end

"""
    plot_calibrated_coeffs(batch_name, EKI_n_iterations, calibrated_ensemble_means)

Plots evolution of calibrated coefficients over the EKI iterations. Plot is saved and returned.
"""
function plot_calibrated_coeffs(batch_name, EKI_n_iterations, calibrated_ensemble_means)

    calibrated_coeffs_fig = MK.Figure(size = (1100, 900), fontsize = 24)
    m_coeff_ax = MK.Axis(calibrated_coeffs_fig[1, 1], ylabel = "m coefficient [-]", title = "$batch_name")
    c_coeff_ax =
        MK.Axis(calibrated_coeffs_fig[1, 2], ylabel = "c coefficient [-]", xlabel = "iteration #", title = "EKI")

    MK.lines!(
        m_coeff_ax,
        collect(1:EKI_n_iterations),
        calibrated_ensemble_means[1],
        label = "ensemble mean",
        color = :orange,
        linewidth = 2.5,
    )
    MK.lines!(
        c_coeff_ax,
        collect(1:EKI_n_iterations),
        calibrated_ensemble_means[2],
        label = "ensemble mean",
        color = :orange,
        linewidth = 2.5,
    )

    MK.save("$batch_name" * "_calibrated_coeffs_fig.svg", calibrated_coeffs_fig)

    return calibrated_coeffs_fig

end

"""
    plot_calibrated_parcels(
        exp_name, Nₜ,
        EKI_parcel, UKI_parcel,
        t_profile, T_profile, P_profile, ICNC_profile,
    )

Plots parcel simulations ran with calibrated parameters. Plot is saved and returned.
"""
function plot_calibrated_parcels(
    exp_name, Nₜ,
    EKI_parcel, UKI_parcel,
    t_profile, T_profile, P_profile, ICNC_profile,
)

    calibrated_parcel_fig = MK.Figure(size = (1500, 800), fontsize = 20)
    ax_parcel_1 =
        MK.Axis(calibrated_parcel_fig[1, 1], ylabel = "saturation [-]", xlabel = "time [s]", title = "$exp_name")
    ax_parcel_2 = MK.Axis(calibrated_parcel_fig[2, 1], ylabel = "liq mixing ratio [g/kg]", xlabel = "time [s]")
    ax_parcel_3 = MK.Axis(calibrated_parcel_fig[1, 2], ylabel = "temperature [K]", xlabel = "time [s]")
    ax_parcel_4 = MK.Axis(calibrated_parcel_fig[2, 2], ylabel = "qᵢ [g/kg]", xlabel = "time [s]")
    ax_parcel_5 = MK.Axis(calibrated_parcel_fig[3, 1], ylabel = "Nₗ [m^-3]", xlabel = "time [s]")
    ax_parcel_6 = MK.Axis(calibrated_parcel_fig[3, 2], ylabel = "Nᵢ [m^-3]", xlabel = "time [s]")
    ax_parcel_7 = MK.Axis(calibrated_parcel_fig[1, 3], ylabel = "pressure [Pa]", xlabel = "time [s]")
    ax_parcel_8 = MK.Axis(calibrated_parcel_fig[2, 3], ylabel = "Nₐ [m^-3]", xlabel = "time [s]")

    MK.lines!(ax_parcel_1, EKI_parcel.t, EKI_parcel[1, :], label = "EKI Calib Liq", color = :orange) # label = "liquid"
    MK.lines!(ax_parcel_1, UKI_parcel.t, UKI_parcel[1, :], label = "UKI Calib Liq", color = :fuchsia) # label = "liquid"
    MK.lines!(
        ax_parcel_1,
        EKI_parcel.t,
        S_i.(tps, EKI_parcel[3, :], EKI_parcel[1, :]),
        label = "EKI Calib Ice",
        color = :orange,
        linestyle = :dash,
    )
    MK.lines!(
        ax_parcel_1,
        UKI_parcel.t,
        S_i.(tps, UKI_parcel[3, :], UKI_parcel[1, :]),
        label = "UKI Calib Ice",
        color = :fuchsia,
        linestyle = :dash,
    )
    # MK.lines!(ax_parcel_1, t_profile, S_l_profile, label = "chamber", color = :blue)

    MK.lines!(ax_parcel_2, EKI_parcel.t, EKI_parcel[5, :], color = :orange)
    MK.lines!(ax_parcel_2, UKI_parcel.t, UKI_parcel[5, :], color = :fuchsia)

    MK.lines!(ax_parcel_3, EKI_parcel.t, EKI_parcel[3, :], color = :orange)
    MK.lines!(ax_parcel_3, UKI_parcel.t, UKI_parcel[3, :], color = :fuchsia)
    MK.lines!(ax_parcel_3, t_profile, T_profile, color = :blue, linestyle = :dash)

    MK.lines!(ax_parcel_4, EKI_parcel.t, EKI_parcel[6, :], color = :orange)
    MK.lines!(ax_parcel_4, UKI_parcel.t, UKI_parcel[6, :], color = :fuchsia)

    MK.lines!(ax_parcel_5, EKI_parcel.t, EKI_parcel[8, :], color = :orange)
    MK.lines!(ax_parcel_5, UKI_parcel.t, UKI_parcel[8, :], color = :fuchsia)

    MK.lines!(ax_parcel_6, EKI_parcel.t, EKI_parcel[9, :], color = :orange, label = "EKI")
    MK.lines!(ax_parcel_6, UKI_parcel.t, UKI_parcel[9, :], color = :fuchsia, label = "UKI")
    MK.lines!(ax_parcel_6, t_profile, ICNC_profile, color = :blue, label = "AIDA")

    error = (ICNC_profile ./ Nₜ) .* 0.1
    MK.errorbars!(ax_parcel_6, t_profile, ICNC_profile ./ Nₜ, error, color = (:blue, 0.3))

    MK.lines!(ax_parcel_7, EKI_parcel.t, EKI_parcel[2, :], color = :orange)
    MK.lines!(ax_parcel_7, UKI_parcel.t, UKI_parcel[2, :], color = :fuchsia)
    MK.lines!(ax_parcel_7, t_profile, P_profile, color = :blue, linestyle = :dash)

    MK.lines!(ax_parcel_8, EKI_parcel.t, EKI_parcel[7, :], color = :orange)
    MK.lines!(ax_parcel_8, UKI_parcel.t, UKI_parcel[7, :], color = :fuchsia)

    MK.axislegend(ax_parcel_6, framevisible = false, labelsize = 16, position = :rb)
    MK.save("$exp_name" * "_calibrated_parcel_fig.svg", calibrated_parcel_fig)

    return calibrated_parcel_fig

end

"""
    plot_compare_ICNC(
        exp_name,
        Nₜ, EKI_parcel, UKI_parcel,
        t_profile, frozen_frac_moving_mean, frozen_frac,
    )

Plots frozen fraction evolution for EKI and UKI calibrated parcel simulations and AIDA data.
Plot is saved and returned.
"""
function plot_compare_ICNC(
    exp_name,
    Nₜ, EKI_parcel, UKI_parcel,
    t_profile, frozen_frac_moving_mean, frozen_frac,
)

    ICNC_comparison_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_compare =
        MK.Axis(ICNC_comparison_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")
    MK.lines!(
        ax_compare,
        EKI_parcel.t,
        EKI_parcel[9, :] ./ Nₜ,
        label = "CM.jl Parcel (EKI Calibrated)",
        linewidth = 2.5,
        color = :orange,
    )
    MK.lines!(
        ax_compare,
        UKI_parcel.t,
        UKI_parcel[9, :] ./ Nₜ,
        label = "CM.jl Parcel (UKI Calibrated)",
        linewidth = 2.5,
        color = :fuchsia,
        linestyle = :dash,
    )
    error = frozen_frac_moving_mean .* 0.1
    MK.errorbars!(ax_compare, t_profile, frozen_frac_moving_mean, error, color = (:blue, 0.3))
    MK.lines!(
        ax_compare,
        t_profile,
        frozen_frac,
        label = "Raw AIDA",
        linewidth = 2,
        color = :blue,
        linestyle = :dash,
    )
    MK.lines!(
        ax_compare,
        t_profile,
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color = :blue,
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

    return plot_UKI_spread

end

"""
    plot_ICNC_overview(
        EKI_calibrated_parcel, UKI_parcel,
    )

Plots frozen fraction evolution for EKI and UKI calibrated parcel simulations
and AIDA data for all batch experiments. Plot is saved and returned.
"""
function plot_ICNC_overview(overview_data)

    EKI_calibrated_parcel = overview_data.EKI_calibrated_parcel
    UKI_calibrated_parcel = overview_data.UKI_calibrated_parcel
    Nₜ_list = overview_data.Nₜ_list
    t_profile_list = overview_data.t_profile_list
    frozen_frac_moving_mean_list = overview_data.frozen_frac_moving_mean_list
    frozen_frac_list = overview_data.frozen_frac_list

    overview_fig = MK.Figure(size = (1000, 800), fontsize = 24)

    exp_name_list = [
        "IN05_17", "IN05_18", "TROPIC04",
        "IN07_01", "IN07_19", "EXP45",
        "ACI04_22", "EXP19",
    ]
    position = [
        overview_fig[1, 1], overview_fig[2, 1], overview_fig[3, 1],
        overview_fig[1, 2], overview_fig[2, 2], overview_fig[3, 2],
        overview_fig[1, 3], overview_fig[2, 3],
    ]

    for (i, exp_name) in enumerate(exp_name_list)

        EKI_parcel = EKI_calibrated_parcel[i]
        UKI_parcel = UKI_calibrated_parcel[i]
        Nₜ = Nₜ_list[i]
        frozen_frac_moving_mean = frozen_frac_moving_mean_list[i]
        t_profile = t_profile_list[i]
        frozen_frac = frozen_frac_list[i]

        ax = MK.Axis(position[i], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")

        MK.lines!(
            ax,
            EKI_parcel.t,
            EKI_parcel[9, :] ./ Nₜ,
            label = "CM.jl Parcel (EKI Calibrated)",
            linewidth = 2.5,
            color = :orange,
        )
        MK.lines!(
            ax,
            UKI_parcel.t,
            UKI_parcel[9, :] ./ Nₜ,
            label = "CM.jl Parcel (UKI Calibrated)",
            linewidth = 2.5,
            color = :fuchsia,
            linestyle = :dash,
        )
        error = frozen_frac_moving_mean .* 0.1
        MK.errorbars!(ax, t_profile, frozen_frac_moving_mean, error, color = (:blue, 0.3))
        MK.lines!(
            ax,
            t_profile,
            frozen_frac,
            label = "Raw AIDA",
            linewidth = 2,
            color = :blue,
            linestyle = :dash,
        )
        MK.lines!(
            ax,
            t_profile,
            frozen_frac_moving_mean,
            label = "AIDA Moving Avg",
            linewidth = 2.5,
            color = :blue,
        )
    end

    legend_axis = MK.Axis(overview_fig[3, 3], title = "Legend")
    MK.lines!(legend_axis, 1, 1, label = "CM.jl Parcel (EKI Calibrated)", color = :orange)
    MK.lines!(legend_axis, 1, 1, label = "CM.jl Parcel (UKI Calibrated)", color = :fuchsia, linestyle = :dash)
    MK.lines!(legend_axis, 1, 1, label = "Raw AIDA", color = :blue, linestyle = :dash)
    MK.lines!(legend_axis, 1, 1, label = "AIDA Moving Avg", color = :blue)
    MK.axislegend(legend_axis, framevisible = false, labelsize = 18, position = :rb)

    MK.hidespines!(legend_axis)
    MK.hidedecorations!(legend_axis)

    MK.save("Overview_ICNC_fig.svg", overview_fig)

    return overview_fig
end
