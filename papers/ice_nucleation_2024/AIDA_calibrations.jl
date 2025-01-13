import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD
import CairoMakie as MK

# To grab data
import DelimitedFiles
using LazyArtifacts
using ClimaUtilities.ClimaArtifacts

FT = Float64
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"))

# Helper functions
function unpack_data(data_file_name)

    file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name)
    return DelimitedFiles.readdlm(file_path, skipstart = 125)
end
function grab_data(unpacked_data)

    AIDA_t_profile = unpacked_data[:, 1]
    AIDA_T_profile = unpacked_data[:, 3]
    AIDA_P_profile = unpacked_data[:, 2] * 1e2  # hPa to Pa
    AIDA_ICNC = unpacked_data[:, 6] .* 1e6      # Nᵢ [m^-3]
    AIDA_e = unpacked_data[:, 4]

    data = (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC, AIDA_e)
    return data
end
moving_average(data, n) =
    [sum(@view data[i:(i + n)]) / n for i in 1:(length(data) - n)]

# Defining data names, start/end times, etc.
data_file_names = [
    "in05_17_aida.edf",
    "in05_18_aida.edf",
    "in05_19_aida.edf",
    "in07_01_aida.edf",
    "in07_19_aida.edf",
]
plot_names = ["IN0517", "IN0518", "IN0519", "IN0701", "IN0719"]
end_sim = 25                                            # Loss func looks at last end_sim timesteps only
start_time_list = [150, 180, 80, 50, 35]                # freezing onset
end_time_list = [290, 290, 170, 400, 400]               # approximate time freezing stops
moving_average_n = 20                                   # average every n points
updrafts = [FT(1.5), FT(1.4), FT(5), FT(1.5), FT(1.5)]  # updrafts matching AIDA cooling rate

# Additional definitions
tps = TD.Parameters.ThermodynamicsParameters(FT)
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v


for (exp_index, data_file_name) in enumerate(data_file_names)
    @info(data_file_name)
    #! format: off
    ### Unpacking experiment-specific variables.
    plot_name = plot_names[exp_index]
    w = updrafts[exp_index]
    start_time = start_time_list[exp_index]
    start_time_index = start_time .+ 100
    end_time = end_time_list[exp_index]
    end_time_index = end_time .+ 100

    nuc_mode = plot_name == "IN0701" || plot_name == "IN0719" ? "ABDINM" : "ABHOM"

    ## Moving average to smooth data.
    moving_average_start_index = Int32(start_time_index + (moving_average_n / 2))
    moving_average_end_index = Int32(end_time_index - (moving_average_n / 2))
    t_max = moving_average_end_index - moving_average_start_index

    ### Check for and grab data in AIDA_data folder.
    chamber_data = unpack_data(data_file_name)
    AIDA_data = grab_data(chamber_data)
    (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC, AIDA_e) = AIDA_data
    t_profile = AIDA_t_profile[moving_average_start_index:moving_average_end_index] .- (moving_average_start_index - 101)
    T_profile = AIDA_T_profile[moving_average_start_index:moving_average_end_index]
    P_profile = AIDA_P_profile[moving_average_start_index:moving_average_end_index]
    ICNC_profile = AIDA_ICNC[moving_average_start_index:moving_average_end_index]
    e_profile = AIDA_e[moving_average_start_index:moving_average_end_index]
    S_l_profile = e_profile ./ (TD.saturation_vapor_pressure.(tps, T_profile, TD.Liquid()))

    params =
        nuc_mode == "ABHOM" ?
        AIDA_IN05_params(FT, w, t_max, t_profile, T_profile, P_profile) :
        AIDA_IN07_params(FT, w, t_max, t_profile, T_profile, P_profile, plot_name)
    IC =
        nuc_mode == "ABHOM" ?
        AIDA_IN05_IC(FT, data_file_name) :
        AIDA_IN07_IC(FT, data_file_name)

    Nₜ = IC[7] + IC[8] + IC[9]
    frozen_frac = AIDA_ICNC[start_time_index:end_time_index] ./ Nₜ
    frozen_frac_moving_mean = moving_average(frozen_frac, moving_average_n)
    ICNC_moving_avg = moving_average(AIDA_ICNC[start_time_index:end_time_index], moving_average_n)
    Γ = 0.1 * LinearAlgebra.I * (maximum(frozen_frac_moving_mean) - minimum(frozen_frac_moving_mean)) # coeff is an estimated of the noise

    ### Calibration.
    EKI_output = calibrate_J_parameters_EKI(FT, nuc_mode, params, IC, frozen_frac_moving_mean, end_sim, Γ)
    UKI_output = calibrate_J_parameters_UKI(FT, nuc_mode, params, IC, frozen_frac_moving_mean, end_sim, Γ)
    
    EKI_n_iterations = size(EKI_output[2])[1]
    EKI_n_ensembles = size(EKI_output[2][1])[2]
    
    EKI_calibrated_parameters = EKI_output[1]
    UKI_calibrated_parameters = UKI_output[1]
    calibrated_ensemble_means = ensemble_means(EKI_output[2], EKI_n_iterations, EKI_n_ensembles)

    ## Calibrated parcel.
    EKI_parcel = run_model(params, nuc_mode, EKI_calibrated_parameters, FT, IC, end_sim)
    UKI_parcel = run_model(params, nuc_mode, UKI_calibrated_parameters, FT, IC, end_sim)

    ### Plots.
    ## Plotting AIDA data.
    AIDA_data_fig = MK.Figure(size = (800, 600), fontsize = 24)
    ax1 = MK.Axis(AIDA_data_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "AIDA data $plot_name")
    MK.lines!(ax1, AIDA_t_profile, AIDA_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
    MK.lines!(ax1, t_profile .+ moving_average_start_index .- 100, ICNC_moving_avg, label = "AIDA moving mean", linewidth = 2.5, color =:blue)
    MK.axislegend(ax1, framevisible = true, labelsize = 12, position = :rc)
    MK.save("$plot_name"*"_ICNC.svg", AIDA_data_fig)

    ## Calibrated coefficients.
    #  Did they converge?
    calibrated_coeffs_fig = MK.Figure(size = (600, 600), fontsize = 24)
    ax3 = MK.Axis(calibrated_coeffs_fig[1, 1], ylabel = "m coefficient [-]", title = "$plot_name")
    ax4 = MK.Axis(calibrated_coeffs_fig[2, 1], ylabel = "c coefficient [-]", xlabel = "iteration #")
    MK.lines!(ax3, collect(1:EKI_n_iterations), calibrated_ensemble_means[1], label = "ensemble mean", color = :orange, linewidth = 2.5)
    MK.lines!(ax4, collect(1:EKI_n_iterations), calibrated_ensemble_means[2], label = "ensemble mean", color = :orange, linewidth = 2.5)
    MK.save("$plot_name"*"_calibrated_coeffs_fig.svg", calibrated_coeffs_fig)

    ## Calibrated parcel simulations.
    #  Does the calibrated parcel give reasonable outputs?
    calibrated_parcel_fig = MK.Figure(size = (1500, 800), fontsize = 20)
    ax_parcel_1 = MK.Axis(calibrated_parcel_fig[1, 1], ylabel = "saturation [-]", xlabel = "time [s]", title = "$plot_name")
    ax_parcel_2 = MK.Axis(calibrated_parcel_fig[2, 1], ylabel = "liq mixing ratio [g/kg]", xlabel = "time [s]")
    ax_parcel_3 = MK.Axis(calibrated_parcel_fig[1, 2], ylabel = "temperature [K]", xlabel = "time [s]")
    ax_parcel_4 = MK.Axis(calibrated_parcel_fig[2, 2], ylabel = "qᵢ [g/kg]", xlabel = "time [s]")
    ax_parcel_5 = MK.Axis(calibrated_parcel_fig[3, 1], ylabel = "Nₗ [m^-3]", xlabel = "time [s]")
    ax_parcel_6 = MK.Axis(calibrated_parcel_fig[3, 2], ylabel = "Nᵢ [m^-3]", xlabel = "time [s]")
    ax_parcel_7 = MK.Axis(calibrated_parcel_fig[1, 3], ylabel = "pressure [Pa]", xlabel = "time [s]")
    ax_parcel_8 = MK.Axis(calibrated_parcel_fig[2, 3], ylabel = "Nₐ [m^-3]", xlabel = "time [s]")
    
    MK.lines!(ax_parcel_1, EKI_parcel.t, EKI_parcel[1, :], label = "EKI Calib Liq", color = :orange) # label = "liquid"
    MK.lines!(ax_parcel_1, UKI_parcel.t, UKI_parcel[1, :], label = "UKI Calib Liq", color = :fuchsia) # label = "liquid"
    MK.lines!(ax_parcel_1, EKI_parcel.t, S_i.(tps, EKI_parcel[3, :], EKI_parcel[1, :]), label = "EKI Calib Ice", color = :green)
    MK.lines!(ax_parcel_1, t_profile, S_l_profile, label = "chamber", color = :blue)
    
    MK.lines!(ax_parcel_2, EKI_parcel.t, EKI_parcel[5, :], color = :orange)
    MK.lines!(ax_parcel_2, UKI_parcel.t, UKI_parcel[5, :], color = :fuchsia)

    MK.lines!(ax_parcel_3, EKI_parcel.t, EKI_parcel[3, :], color = :orange)
    MK.lines!(ax_parcel_3, UKI_parcel.t, UKI_parcel[3, :], color = :fuchsia)
    MK.lines!(ax_parcel_3, t_profile, T_profile, color = :blue, linestyle =:dash)

    MK.lines!(ax_parcel_4, EKI_parcel.t, EKI_parcel[6, :], color = :orange)
    MK.lines!(ax_parcel_4, UKI_parcel.t, UKI_parcel[6, :], color = :fuchsia)

    MK.lines!(ax_parcel_5, EKI_parcel.t, EKI_parcel[8, :], color = :orange)
    MK.lines!(ax_parcel_5, UKI_parcel.t, UKI_parcel[8, :], color = :fuchsia)    

    MK.lines!(ax_parcel_6, EKI_parcel.t, EKI_parcel[9, :], color = :orange, label = "EKI")
    MK.lines!(ax_parcel_6, UKI_parcel.t, UKI_parcel[9, :], color = :fuchsia, label = "UKI")
    MK.lines!(ax_parcel_6, t_profile, ICNC_profile, color = :blue, label = "AIDA",)
    
    error = fill(sqrt(Γ[1,1]) * 2, length(AIDA_t_profile[start_time_index:end_time_index] .- start_time))
    MK.errorbars!(ax_parcel_6, AIDA_t_profile[start_time_index:end_time_index] .- start_time, AIDA_ICNC[start_time_index:end_time_index] ./ Nₜ, error)

    MK.lines!(ax_parcel_7, EKI_parcel.t, EKI_parcel[2, :], color = :orange)
    MK.lines!(ax_parcel_7, UKI_parcel.t, UKI_parcel[2, :], color = :fuchsia)
    MK.lines!(ax_parcel_7, t_profile, P_profile, color = :blue) #, linestyle =:dash

    MK.lines!(ax_parcel_8, EKI_parcel.t, EKI_parcel[7, :], color = :orange)
    MK.lines!(ax_parcel_8, UKI_parcel.t, UKI_parcel[7, :], color = :fuchsia)
    
    MK.axislegend(ax_parcel_6, framevisible = false, labelsize = 16, position = :rb)
    MK.save("$plot_name"*"_calibrated_parcel_fig.svg", calibrated_parcel_fig)

    ## Comparing AIDA data and calibrated parcel.
    #  Does calibrated parcel look like observations?
    ICNC_comparison_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_compare = MK.Axis(ICNC_comparison_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$plot_name")
    MK.lines!(
        ax_compare,
        EKI_parcel.t,
        EKI_parcel[9, :]./ Nₜ,
        label = "CM.jl Parcel (EKI Calibrated)",
        linewidth = 2.5,
        color =:orange,
    )
    MK.lines!(
        ax_compare,
        UKI_parcel.t,
        UKI_parcel[9, :]./ Nₜ,
        label = "CM.jl Parcel (UKI Calibrated)",
        linewidth = 2.5,
        color =:fuchsia,
    )
    MK.lines!(
        ax_compare,
        AIDA_t_profile[start_time_index:end_time_index] .- start_time,
        frozen_frac,
        label = "Raw AIDA",
        linewidth = 2,
        color =:blue,
        linestyle =:dash,
    )
    MK.lines!(
        ax_compare,
        t_profile,
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color =:blue
    )
    error = fill(sqrt(Γ[1,1]) * 2, length(frozen_frac_moving_mean))
    MK.errorbars!(ax_compare, t_profile, frozen_frac_moving_mean, error)

    MK.axislegend(ax_compare, framevisible = false, labelsize = 20, position = :rb)
    MK.save("$plot_name"*"_ICNC_comparison_fig.svg", ICNC_comparison_fig)

    ## Looking at spread in UKI calibrated parameters
    ϕ_UKI = UKI_output[2]
    UKI_parcel_1 = run_model(params, nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, IC, end_sim)
    UKI_parcel_2 = run_model(params, nuc_mode, [ϕ_UKI[1,2], ϕ_UKI[2,2]], FT, IC, end_sim)
    UKI_parcel_3 = run_model(params, nuc_mode, [ϕ_UKI[1,3], ϕ_UKI[2,3]], FT, IC, end_sim)
    UKI_parcel_4 = run_model(params, nuc_mode, [ϕ_UKI[1,4], ϕ_UKI[2,4]], FT, IC, end_sim)
    UKI_parcel_5 = run_model(params, nuc_mode, [ϕ_UKI[1,5], ϕ_UKI[2,5]], FT, IC, end_sim)

    UKI_spread_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_spread = MK.Axis(UKI_spread_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$plot_name")
    MK.lines!(
        ax_spread,
        t_profile,
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color =:blue
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_1.t,
        UKI_parcel_1[9, :]./ Nₜ,
        linewidth = 2.5,
        color =:grey80,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_2.t,
        UKI_parcel_2[9, :]./ Nₜ,
        linewidth = 2.5,
        color =:grey60,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_3.t,
        UKI_parcel_3[9, :]./ Nₜ,
        linewidth = 2.5,
        color =:grey40,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_4.t,
        UKI_parcel_4[9, :]./ Nₜ,
        linewidth = 2.5,
        color =:grey25,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_5.t,
        UKI_parcel_5[9, :]./ Nₜ,
        linewidth = 2.5,
        color =:grey12,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel.t,
        UKI_parcel[9, :]./ Nₜ,
        label = "CM.jl Parcel (UKI Calibrated)",
        linewidth = 2.5,
        color =:fuchsia,
        linestyle = :dash,
    )
    error = fill(sqrt(Γ[1,1]) * 2, length(frozen_frac_moving_mean))
    MK.errorbars!(ax_spread, t_profile, frozen_frac_moving_mean, error)
    MK.save("$plot_name"*"_UKI_spread_fig.svg", UKI_spread_fig)

    ## Freezing rates
    freezing_rates_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_ABDINM = MK.Axis(freezing_rates_fig[1, 1], ylabel = "ABDINM", xlabel = "time [s]", title = "$plot_name")
    ax_ABIFM = MK.Axis(freezing_rates_fig[1, 2], ylabel = "ABIFM", xlabel = "time [s]", title = "dNᵢ/dt [#/m^3/s]")
    ax_ABHOM = MK.Axis(freezing_rates_fig[2, 1], ylabel = "ABHOM", xlabel = "time [s]")
    MK.lines!(
        ax_ABDINM,
        EKI_parcel.t,
        EKI_parcel[11, :],
        label = "EKI ABDINM",
        linewidth = 2.5,
        color =:orange
    )
    MK.lines!(
        ax_ABDINM,
        UKI_parcel.t,
        UKI_parcel[11, :],
        label = "UKI ABDINM",
        linewidth = 2.5,
        color =:fuchsia
    )
    MK.lines!(
        ax_ABIFM,
        EKI_parcel.t,
        EKI_parcel[12, :],
        label = "EKI ABIFM",
        linewidth = 2.5,
        color =:orange
    )
    MK.lines!(
        ax_ABIFM,
        UKI_parcel.t,
        UKI_parcel[12, :],
        label = "UKI ABIFM",
        linewidth = 2.5,
        color =:fuchsia
    )
    MK.lines!(
        ax_ABHOM,
        EKI_parcel.t,
        EKI_parcel[13, :],
        label = "EKI ABHOM",
        linewidth = 2.5,
        color =:orange
    )
    MK.lines!(
        ax_ABHOM,
        UKI_parcel.t,
        UKI_parcel[13, :],
        label = "UKI ABHOM",
        linewidth = 2.5,
        color =:fuchsia
    )
    MK.save("$plot_name"*"_freezing_rates_fig.svg", freezing_rates_fig)


    #! format: on
end
