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
moving_average(data, n) =
    [sum(@view data[i:(i + n - 1)]) / n for i in 1:(length(data) - (n - 1))]

# Defining data names, start/end times, etc.
data_file_names = ["in05_17_aida.edf", "in05_18_aida.edf", "in05_19_aida.edf"]
plot_names = ["IN0517", "IN0518", "IN0519"]
end_sim = 25                               # Loss func looks at last end_sim timesteps only
start_time_list = [195, 180, 80] .+ 100    # AIDA time starts at -100 seconds. Add freezing onset as needed.
end_time_list = [290, 290, 170] .+ 100     # approximate time freezing stops
moving_average_n = 20                      # average every n points
updrafts = [FT(1.5), FT(1.4), FT(5)]       # updrafts matching AIDA cooling rate

# Additional definitions
ips = CMP.IceNucleationParameters(FT)
tps = TD.Parameters.ThermodynamicsParameters(FT)
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v


for (exp_index, data_file_name) in enumerate(data_file_names)
    #! format: off
    ### Unpacking experiment-specific variables.
    plot_name = plot_names[exp_index]
    w = updrafts[exp_index]
    start_time = start_time_list[exp_index]
    end_time = end_time_list[exp_index]

    ## Moving average to smooth data.
    moving_average_start = Int32(start_time + (moving_average_n / 2 - 1))
    moving_average_end = Int32(end_time - (moving_average_n / 2))
    t_max =
        (end_time - 100 - (moving_average_n / 2)) -     # AIDA time starts at -100 seconds (t = 0 at index 101)
        (start_time - 100 + (moving_average_n / 2 - 1)) # duration of simulation

    params = AIDA_IN05_params(FT, w, t_max)
    IC = AIDA_IN05_IC(FT, data_file_name)

    ### Check for data file in AIDA_data folder.
    chamber_data = unpack_data(data_file_name)
    ICNC = chamber_data[start_time:end_time, 6] .* 1e6 # Nᵢ [m^-3]

    frozen_frac = ICNC ./ (IC[7] + IC[8] + IC[9])
    frozen_frac_moving_mean = moving_average(frozen_frac, moving_average_n)
    ICNC_moving_avg = moving_average(ICNC, moving_average_n)
    Γ = 0.1 * LinearAlgebra.I * (maximum(frozen_frac_moving_mean) - minimum(frozen_frac_moving_mean)) # coeff is an estimated of the noise

    ## Calculating some variables from AIDA data.
    chamber_S_l = zeros(FT, size(chamber_data[start_time:end_time, 4], 1), size(chamber_data[start_time:end_time, 4], 2))
    for (i, e) in enumerate(chamber_data[start_time:end_time, 4])
        temp = chamber_data[start_time:end_time, 3][i]
        eₛ = TD.saturation_vapor_pressure(tps, temp, TD.Liquid())
        chamber_S_l[i] = e / eₛ
    end
    chamber_r_l = zeros(FT, size(chamber_data[start_time:end_time, 8], 1), size(chamber_data[start_time:end_time, 8], 2))
    for (i, r_l) in enumerate(chamber_data[start_time:end_time, 8] .* 1e-6)
        chamber_r_l[i] = r_l < 0 ? FT(0) : r_l
    end

    ### Calibration.
    EKI_output = calibrate_J_parameters_EKI(FT, "ABHOM", params, IC, frozen_frac_moving_mean, end_sim, Γ)
    UKI_output = calibrate_J_parameters_UKI(FT, "ABHOM", params, IC, frozen_frac_moving_mean, end_sim, Γ)
    
    EKI_n_iterations = size(EKI_output[2])[1]
    EKI_n_ensembles = size(EKI_output[2][1])[2]
    
    EKI_calibrated_parameters = EKI_output[1]
    UKI_calibrated_parameters = UKI_output[1]
    calibrated_ensemble_means = ensemble_means(EKI_output[2], EKI_n_iterations, EKI_n_ensembles)

    ## Calibrated parcel.
    EKI_parcel = run_calibrated_model(FT, "ABHOM", EKI_calibrated_parameters, params, IC)
    UKI_parcel = run_calibrated_model(FT, "ABHOM", UKI_calibrated_parameters, params, IC)
    parcel_default = run_calibrated_model(FT, "ABHOM", [FT(255.927125), FT(-68.553283)], params, IC)



    ### Plots.
    ## Plotting AIDA data.
    AIDA_data_fig = MK.Figure(size = (800, 600), fontsize = 24)
    ax1 = MK.Axis(AIDA_data_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "AIDA data $plot_name")
    MK.lines!(
        ax1,
        chamber_data[start_time:end_time, 1],
        ICNC,
        label = "Raw AIDA",
        color =:blue,
        linestyle =:dash,
        linewidth = 2,
    )
    MK.lines!(
        ax1,
        chamber_data[moving_average_start:moving_average_end, 1],
        ICNC_moving_avg,
        label = "AIDA moving mean",
        linewidth = 2.5,
        color =:blue,
    )
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
    
    MK.lines!(ax_parcel_1, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[1, :], label = "EKI Calib Liq", color = :orange) # label = "liquid"
    MK.lines!(ax_parcel_1, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[1, :], label = "UKI Calib Liq", color = :fuchsia) # label = "liquid"
    MK.lines!(ax_parcel_1, parcel_default.t .+ (moving_average_n / 2), parcel_default[1, :], label = "default", color = :darkorange2)
    MK.lines!(ax_parcel_1, EKI_parcel.t .+ (moving_average_n / 2), S_i.(tps, EKI_parcel[3, :], EKI_parcel[1, :]), label = "EKI Calib Ice", color = :green)
    # MK.lines!(
    #     ax_parcel_1,
    #     chamber_data[start_time:end_time, 1] .- (start_time - 100),
    #     vec(chamber_S_l),
    #     label = "chamber",
    #     color = :blue
    # )
    
    MK.lines!(ax_parcel_2, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[5, :], color = :orange)
    MK.lines!(ax_parcel_2, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[5, :], color = :fuchsia)
    MK.lines!(ax_parcel_2, parcel_default.t .+ (moving_average_n / 2), parcel_default[5,:], color = :darkorange2)
    
    MK.lines!(ax_parcel_3, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[3, :], color = :orange)
    MK.lines!(ax_parcel_3, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[3, :], color = :fuchsia)
    MK.lines!(ax_parcel_3, parcel_default.t .+ (moving_average_n / 2), parcel_default[3, :], color = :darkorange2)
    MK.lines!(
        ax_parcel_3,
        chamber_data[start_time:end_time, 1] .- (start_time - 100),
        chamber_data[start_time:end_time, 3],
        color = :blue
    )

    MK.lines!(ax_parcel_4, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[6, :], color = :orange)
    MK.lines!(ax_parcel_4, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[6, :], color = :fuchsia)
    MK.lines!(ax_parcel_4, parcel_default.t .+ (moving_average_n / 2), parcel_default[6, :], color = :darkorange2)

    MK.lines!(ax_parcel_5, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[8, :], color = :orange)
    MK.lines!(ax_parcel_5, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[8, :], color = :fuchsia)    
    MK.lines!(ax_parcel_5, parcel_default.t .+ (moving_average_n / 2), parcel_default[8, :], color = :darkorange2)

    MK.lines!(ax_parcel_6, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[9, :], color = :orange, label = "EKI")
    MK.lines!(ax_parcel_6, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[9, :], color = :fuchsia, label = "UKI")
    MK.lines!(ax_parcel_6, parcel_default.t .+ (moving_average_n / 2), parcel_default[9, :], color = :darkorange2, label = "Pre-Calib")
    MK.lines!(ax_parcel_6, 
        chamber_data[start_time:end_time, 1] .- (start_time - 100),
        ICNC,
        color = :blue,
        label = "AIDA",
    )
    error = fill(sqrt(Γ[1,1]) * 2, length(chamber_data[start_time:end_time, 1] .- (start_time - 100)))
    MK.errorbars!(ax_parcel_6, chamber_data[start_time:end_time, 1] .- (start_time - 100), ICNC ./ (IC[7] + IC[8] + IC[9]), error)

    MK.lines!(ax_parcel_7, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[2, :], color = :orange)
    MK.lines!(ax_parcel_7, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[2, :], color = :fuchsia)
    MK.lines!(
        ax_parcel_7,
        chamber_data[start_time:end_time, 1] .- (start_time - 100),
        chamber_data[start_time:end_time, 2] .* 1e2,
        color = :blue
    )

    MK.lines!(ax_parcel_8, EKI_parcel.t .+ (moving_average_n / 2), EKI_parcel[7, :], color = :orange)
    MK.lines!(ax_parcel_8, UKI_parcel.t .+ (moving_average_n / 2), UKI_parcel[7, :], color = :fuchsia)
    
    MK.axislegend(ax_parcel_6, framevisible = false, labelsize = 16, position = :rb)
    MK.save("$plot_name"*"_calibrated_parcel_fig.svg", calibrated_parcel_fig)

    ## Comparing AIDA data and calibrated parcel.
    #  Does calibrated parcel look like observations?
    ICNC_comparison_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_compare = MK.Axis(ICNC_comparison_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$plot_name")
    # MK.lines!(
    #     ax_compare,
    #     parcel_default.t .+ (moving_average_n / 2),
    #     parcel_default[9, :]./  (IC[7] + IC[8] + IC[9]),
    #     label = "CM.jl Parcel (Pre-Calib)",
    #     linewidth = 2.5,
    #     color =:orangered,
    # )
    MK.lines!(
        ax_compare,
        EKI_parcel.t .+ (moving_average_n / 2),
        EKI_parcel[9, :]./ (IC[7] + IC[8] + IC[9]),
        label = "CM.jl Parcel (EKI Calibrated)",
        linewidth = 2.5,
        color =:orange,
    )
    MK.lines!(
        ax_compare,
        UKI_parcel.t .+ (moving_average_n / 2),
        UKI_parcel[9, :]./ (IC[7] + IC[8] + IC[9]),
        label = "CM.jl Parcel (UKI Calibrated)",
        linewidth = 2.5,
        color =:fuchsia,
    )
    MK.lines!(
        ax_compare,
        chamber_data[start_time:end_time, 1] .- (start_time - 100),
        frozen_frac,
        label = "Raw AIDA",
        linewidth = 2,
        color =:blue,
        linestyle =:dash,
    )
    MK.lines!(
        ax_compare,
        chamber_data[moving_average_start:moving_average_end, 1] .- (start_time - 100),
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color =:blue
    )
    error = fill(sqrt(Γ[1,1]) * 2, length(frozen_frac_moving_mean))
    MK.errorbars!(ax_compare, chamber_data[moving_average_start:moving_average_end, 1] .- (start_time - 100), frozen_frac_moving_mean, error)

    MK.axislegend(ax_compare, framevisible = false, labelsize = 20, position = :rb)
    MK.save("$plot_name"*"_ICNC_comparison_fig.svg", ICNC_comparison_fig)

    ## Looking at spread in UKI calibrated parameters
    ϕ_UKI = UKI_output[2]
    UKI_parcel_1 = run_calibrated_model(FT, "ABHOM", [ϕ_UKI[1,1], ϕ_UKI[2,1]], params, IC)
    UKI_parcel_2 = run_calibrated_model(FT, "ABHOM", [ϕ_UKI[1,2], ϕ_UKI[2,2]], params, IC)
    UKI_parcel_3 = run_calibrated_model(FT, "ABHOM", [ϕ_UKI[1,3], ϕ_UKI[2,3]], params, IC)
    UKI_parcel_4 = run_calibrated_model(FT, "ABHOM", [ϕ_UKI[1,4], ϕ_UKI[2,4]], params, IC)
    UKI_parcel_5 = run_calibrated_model(FT, "ABHOM", [ϕ_UKI[1,5], ϕ_UKI[2,5]], params, IC)

    UKI_spread_fig = MK.Figure(size = (700, 600), fontsize = 24)
    ax_spread = MK.Axis(UKI_spread_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$plot_name")
    MK.lines!(
        ax_spread,
        chamber_data[moving_average_start:moving_average_end, 1] .- (start_time - 100),
        frozen_frac_moving_mean,
        label = "AIDA Moving Avg",
        linewidth = 2.5,
        color =:blue
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_1.t .+ (moving_average_n / 2),
        UKI_parcel_1[9, :]./ (IC[7] + IC[8] + IC[9]),
        linewidth = 2.5,
        color =:grey80,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_2.t .+ (moving_average_n / 2),
        UKI_parcel_2[9, :]./ (IC[7] + IC[8] + IC[9]),
        linewidth = 2.5,
        color =:grey60,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_3.t .+ (moving_average_n / 2),
        UKI_parcel_3[9, :]./ (IC[7] + IC[8] + IC[9]),
        linewidth = 2.5,
        color =:grey40,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_4.t .+ (moving_average_n / 2),
        UKI_parcel_4[9, :]./ (IC[7] + IC[8] + IC[9]),
        linewidth = 2.5,
        color =:grey25,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel_5.t .+ (moving_average_n / 2),
        UKI_parcel_5[9, :]./ (IC[7] + IC[8] + IC[9]),
        linewidth = 2.5,
        color =:grey12,
    )
    MK.lines!(
        ax_spread,
        UKI_parcel.t .+ (moving_average_n / 2),
        UKI_parcel[9, :]./ (IC[7] + IC[8] + IC[9]),
        label = "CM.jl Parcel (UKI Calibrated)",
        linewidth = 2.5,
        color =:fuchsia,
        linestyle = :dash,
    )
    error = fill(sqrt(Γ[1,1]) * 2, length(frozen_frac_moving_mean))
    MK.errorbars!(ax_spread, chamber_data[moving_average_start:moving_average_end, 1] .- (start_time - 100), frozen_frac_moving_mean, error)
    MK.save("$plot_name"*"_UKI_spread_fig.svg", UKI_spread_fig)

    #! format: on
end
