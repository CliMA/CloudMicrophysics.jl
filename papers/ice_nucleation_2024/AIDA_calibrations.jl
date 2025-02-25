import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD
import CairoMakie as MK

# To grab data
import DelimitedFiles
using LazyArtifacts
using ClimaUtilities.ClimaArtifacts
import NaNStatistics

FT = Float64
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "unpack_AIDA.jl"))

# Defining data names, start/end times, etc.
global edf_data_names = [
    "in05_17_aida.edf", "in05_18_aida.edf",
    "in07_01_aida.edf", "in07_19_aida.edf",
]
data_file_names = [
    ["in05_17_aida.edf", "in05_18_aida.edf", "TROPIC04"],
    ["in07_01_aida.edf"],
    ["in07_19_aida.edf"],
    ["ACI04_22"],
    ["EXP19"],
    ["EXP45"],
    # ["EXP28"],
]
batch_names = ["HOM", "IN0701", "IN0719", "ACI04_22", "EXP19", "EXP45"]
end_sim = 25            # Loss func looks at last end_sim timesteps only
start_time_list = [     # freezing onset
    [Int32(150), Int32(180), Int32(0)],
    [Int32(50)],
    [Int32(35)],
    [Int32(0)],
    [Int32(0)],
    [Int32(0)],
    # [Int32(0)],
]
end_time_list = [       # approximate time freezing stops
    [Int32(295), Int32(290), Int32(700)],
    [Int32(375)],
    [Int32(375)],
    [Int32(300)],
    [Int32(200)],
    [Int32(200)],
    # [Int32(300)]
]
moving_average_n = 5    # average every length(data) / n points
updrafts = [            # updrafts matching AIDA cooling rate
    [FT(1.5), FT(1.4), FT(1.4)],
    [FT(1.5)],
    [FT(1.5)],
    [FT(1.5)],
    [FT(9)],
    [FT(1.5)],
    # [FT(1.5)],
]

# Additional definitions
tps = TD.Parameters.ThermodynamicsParameters(FT)
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v

global EKI_calibratated_coeff_dict = Dict()
global UKI_calibratated_coeff_dict = Dict()

for (batch_index, batch_name) in enumerate(batch_names)
    @info(batch_name)
    #! format: off
    ### Unpacking experiment-specific variables.
    data_file_name_list = data_file_names[batch_index]
    w = updrafts[batch_index]
    start_time = start_time_list[batch_index]
    end_time = end_time_list[batch_index]
    t_max = end_time .- start_time
    start_time_index = start_time .+ 1
    end_time_index = end_time .+ 1

    if batch_name == "HOM"
        nuc_mode = "ABHOM"
    elseif batch_name == "IN0701" || batch_name == "IN0719" || batch_name == "EXP45"
        nuc_mode = "ABDINM"
    elseif batch_name == "ACI04_22" ||batch_name == "EXP19"
        nuc_mode = "ABIFM"
    end

    ### Check for and grab data in AIDA_data folder.
    calib_variables = data_to_calib_inputs(
        FT, batch_name, data_file_name_list, nuc_mode,
        start_time, end_time,
        w, t_max, end_sim,
    )
    (;
        t_profile, T_profile, P_profile, ICNC_profile, e_profile, S_l_profile,
        params_list, IC_list,
        frozen_frac, frozen_frac_moving_mean, ICNC_moving_avg,
        Nₜ, y_truth,
    ) = calib_variables

    ff_max = maximum.(frozen_frac_moving_mean)
    ff_min = minimum.(frozen_frac_moving_mean)
    Γ = 0.01 / 9 * LinearAlgebra.I * (maximum(ff_max) - minimum(ff_min))

    ### Calibration.
    EKI_output = calibrate_J_parameters_EKI(FT, nuc_mode, params_list, IC_list, y_truth, end_sim, Γ)
    UKI_output = calibrate_J_parameters_UKI(FT, nuc_mode, params_list, IC_list, y_truth, end_sim, Γ)

    EKI_n_iterations = size(EKI_output[2])[1]
    EKI_n_ensembles = size(EKI_output[2][1])[2]

    EKI_calibrated_parameters = EKI_output[1]
    UKI_calibrated_parameters = UKI_output[1]
    calibrated_ensemble_means = ensemble_means(EKI_output[2], EKI_n_iterations, EKI_n_ensembles)
    merge!(EKI_calibratated_coeff_dict, Dict(batch_name => EKI_calibrated_parameters))
    merge!(UKI_calibratated_coeff_dict, Dict(batch_name => UKI_calibrated_parameters))

    ### Plot.
    for (exp_index, data_file) in enumerate(data_file_name_list)
        if nuc_mode == "ABHOM"
            if exp_index == 1
                exp_name = "IN0517"
            elseif exp_index == 2
                exp_name = "IN0518"
            elseif exp_index == 3
                exp_name = "TROPIC04"
            end
        else
            exp_name = batch_name
        end

        ## Calibrated parcel.
        EKI_parcel = run_model([params_list[exp_index]], nuc_mode, EKI_calibrated_parameters, FT, [IC_list[exp_index]], end_sim)
        UKI_parcel = run_model([params_list[exp_index]], nuc_mode, UKI_calibrated_parameters, FT, [IC_list[exp_index]], end_sim)

        ## Plots.
        ## Plotting AIDA data.
        AIDA_data = data_file in edf_data_names ? unpack_data(data_file) : unpack_data(data_file, total_t = t_max[exp_index])
        (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC_profile, AIDA_e_profile) = AIDA_data

        AIDA_data_fig = MK.Figure(size = (1000, 600), fontsize = 24)
        data_ax1 = MK.Axis(AIDA_data_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "AIDA data $exp_name")
        data_ax2 = MK.Axis(AIDA_data_fig[1, 2], ylabel = "Frozen Frac Moving Mean [-]", xlabel = "time [s]", title = "AIDA data $exp_name")
        MK.lines!(data_ax1, AIDA_t_profile, AIDA_ICNC_profile, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
        MK.lines!(data_ax1, t_profile[exp_index], ICNC_moving_avg[exp_index], label = "AIDA ICNC moving mean", linewidth = 2.5, color =:blue)
        MK.lines!(data_ax2, t_profile[exp_index], frozen_frac_moving_mean[exp_index], linewidth = 2.5, color =:blue)
        MK.axislegend(data_ax1, framevisible = true, labelsize = 12, position = :rc)
        MK.save("$exp_name"*"_ICNC.svg", AIDA_data_fig)

        ## Calibrated coefficients.
        #  Did they converge?
        calibrated_coeffs_fig = MK.Figure(size = (1100, 900), fontsize = 24)
        ax3 = MK.Axis(calibrated_coeffs_fig[1, 1], ylabel = "m coefficient [-]", title = "$batch_name")
        ax4 = MK.Axis(calibrated_coeffs_fig[1, 2], ylabel = "c coefficient [-]", xlabel = "iteration #", title = "EKI")
    
        MK.lines!(ax3, collect(1:EKI_n_iterations), calibrated_ensemble_means[1], label = "ensemble mean", color = :orange, linewidth = 2.5)
        MK.lines!(ax4, collect(1:EKI_n_iterations), calibrated_ensemble_means[2], label = "ensemble mean", color = :orange, linewidth = 2.5)
    
        MK.save("$batch_name"*"_calibrated_coeffs_fig.svg", calibrated_coeffs_fig)

        ## Calibrated parcel simulations.
        #  Does the calibrated parcel give reasonable outputs?
        calibrated_parcel_fig = MK.Figure(size = (1500, 800), fontsize = 20)
        ax_parcel_1 = MK.Axis(calibrated_parcel_fig[1, 1], ylabel = "saturation [-]", xlabel = "time [s]", title = "$exp_name")
        ax_parcel_2 = MK.Axis(calibrated_parcel_fig[2, 1], ylabel = "liq mixing ratio [g/kg]", xlabel = "time [s]")
        ax_parcel_3 = MK.Axis(calibrated_parcel_fig[1, 2], ylabel = "temperature [K]", xlabel = "time [s]")
        ax_parcel_4 = MK.Axis(calibrated_parcel_fig[2, 2], ylabel = "qᵢ [g/kg]", xlabel = "time [s]")
        ax_parcel_5 = MK.Axis(calibrated_parcel_fig[3, 1], ylabel = "Nₗ [m^-3]", xlabel = "time [s]")
        ax_parcel_6 = MK.Axis(calibrated_parcel_fig[3, 2], ylabel = "Nᵢ [m^-3]", xlabel = "time [s]")
        ax_parcel_7 = MK.Axis(calibrated_parcel_fig[1, 3], ylabel = "pressure [Pa]", xlabel = "time [s]")
        ax_parcel_8 = MK.Axis(calibrated_parcel_fig[2, 3], ylabel = "Nₐ [m^-3]", xlabel = "time [s]")
        
        MK.lines!(ax_parcel_1, EKI_parcel.t, EKI_parcel[1, :], label = "EKI Calib Liq", color = :orange) # label = "liquid"
        MK.lines!(ax_parcel_1, UKI_parcel.t, UKI_parcel[1, :], label = "UKI Calib Liq", color = :fuchsia) # label = "liquid"
        MK.lines!(ax_parcel_1, EKI_parcel.t, S_i.(tps, EKI_parcel[3, :], EKI_parcel[1, :]), label = "EKI Calib Ice", color = :orange, linestyle = :dash)
        MK.lines!(ax_parcel_1, UKI_parcel.t, S_i.(tps, UKI_parcel[3, :], UKI_parcel[1, :]), label = "UKI Calib Ice", color = :fuchsia, linestyle = :dash)
        # MK.lines!(ax_parcel_1, t_profile[exp_index], S_l_profile[exp_index], label = "chamber", color = :blue)
        
        MK.lines!(ax_parcel_2, EKI_parcel.t, EKI_parcel[5, :], color = :orange)
        MK.lines!(ax_parcel_2, UKI_parcel.t, UKI_parcel[5, :], color = :fuchsia)

        MK.lines!(ax_parcel_3, EKI_parcel.t, EKI_parcel[3, :], color = :orange)
        MK.lines!(ax_parcel_3, UKI_parcel.t, UKI_parcel[3, :], color = :fuchsia)
        MK.lines!(ax_parcel_3, t_profile[exp_index], T_profile[exp_index], color = :blue, linestyle =:dash)

        MK.lines!(ax_parcel_4, EKI_parcel.t, EKI_parcel[6, :], color = :orange)
        MK.lines!(ax_parcel_4, UKI_parcel.t, UKI_parcel[6, :], color = :fuchsia)

        MK.lines!(ax_parcel_5, EKI_parcel.t, EKI_parcel[8, :], color = :orange)
        MK.lines!(ax_parcel_5, UKI_parcel.t, UKI_parcel[8, :], color = :fuchsia)    

        MK.lines!(ax_parcel_6, EKI_parcel.t, EKI_parcel[9, :], color = :orange, label = "EKI")
        MK.lines!(ax_parcel_6, UKI_parcel.t, UKI_parcel[9, :], color = :fuchsia, label = "UKI")
        MK.lines!(ax_parcel_6, t_profile[exp_index], ICNC_profile[exp_index], color = :blue, label = "AIDA",)

        error = (AIDA_ICNC_profile[start_time_index[exp_index]:end_time_index[exp_index]] ./ Nₜ[exp_index]) .* 0.1
        MK.errorbars!(ax_parcel_6, AIDA_t_profile[start_time_index[exp_index]:end_time_index[exp_index]] .- start_time[exp_index], AIDA_ICNC_profile[start_time_index[exp_index]:end_time_index[exp_index]] ./ Nₜ[exp_index], error, color = (:blue, 0.3))

        MK.lines!(ax_parcel_7, EKI_parcel.t, EKI_parcel[2, :], color = :orange)
        MK.lines!(ax_parcel_7, UKI_parcel.t, UKI_parcel[2, :], color = :fuchsia)
        MK.lines!(ax_parcel_7, t_profile[exp_index], P_profile[exp_index], color = :blue, linestyle =:dash)

        MK.lines!(ax_parcel_8, EKI_parcel.t, EKI_parcel[7, :], color = :orange)
        MK.lines!(ax_parcel_8, UKI_parcel.t, UKI_parcel[7, :], color = :fuchsia)
        
        MK.axislegend(ax_parcel_6, framevisible = false, labelsize = 16, position = :rb)            
        MK.save("$exp_name"*"_calibrated_parcel_fig.svg", calibrated_parcel_fig)

        ## Comparing AIDA data and calibrated parcel.
        #  Does calibrated parcel look like observations?
        ICNC_comparison_fig = MK.Figure(size = (700, 600), fontsize = 24)
        ax_compare = MK.Axis(ICNC_comparison_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")
        MK.lines!(
            ax_compare,
            EKI_parcel.t,
            EKI_parcel[9, :]./ Nₜ[exp_index],
            label = "CM.jl Parcel (EKI Calibrated)",
            linewidth = 2.5,
            color =:orange,
        )
        MK.lines!(
            ax_compare,
            UKI_parcel.t,
            UKI_parcel[9, :]./ Nₜ[exp_index],
            label = "CM.jl Parcel (UKI Calibrated)",
            linewidth = 2.5,
            color =:fuchsia,
        )
        error = frozen_frac_moving_mean[exp_index] .* 0.1
        MK.errorbars!(ax_compare, t_profile[exp_index], frozen_frac_moving_mean[exp_index], error, color = (:blue, 0.3))
        MK.lines!(
            ax_compare,
            AIDA_t_profile[start_time_index[exp_index]:end_time_index[exp_index]] .- start_time[exp_index],
            frozen_frac[exp_index],
            label = "Raw AIDA",
            linewidth = 2,
            color =:blue,
            linestyle =:dash,
        )
        MK.lines!(
            ax_compare,
            t_profile[exp_index],
            frozen_frac_moving_mean[exp_index],
            label = "AIDA Moving Avg",
            linewidth = 2.5,
            color =:blue
        )
        
        MK.axislegend(ax_compare, framevisible = false, labelsize = 20, position = :rb)
        MK.save("$exp_name"*"_ICNC_comparison_fig.svg", ICNC_comparison_fig)

        ## Looking at spread in UKI calibrated parameters
        ϕ_UKI = UKI_output[2]
        UKI_parcel_1 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_2 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_3 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_4 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_5 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)

        UKI_spread_fig = MK.Figure(size = (700, 600), fontsize = 24)
        ax_spread = MK.Axis(UKI_spread_fig[1, 1], ylabel = "Frozen Fraction [-]", xlabel = "time [s]", title = "$exp_name")
        MK.lines!(
            ax_spread,
            t_profile[exp_index],
            frozen_frac_moving_mean[exp_index],
            label = "AIDA Moving Avg",
            linewidth = 2.5,
            color =:blue
        )
        MK.lines!(
            ax_spread,
            UKI_parcel_1.t,
            UKI_parcel_1[9, :]./ Nₜ[exp_index],
            linewidth = 2.5,
            color =:grey80,
        )
        MK.lines!(
            ax_spread,
            UKI_parcel_2.t,
            UKI_parcel_2[9, :]./ Nₜ[exp_index],
            linewidth = 2.5,
            color =:grey60,
        )
        MK.lines!(
            ax_spread,
            UKI_parcel_3.t,
            UKI_parcel_3[9, :]./ Nₜ[exp_index],
            linewidth = 2.5,
            color =:grey40,
        )
        MK.lines!(
            ax_spread,
            UKI_parcel_4.t,
            UKI_parcel_4[9, :]./ Nₜ[exp_index],
            linewidth = 2.5,
            color =:grey25,
        )
        MK.lines!(
            ax_spread,
            UKI_parcel_5.t,
            UKI_parcel_5[9, :]./ Nₜ[exp_index],
            linewidth = 2.5,
            color =:grey12,
        )
        MK.lines!(
            ax_spread,
            UKI_parcel.t,
            UKI_parcel[9, :]./ Nₜ[exp_index],
            label = "CM.jl Parcel (UKI Calibrated)",
            linewidth = 2.5,
            color =:fuchsia,
            linestyle = :dash,
        )
        error = frozen_frac_moving_mean[exp_index] .* 0.1
        MK.errorbars!(ax_spread, t_profile[exp_index], frozen_frac_moving_mean[exp_index], error, color = (:blue, 0.3))
        MK.save("$exp_name"*"_UKI_spread_fig.svg", UKI_spread_fig)
    end # plotting
    #! format: on
end
