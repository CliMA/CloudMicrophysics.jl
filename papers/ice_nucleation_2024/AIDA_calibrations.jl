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
function moving_average(data, n)
    window_size = length(data) / n
    moving_avg = NaNStatistics.movmean(data, window_size)
    #[sum(@view data[i:(i + n)]) / n for i in 1:(length(data) - n)]
    return moving_avg
end

# Defining data names, start/end times, etc.
data_file_names = [
    ["in05_17_aida.edf", "in05_18_aida.edf"],
    ["in07_01_aida.edf"],
    ["in07_19_aida.edf"],
]
batch_names = ["IN05", "IN0701", "IN0719"]
end_sim = 25               # Loss func looks at last end_sim timesteps only
start_time_list =          # freezing onset
    [[Int32(150), Int32(180)], [Int32(50)], [Int32(35)]]
end_time_list =            # approximate time freezing stops
    [[Int32(315), Int32(290)], [Int32(375)], [Int32(375)]]
moving_average_n = 3      # average every length(data) / n points
updrafts = [[FT(1.5), FT(1.4)], [FT(1.5)], [FT(1.5)]]  # updrafts matching AIDA cooling rate

# Additional definitions
tps = TD.Parameters.ThermodynamicsParameters(FT)
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v

global EKI_calibratated_coeff_dict = Dict()
global UKI_calibratated_coeff_dict = Dict()

for (calib_index, batch_name) in enumerate(batch_names)
    @info(batch_name)
    #! format: off
    ### Unpacking experiment-specific variables.
    data_file_name_list = data_file_names[calib_index]
    w = updrafts[calib_index]
    start_time = start_time_list[calib_index]
    start_time_index = start_time .+ Int32(100)
    end_time = end_time_list[calib_index]
    end_time_index = end_time .+ Int32(100)
    t_max = end_time_index .- start_time_index

    nuc_mode = batch_name == "IN05" ? "ABHOM" : "ABDINM"

    ### Check for and grab data in AIDA_data folder.
    t_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    T_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    P_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    ICNC_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    e_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    S_l_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)

    params_list = Vector{NamedTuple{
                    (:const_dt, :w, :t_max, :ips,
                     :prescribed_thermodynamics, :t_profile, :T_profile, :P_profile,
                     :aerosol_act, :aerosol, :r_nuc, :aero_σ_g,
                     :condensation_growth, :deposition_growth,
                     :liq_size_distribution, :ice_size_distribution,
                     :dep_nucleation, :heterogeneous, :homogeneous
                    ),
                    Tuple{
                        Float64, Float64, Int32, CM.Parameters.IceNucleationParameters{Float64, CM.Parameters.Mohler2006{Float64}, CM.Parameters.Koop2000{Float64}, CM.Parameters.MorrisonMilbrandt2014{Float64}},
                        Bool, Vector{Float64}, Vector{Float64}, Vector{Float64}, 
                        String, CM.Parameters.ParametersType{Float64}, Float64, Float64,
                        Vararg{String, 7}
                    }
                    }}(undef, length(data_file_name_list))
    IC_list = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)

    frozen_frac = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    frozen_frac_moving_mean = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    ICNC_moving_avg = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    Nₜ = Array{Float64}(undef, length(data_file_name_list), 1)
    y_truth = []

    for (exp_index, data_file_name) in enumerate(data_file_name_list)
        AIDA_data = grab_data(unpack_data(data_file_name))
        (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC, AIDA_e) = AIDA_data
        
        t_profile[exp_index] = AIDA_t_profile[start_time_index[exp_index]:end_time_index[exp_index]] .- (start_time_index[exp_index] - 101)
        T_profile[exp_index] = AIDA_T_profile[start_time_index[exp_index]:end_time_index[exp_index]]
        P_profile[exp_index] = AIDA_P_profile[start_time_index[exp_index]:end_time_index[exp_index]]
        ICNC_profile[exp_index] = AIDA_ICNC[start_time_index[exp_index]:end_time_index[exp_index]]
        e_profile[exp_index] = AIDA_e[start_time_index[exp_index]:end_time_index[exp_index]]
        S_l_profile[exp_index] = e_profile[exp_index] ./ (TD.saturation_vapor_pressure.(tps, T_profile[exp_index], TD.Liquid()))

        params_list[exp_index] =
            nuc_mode == "ABHOM" ?
            AIDA_IN05_params(FT, w[exp_index], t_max[exp_index], t_profile[exp_index], T_profile[exp_index], P_profile[exp_index]) :
            AIDA_IN07_params(FT, w[exp_index], t_max[exp_index], t_profile[exp_index], T_profile[exp_index], P_profile[exp_index], batch_name)
        IC_list[exp_index] = nuc_mode == "ABHOM" ?
            AIDA_IN05_IC(FT, data_file_name) :
            AIDA_IN07_IC(FT, data_file_name)

        Nₜ[exp_index] = IC_list[exp_index][7] + IC_list[exp_index][8] + IC_list[exp_index][9]
        frozen_frac[exp_index] = AIDA_ICNC[start_time_index[exp_index]:end_time_index[exp_index]] ./ Nₜ[exp_index]
        frozen_frac_moving_mean[exp_index] = moving_average(frozen_frac[exp_index], moving_average_n)
        ICNC_moving_avg[exp_index] = moving_average(AIDA_ICNC[start_time_index[exp_index]:end_time_index[exp_index]], moving_average_n)
        
        pseudo_ss_ff = NaNStatistics.nanmean(frozen_frac_moving_mean[exp_index][end - end_sim: end])
        pseudo_ss_ff_array = zeros(length(frozen_frac_moving_mean[exp_index][end - end_sim: end])) .+ pseudo_ss_ff
        append!(y_truth, pseudo_ss_ff_array)
    end

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
            exp_name = exp_index == 1 ? "IN0517" : "IN0518"
        else
            exp_name = batch_name
        end

        ## Calibrated parcel.
        EKI_parcel = run_model([params_list[exp_index]], nuc_mode, EKI_calibrated_parameters, FT, [IC_list[exp_index]], end_sim)
        UKI_parcel = run_model([params_list[exp_index]], nuc_mode, UKI_calibrated_parameters, FT, [IC_list[exp_index]], end_sim)

        ## Plots.
        ## Plotting AIDA data.
        AIDA_data = grab_data(unpack_data(data_file))
        (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC, AIDA_e) = AIDA_data

        AIDA_data_fig = MK.Figure(size = (1000, 600), fontsize = 24)
        data_ax1 = MK.Axis(AIDA_data_fig[1, 1], ylabel = "ICNC [m^-3]", xlabel = "time [s]", title = "AIDA data $exp_name")
        data_ax2 = MK.Axis(AIDA_data_fig[1, 2], ylabel = "Frozen Frac Moving Mean [-]", xlabel = "time [s]", title = "AIDA data $exp_name")
        MK.lines!(data_ax1, AIDA_t_profile, AIDA_ICNC, label = "Raw AIDA", color =:blue, linestyle =:dash, linewidth = 2)
        MK.lines!(data_ax1, t_profile[exp_index] .+ start_time_index[exp_index] .- 100, ICNC_moving_avg[exp_index], label = "AIDA ICNC moving mean", linewidth = 2.5, color =:blue)
        MK.lines!(data_ax2, t_profile[exp_index] .+ start_time_index[exp_index] .- 100, frozen_frac_moving_mean[exp_index], linewidth = 2.5, color =:blue)
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

        error = (AIDA_ICNC[start_time_index[exp_index]:end_time_index[exp_index]] ./ Nₜ[exp_index]) .* 0.1
        MK.errorbars!(ax_parcel_6, AIDA_t_profile[start_time_index[exp_index]:end_time_index[exp_index]] .- start_time[exp_index], AIDA_ICNC[start_time_index[exp_index]:end_time_index[exp_index]] ./ Nₜ[exp_index], error, color = (:blue, 0.3))

        MK.lines!(ax_parcel_7, EKI_parcel.t, EKI_parcel[2, :], color = :orange)
        MK.lines!(ax_parcel_7, UKI_parcel.t, UKI_parcel[2, :], color = :fuchsia)
        MK.lines!(ax_parcel_7, t_profile[exp_index], P_profile[exp_index], color = :blue) #, linestyle =:dash

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
