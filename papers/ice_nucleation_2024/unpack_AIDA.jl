import Thermodynamics as TD
import DelimitedFiles
using LazyArtifacts
using ClimaUtilities.ClimaArtifacts
import Interpolations as Intp

include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration_setup.jl"))
FT = Float64

function interpolated_data(raw_data, time_array)
    data_at_t = Intp.linear_interpolation(raw_data[:, 1], raw_data[:, 2])
    return data_at_t.(time_array)
end

function moving_average(data, n)
    window_size = length(data) / n
    moving_avg = NaNStatistics.movmean(data, window_size)
    return moving_avg
end

function unpack_data(data_file_name, ; total_t = 0)
    edf_data_names = [
        "in05_17_aida.edf", "in05_18_aida.edf",
        "in07_01_aida.edf", "in07_19_aida.edf",
    ]

    if data_file_name in edf_data_names

        file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name)
        unpacked_data = DelimitedFiles.readdlm(file_path, skipstart = 225)

        AIDA_t_profile = unpacked_data[:, 1]
        AIDA_T_profile = unpacked_data[:, 3]
        AIDA_P_profile = unpacked_data[:, 2] .* 1e2          # hPa to Pa
        AIDA_ICNC_profile = unpacked_data[:, 6] .* 1e6       # Nᵢ [m^-3]
        AIDA_e_profile = unpacked_data[:, 4]

    elseif data_file_name == "TROPIC04"
        AIDA_t_profile = collect(0:1:total_t)

        AIDA_T_profile = zeros(total_t)
        AIDA_P_profile = zeros(total_t)
        AIDA_ICNC_profile = zeros(total_t)
        AIDA_e_profile = zeros(total_t)

        T_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_T.csv")
        P_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_P.csv")
        act_frac_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_act_frac.csv")
        S_ice_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_S_ice.csv")

        raw_T_data = DelimitedFiles.readdlm(T_file_path, ',', header = false)
        raw_P_data = DelimitedFiles.readdlm(P_file_path, ',', header = false)  # already in Pa
        raw_act_frac_data = DelimitedFiles.readdlm(act_frac_file_path, ',', header = false)
        raw_S_ice_data = DelimitedFiles.readdlm(S_ice_file_path, ',', header = false)

        AIDA_T_profile = interpolated_data(raw_T_data, AIDA_t_profile)
        AIDA_P_profile = interpolated_data(raw_P_data, AIDA_t_profile)
        AIDA_act_frac_profile = interpolated_data(raw_act_frac_data, AIDA_t_profile)
        AIDA_S_ice_profile = interpolated_data(raw_S_ice_data, AIDA_t_profile)

        # converting activated frac and S_ice to ICNC and e respectively
        N_total = 66.19233e6  # m^-3; should match initial conditions; Nₐ + Nₗ + Nᵢ
        AIDA_ICNC_profile = AIDA_act_frac_profile .* N_total

        eₛ = [TD.saturation_vapor_pressure(tps, T₀, TD.Ice()) for T₀ in AIDA_T_profile]
        AIDA_e_profile = AIDA_S_ice_profile .* eₛ

    else # all other .csv files

        AIDA_t_profile = collect(0:1:total_t)

        AIDA_T_profile = zeros(total_t)
        AIDA_P_profile = zeros(total_t)
        AIDA_ICNC_profile = zeros(total_t)
        AIDA_e_profile = zeros(total_t)

        T_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_T.csv")
        P_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_P.csv")
        ICNC_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_N_ice.csv")
        RH_file_path = CM.ArtifactCalling.AIDA_ice_nucleation(data_file_name * "_RH_water.csv")

        raw_T_data = DelimitedFiles.readdlm(T_file_path, ',', header = false)
        raw_P_data = DelimitedFiles.readdlm(P_file_path, ',', header = false)
        raw_ICNC_data = DelimitedFiles.readdlm(ICNC_file_path, ',', header = false)
        raw_RH_data = DelimitedFiles.readdlm(RH_file_path, ',', header = false)

        AIDA_T_profile = interpolated_data(raw_T_data, AIDA_t_profile)
        AIDA_P_profile = interpolated_data(raw_P_data, AIDA_t_profile) .* 100
        AIDA_ICNC_profile = interpolated_data(raw_ICNC_data, AIDA_t_profile) .* 1e6
        AIDA_RH_profile =
            data_file_name == "ACI04_22" ?
            interpolated_data(raw_RH_data, AIDA_t_profile) ./ 100 :
            interpolated_data(raw_RH_data, AIDA_t_profile)

        eₛ = [TD.saturation_vapor_pressure(tps, T, TD.Liquid()) for T in AIDA_T_profile]
        AIDA_e_profile = AIDA_RH_profile .* eₛ

    end

    data = (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC_profile, AIDA_e_profile)
    return data

end

function data_to_calib_inputs(
    FT, batch_name, data_file_name_list, nuc_mode,
    batch_start_times, batch_end_times,
    w, t_max, end_sim,
)
    edf_data_names = [
        "in05_17_aida.edf", "in05_18_aida.edf",
        "in07_01_aida.edf", "in07_19_aida.edf",
    ]

    # pre-defined variable arrays
    global t_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    global T_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    global P_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    global ICNC_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    global e_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    global S_l_profile = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)

    params_list = Vector{
        NamedTuple{
            (:const_dt, :w, :t_max, :ips,
                :prescribed_thermodynamics, :t_profile, :T_profile, :P_profile,
                :aerosol_act, :aerosol, :r_nuc, :aero_σ_g, :A_aer,
                :condensation_growth, :deposition_growth,
                :liq_size_distribution, :ice_size_distribution,
                :dep_nucleation, :heterogeneous, :homogeneous,
            ),
            Tuple{
                Float64, Float64, Int32,
                CM.Parameters.IceNucleationParameters{
                    Float64,
                    CM.Parameters.Mohler2006{Float64},
                    CM.Parameters.Koop2000{Float64},
                    CM.Parameters.MorrisonMilbrandt2014{Float64},
                },
                Bool, Vector{Float64}, Vector{Float64}, Vector{Float64},
                String, CM.Parameters.ParametersType{Float64}, Float64, Float64, Float64,
                Vararg{String, 7},
            },
        },
    }(
        undef,
        length(data_file_name_list),
    )
    IC_list = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)

    frozen_frac = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    frozen_frac_moving_mean = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    ICNC_moving_avg = Array{Vector{Float64}}(undef, length(data_file_name_list), 1)
    Nₜ = Array{Float64}(undef, length(data_file_name_list), 1)
    y_truth = []

    # filling in empty variable arrays
    for (exp_index, data_file_name) in enumerate(data_file_name_list)
        # Takes cropped version of data in the window of time that aligns with our simulation
        # and creates y_truth to be used for calibration
        if data_file_name in edf_data_names
            AIDA_data = unpack_data(data_file_name)
            (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC_profile, AIDA_e_profile) = AIDA_data

            t_profile[exp_index] =
                AIDA_t_profile[batch_start_times[exp_index]:batch_end_times[exp_index]] .- batch_start_times[exp_index]
            T_profile[exp_index] = AIDA_T_profile[batch_start_times[exp_index]:batch_end_times[exp_index]]
            P_profile[exp_index] = AIDA_P_profile[batch_start_times[exp_index]:batch_end_times[exp_index]]
            ICNC_profile[exp_index] = AIDA_ICNC_profile[batch_start_times[exp_index]:batch_end_times[exp_index]]
            e_profile[exp_index] = AIDA_e_profile[batch_start_times[exp_index]:batch_end_times[exp_index]]
            S_l_profile[exp_index] =
                e_profile[exp_index] ./ (TD.saturation_vapor_pressure.(tps, T_profile[exp_index], TD.Liquid()))

            params_list[exp_index] =
                nuc_mode == "ABHOM" ?
                AIDA_IN05_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                ) :
                AIDA_IN07_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                    batch_name,
                )
            IC_list[exp_index] =
                nuc_mode == "ABHOM" ?
                AIDA_IN05_IC(FT, data_file_name) :
                AIDA_IN07_IC(FT, data_file_name)

        else
            AIDA_data = unpack_data(data_file_name, total_t = t_max[exp_index])
            (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC_profile, AIDA_e_profile) = AIDA_data
            t_profile[exp_index] = AIDA_t_profile
            T_profile[exp_index] = AIDA_T_profile
            P_profile[exp_index] = AIDA_P_profile
            ICNC_profile[exp_index] = AIDA_ICNC_profile
            e_profile[exp_index] = AIDA_e_profile

            S_l_profile[exp_index] =
                e_profile[exp_index] ./ (TD.saturation_vapor_pressure.(tps, T_profile[exp_index], TD.Liquid()))
            if data_file_name == "TROPIC04"
                params_list[exp_index] = TROPIC04_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                )
                IC_list[exp_index] = TROPIC04_IC(FT)
            elseif data_file_name == "ACI04_22"
                params_list[exp_index] = ACI04_22_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                    batch_name,
                )
                IC_list[exp_index] = ACI04_22_IC(FT)
            elseif data_file_name == "EXP19"
                params_list[exp_index] = EXP19_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                    batch_name,
                )
                IC_list[exp_index] = EXP19_IC(FT)
            elseif data_file_name == "EXP45"
                params_list[exp_index] = EXP45_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                    batch_name,
                )
                IC_list[exp_index] = EXP45_IC(FT)
            elseif data_file_name == "EXP28"
                params_list[exp_index] = EXP45_params(
                    FT,
                    w[exp_index],
                    t_max[exp_index],
                    t_profile[exp_index],
                    T_profile[exp_index],
                    P_profile[exp_index],
                    batch_name,
                )
                IC_list[exp_index] = EXP45_IC(FT)
            end
        end

        Nₜ[exp_index] = IC_list[exp_index][7] + IC_list[exp_index][8] + IC_list[exp_index][9]
        frozen_frac[exp_index] = ICNC_profile[exp_index] ./ Nₜ[exp_index]
        frozen_frac_moving_mean[exp_index] = moving_average(frozen_frac[exp_index], moving_average_n)
        ICNC_moving_avg[exp_index] = moving_average(ICNC_profile[exp_index], moving_average_n)

        pseudo_ss_ff = NaNStatistics.nanmean(frozen_frac_moving_mean[exp_index][(end - end_sim):end])
        pseudo_ss_ff_array = zeros(length(frozen_frac_moving_mean[exp_index][(end - end_sim):end])) .+ pseudo_ss_ff
        append!(y_truth, pseudo_ss_ff_array)

    end


    calib_variables = (;
        t_profile, T_profile, P_profile, ICNC_profile, e_profile, S_l_profile,
        params_list, IC_list,
        frozen_frac, frozen_frac_moving_mean, ICNC_moving_avg,
        Nₜ, y_truth,
    )
    return calib_variables

end
