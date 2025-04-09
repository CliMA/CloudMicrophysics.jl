import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD
import CairoMakie as MK
import OrdinaryDiffEq as ODE

# To grab data
import DelimitedFiles
using LazyArtifacts
using ClimaUtilities.ClimaArtifacts
import NaNStatistics

FT = Float64
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "unpack_AIDA.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "plots.jl"))

# Defining data names, start/end times, etc.
global edf_data_names = [
    "in05_17_aida.edf", "in05_18_aida.edf",
    "in07_01_aida.edf", "in07_19_aida.edf",
]
data_file_names = [
    ["in05_17_aida.edf", "in05_18_aida.edf"],
    ["in07_01_aida.edf"],
    ["in07_19_aida.edf"],
    ["ACI04_22"],
    ["EXP19"],
    ["EXP45"],
    ["in07_01_aida.edf", "in07_19_aida.edf", "EXP45"],
    ["ACI04_22", "EXP19"],
]
batch_names = ["HOM", "IN0701", "IN0719", "ACI04_22", "EXP19", "EXP45", "DEP", "IMM"]
end_sim = 25            # Loss func looks at last end_sim timesteps only
start_time_list = [     # freezing onset
    [Int32(150), Int32(180)],
    [Int32(50)],
    [Int32(35)],
    [Int32(0)],
    [Int32(0)],
    [Int32(0)],
    [Int32(50), Int32(35), Int32(0)],
    [Int32(0), Int32(0)],
]
end_time_list = [       # approximate time freezing stops
    [Int32(295), Int32(290)],
    [Int32(375)],
    [Int32(375)],
    [Int32(300)],
    [Int32(200)],
    [Int32(200)],
    [Int32(375), Int32(375), Int32(200)],
    [Int32(300), Int32(200)],
]
moving_average_n = 5    # average every length(data) / n points
updrafts = [            # updrafts matching AIDA cooling rate
    [FT(1.5), FT(1.4)],
    [FT(1.5)],
    [FT(1.5)],
    [FT(1.5)],
    [FT(9)],
    [FT(1.5)],
    [FT(1.5), FT(1.5), FT(1.5)],
    [FT(1.5), FT(9)],
]

# Supercooled water thermodynamic properties
override_file = Dict( 
    "pressure_triple_point" =>
        Dict("value" => FT(25.16299251 * 1.0047), "type" => "float"),
    "temperature_triple_point" =>
        Dict("value" => FT(236), "type" => "float"),
    "thermodynamics_temperature_reference" =>
        Dict("value" => FT(236), "type" => "float"),
    "latent_heat_vaporization_at_reference" =>
        Dict("value" => FT(2600307.36), "type" => "float"),
    "latent_heat_sublimation_at_reference" =>
        Dict("value" => FT(2841004.893), "type" => "float"),
    "isobaric_specific_heat_vapor" =>
        Dict("value" => FT(1853.675619), "type" => "float"),
    "isobaric_specific_heat_liquid" =>
        Dict("value" => FT(5396.1059), "type" => "float"),
    "isobaric_specific_heat_ice" =>
        Dict("value" => FT(1823.342958), "type" => "float"),
)
toml_dict = CP.create_toml_dict(FT; override_file)
tps = TD.Parameters.ThermodynamicsParameters(toml_dict)

# Additional definitions
R_v = TD.Parameters.R_v(tps)
R_d = TD.Parameters.R_d(tps)
ϵₘ = R_d / R_v

global EKI_calibrated_coeff_dict = Dict()
global UKI_calibrated_coeff_dict = Dict()
mutable struct all_calibration_data
    EKI_calibrated_parcel::Vector{ODE.ODESolution}
    UKI_calibrated_parcel::Vector{ODE.ODESolution}
    Nₜ_list::Vector{Any}
    t_profile_list::Vector{Any}
    frozen_frac_moving_mean_list::Vector{Any}
    frozen_frac_list::Vector{Any}
end
global overview_data = all_calibration_data(Vector{ODE.ODESolution}(), Vector{ODE.ODESolution}(), [], [], [], [])

for (batch_index, batch_name) in enumerate(batch_names)
    @info(batch_name)
    #! format: off
    ### Unpacking experiment-specific variables.
    data_file_name_list = data_file_names[batch_index]
    w = updrafts[batch_index]
    start_time = start_time_list[batch_index]
    end_time = end_time_list[batch_index]
    t_max = end_time .- start_time

    if batch_name == "HOM"
        nuc_mode = "ABHOM"
    elseif batch_name in ["IN0701", "IN0719", "EXP45", "DEP"]
        nuc_mode = "ABDINM"
    elseif batch_name in ["ACI04_22", "EXP19", "IMM"]
        nuc_mode = "ABIFM"
    end

    ### Check for and grab data in AIDA_data folder.
    calib_variables = data_to_calib_inputs(
        FT, tps,
        batch_name, data_file_name_list, nuc_mode,
        start_time, end_time,
        w, t_max, end_sim,
    )
    (;
        t_profile, T_profile, P_profile, ICNC_profile, e_profile, S_l_profile,
        params_list, IC_list,
        frozen_frac, frozen_frac_moving_mean, ICNC_moving_avg,
        Nₜ, y_truth,
    ) = calib_variables

    if batch_name in ["HOM", "DEP", "IMM"]
        append!(overview_data.Nₜ_list, Nₜ)
        append!(overview_data.t_profile_list, t_profile)
        append!(overview_data.frozen_frac_moving_mean_list, frozen_frac_moving_mean)
        append!(overview_data.frozen_frac_list, frozen_frac)
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
    merge!(EKI_calibrated_coeff_dict, Dict(batch_name => EKI_calibrated_parameters))
    merge!(UKI_calibrated_coeff_dict, Dict(batch_name => UKI_calibrated_parameters))

    ### Plot for individual experiments.
    for (exp_index, data_file) in enumerate(data_file_name_list)
        if batch_name == "HOM"
            exp_names = ["HOM_Batch_IN0517", "HOM_Batch_IN0518", "HOM_Batch_TROPIC04"]
        elseif batch_name == "DEP"
            exp_names = ["DEP_Batch_IN0701", "DEP_Batch_IN0719", "DEP_Batch_EXP45"]
        elseif batch_name == "IMM"
            exp_names = ["IMM_Batch_ACI04_22", "IMM_Batch_EXP19"]
        else
            exp_names = [batch_name]
        end

        ## Calibrated parcel.
        EKI_parcel = run_model([params_list[exp_index]], nuc_mode, EKI_calibrated_parameters, FT, [IC_list[exp_index]], end_sim)
        UKI_parcel = run_model([params_list[exp_index]], nuc_mode, UKI_calibrated_parameters, FT, [IC_list[exp_index]], end_sim)
        if batch_name in ["HOM", "DEP", "IMM"]
            # the order it appends is "in05_17_aida.edf", "in05_18_aida.edf", "TROPIC04",
            # "in07_01_aida.edf", "in07_19_aida.edf", "EXP45", "ACI04_22", "EXP19".
            push!(overview_data.EKI_calibrated_parcel, EKI_parcel)
            push!(overview_data.UKI_calibrated_parcel, UKI_parcel)
        end

        ## Plots.
        ## Plotting AIDA data.
        AIDA_data = data_file in edf_data_names ? unpack_data(data_file) : unpack_data(data_file, total_t = t_max[exp_index])
        (; AIDA_t_profile, AIDA_T_profile, AIDA_P_profile, AIDA_ICNC_profile, AIDA_e_profile) = AIDA_data

        AIDA_ICNC_data_fig = plot_AIDA_ICNC_data(
            exp_names[exp_index],
            AIDA_t_profile, AIDA_ICNC_profile,
            t_profile[exp_index], ICNC_moving_avg[exp_index], frozen_frac_moving_mean[exp_index],
        )

        ## Calibrated coefficients.
        #  Did they converge?
        calibrated_coeffs_fig = plot_calibrated_coeffs(batch_name, EKI_n_iterations, calibrated_ensemble_means)

        ## Calibrated parcel simulations.
        #  Does the calibrated parcel give reasonable outputs?
        calibrated_parcels_fig = plot_calibrated_parcels(
            exp_names[exp_index], Nₜ[exp_index],
            EKI_parcel, UKI_parcel,
            t_profile[exp_index], T_profile[exp_index], P_profile[exp_index], ICNC_profile[exp_index],
        )

        ## Comparing AIDA data and calibrated parcel.
        #  Does calibrated parcel look like observations?
        compare_ICNC_fig = plot_compare_ICNC(
            exp_names[exp_index],
            Nₜ[exp_index], EKI_parcel, UKI_parcel,
            t_profile[exp_index], frozen_frac_moving_mean[exp_index], frozen_frac[exp_index],
        )

        ## Looking at spread in UKI calibrated parameters
        ϕ_UKI = UKI_output[2]
        UKI_parcel_1 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_2 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_3 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_4 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)
        UKI_parcel_5 = run_model([params_list[exp_index]], nuc_mode, [ϕ_UKI[1,1], ϕ_UKI[2,1]], FT, [IC_list[exp_index]], end_sim)

        UKI_spread_fig = plot_UKI_spread(
            exp_names[exp_index],
            t_profile[exp_index], frozen_frac_moving_mean[exp_index], Nₜ[exp_index],
            UKI_parcel_1, UKI_parcel_2, UKI_parcel_3, UKI_parcel_4, UKI_parcel_5, UKI_parcel,
        )

    end # plotting individual experiments
end

# Plotting overview
plot_ICNC_overview(overview_data)

#! format: on
