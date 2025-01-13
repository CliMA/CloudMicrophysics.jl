import EnsembleKalmanProcesses as EKP
import EnsembleKalmanProcesses.ParameterDistributions
import Random
import Distributions
import LinearAlgebra
import OrdinaryDiffEq as ODE

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import Thermodynamics as TD
import StatsBase as SB

#! format: off
# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration_setup.jl"))

# Define model which wraps around parcel and overwrites calibrated parameters
function run_model(p, IN_mode, coefficients, FT, IC, end_sim; calibration = false)
    # grabbing parameters
    m_calibrated, c_calibrated = coefficients
    (; prescribed_thermodynamics, t_profile, T_profile, P_profile) = p
    (; const_dt, w, t_max, ips) = p
    (; aerosol_act, aerosol, r_nuc) = p
    (; deposition_growth, condensation_growth) = p
    (; dep_nucleation, heterogeneous, homogeneous) = p
    (; liq_size_distribution, ice_size_distribution, aero_σ_g) = p

    # overwriting
    if aerosol == CMP.Kaolinite(FT)
        if IN_mode == "ABDINM"
            override_file = Dict(
            "China2017_J_deposition_m_Kaolinite" =>
                Dict("value" => m_calibrated, "type" => "float"),
            "China2017_J_deposition_c_Kaolinite" =>
                Dict("value" => c_calibrated, "type" => "float"),
            )
        elseif IN_mode == "ABIFM"
            override_file = Dict(
                "KnopfAlpert2013_J_ABIFM_m_Kaolinite" =>
                    Dict("value" => m_calibrated, "type" => "float"),
                "KnopfAlpert2013_J_ABIFM_c_Kaolinite" =>
                    Dict("value" => c_calibrated, "type" => "float"),
            )
        end
        kaolinite_calibrated = CP.create_toml_dict(FT; override_file)
        aerosol = CMP.Kaolinite(kaolinite_calibrated)
    elseif aerosol == CMP.ArizonaTestDust(FT)
        if IN_mode == "ABDINM"
            override_file = Dict(
                "J_ABDINM_m_ArizonaTestDust" =>
                    Dict("value" => m_calibrated, "type" => "float"),
                "J_ABDINM_c_ArizonaTestDust" =>
                    Dict("value" => c_calibrated, "type" => "float"),
            )
        elseif IN_mode == "ABIFM"
            override_file = Dict(
                "J_ABIFM_m_ArizonaTestDust" =>
                    Dict("value" => m_calibrated, "type" => "float"),
                "J_ABIFM_c_ArizonaTestDust" =>
                    Dict("value" => c_calibrated, "type" => "float"),
            )
        end
        ATD_calibrated = CP.create_toml_dict(FT; override_file)
        aerosol = CMP.ArizonaTestDust(ATD_calibrated)
    elseif aerosol == CMP.Illite(FT)
        if IN_mode == "ABDINM"
            override_file = Dict(
                "J_ABDINM_m_Illite" =>
                    Dict("value" => m_calibrated, "type" => "float"),
                "J_ABDINM_c_Illite" =>
                    Dict("value" => c_calibrated, "type" => "float"),
            )
        elseif IN_mode == "ABIFM"
            override_file = Dict(
                "KnopfAlpert2013_J_ABIFM_m_Illite" =>
                    Dict("value" => m_calibrated, "type" => "float"),
                "KnopfAlpert2013_J_ABIFM_c_Illite" =>
                    Dict("value" => c_calibrated, "type" => "float"),
            )
        end
        illite_calibrated = CP.create_toml_dict(FT; override_file)
        aerosol = CMP.Illite(illite_calibrated)
    elseif aerosol == CMP.Sulfate(FT)
        override_file = Dict(
            "Linear_J_hom_coeff2" =>
                Dict("value" => m_calibrated, "type" => "float"),
            "Linear_J_hom_coeff1" =>
                Dict("value" => c_calibrated, "type" => "float"),
        )
        ip_calibrated = CP.create_toml_dict(FT; override_file)
        ips = CMP.IceNucleationParameters(ip_calibrated)
    end

    # run parcel with new coefficients
    local params = parcel_params{FT}(
        const_dt = const_dt,
        w = w,
        prescribed_thermodynamics = prescribed_thermodynamics,
        t_profile = t_profile,
        T_profile = T_profile,
        P_profile = P_profile,
        aerosol_act = aerosol_act,
        aerosol = aerosol,
        aero_σ_g = aero_σ_g,
        homogeneous = homogeneous,
        deposition = dep_nucleation,
        heterogeneous = heterogeneous,
        condensation_growth = condensation_growth,
        deposition_growth = deposition_growth,
        liq_size_distribution = liq_size_distribution,
        ice_size_distribution = ice_size_distribution,
        ips = ips,
        r_nuc = r_nuc,
    )

    # solve ODE
    local sol = run_parcel(IC, FT(0), t_max, params)
    if calibration == true
        return sol[9, (end - end_sim):end] ./ (IC[7] + IC[8] + IC[9])
    elseif calibration == false
        return sol
    end
end

function create_prior(FT, IN_mode, ; perfect_model = false, aerosol_type = nothing)
    # TODO - add perfect_model flag to plot_ensemble_mean.jl
    observation_data_names = ["m_coeff", "c_coeff"]

    # stats = [mean, std dev, lower bound, upper bound]
    # for perfect model calibration
    if perfect_model == true
        if IN_mode == "ABDINM"
            m_stats = [FT(20), FT(8), FT(0), Inf]
            c_stats = [FT(-1), FT(2), -Inf, Inf]
        elseif IN_mode == "ABIFM"
            m_stats = [FT(50), FT(7), FT(0), Inf]
            c_stats = [FT(-7), FT(4), -Inf, Inf]
        elseif IN_mode == "ABHOM"
            m_stats = [FT(251), FT(6), FT(0), Inf]
            c_stats = [FT(-66), FT(3), -Inf, Inf]
        end
    elseif perfect_model == false
        if IN_mode == "ABDINM"
            if aerosol_type == CMP.ArizonaTestDust(FT)
                m_stats = [FT(40), FT(20), FT(0), Inf]
                c_stats = [FT(0.5), FT(20), -Inf, Inf]
            elseif aerosol_type == CMP.Illite(FT)
                m_stats = [FT(30), FT(20), FT(0), Inf]
                c_stats = [FT(0.7), FT(7), -Inf, Inf]
            end
        elseif IN_mode == "ABIFM"
            # m_stats = [FT(50), FT(1), FT(0), Inf]
            # c_stats = [FT(-7), FT(1), -Inf, Inf]
            error("ABIFM calibrations not yet supported.")
        elseif IN_mode == "ABHOM"
            m_stats = [FT(260.927125), FT(25), FT(0), Inf]
            c_stats = [FT(-68.553283), FT(10), -Inf, Inf]
        end
    end

    m_prior = EKP.constrained_gaussian(
        observation_data_names[1],
        m_stats[1],
        m_stats[2],
        m_stats[3],
        m_stats[4],
    )
    c_prior = EKP.constrained_gaussian(
        observation_data_names[2],
        c_stats[1],
        c_stats[2],
        c_stats[3],
        c_stats[4],
    )
    prior = EKP.combine_distributions([m_prior, c_prior])
    return prior
end

function calibrate_J_parameters_EKI(FT, IN_mode, params, IC, y_truth, end_sim, Γ,; perfect_model = false)
    @info("Starting EKI calibration")
    # Random number generator
    rng_seed = 24
    rng = Random.seed!(Random.GLOBAL_RNG, rng_seed)

    (; aerosol) = params

    prior = create_prior(FT, IN_mode, perfect_model = perfect_model, aerosol_type = aerosol)
    N_ensemble = 25       # runs N_ensemble trials per iteration
    N_iterations = 50     # number of iterations the inverse problem goes through

    # Generate initial ensemble and set up EKI
    initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ensemble)
    EKI_obj = EKP.EnsembleKalmanProcess(
        initial_ensemble,
        y_truth[(end - end_sim):end],
        Γ,
        EKP.Inversion();
        rng = rng,
        verbose = true,
        localization_method = EKP.Localizers.NoLocalization(), # no localization
    )

    # Carry out the EKI calibration
    # ϕ_n_values[iteration] stores ensembles of calibrated coeffs in that iteration
    global ϕ_n_values = []
    global final_iter = N_iterations
    for n in 1:N_iterations
        ϕ_n = EKP.get_ϕ_final(prior, EKI_obj)
        G_ens = hcat(
            [
                run_model(params, IN_mode, ϕ_n[:, i], FT, IC, end_sim, calibration = true) for
                i in 1:N_ensemble
            ]...,
        )
        # Update ensemble
        terminated = EKP.update_ensemble!(EKI_obj, G_ens)
        # if termination flagged, can stop at earlier iteration
        if !isnothing(terminated)
            final_iter = n - 1
            break
        end

        ϕ_n_values = vcat(ϕ_n_values, [ϕ_n])
    end

    # Mean coefficients of all ensembles in the final iteration
    m_coeff_ekp = round(
        Distributions.mean(ϕ_n_values[final_iter][1, 1:N_ensemble]),
        digits = 6,
    )
    c_coeff_ekp = round(
        Distributions.mean(ϕ_n_values[final_iter][2, 1:N_ensemble]),
        digits = 6,
    )

    calibrated_coeffs = [m_coeff_ekp, c_coeff_ekp]

    return [calibrated_coeffs, ϕ_n_values]
end

function calibrate_J_parameters_UKI(FT, IN_mode, params, IC, y_truth, end_sim, Γ,; perfect_model = false)
    @info("Starting UKI calibration")
    (; aerosol) = params

    prior = create_prior(FT, IN_mode, perfect_model = perfect_model, aerosol_type = aerosol)
    N_iterations = 25
    α_reg = 1.0
    update_freq = 1

    # truth = EKP.Observations.Observation(y_truth, Γ, "y_truth")
    truth = EKP.Observation(
        Dict("samples" => vec(SB.mean(y_truth[end - end_sim: end], dims = 2)), "covariances" => Γ, "names" => "y_truth")
    )

    # Generate initial ensemble and set up UKI
    process = EKP.Unscented(
        SB.mean(prior),
        SB.cov(prior);
        α_reg = α_reg,
        update_freq = update_freq,
        impose_prior = false,
    )
    UKI_obj = EKP.EnsembleKalmanProcess(truth, process; verbose = true)

    global err = []
    global final_iter =[N_iterations]
    for n in 1:N_iterations
        # Return transformed parameters in physical/constrained space
        ϕ_n = EKP.get_ϕ_final(prior, UKI_obj)
        # Evaluate forward map
        G_n = [
            run_model(params, IN_mode, ϕ_n[:, i], FT, IC, end_sim, calibration = true) for
            i in 1:size(ϕ_n)[2]  #i in 1:N_ensemble
        ]
        # Reformat into `d x N_ens` matrix
        G_ens = hcat(G_n...)
        # Update ensemble
        terminate = EKP.EnsembleKalmanProcesses.update_ensemble!(UKI_obj, G_ens)
        push!(err, EKP.get_error(UKI_obj)[end])
        # println(
        #     "Iteration: " *
        #     string(n) *
        #     ", Error: " *
        #     string(err[n]) *
        #     " norm(Cov):" *
        #     string(Distributions.norm(EKP.get_process(UKI_obj).uu_cov[n]))
        # )
        if !isnothing(terminate)
            final_iter[1] = n - 1
            break
        end
    end

    UKI_mean_u_space = EKP.get_u_mean_final(UKI_obj)
    UKI_mean = EKP.transform_unconstrained_to_constrained(prior, UKI_mean_u_space)

    ϕ_n = EKP.get_ϕ_final(prior, UKI_obj)

    return [UKI_mean, ϕ_n, final_iter]
end

function ensemble_means(ϕ_n_values, N_iterations, N_ensemble)
    iterations = collect(1:N_iterations)
    m_mean = zeros(length(iterations))
    c_mean = zeros(length(iterations))

    for iter in iterations
        m_mean[iter] =
            Distributions.mean(ϕ_n_values[iter][1, i] for i in 1:N_ensemble)
        c_mean[iter] =
            Distributions.mean(ϕ_n_values[iter][2, i] for i in 1:N_ensemble)
    end

    return [m_mean, c_mean]
end
