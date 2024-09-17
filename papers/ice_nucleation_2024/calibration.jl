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

#! format: off
# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))
include(joinpath(pkgdir(CM), "papers", "ice_nucleation_2024", "calibration_setup.jl"))

# Define model which wraps around parcel and overwrites calibrated parameters
function run_model(p, coefficients, IN_mode, FT, IC)
    # grabbing parameters
    m_calibrated, c_calibrated = coefficients
    (; const_dt, w, t_max, aerosol_act, aerosol, r_nuc) = p
    (; deposition_growth, condensation_growth) = p
    (; dep_nucleation, heterogeneous, homogeneous) = p
    (; liq_size_distribution, ice_size_distribution, aero_σ_g) = p

    if IN_mode == "ABDINM"
        # overwriting
        override_file = Dict(
            "China2017_J_deposition_m_Kaolinite" =>
                Dict("value" => m_calibrated, "type" => "float"),
            "China2017_J_deposition_c_Kaolinite" =>
                Dict("value" => c_calibrated, "type" => "float"),
        )
        kaolinite_calibrated = CP.create_toml_dict(FT; override_file)
        overwrite = CMP.Kaolinite(kaolinite_calibrated)

        # run parcel with new coefficients
        local params = parcel_params{FT}(
            const_dt = const_dt,
            w = w,
            aerosol = overwrite,
            deposition = dep_nucleation,
            condensation_growth = condensation_growth,
            deposition_growth = deposition_growth,
            liq_size_distribution = liq_size_distribution,
            ice_size_distribution = ice_size_distribution,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol[9, :] .* 1e10 # ICNC magnified

    elseif IN_mode == "ABIFM"
        # overwriting
        override_file = Dict(
            "KnopfAlpert2013_J_ABIFM_m_Kaolinite" =>
                Dict("value" => m_calibrated, "type" => "float"),
            "KnopfAlpert2013_J_ABIFM_c_Kaolinite" =>
                Dict("value" => c_calibrated, "type" => "float"),
        )
        kaolinite_calibrated = CP.create_toml_dict(FT; override_file)
        overwrite = CMP.Kaolinite(kaolinite_calibrated)

        # run parcel with new coefficients
        local params = parcel_params{FT}(
            const_dt = const_dt,
            w = w,
            aerosol = overwrite,
            heterogeneous = heterogeneous,
            condensation_growth = condensation_growth,
            deposition_growth = deposition_growth,
            liq_size_distribution = liq_size_distribution,
            ice_size_distribution = ice_size_distribution,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol[9, :] # ICNC

    elseif IN_mode == "ABHOM"
        # overwriting
        override_file = Dict(
            "Linear_J_hom_coeff2" =>
                Dict("value" => m_calibrated, "type" => "float"),
            "Linear_J_hom_coeff1" =>
                Dict("value" => c_calibrated, "type" => "float"),
        )
        ip_calibrated = CP.create_toml_dict(FT; override_file)
        overwrite = CMP.IceNucleationParameters(ip_calibrated)

        # run parcel with new coefficients
        local params = parcel_params{FT}(
            const_dt = const_dt,
            w = w,
            aerosol_act = aerosol_act,
            aerosol = aerosol,
            aero_σ_g = aero_σ_g,
            homogeneous = homogeneous,
            condensation_growth = condensation_growth,
            deposition_growth = deposition_growth,
            liq_size_distribution = liq_size_distribution,
            ice_size_distribution = ice_size_distribution,
            ips = overwrite,
            r_nuc = r_nuc,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol[9, :] ./  (IC[7] + IC[8] + IC[9])  # frozen fraction
    end
end

function run_calibrated_model(FT, IN_mode, coefficients, p, IC)
    # grabbing parameters
    m_calibrated, c_calibrated = coefficients
    (; const_dt, w, t_max, aerosol_act, aerosol) = p
    (; deposition_growth, condensation_growth) = p
    (; homogeneous) = p
    (; liq_size_distribution, ice_size_distribution, aero_σ_g) = p
    
    if IN_mode == "ABHOM"
        # overwriting
        override_file = Dict(
            "Linear_J_hom_coeff2" =>
                Dict("value" => m_calibrated, "type" => "float"),
            "Linear_J_hom_coeff1" =>
                Dict("value" => c_calibrated, "type" => "float"),
        )
        ip_calibrated = CP.create_toml_dict(FT; override_file)
        overwrite = CMP.IceNucleationParameters(ip_calibrated)

        # run parcel with new coefficients
        local params = parcel_params{FT}(
            const_dt = const_dt,
            w = w,
            aerosol_act = aerosol_act,
            aerosol = aerosol,
            aero_σ_g = aero_σ_g,
            homogeneous = homogeneous,
            condensation_growth = condensation_growth,
            deposition_growth = deposition_growth,
            liq_size_distribution = liq_size_distribution,
            ice_size_distribution = ice_size_distribution,
            ips = overwrite,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol
    end
end

function create_prior(FT, IN_mode, ; perfect_model = false)
    # TODO - add perfect_model flag to plot_ensemble_mean.jl
    observation_data_names = ["m_coeff", "c_coeff"]

    # stats = [mean, std dev, lower bound, upper bound]
    # for perfect model calibration
    if perfect_model == true
        if IN_mode == "ABDINM"
            m_stats = [FT(20), FT(1), FT(0), Inf]
            c_stats = [FT(-1), FT(1), -Inf, Inf]
        elseif IN_mode == "ABIFM"
            m_stats = [FT(50), FT(1), FT(0), Inf]
            c_stats = [FT(-7), FT(1), -Inf, Inf]
        elseif IN_mode == "ABHOM"
            m_stats = [FT(260), FT(1), FT(0), Inf]
            c_stats = [FT(-70), FT(1), -Inf, Inf]
        end
    elseif perfect_model == false
        if IN_mode == "ABDINM"
            # m_stats = [FT(20), FT(1), FT(0), Inf]
            # c_stats = [FT(-1), FT(1), -Inf, Inf]
            println("Calibration for ABDINM with AIDA not yet implemented.")
        elseif IN_mode == "ABIFM"
            # m_stats = [FT(50), FT(1), FT(0), Inf]
            # c_stats = [FT(-7), FT(1), -Inf, Inf]
            println("Calibration for ABIFM with AIDA not yet implemented.")
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

function calibrate_J_parameters(FT, IN_mode, params, IC, y_truth, Γ,; perfect_model = false)
    # Random number generator
    rng_seed = 24
    rng = Random.seed!(Random.GLOBAL_RNG, rng_seed)

    prior = create_prior(FT, IN_mode, perfect_model = perfect_model)
    N_ensemble = 10      # runs N_ensemble trials per iteration
    N_iterations = 150    # number of iterations the inverse problem goes through

    # Generate initial ensemble and set up EKI
    initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ensemble)
    EKI_obj = EKP.EnsembleKalmanProcess(
        initial_ensemble,
        y_truth,
        Γ,
        EKP.Inversion();
        rng = rng,
    )

    # Carry out the EKI calibration
    # ϕ_n_values[iteration] stores ensembles of calibrated coeffs in that iteration
    global ϕ_n_values = []
    for n in 1:N_iterations
        ϕ_n = EKP.get_ϕ_final(prior, EKI_obj)
        G_ens = hcat(
            [
                run_model(params, ϕ_n[:, i], IN_mode, FT, IC) for
                i in 1:N_ensemble
            ]...,
        )
        EKP.update_ensemble!(EKI_obj, G_ens)

        ϕ_n_values = vcat(ϕ_n_values, [ϕ_n])
    end

    # Mean coefficients of all ensembles in the final iteration
    m_coeff_ekp = round(
        Distributions.mean(ϕ_n_values[N_iterations][1, 1:N_ensemble]),
        digits = 6,
    )
    c_coeff_ekp = round(
        Distributions.mean(ϕ_n_values[N_iterations][2, 1:N_ensemble]),
        digits = 6,
    )

    return [m_coeff_ekp, c_coeff_ekp, ϕ_n_values]
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
