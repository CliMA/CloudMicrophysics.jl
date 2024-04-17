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

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

# Define model which wraps around parcel and overwrites calibrated parameters
function run_model(p, coefficients, IN_mode, FT, IC)
    # grabbing parameters
    m_calibrated, c_calibrated = coefficients
    (; const_dt, w, deposition_growth, size_distribution) = p

    t_max = FT(100)

    if IN_mode == "ABDINM"
        (; dep_nucleation) = p

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
            deposition_growth = deposition_growth,
            size_distribution = size_distribution,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol[9, :] .* 1e10 # ICNC magnified

    elseif IN_mode == "ABIFM"
        (; heterogeneous, condensation_growth) = p

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
            size_distribution = size_distribution,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol[9, :] # ICNC

    elseif IN_mode == "ABHOM"
        (; homogeneous) = p

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
            homogeneous = homogeneous,
            deposition_growth = deposition_growth,
            size_distribution = size_distribution,
            ips = overwrite,
        )

        # solve ODE
        local sol = run_parcel(IC, FT(0), t_max, params)
        return sol[9, :] # ICNC
    end
end

# Creating noisy pseudo-observations
function create_prior(FT, IN_mode)
    observation_data_names = ["m_coeff", "c_coeff"]

    # Define prior distributions for the coefficients
    # stats = [mean, std dev, lower bound, upper bound]
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

function calibrate_J_parameters(FT, IN_mode)
    # Random number generator
    rng_seed = 24
    rng = Random.seed!(Random.GLOBAL_RNG, rng_seed)

    # Get free parameters
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    wps = CMP.WaterProperties(FT)

    # Define other parameters and initial conditions
    ρₗ = wps.ρw
    R_d = TD.Parameters.R_d(tps)
    R_v = TD.Parameters.R_v(tps)
    ϵₘ = R_d / R_v

    const_dt = FT(1)

    # Create prior
    prior = create_prior(FT, IN_mode)

    if IN_mode == "ABDINM"
        Nₐ = FT(8e6) # according to fig 2 panel b deposition nucleation experiment
        Nₗ = FT(0)
        Nᵢ = FT(0)
        r₀ = FT(0.5e-6)
        p₀ = FT(987.018 * 1e2)
        T₀ = FT(212.978)

        C_v = FT(10.8509 * 1e-6)
        qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
        qₗ = FT(0)
        qᵢ = FT(0)

        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
        w = FT(0.4)
        IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]

        dep_nucleation = "ABDINM"
        deposition_growth = "Deposition"
        size_distribution = "Monodisperse"

        params = (;
            const_dt,
            w,
            dep_nucleation,
            deposition_growth,
            size_distribution,
        )

        coeff_true = [FT(27.551), FT(-2.2209)]

    elseif IN_mode == "ABIFM"
        Nₐ = FT(0)
        Nₗ = FT(8e6)
        Nᵢ = FT(0)
        r₀ = FT(0.5e-6)                 # droplets maybe bigger?
        p₀ = FT(987.018 * 1e2)
        T₀ = FT(212.978)

        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        C_l = FT(qₗ / ((1 - qₗ) * ϵₘ + qₗ))  # concentration/mol fraction of liquid
        C_v = FT(10.8509 * 1e-6 - C_l)     # concentration/mol fraction of vapor
        qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
        qᵢ = FT(0)

        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
        w = FT(0.4)
        IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
        Sᵢ = ξ(tps, T₀) * Sₗ

        heterogeneous = "ABIFM"
        condensation_growth = "Condensation"
        deposition_growth = "Deposition"
        size_distribution = "Monodisperse"

        params = (;
            const_dt,
            w,
            heterogeneous,
            condensation_growth,
            deposition_growth,
            size_distribution,
        )

        coeff_true = [FT(54.58834), FT(-10.54758)]

    elseif IN_mode == "ABHOM"
        Nₐ = FT(0)
        # Nₗ = FT(360 * 1e6)
        Nₗ = FT(300 * 1e6)
        Nᵢ = FT(0)
        #r₀ = FT(2e-7)                 # droplets maybe bigger?
        r₀ = FT(25 * 1e-9)
        #p₀ = FT(987.345 * 1e2)
        p₀ = FT(9712.183)
        # T₀ = FT(243.134)
        T₀ = FT(190)

        qₗ = FT(Nₗ * 4 / 3 * FT(π) * r₀^3 * ρₗ / FT(1.2)) # 1.2 should be ρₐ
        C_l = FT(qₗ / ((1 - qₗ) * ϵₘ + qₗ))  # concentration/mol fraction of liquid
        #C_v = FT(357.096 * 1e-6 - C_l)     # concentration/mol fraction of vapor
        C_v = FT(5 * 1e-6)
        qᵥ = ϵₘ / (ϵₘ - 1 + 1 / C_v)
        qᵢ = FT(0)

        q = TD.PhasePartition.(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
        Rₐ = TD.gas_constant_air(tps, q)
        eₛ = TD.saturation_vapor_pressure(tps, T₀, TD.Liquid())
        e = eᵥ(qᵥ, p₀, Rₐ, R_v)
        Sₗ = FT(e / eₛ)
        w = FT(1)
        IC = [Sₗ, p₀, T₀, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, FT(0)]
        Sᵢ = ξ(tps, T₀) * Sₗ

        homogeneous = "ABHOM"
        deposition_growth = "Deposition"
        size_distribution = "Monodisperse"

        params =
            (; const_dt, w, homogeneous, deposition_growth, size_distribution)

        coeff_true = [FT(255.927125), FT(-68.553283)]
    end

    n_samples = 10
    G_truth = run_model(params, coeff_true, IN_mode, FT, IC)     # ICNC from running parcel w default values
    y_truth = zeros(length(G_truth), n_samples) # where noisy ICNC will be stored

    dim_output = length(G_truth)
    Γ = 0.005 * LinearAlgebra.I * (maximum(G_truth) - minimum(G_truth))
    noise_dist = Distributions.MvNormal(zeros(dim_output), Γ)
    y_truth = G_truth .+ rand(noise_dist)

    # Generate initial ensemble and set up EKI
    # TODO - make compatible with Float32. Only works with Float64 for now.
    N_ensemble = 10      # runs N_ensemble trials per iteration
    N_iterations = 15    # number of iterations the inverse problem goes through
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

    return [m_coeff_ekp, c_coeff_ekp, coeff_true[1], coeff_true[2], ϕ_n_values]
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
