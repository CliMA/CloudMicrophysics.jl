import MLJ
import StatsBase
import Random
import EnsembleKalmanProcesses as EKP
using EnsembleKalmanProcesses.ParameterDistributions

# Get the testing package
import Test as TT

# Get the CliMA packages
import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP

# Load aerosol data reading and preprocessing functions
include(joinpath(pkgdir(CM), "ext", "Common.jl"))

function calibrate_ARG(
    FT;
    fname = "2modal_dataset1_train.csv",
    sample_size = 500,
    N_samples = 1,
    N_ensemble = 60,
    N_iterations_per_sample = 20,
)
    aip = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    # download (if necessary) and load the data
    fpath = joinpath(pkgdir(CM), "test", "data", fname)
    if !isfile(fpath)
        # see above note on downloading from dropbox
        url = "https://www.dropbox.com/scl/fi/qgq6ujvqenebjkskqvht5/2modal_dataset1_train.csv?rlkey=53qtqz0mtce993gy5jtnpdfz5&dl=0"
        download(url, fname)
        mkpath(dirname(fpath))
        mv(fname, fpath, force = true)
    end
    X_train, Y_train, initial_data = read_aerosol_dataset(fpath)

    rng = Random.MersenneTwister(1)

    prior = combine_distributions([
        constrained_gaussian("f_coeff_1_ARG2000", 0.5, 0.5, 0, Inf),
        constrained_gaussian("f_coeff_2_ARG2000", 2.5, 0.5, 0, Inf),
        constrained_gaussian("g_coeff_1_ARG2000", 1.0, 0.5, 0, Inf),
        constrained_gaussian("g_coeff_2_ARG2000", 0.25, 0.5, 0, Inf),
        constrained_gaussian("pow_1_ARG2000", 1.5, 0.5, 0, Inf),
        constrained_gaussian("pow_2_ARG2000", 0.75, 0.5, 0, Inf),
    ])

    all_params_means = []
    all_mean_error_metrics = []

    for i in 1:N_samples
        sample_inds = StatsBase.sample(1:DF.nrow(X_train), sample_size, replace = false)
        X_sample = X_train[sample_inds, :]
        Y_sample = Y_train[sample_inds]

        initial_ensemble = EKP.construct_initial_ensemble(rng, prior, N_ensemble)
        ensemble_kalman_process = EKP.EnsembleKalmanProcess(
            initial_ensemble, Y_sample, 1.0 * EKP.I, EKP.Inversion();
            localization_method = EKP.Localizers.NoLocalization(),
        )

        mean_error_metrics = []

        for j in 1:N_iterations_per_sample
            params_cur = EKP.get_ϕ_final(prior, ensemble_kalman_process)
            errs = calibration_error_metrics(X_sample, Y_sample, params_cur, aip, tps, FT)
            push!(mean_error_metrics, StatsBase.mean(errs))

            get_pred(p) = calibrated_prediction(X_sample, p, aip, tps, FT)
            pred_ens = hcat(map(get_pred, eachcol(params_cur))...)

            terminate = EKP.update_ensemble!(ensemble_kalman_process, pred_ens)
            terminate == true && break
        end

        params_final = EKP.get_ϕ_final(prior, ensemble_kalman_process)
        errs = calibration_error_metrics(X_sample, Y_sample, params_final, aip, tps, FT)
        push!(mean_error_metrics, StatsBase.mean(errs))

        params_mean_final = EKP.get_ϕ_mean_final(prior, ensemble_kalman_process)
        push!(all_params_means, params_mean_final)
        push!(all_mean_error_metrics, mean_error_metrics)
    end
    final_params_means = StatsBase.mean(all_params_means)

    return final_params_means, all_mean_error_metrics
end

function test_emulator(FT; rtols = [1e-4, 1e-3, 0.26], N_samples_calib = 2)

    aip = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    ap = CMP.AerosolActivationParameters(FT)

    # Atmospheric conditions
    T = FT(294)    # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)    # vertical velocity m/s
    p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    q_vs = 1 / (1 - 1 / TDI.Rd_over_Rv(tps) * (p_vs - p) / p_vs)

    # Aerosol size distribution
    salt = CMP.Seasalt(FT)
    # Accumulation mode
    r1 = FT(0.243 * 1e-6) # m
    σ1 = FT(1.4)          # -
    N1 = FT(100 * 1e6)    # 1/m3
    # Coarse Mode
    r2 = FT(1.5 * 1e-6)   # m
    σ2 = FT(2.1)          # -
    N2 = FT(1e6)          # 1/m3
    acc = AM.Mode_κ(r1, σ1, N1, (FT(1.0),), (FT(1.0),), (salt.M,), (salt.κ,))
    crs = AM.Mode_κ(r2, σ2, N2, (FT(1.0),), (FT(1.0),), (salt.M,), (salt.κ,))
    ad = AM.AerosolDistribution(crs, acc)

    calib_params, errs = calibrate_ARG(FT, N_samples = N_samples_calib)
    ap_calib = CMP.AerosolActivationParameters(calib_params)

    TT.@test AA.N_activated_per_mode(ap_calib, ad, aip, tps, T, p, w, q_vs, FT(0), FT(0))[1] ≈
             AA.N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_vs, FT(0), FT(0))[1] rtol =
        rtols[1]
    TT.@test AA.N_activated_per_mode(ap_calib, ad, aip, tps, T, p, w, q_vs, FT(0), FT(0))[2] ≈
             AA.N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_vs, FT(0), FT(0))[2] rtol =
        rtols[2]
    for n in 1:N_samples_calib
        TT.@test errs[n][end] < rtols[3]
    end
end

TT.@testset "Aerosol activation test: Calibrated ARG" begin
    test_emulator(Float32, N_samples_calib = 5, rtols = [1e-4, 1e-3, 0.26])
end
nothing
