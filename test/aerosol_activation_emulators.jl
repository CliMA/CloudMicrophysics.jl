# Get ML packages
import MLJ
import Flux
import MLJModels
import MLJFlux
import GaussianProcesses
import StatsBase

Standardizer = MLJ.@load Standardizer pkg = MLJModels
NeuralNetworkRegressor = MLJ.@load NeuralNetworkRegressor pkg = MLJFlux
EvoTreeRegressor = MLJ.@load EvoTreeRegressor pkg = EvoTrees

# Get the testing package
import Test as TT

# Get the CliMA packages
import CloudMicrophysics as CM
import ClimaParams as CP
import Thermodynamics as TD
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP

# NN helpers
# Container for the NN we are about to build
struct NNBuilder <: MLJFlux.Builder
    layer_sizes::Vector{Integer}
    dropout::Vector
end

# Define the NN structure
function MLJFlux.build(builder::NNBuilder, rng, n_in, n_out)
    @assert length(builder.layer_sizes) == length(builder.dropout)
    num_hidden_layers = length(builder.layer_sizes)
    init = Flux.glorot_uniform(rng)
    layers::Vector{Any} = []
    if num_hidden_layers == 0
        push!(layers, Flux.Dense(n_in => n_out, init = init))
    else
        push!(
            layers,
            Flux.Dense(
                n_in => builder.layer_sizes[1],
                Flux.sigmoid_fast,
                init = init,
            ),
        )
    end
    for i in 1:num_hidden_layers
        push!(layers, Flux.Dropout(builder.dropout[i]))
        if i == num_hidden_layers
            push!(
                layers,
                Flux.Dense(builder.layer_sizes[i] => n_out, init = init),
            )
        else
            push!(
                layers,
                Flux.Dense(
                    builder.layer_sizes[i] => builder.layer_sizes[i + 1],
                    Flux.sigmoid_fast,
                    init = init,
                ),
            )
        end
    end
    return Flux.Chain(layers...)
end

# GP helpers
# Container for the GP we want to build
mutable struct GPRegressor <: MLJ.Deterministic
    num_gps::Integer
    sample_size::Integer
    use_ARG_weights::Bool
    use_DTC::Bool
    sample_size_inducing::Integer
end

# Define the function that fits the GP model
function MLJ.fit(model::GPRegressor, verbosity, X, y)
    gps = []
    for i in 1:(model.num_gps)
        if model.use_ARG_weights
            weights = StatsBase.Weights([
                Distributions.pdf(Distributions.Normal(0.0, 0.5), x) for
                x in X.mode_1_ARG_act_frac
            ])
            inds = StatsBase.sample(
                1:DF.nrow(X),
                weights,
                model.sample_size,
                replace = false,
            )
            inds_inducing = StatsBase.sample(
                1:DF.nrow(X),
                weights,
                model.sample_size_inducing,
                replace = false,
            )
        else
            inds = StatsBase.sample(
                1:DF.nrow(X),
                model.sample_size,
                replace = false,
            )
            inds_inducing = StatsBase.sample(
                1:DF.nrow(X),
                model.sample_size_inducing,
                replace = false,
            )
        end
        if model.use_DTC
            gp = GaussianProcesses.DTC(
                Matrix(X[inds, :])',
                Matrix(X[inds_inducing, :])',
                y[inds],
                GaussianProcesses.MeanConst(StatsBase.mean(y)),
                GaussianProcesses.SEArd(fill(4.0, DF.ncol(X)), 0.0),
                2.0,
            )
        else
            gp = GaussianProcesses.GPA(
                Matrix(X[inds, :])',
                y[inds],
                GaussianProcesses.MeanConst(StatsBase.mean(y)),
                GaussianProcesses.SEArd(fill(4.0, DF.ncol(X)), 0.0),
                GaussianProcesses.GaussLik(2.0),
            )
        end
        GaussianProcesses.optimize!(gp)
        push!(gps, gp)
    end
    return gps, nothing, nothing
end

# Define the function that predicts by using the GP model
function MLJ.predict(::GPRegressor, fitresult, Xnew)
    means = reduce(
        hcat,
        [GaussianProcesses.predict_f(gp, Matrix(Xnew)')[1] for gp in fitresult],
    )
    variances = reduce(
        hcat,
        [GaussianProcesses.predict_f(gp, Matrix(Xnew)')[2] for gp in fitresult],
    )
    return (sum(
        means ./ variances,
        dims = 2,
    ) ./ sum(1.0 ./ variances, dims = 2))[
        :,
        1,
    ]
end

# Load aerosol data reading and preprocessing functions
include(joinpath(pkgdir(CM), "ext", "Common.jl"))

function preprocess_aerosol_data_with_ARG_act_frac_FT32(x::DataFrame)
    aip = CMP.AirProperties(Float32)
    tps = TD.Parameters.ThermodynamicsParameters(Float32)
    ap = CMP.AerosolActivationParameters(Float32)
    return preprocess_aerosol_data_with_ARG_act_frac(x, ap, aip, tps, Float32)
end

# Get the ML model
function get_2modal_model_FT32(;
    ml_model = "NN",
    with_ARG_act_frac = false,
    machine_name = "2modal_nn_machine_naive.jls",
)
    FT = Float32

    # If the ML model already exists load it in.
    # If it does not exist, train
    fpath = joinpath(pkgdir(CM), "test")
    if isfile(joinpath(fpath, machine_name))
        # Read-in the saved ML model
        emulator_filepath = joinpath(fpath, machine_name)
        return MLJ.machine(emulator_filepath)
    else
        # Download data files if you don't have them
        # (not the most elegant way...)
        fname = "2modal_dataset1_train.csv"
        if !isfile(joinpath(fpath, "data", fname))
            # For reasons unknown to humankind uploading/downloading the file
            # from Caltech box changes the file. We are using dropbox for now.
            # In the end we should upload the files to Caltech Data
            # (and pray they are not changed there...)
            url = "https://www.dropbox.com/scl/fi/qgq6ujvqenebjkskqvht5/2modal_dataset1_train.csv?rlkey=53qtqz0mtce993gy5jtnpdfz5&dl=0"
            download(url, fname)
            mkdir(joinpath(fpath, "data"))
            mv(fname, joinpath(fpath, "data", fname), force = true)
        end

        # Read the aerosol data
        X_train, Y_train, initial_data =
            read_aerosol_dataset(joinpath(fpath, "data", fname))

        # Define the training pipeline
        if with_ARG_act_frac
            _preprocess_aerosol_data =
                preprocess_aerosol_data_with_ARG_act_frac_FT32
        else
            _preprocess_aerosol_data = preprocess_aerosol_data
        end
        if ml_model == "NN"
            pipeline =
                _preprocess_aerosol_data |>
                Standardizer() |>
                MLJ.TransformedTargetModel(
                    NeuralNetworkRegressor(
                        builder = NNBuilder(
                            [250, 50, 5],
                            [FT(0.3), FT(0), FT(0)],
                        ),
                        optimiser = Flux.Optimise.Adam(
                            0.001,
                            (0.9, 0.999),
                            1.0e-8,
                        ),
                        epochs = 2000,
                        loss = Flux.mse,
                        batch_size = 1000,
                    ),
                    transformer = target_transform,
                    inverse = inverse_target_transform,
                )
        elseif ml_model == "GP"
            pipeline =
                _preprocess_aerosol_data |>
                Standardizer() |>
                MLJ.TransformedTargetModel(
                    GPRegressor(5, 2000, false, true, 20),
                    transformer = target_transform,
                    inverse = inverse_target_transform,
                )
        elseif ml_model == "ET"
            pipeline =
                _preprocess_aerosol_data |>
                Standardizer() |>
                MLJ.TransformedTargetModel(
                    MLJ.TunedModel(
                        tuning = MLJ.Grid(goal = 30),
                        model = EvoTreeRegressor(),
                        resampling = MLJ.CV(nfolds = 5),
                        range = [
                            range(
                                EvoTreeRegressor(),
                                :eta,
                                lower = 0.05,
                                upper = 1,
                                scale = :log,
                            ),
                            range(
                                EvoTreeRegressor(),
                                :max_depth,
                                lower = 3,
                                upper = 15,
                            ),
                        ],
                    ),
                    transformer = target_transform,
                    inverse = inverse_target_transform,
                )
        else
            error("Unknown ML method!")
        end

        # Create the untrained ML model
        mach = MLJ.machine(pipeline, X_train, Y_train)
        # Train a new ML model
        @info(
            "Training a new ML model ($(ml_model)). This may take a minute or two."
        )
        MLJ.fit!(mach, verbosity = 0)
        # Save the ML model to a binary file that can be re-used by CloudMicrophysics
        MLJ.save(joinpath(fpath, machine_name), mach)
        return mach
    end
end

function test_emulator(
    FT;
    ml_model = "NN",
    with_ARG_act_frac = false,
    machine_name = "2modal_nn_machine_naive.jls",
    rtols = [1e-4, 1e-3],
)

    aip = CMP.AirProperties(FT)
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    ap = CMP.AerosolActivationParameters(FT)

    # Atmospheric conditions
    T = FT(294)    # air temperature K
    p = FT(1e5)    # air pressure Pa
    w = FT(0.5)    # vertical velocity m/s
    p_vs = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    q_vs = 1 / (1 - TD.Parameters.molmass_ratio(tps) * (p_vs - p) / p_vs)
    q = TD.PhasePartition(q_vs)

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
    acc = AM.Mode_κ(r1, σ1, N1, (FT(1.0),), (FT(1.0),), (salt.M,), (salt.κ,), 1)
    crs = AM.Mode_κ(r2, σ2, N2, (FT(1.0),), (FT(1.0),), (salt.M,), (salt.κ,), 1)
    ad = AM.AerosolDistribution((crs, acc))

    # Get the ML model
    mach = get_2modal_model_FT32(
        ml_model = ml_model,
        with_ARG_act_frac = with_ARG_act_frac,
        machine_name = machine_name,
    )

    TT.@test AA.N_activated_per_mode(mach, ap, ad, aip, tps, T, p, w, q)[1] ≈
             AA.N_activated_per_mode(ap, ad, aip, tps, T, p, w, q)[1] rtol =
        rtols[1]
    TT.@test AA.N_activated_per_mode(mach, ap, ad, aip, tps, T, p, w, q)[2] ≈
             AA.N_activated_per_mode(ap, ad, aip, tps, T, p, w, q)[2] rtol =
        rtols[2]
end

@info "Aerosol activation test"

TT.@testset "Neural Network" begin
    test_emulator(
        Float32,
        ml_model = "NN",
        with_ARG_act_frac = false,
        machine_name = "2modal_nn_machine_naive.jls",
        rtols = [1e-4, 1e-3],
    )
end

TT.@testset "Gaussian Processes" begin
    test_emulator(
        Float32,
        ml_model = "GP",
        with_ARG_act_frac = false,
        machine_name = "2modal_gp_machine_naive.jls",
        rtols = [1e-2, 5e-2],
    )
end

TT.@testset "Evo Trees" begin
    test_emulator(
        Float32,
        ml_model = "ET",
        with_ARG_act_frac = false,
        machine_name = "2modal_et_machine_naive.jls",
        rtols = [5e-3, 5e-3],
    )
end

TT.@testset "Neural Network with ARG-informed data" begin
    test_emulator(
        Float32,
        ml_model = "NN",
        with_ARG_act_frac = true,
        machine_name = "2modal_nn_machine_naive_with_ARG_act_frac.jls",
        rtols = [1e-4, 1e-3],
    )
end
