# Get ML packages
import MLJ
import Flux
import MLJModels
import MLJFlux
Standardizer = MLJ.@load Standardizer pkg = MLJModels
NeuralNetworkRegressor = MLJ.@load NeuralNetworkRegressor pkg = MLJFlux

# Get the testing package
import Test as TT

# Get the CliMA packages
import CloudMicrophysics as CM
import ClimaParams as CP
import Thermodynamics as TD
import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.Parameters as CMP

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

# Load aerosol data reading and preprocessing functions
include(joinpath(pkgdir(CM), "ext", "Common.jl"))

# Get the ML model
function get_2modal_NN_model_FT32()
    FT = Float32
    machine_name = "2modal_nn_machine_naive.jls"

    # If the ML model already exists load it in.
    # If it does not exist, train a NN
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
        pipeline =
            preprocess_aerosol_data |>
            Standardizer() |>
            MLJ.TransformedTargetModel(
                NeuralNetworkRegressor(
                    builder = NNBuilder([250, 50, 5], [FT(0.3), FT(0), FT(0)]),
                    optimiser = Flux.Optimise.Adam(0.001, (0.9, 0.999), 1.0e-8),
                    epochs = 2000,
                    loss = Flux.mse,
                    batch_size = 1000,
                ),
                transformer = target_transform,
                inverse = inverse_target_transform,
            )
        # Create the untrained ML model
        mach = MLJ.machine(pipeline, X_train, Y_train)
        # Train a new ML model
        @info("Training a new ML model. This may take a minute or two.")
        MLJ.fit!(mach, verbosity = 0)
        # Save the ML model to a binary file that can be re-used by CloudMicrophysics
        MLJ.save(joinpath(fpath, machine_name), mach)
        return mach
    end
end

function test_emulator_NN(FT)

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
    mach = get_2modal_NN_model_FT32()

    TT.@test AA.N_activated_per_mode(mach, ap, ad, aip, tps, T, p, w, q)[1] ≈
             AA.N_activated_per_mode(ap, ad, aip, tps, T, p, w, q)[1] rtol =
        1e-5
    TT.@test AA.N_activated_per_mode(mach, ap, ad, aip, tps, T, p, w, q)[2] ≈
             AA.N_activated_per_mode(ap, ad, aip, tps, T, p, w, q)[2] rtol =
        1e-3
end

@info "Aerosol activation NN test"
test_emulator_NN(Float32)
