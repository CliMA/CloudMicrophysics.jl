module ActivationEmulatorModels

import MLJ
import MLJFlux
import Flux
import DataFrames as DF
import GaussianProcesses
import StatsBase
import Distributions

struct NNBuilder <: MLJFlux.Builder
    layer_sizes::Vector{Integer}
    dropout::Vector{Float64}
end

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

mutable struct GPRegressor <: MLJ.Deterministic
    num_gps::Integer
    sample_size::Integer
    use_ARG_weights::Bool
    use_DTC::Bool
    sample_size_inducing::Integer
end

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

end # module
