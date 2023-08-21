module ActivationEmulatorModels

import MLJ
import MLJFlux
import Flux
import DataFrames as DF
import GaussianProcesses
import StatsBase
import Distributions

struct NNBuilder <: MLJFlux.Builder
    n1::Int
    n2::Int
    n3::Int
    dropout1::Float64
    dropout2::Float64
    dropout3::Float64
end

function MLJFlux.build(builder::NNBuilder, rng, n_in, n_out)
    init = Flux.glorot_uniform(rng)
    return Flux.Chain(
        Flux.Dense(n_in => builder.n1, Flux.relu, init = init),
        Flux.Dropout(builder.dropout1),
        Flux.Dense(builder.n1 => builder.n2, Flux.relu, init = init),
        Flux.Dropout(builder.dropout2),
        Flux.Dense(builder.n2 => builder.n3, Flux.relu, init = init),
        Flux.Dropout(builder.dropout3),
        Flux.Dense(builder.n3 => n_out, init = init),
    )
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
