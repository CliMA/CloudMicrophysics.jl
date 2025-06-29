import CSV
import DataFrames as DF
using DataFramesMeta
using StatsBase

import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.AerosolActivation as AA
import CloudMicrophysics.ThermodynamicsInterface.TD.Parameters as TDP

function get_num_modes(df::DataFrame)
    i = 1
    while true
        if !("mode_$(i)_N" in names(df))
            return i - 1
        end
        i += 1
    end
end

function get_num_modes(data_row::NamedTuple)
    i = 1
    while true
        if !(Symbol("mode_$(i)_N") in keys(data_row))
            return i - 1
        end
        i += 1
    end
end

function read_aerosol_dataset(
    dataset_filename::String,
    Y_name = :mode_1_act_frac_S_interp,
)
    initial_data = DF.DataFrame(CSV.File(dataset_filename))
    df = filter(row -> row.S_max > 0 && row.S_max < 0.2, initial_data)
    selected_columns_X = []
    num_modes = get_num_modes(df)
    for i in 1:num_modes
        append!(
            selected_columns_X,
            Symbol.([
                "mode_$(i)_N",
                "mode_$(i)_mean",
                "mode_$(i)_stdev",
                "mode_$(i)_kappa",
            ]),
        )
    end
    append!(
        selected_columns_X,
        [:velocity, :initial_temperature, :initial_pressure],
    )
    X = df[:, selected_columns_X]
    Y = df[:, Y_name]
    return (X, Y, initial_data)
end

function preprocess_aerosol_data(X::DataFrame)
    num_modes = get_num_modes(X)
    for i in 1:num_modes
        X = DF.transform(
            X,
            Symbol("mode_$(i)_N") => ByRow(log) => Symbol("mode_$(i)_N"),
        )
        X = DF.transform(
            X,
            Symbol("mode_$(i)_mean") =>
                ByRow(log) => Symbol("mode_$(i)_mean"),
        )
    end
    X = DF.transform(X, :velocity => ByRow(log) => :velocity)
    return X
end

function get_ARG_act_frac(
    data_row::NamedTuple,
    ap::CMP.AerosolActivationParameters,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    FT::DataType,
)

    num_modes = get_num_modes(data_row)
    @assert num_modes > 0
    mode_Ns = []
    mode_means = []
    mode_stdevs = []
    mode_kappas = []
    w = data_row.velocity
    T = data_row.initial_temperature
    p = data_row.initial_pressure
    for i in 1:num_modes
        push!(mode_Ns, data_row[Symbol("mode_$(i)_N")])
        push!(mode_means, data_row[Symbol("mode_$(i)_mean")])
        push!(mode_stdevs, data_row[Symbol("mode_$(i)_stdev")])
        push!(mode_kappas, data_row[Symbol("mode_$(i)_kappa")])
    end
    ad = AM.AerosolDistribution(
        (
            AM.Mode_κ(
                mode_means[i],
                mode_stdevs[i],
                mode_Ns[i],
                FT(1),
                FT(1),
                FT(0),
                FT(mode_kappas[i]),
            ) for i in 1:num_modes
        )...,
    )
    pv0 = TDI.saturation_vapor_pressure_over_liquid(tps, FT(T))
    vapor_mix_ratio = pv0 * TDI.Rd_over_Rv(tps) / (p - pv0)
    q_vap = vapor_mix_ratio / (vapor_mix_ratio + 1)

    return collect(
        AA.N_activated_per_mode(ap, ad, aip, tps, FT(T), FT(p), FT(w), FT(q_vap), FT(0), FT(0)),
    ) ./ mode_Ns
end

function get_ARG_act_frac(
    X::DataFrame,
    ap::CMP.AerosolActivationParameters,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    FT::DataType,
)
    return transpose(
        hcat(get_ARG_act_frac.(NamedTuple.(eachrow(X)), ap, aip, tps, FT)...),
    )
end

function preprocess_aerosol_data_with_ARG_act_frac(
    X::DataFrame,
    ap::CMP.AerosolActivationParameters,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    FT::DataType,
)
    f(y) = get_ARG_act_frac(y, ap, aip, tps, FT)
    num_modes = get_num_modes(X)
    X = DF.transform(
        X,
        AsTable(All()) =>
            ByRow(x -> f(x)) =>
                [Symbol("mode_$(i)_ARG_act_frac") for i in 1:num_modes],
    )
    return preprocess_aerosol_data(X)
end

function target_transform(act_frac)
    return @. atanh(2.0 * 0.99 * (act_frac - 0.5))
end

function inverse_target_transform(transformed_act_frac)
    return @. (1.0 / (2.0 * 0.99)) * tanh(transformed_act_frac) + 0.5
end

function calibration_error_metrics(X, Y, ensemble, aip, tps, FT)
    N_ensemble = length(ensemble[1, :])
    rmse = zeros(N_ensemble)
    for i in 1:N_ensemble
        pred = calibrated_prediction(X, ensemble[:, i], aip, tps, FT)
        rmse[i] = StatsBase.mean(sqrt.((pred .- Y) .^ 2))
    end
    return rmse
end

function calibrated_prediction(X, ensemble, aip, tps, FT)
    param_set = CMP.AerosolActivationParameters(ensemble)
    return get_ARG_act_frac(X, param_set, aip, tps, FT)[:, 1]
end
