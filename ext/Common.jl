import CSV
import DataFrames as DF
using DataFramesMeta

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
    @info(num_modes)
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

function target_transform(act_frac)
    return @. atanh(2.0 * 0.99 * (act_frac - 0.5))
end

function inverse_target_transform(transformed_act_frac)
    return @. (1.0 / (2.0 * 0.99)) * tanh(transformed_act_frac) + 0.5
end
