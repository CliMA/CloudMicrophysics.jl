module PreprocessAerosolData

import CloudMicrophysics as CM
import ..AerosolActivation as AA
import ..AerosolModel as AM
import ..Parameters as CMP
import ..CommonTypes as CT
import ..Parameters as CMP
import CLIMAParameters as CP
import Thermodynamics as TD
import DataFrames as DF
using DataFramesMeta

const FT = Float64
const APS = CMP.AbstractCloudMicrophysicsParameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
default_param_set = cloud_microphysics_parameters(toml_dict)

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

function convert_to_ARG_params(data_row::NamedTuple, param_set::APS)
    num_modes = get_num_modes(data_row)
    @assert num_modes > 0
    mode_Ns = []
    mode_means = []
    mode_stdevs = []
    mode_kappas = []
    velocity = data_row.velocity
    temperature = data_row.initial_temperature
    pressure = data_row.initial_pressure
    for i in 1:num_modes
        push!(mode_Ns, data_row[Symbol("mode_$(i)_N")])
        push!(mode_means, data_row[Symbol("mode_$(i)_mean")])
        push!(mode_stdevs, data_row[Symbol("mode_$(i)_stdev")])
        push!(mode_kappas, data_row[Symbol("mode_$(i)_kappa")])
    end
    ad = AM.AerosolDistribution(
        Tuple(
            AM.Mode_κ(
                mode_means[i],
                mode_stdevs[i],
                mode_Ns[i],
                FT(1),
                FT(1),
                FT(0),
                mode_kappas[i],
                1,
            ) for i in 1:num_modes
        ),
    )
    thermo_params = CMP.thermodynamics_params(param_set)
    pv0 = TD.saturation_vapor_pressure(thermo_params, temperature, TD.Liquid())
    vapor_mix_ratio =
        pv0 / TD.Parameters.molmass_ratio(thermo_params) / (pressure - pv0)
    q_vap = vapor_mix_ratio / (vapor_mix_ratio + 1)
    q = TD.PhasePartition(q_vap, FT(0), FT(0))
    return (; ad, temperature, pressure, velocity, q, mode_Ns)
end

function convert_to_ARG_params(data_row::NamedTuple)
    return convert_to_ARG_params(data_row, default_param_set)
end

function get_ARG_S_max(data_row::NamedTuple, param_set::APS)
    (; ad, temperature, pressure, velocity, q) =
        convert_to_ARG_params(data_row, param_set)
    max_supersaturation = AA.max_supersaturation(
        param_set,
        CT.ARG2000Type(),
        ad,
        temperature,
        pressure,
        velocity,
        q,
    )
    return max_supersaturation
end

function get_ARG_S_max(data_row::NamedTuple)
    return get_ARG_S_max(data_row, default_param_set)
end

function get_ARG_S_max(X::DataFrame, param_set::APS)
    return get_ARG_S_max.(NamedTuple.(eachrow(X)), param_set)
end

function get_ARG_S_max(X::DataFrame)
    return get_ARG_S_max(X, default_param_set)
end

function get_ARG_S_crit(data_row::NamedTuple, param_set::APS)
    (; ad, temperature) = convert_to_ARG_params(data_row, param_set)
    return AA.critical_supersaturation(param_set, ad, temperature)
end

function get_ARG_S_crit(data_row::NamedTuple)
    return get_ARG_S_crit(data_row, default_param_set)
end

function get_ARG_S_crit(X::DataFrame, param_set::APS)
    return get_ARG_S_crit.(NamedTuple.(eachrow(X)), param_set)
end

function get_ARG_S_crit(X::DataFrame)
    return get_ARG_S_crit(X, default_param_set)
end

function get_ARG_act_N(
    data_row::NamedTuple,
    param_set::APS,
    S_max = nothing,
)
    (; ad, temperature, pressure, velocity, q) =
        convert_to_ARG_params(data_row, param_set)
    if S_max === nothing
        return collect(
            AA.N_activated_per_mode(
                param_set,
                CT.ARG2000Type(),
                ad,
                temperature,
                pressure,
                velocity,
                q,
            ),
        )
    else
        critical_supersaturation =
            AA.critical_supersaturation(param_set, ad, temperature)
        return collect(
            AA.N_activated_per_mode(
                param_set,
                CT.ARG2000Type(),
                ad,
                temperature,
                pressure,
                velocity,
                q,
                S_max,
                critical_supersaturation,
            ),
        )
    end
end

function get_ARG_act_N(data_row::NamedTuple, S_max = nothing)
    return get_ARG_act_N(data_row, default_param_set, S_max)
end

function get_ARG_act_N(X::DataFrame, param_set::APS, S_max = nothing)
    return transpose(
        hcat(get_ARG_act_N.(NamedTuple.(eachrow(X)), param_set, S_max)...),
    )
end

function get_ARG_act_N(X::DataFrame, S_max = nothing)
    return get_ARG_act_N(X, default_param_set, S_max)
end

function get_ARG_act_frac(
    data_row::NamedTuple,
    param_set::APS,
    S_max = nothing,
)
    (; mode_Ns) = convert_to_ARG_params(data_row, param_set)
    return get_ARG_act_N(data_row, param_set, S_max) ./ mode_Ns
end

function get_ARG_act_frac(data_row::NamedTuple, S_max = nothing)
    return get_ARG_act_frac(data_row, default_param_set, S_max)
end

function get_ARG_act_frac(X::DataFrame, param_set::APS, S_max = nothing)
    return transpose(
        hcat(get_ARG_act_frac.(NamedTuple.(eachrow(X)), param_set, S_max)...),
    )
end

function get_ARG_act_frac(X::DataFrame, S_max = nothing)
    return get_ARG_act_frac(X, default_param_set, S_max)
end

function preprocess_aerosol_data(X::DataFrame, add_ARG_act_frac::Bool)
    num_modes = get_num_modes(X)
    if add_ARG_act_frac
        X = DF.transform(
            X,
            AsTable(All()) =>
                ByRow(x -> get_ARG_act_frac(x)) =>
                    [Symbol("mode_$(i)_ARG_act_frac") for i in 1:num_modes],
        )
    end
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

function preprocess_aerosol_data_standard(X::DataFrame)
    return preprocess_aerosol_data(X, false)
end

function preprocess_aerosol_data_with_ARG_act_frac(X::DataFrame)
    return preprocess_aerosol_data(X, true)
end

function preprocess_aerosol_data_with_ARG_intermediates(X::DataFrame)
    
end

function target_transform(act_frac)
    return @. atanh(2.0 * 0.99 * (act_frac - 0.5))
end

function inverse_target_transform(transformed_act_frac)
    return @. (1.0 / (2.0 * 0.99)) * tanh(transformed_act_frac) + 0.5
end

end # module
