module PreprocessAerosolData

import CloudMicrophysics as CM
import ..AerosolActivation as AA
import ..AerosolModel as AM
import ..Parameters as CMP
import ..CommonTypes as CT
import ..Parameters as CMP
import ..Common as CO
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
    pv0 = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    vapor_mix_ratio =
        pv0 / TD.Parameters.molmass_ratio(thermo_params) / (p - pv0)
    q_vap = vapor_mix_ratio / (vapor_mix_ratio + 1)
    q = TD.PhasePartition(q_vap, FT(0), FT(0))
    return (; ad, T, p, w, q, mode_Ns)
end

function convert_to_ARG_params(data_row::NamedTuple)
    return convert_to_ARG_params(data_row, default_param_set)
end

function convert_to_ARG_intermediates(data_row::NamedTuple, param_set::APS)
    num_modes = get_num_modes(data_row)
    @assert num_modes > 0

    (; ad, T, p, w, q) = convert_to_ARG_params(data_row, param_set)

    thermo_params = CMP.thermodynamics_params(param_set)
    _grav::FT = CMP.grav(param_set)
    _ρ_cloud_liq::FT = CMP.ρ_cloud_liq(param_set)

    _ϵ::FT = 1 / CMP.molmass_ratio(param_set)
    R_m::FT = TD.gas_constant_air(thermo_params, q)
    cp_m::FT = TD.cp_m(thermo_params, q)

    L::FT = TD.latent_heat_vapor(thermo_params, T)
    p_vs::FT = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    G::FT = CO.G_func(param_set, T, TD.Liquid()) / _ρ_cloud_liq

    α::FT = L * _grav * _ϵ / R_m / cp_m / T^2 - _grav / R_m / T
    γ::FT = R_m * T / _ϵ / p_vs + _ϵ * L^2 / cp_m / T / p

    A::FT = AA.coeff_of_curvature(param_set, T)
    ζ::FT = 2 * A / 3 * sqrt(α * w / G)

    Sm = AA.critical_supersaturation(param_set, ad, T)
    η = [
        (α * w / G)^FT(3 / 2) / (FT(2 * pi) * _ρ_cloud_liq * γ * ad.Modes[i].N) for i in 1:num_modes
    ]

    per_mode_intermediates = [
        (;
            Symbol("mode_$(i)_log_stdev") => log(ad.Modes[i].stdev),
            Symbol("mode_$(i)_η") => η[i],
            Symbol("mode_$(i)_Sm") => Sm[i],
            Symbol("mode_$(i)_term1") => (ζ / η[i])^(3 / 2),
            Symbol("mode_$(i)_term2") => (Sm[i]^2 / (η[i] + 3 * ζ))^(3 / 4),
        ) for i in 1:num_modes
    ]

    return merge(reduce(merge, per_mode_intermediates), (; ζ))
end

function convert_to_ARG_intermediates(data_row::NamedTuple)
    return convert_to_ARG_intermediates(data_row, default_param_set)
end

function get_ARG_S_max(data_row::NamedTuple, param_set::APS)
    (; ad, T, p, w, q) = convert_to_ARG_params(data_row, param_set)
    max_supersaturation =
        AA.max_supersaturation(param_set, CT.ARG2000Type(), ad, T, p, w, q)
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
    (; ad, T) = convert_to_ARG_params(data_row, param_set)
    return AA.critical_supersaturation(param_set, ad, T)
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

function get_ARG_act_N(data_row::NamedTuple, param_set::APS, S_max = nothing)
    (; ad, T, p, w, q) = convert_to_ARG_params(data_row, param_set)
    if S_max === nothing
        return collect(
            AA.N_activated_per_mode(
                param_set,
                CT.ARG2000Type(),
                ad,
                T,
                p,
                w,
                q,
            ),
        )
    else
        critical_supersaturation = AA.critical_supersaturation(param_set, ad, T)
        return collect(
            AA.N_activated_per_mode(
                param_set,
                CT.ARG2000Type(),
                ad,
                T,
                p,
                w,
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

function get_ARG_act_frac(data_row::NamedTuple, param_set::APS, S_max = nothing)
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
    return DF.DataFrame(convert_to_ARG_intermediates.(NamedTuple.(eachrow(X))))
end

function target_transform(act_frac)
    return @. atanh(2.0 * 0.99 * (act_frac - 0.5))
end

function inverse_target_transform(transformed_act_frac)
    return @. (1.0 / (2.0 * 0.99)) * tanh(transformed_act_frac) + 0.5
end

end # module
