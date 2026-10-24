using Test: @testset, @test, @test_throws, @inferred
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import ClimaParams as CP
import ForwardDiff as FD
import BenchmarkTools as BT
import Adapt
import Random

# Adaptor that wraps every array in a view, for the Adapt path of the tables
struct ViewAdaptor end
Adapt.adapt_storage(::ViewAdaptor, a::Array) = view(a, axes(a)...)

# Point counts of the axes for the tests, coarser than the default grid
const TEST_ICE_AXES = (; logλ = 64, w = 9, ρ_rim = 5, logρₐ = 8)
const TEST_GRID = (;
    ice = TEST_ICE_AXES,
    liquid_ice = (; logλ = 48, w = 9, ρ_rim = 4, logρₐ = 5, cloud_logx̄ = 10, rain_logx̄ = 6),
    rime_density = (; logλ = 12, w = 4, ρ_rim = 3, logρₐ = 3, cloud_logx̄ = 5, rain_logx̄ = 4, logR = 24),
    wet_set = (;
        logλ = 16, w = 5, ρ_rim = 3, logρₐ = 3, cloud_logx̄ = 5, rain_logx̄ = 4,
        logD = 48, number_logit = 17, μ = 32,
    ),
    particle_rime_density = (; logD = 16, w = 3, ρ_rim = 3, logρₐ = 3, cloud_logx̄ = 4, rain_logx̄ = 4, logR = 12),
)

# A terminal velocity that throws when evaluated, so a test detects any quadrature node
struct ThrowingVelocity{FT}
    cutoff::FT
end
CO.particle_terminal_velocity(::ThrowingVelocity, ρₐ, ρ_i) = D -> error("a quadrature node was evaluated")
CO.particle_terminal_velocity(::ThrowingVelocity, ρₐ) = D -> error("a quadrature node was evaluated")

# Seeded states `(log x̄, F_rim, ρ_rim)` between the grid points: a fifth unrimed, a fifth with
# `F_rim` log-uniform in [1e-8, 1e-2], two fifths with `F_rim` uniform in [0.01, 0.9], and a fifth
# with `1 - F_rim` log-uniform in [1e-3, 1e-1]
function sample_ice_states(params, n)
    rng = Random.Xoshiro(1)
    logx̄_lo, logx̄_hi = log(P3.ice_mean_particle_mass_min(params)), log(P3.ice_mean_particle_mass_max(params))
    ρ_rim_lo, ρ_rim_hi = P3.rime_density_bounds(params)
    return map(1:n) do i
        logx̄ = logx̄_lo + (logx̄_hi - logx̄_lo) * rand(rng)
        ρ_rim = ρ_rim_lo + (ρ_rim_hi - ρ_rim_lo) * rand(rng)
        F_rim =
            i % 5 == 0 ? 0.0 :
            i % 5 == 1 ? exp10(-8 + 6 * rand(rng)) :
            i % 5 == 4 ? 1 - exp10(-1 - 2 * rand(rng)) : 0.01 + 0.89 * rand(rng)
        (logx̄, F_rim, ρ_rim)
    end
end

# The states of `sample_ice_states` with seeded air densities, log-uniform in [0.005, 1.7] kg/m³
function sample_ice_states_in_air(params, n)
    rng = Random.Xoshiro(2)
    return [(s..., exp(log(0.005) + log(1.7 / 0.005) * rand(rng))) for s in sample_ice_states(params, n)]
end

# `logλ` from a Float64 Brent solve to the rounding floor
function reference_logλ(params, logx̄, F_rim, ρ_rim)
    particles = P3.P3Particles(params, F_rim, ρ_rim)
    lo, hi = P3.logλ_bounds(params)
    f(logλ) = P3.logLdivN(particles, logλ) - logx̄
    method = P3.RS.BrentsMethod(lo, hi)
    return P3.RS.find_zero(f, method, P3.RS.CompactSolution(), P3.RS.SolutionTolerance(1e-14), 200).root
end

# The ice integrals of a state under the integration `rule`, at 2 °C for the melting
function ice_integrals(ice, rule, (logx̄, F_rim, ρ_rim, ρₐ), aps, tps)
    FT = eltype(ice.scheme)
    ρn_ice = FT(1e4)
    ρq_ice = FT(exp(logx̄)) * ρn_ice
    ρq_rim = FT(F_rim) * ρq_ice
    st = P3.P3State(ice.scheme, rule, ρq_ice, ρn_ice, ρq_rim, ρq_rim / FT(ρ_rim))
    vel = ice.terminal_velocity
    ρₐ = FT(ρₐ)
    (; M_1, M_v) = P3.ice_ventilation_moments(vel, st, ρₐ)
    (; v_n, v_m) = P3.ice_terminal_velocities(vel, st, ρₐ)
    return (;
        M_1,
        M_v,
        ventilation = P3.ice_ventilation_integral(vel, aps, st, ρₐ),
        self_collection = P3.ice_self_collection(vel, st, ρₐ).dNdt,
        v_n,
        v_m,
        melt = P3.ice_melt(vel, aps, tps, st, ρₐ, FT(275.15)).dLdt,
    )
end

# The median, 99th percentile and maximum of `errors`
function quantiles(errors)
    e = sort(collect(errors))
    return (; median = e[cld(length(e), 2)], p99 = e[ceil(Int, 0.99 * length(e))], max = e[end])
end

# The quantiles of the relative errors of `tabulated` against the nonzero values of `reference`,
# and the number of states where only one of the two is zero
function error_quantiles(reference, tabulated)
    nonzero = [(r, t) for (r, t) in zip(reference, tabulated) if r != 0]
    mismatch = count(((r, t),) -> iszero(r) != iszero(t), zip(reference, tabulated))
    return (; quantiles(abs(t / r - 1) for (r, t) in nonzero)..., mismatch)
end

function test_ice_table_accuracy(ice, tables; n_states, median, p99)
    FT = eltype(ice.scheme)
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    states = sample_ice_states_in_air(CMP.ParametersP3(CP.create_toml_dict(Float64)), n_states)
    # Dense cold surface air and thin air above 20 km
    states = [states; [(s[1:3]..., ρₐ) for s in states[1:10] for ρₐ in (0.02, 1.5)]]
    reference = [ice_integrals(ice, ice.integration, s, aps, tps) for s in states]
    tabulated = [ice_integrals(ice, tables, s, aps, tps) for s in states]
    for name in keys(first(reference))
        q = error_quantiles(getproperty.(reference, name), getproperty.(tabulated, name))
        @test q.median < median
        @test q.p99 < p99
        @test q.mismatch == 0
    end
end

# The tables of the point functions `points`, a `NamedTuple` with keys of `P3.table_axes`, on `grid`
function tables_of(ice, grid, points, partition_rule = nothing)
    FT = eltype(ice.scheme)
    axes = P3.table_axes(grid)
    groups =
        map((name, f) -> P3.group_table(Val(name), P3.tabulate(f, axes[name]), grid, FT), keys(points), values(points))
    return P3.P3Tables(NamedTuple{keys(points)}(groups), partition_rule, P3.table_fingerprint(ice))
end
