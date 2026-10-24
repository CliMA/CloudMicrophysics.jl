using Test: @testset, @test, @test_throws, @inferred
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
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

# `logλ` from a Float64 Brent solve to the rounding floor
function reference_logλ(params, logx̄, F_rim, ρ_rim)
    particles = P3.P3Particles(params, F_rim, ρ_rim)
    lo, hi = P3.logλ_bounds(params)
    f(logλ) = P3.logLdivN(particles, logλ) - logx̄
    method = P3.RS.BrentsMethod(lo, hi)
    return P3.RS.find_zero(f, method, P3.RS.CompactSolution(), P3.RS.SolutionTolerance(1e-14), 200).root
end

# The tables of the point functions `points`, a `NamedTuple` with keys of `P3.table_axes`, on `grid`
function tables_of(ice, grid, points, partition_rule = nothing)
    FT = eltype(ice.scheme)
    axes = P3.table_axes(grid)
    groups =
        map((name, f) -> P3.group_table(Val(name), P3.tabulate(f, axes[name]), grid, FT), keys(points), values(points))
    return P3.P3Tables(NamedTuple{keys(points)}(groups), partition_rule, P3.table_fingerprint(ice))
end
