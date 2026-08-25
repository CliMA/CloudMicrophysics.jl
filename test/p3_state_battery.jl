"""
Comprehensive, physically-realizable P3 state battery for quadrature/accuracy studies.

Supplements `generate_column_states`/`hail_core_states` (`p3_quadrature_error_study.jl`)
with denser systematic coverage of the full P3 prognostic and environment space, plus
explicit physically-important corners (multi-regime, active wet growth, velocity
crossover in the liquid PSD bulk, extreme rime-density contrast, near-degenerate
populations).  Every returned state is realizable: `P3.state_from_prognostic` succeeds
and its shape parameters, bounds, and mean diameter are finite.

Entry point: `generate_comprehensive_states(FT; n_lhs, seed)`.
"""

import Random
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI

"""
    latin_hypercube(n, d; rng)

`n x d` Latin hypercube sample in `[0, 1)^d`: each column stratified into `n` equal
bins (one sample per bin, uniform offset within the bin), independently permuted
across rows per column.
"""
function latin_hypercube(n::Int, d::Int; rng)
    X = zeros(n, d)
    for j in 1:d
        perm = Random.shuffle(rng, 0:(n - 1))
        X[:, j] = (perm .+ rand(rng, n)) ./ n
    end
    return X
end

log_range(u, lo, hi) = exp(log(lo) + u * (log(hi) - log(lo)))
lin_range(u, lo, hi) = lo + u * (hi - lo)

"""
    build_state(FT; ρ, T, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, F_rim, ρ_rim)

Assemble a `generate_column_states`-format `NamedTuple` from physical prognostic
values, computing `q_rim = F_rim * q_ice` and `b_rim = q_rim / ρ_rim` (zero if
`F_rim` is zero), and a saturated-boundary `q_tot`.
"""
function build_state(::Type{FT}; ρ, T, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, F_rim, ρ_rim) where {FT}
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    q_vs =
        T > TDI.TD.Parameters.T_freeze(tps) ?
        TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ) :
        TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    q_rim = F_rim * q_ice
    b_rim = F_rim > 0 ? q_rim / ρ_rim : zero(FT)
    return (;
        ρ = FT(ρ), T = FT(T),
        q_tot = FT(q_vs + q_lcl + q_rai + (q_ice > 0 ? 1e-6 : 0)),
        q_lcl = FT(q_lcl), n_lcl = FT(n_lcl),
        q_rai = FT(q_rai), n_rai = FT(n_rai),
        q_ice = FT(q_ice), n_ice = FT(n_ice),
        q_rim = FT(q_rim), b_rim = FT(b_rim),
    )
end

"""
    realizable(mp, s)

`true` if `s` (a `build_state` `NamedTuple`) produces a well-posed P3 state under
`mp`: `state_from_prognostic` succeeds, and its shape parameter, bounds, and mean
diameter are finite and correctly ordered.
"""
function realizable(mp, s)
    s.q_ice > 0 && s.n_ice > 0 || return true  # no ice: nothing to validate
    try
        state = P3.state_from_prognostic(mp.ice.scheme, s.ρ * s.q_ice, s.ρ * s.n_ice, s.ρ * s.q_rim, s.ρ * s.b_rim)
        logλ = P3.get_distribution_logλ(state)
        isfinite(logλ) || return false
        μ = P3.get_μ(state, logλ)
        isfinite(μ) || return false
        v_i = P3.ice_particle_terminal_velocity(mp.ice.terminal_velocity, s.ρ, state)
        bnds = P3.velocity_integral_bounds(state, logλ, v_i; p = oftype(s.ρ, 1e-5))
        all(isfinite, bnds) && issorted(bnds) || return false
        D_m = P3.D_m(state, logλ)
        isfinite(D_m) && D_m > 0 || return false
        return true
    catch
        return false
    end
end

"""
    lhs_states(FT, mp; n_lhs = 300, seed = 1)

Systematic Latin-hypercube sweep over `log(q_ice), log(n_ice), F_rim, log(ρ_rim),
T, log(ρₐ)`, plus independent Bernoulli presence draws for cloud liquid and rain
(each ~75% present) with log-uniform magnitude when present. Invalid (non-realizable)
draws are redrawn up to a fixed number of times, then dropped.
"""
function lhs_states(::Type{FT}, mp; n_lhs = 300, seed = 1) where {FT}
    rng = Random.MersenneTwister(seed)
    states = NamedTuple[]
    max_attempts = 3 * n_lhs
    attempts = 0
    while length(states) < n_lhs && attempts < max_attempts
        attempts += 1
        u = rand(rng, 6)
        q_ice = log_range(u[1], FT(1e-7), FT(1e-2))
        n_ice = log_range(u[2], FT(1e-1), FT(1e6))
        F_rim = lin_range(u[3], FT(0), FT(0.99))
        ρ_rim = log_range(u[4], FT(50), FT(800))
        T = lin_range(u[5], FT(220), FT(272.5))
        ρ = log_range(u[6], FT(0.2), FT(1.2))
        q_lcl, n_lcl =
            rand(rng) < 0.75 ?
            (log_range(rand(rng), FT(1e-5), FT(3e-3)), log_range(rand(rng), FT(1e6), FT(3e8))) :
            (zero(FT), zero(FT))
        q_rai, n_rai =
            rand(rng) < 0.75 ?
            (log_range(rand(rng), FT(1e-5), FT(3e-3)), log_range(rand(rng), FT(1e2), FT(1e5))) :
            (zero(FT), zero(FT))
        s = build_state(FT; ρ, T, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, F_rim, ρ_rim)
        realizable(mp, s) && push!(states, s)
    end
    return states
end

"""
    corner_states(FT, mp)

Explicit physically-important corners a random sweep may miss or under-sample:
multi-regime (ice mass-regime and velocity-cutoff breakpoints inside the bulk of
the ice PSD), active wet growth (onset diameter inside the bulk), velocity
crossover inside the liquid PSD bulk, extreme rime-density contrast (`ρ_rim` far
from `ρ_i`), near-degenerate populations (trace ice near the `eps(FT)` threshold), and
above-freezing (melting-active).
"""
function corner_states(::Type{FT}, mp) where {FT}
    states = NamedTuple[]

    # Multi-regime: F_rim/ρ_rim pairs that spread D_th, D_gr, D_cr, and the
    # velocity cutoff across the bulk of a moderately broad ice PSD.
    for F_rim in (FT(0.3), FT(0.5), FT(0.7)), ρ_rim in (FT(300), FT(500), FT(700), FT(900))
        push!(
            states,
            build_state(FT;
                ρ = FT(0.85), T = FT(265), q_lcl = FT(3e-4), n_lcl = FT(1e8),
                q_rai = FT(1e-4), n_rai = FT(1e4), q_ice = FT(3e-4 / 0.85), n_ice = FT(1e5 / 0.85),
                F_rim, ρ_rim,
            ),
        )
    end

    # Active wet growth: warm-for-ice T, strong collection, weak freeze capacity,
    # pushing the onset diameter into the bulk of the ice PSD.
    for T in (FT(270), FT(271), FT(272.5)), (Lc, Nc) in ((FT(1e-3), FT(2e8)), (FT(2e-3), FT(3e8)))
        push!(
            states,
            build_state(FT;
                ρ = FT(0.85), T, q_lcl = Lc / FT(0.85), n_lcl = Nc / FT(0.85),
                q_rai = FT(1e-4 / 0.85), n_rai = FT(1e4 / 0.85),
                q_ice = FT(1e-4 / 0.85), n_ice = FT(1e5 / 0.85), F_rim = FT(0.3), ρ_rim = FT(500),
            ),
        )
    end

    # Velocity crossover inside the liquid bulk: rain mean size/fall speed
    # comparable to the ice terminal velocity across a spread of ice sizes.
    for n_ice in (FT(1e2), FT(1e4), FT(1e6)), (Lr, Nr) in ((FT(3e-4), FT(2e3)), (FT(1e-3), FT(1e4)))
        push!(
            states,
            build_state(FT;
                ρ = FT(0.9), T = FT(263), q_lcl = zero(FT), n_lcl = zero(FT),
                q_rai = Lr / FT(0.9), n_rai = Nr / FT(0.9),
                q_ice = FT(5e-4 / 0.9), n_ice = n_ice / FT(0.9), F_rim = FT(0.2), ρ_rim = FT(500),
            ),
        )
    end

    # Extreme rime-density contrast: ρ_rim far below ρ_i so ρ_g deviates
    # strongly from ρ_i at the graupel/partially-rimed thresholds.
    for ρ_rim in (FT(50), FT(100), FT(150)), F_rim in (FT(0.3), FT(0.5), FT(0.7), FT(0.9))
        push!(
            states,
            build_state(FT;
                ρ = FT(0.85), T = FT(265), q_lcl = FT(3e-4 / 0.85), n_lcl = FT(1e8 / 0.85),
                q_rai = FT(1e-4 / 0.85), n_rai = FT(1e4 / 0.85),
                q_ice = FT(3e-4 / 0.85), n_ice = FT(1e5 / 0.85), F_rim, ρ_rim,
            ),
        )
    end

    # Near-degenerate: trace ice just above the `eps(FT)` has-ice threshold, at
    # both very low and very high number concentration (extreme mean-size tails).
    for (q_ice, n_ice) in ((FT(1e-9), FT(1e-3)), (FT(1e-9), FT(1e2)), (FT(1e-6), FT(1e-2)))
        push!(
            states,
            build_state(FT;
                ρ = FT(1.0), T = FT(250), q_lcl = zero(FT), n_lcl = zero(FT),
                q_rai = zero(FT), n_rai = zero(FT), q_ice, n_ice, F_rim = zero(FT), ρ_rim = FT(400),
            ),
        )
    end

    # Melting-active (above freezing): exercises `ice_melt`'s always-on branch.
    for T in (FT(273.5), FT(275), FT(278))
        push!(
            states,
            build_state(FT;
                ρ = FT(1.0), T, q_lcl = zero(FT), n_lcl = zero(FT),
                q_rai = zero(FT), n_rai = zero(FT), q_ice = FT(1e-3), n_ice = FT(5e4),
                F_rim = FT(0.2), ρ_rim = FT(500),
            ),
        )
    end

    # Unrimed end-member and fully-tiny-cloud-ice (single small regime): explicit
    # bookends of the F_rim and mean-size ranges.
    push!(
        states,
        build_state(FT;
            ρ = FT(1.0), T = FT(240), q_lcl = zero(FT), n_lcl = zero(FT),
            q_rai = zero(FT), n_rai = zero(FT), q_ice = FT(1e-6), n_ice = FT(1e5),
            F_rim = zero(FT), ρ_rim = FT(400),
        ),
    )
    push!(
        states,
        build_state(FT;
            ρ = FT(1.1), T = FT(268), q_lcl = FT(1e-4 / 1.1), n_lcl = FT(1e8 / 1.1),
            q_rai = zero(FT), n_rai = zero(FT), q_ice = FT(8e-3), n_ice = FT(50.0),
            F_rim = FT(0.95), ρ_rim = FT(800),
        ),
    )

    return filter(s -> realizable(mp, s), states)
end

"""
    cast_state(FT, s)

Cast a `build_state`-format `NamedTuple`'s fields to `FT`, preserving the exact
physical scenario (same prognostic values, only the numeric type changes).
"""
cast_state(::Type{FT}, s) where {FT} = (;
    ρ = FT(s.ρ), T = FT(s.T), q_tot = FT(s.q_tot),
    q_lcl = FT(s.q_lcl), n_lcl = FT(s.n_lcl),
    q_rai = FT(s.q_rai), n_rai = FT(s.n_rai),
    q_ice = FT(s.q_ice), n_ice = FT(s.n_ice),
    q_rim = FT(s.q_rim), b_rim = FT(s.b_rim),
)

"""
    generate_comprehensive_states(FT; n_lhs = 300, seed = 1)

Full battery: [`lhs_states`](@ref) (systematic sweep) plus [`corner_states`](@ref)
(explicit physically-important corners), all filtered through [`realizable`](@ref).

Always sampled and filtered in `Float64`, then [`cast_state`](@ref) to `FT` — so
`generate_comprehensive_states(Float64; n_lhs, seed)` and
`generate_comprehensive_states(Float32; n_lhs, seed)` return the same physical
scenarios (same state at index `i` in both), differing only in field type.
Sampling and filtering independently per `FT` would let `realizable` reject a
different subset in each precision, misaligning state indices across
precisions.
"""
function generate_comprehensive_states(::Type{FT} = Float64; n_lhs = 300, seed = 1) where {FT}
    probe_mp = CMP.Microphysics2MParams(Float64; with_ice = true, quad = CM.Quadrature.GaussLegendre(Float64, 6))
    states64 = vcat(lhs_states(Float64, probe_mp; n_lhs, seed), corner_states(Float64, probe_mp))
    FT === Float64 && return states64
    return [cast_state(FT, s) for s in states64]
end
