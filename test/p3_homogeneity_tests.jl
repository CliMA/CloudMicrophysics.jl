using Test

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI

include("p3_state_battery.jl")

# A2: the homogeneity-degree battery. For a handful of ice-phase process rates, holding
# every input fixed except the ice number `n_ice`, `rate(s * n_ice) = s^d * rate(n_ice)`
# is an EXACT identity of the code as written, not an approximation - the mechanism is
# that within one `p3_2m_process_rates` evaluation the ice size-distribution shape
# (`logλ`, `μ`) is a FROZEN argument rather than re-derived from `(q_ice, n_ice)`
# (`logλ` is a substep input; `μ` comes from `state`, whose fields - `F_rim = L_rim /
# L_ice`, `ρ_rim = B_rim / L_rim` - do not reference `N_ice` at all: see
# `P3.state_from_prognostic` / `P3.get_μ`'s call sites). A rate formed as a linear
# functional of the ice number-density distribution `N'(D)` at fixed shape is therefore
# exactly linear in `N_ice` (`N'(D) = N_ice * (a shape function of D alone)`), and a rate
# quadratic in `N'(D)` (a self-collision integral) is exactly quadratic in `N_ice`.
#
# This battery does NOT attempt joint (mass, number) proportional-scaling homogeneity
# claims (holding mean particle size fixed) for the warm-rain processes: several of
# those closures gate or weight on the vapor deficit `q_vap = q_tot - q_lcl - q_rai -
# q_ice`, which is not invariant under scaling a single donor species, so "the known
# scaling degree" is not a single clean number there without deeper derivation this file
# has not done. The table below is scoped to what is directly traceable to source: two
# entries (`ice_aggregation.n_ice`, degree 2; the ice-number sublimation pathway inside
# `ice_depsub.n_ice`, degree 2) are quoted almost verbatim from `_jacobian_2mp3_manual`'s
# own docstring in `src/BMT_2mp3_jacobian.jl` (the "Tier 3" and ice-number-sublimation
# commentary); `ice_depsub.q_ice`'s degree 1 follows from that same docstring's
# "1/τ_i ∝ n_ice" statement applied to the bare relaxation rate; the `ice_melting` row is
# this file's own extension of the same linear-functional argument to `P3.ice_melt`, NOT
# quoted anywhere in the source, and is flagged lower confidence below. Per the battery
# spec, this table "encodes claims about the code" and needs Haakon's review before it is
# trusted as a merge gate; until then it documents its own reasoning inline so the review
# can check each row against the claim it makes.
#
# | process        | slot (n_ice degree) | degree | confidence                          |
# |----------------|----------------------|--------|--------------------------------------|
# | ice_aggregation | n_ice                | 2      | documented (Jacobian file, verbatim) |
# | ice_depsub      | q_ice                | 1      | derived from a documented fact       |
# | ice_depsub      | n_ice                | 2      | derived from a documented fact       |
# | ice_melting     | q_rai, n_rai, q_ice, n_ice | 1 | reasoned by analogy, NOT quoted     |
# | ice_melting     | q_rim, b_rim          | 0      | reasoned by analogy, NOT quoted      |
#
# Only `n_ice` scaling is exercised, and only upward (`SCALE > 1`), which is the safe
# direction for both gates this file re-checks per state rather than assumes:
#   - `P3.ice_population_is_present` (guards `ice_aggregation` and `ice_melting`) is
#     `(ρq_ice / m_nuc > ρn_ice) & (ρn_ice > 0)` - NOT monotone in `n_ice`, so this file
#     checks it explicitly at both the base and the scaled `n_ice` rather than assuming
#     "present at the base stays present when scaled up".
#   - `P3.ice_deposition_is_degenerate` (guards `ice_depsub`) tests `τ_dep` against its
#     upper cap, and `1/τ_dep ∝ n_ice`, so scaling `n_ice` up only moves `τ_dep` away
#     from the cap; checking non-degenerate at the base is therefore sufficient.
const SCALE = 3.0

FT_LOCKED = Float64  # the degree table is locked at Float64 only, per the battery spec;
# Float32 is exercised only by the finiteness sweep below, which needs no degree claim.

function context_at(mp, s, ::Type{FT}) where {FT}
    ρ = FT(s.ρ)
    state = P3.state_from_prognostic(
        mp.ice.scheme, ρ * FT(s.q_ice), ρ * FT(s.n_ice), ρ * FT(s.q_rim), ρ * FT(s.b_rim),
    )
    logλ = P3.get_distribution_logλ(state)
    return state, logλ
end

function pp_at(mp, tps, s, ::Type{FT}, logλ; n_ice_scale = one(FT)) where {FT}
    micro = (;
        q_tot = FT(s.q_tot), q_lcl = FT(s.q_lcl), n_lcl = FT(s.n_lcl),
        q_rai = FT(s.q_rai), n_rai = FT(s.n_rai),
        q_ice = FT(s.q_ice), n_ice = FT(s.n_ice) * n_ice_scale,
        q_rim = FT(s.q_rim), b_rim = FT(s.b_rim),
    )
    thermo = (; ρ = FT(s.ρ), T = FT(s.T), w = zero(FT), p = zero(FT), logλ)
    return BMT.p3_2m_process_rates(mp, tps, micro, thermo)
end

"""
    check_degree!(counts, key, base, scaled, s, d; rtol)

Assert `scaled ≈ s^d * base` and bump `counts[key]`, unless both values are negligible
(no meaningful check to make at this state), in which case nothing is asserted.
"""
function check_degree!(counts, key, base, scaled, s, d; rtol)
    FT = typeof(base)
    if abs(base) < eps(FT) && abs(scaled) < eps(FT)
        return nothing
    end
    @test isapprox(scaled, FT(s)^d * base; rtol = FT(rtol))
    counts[key] = get(counts, key, 0) + 1
    return nothing
end

@testset "A2 homogeneity degrees (n_ice, $FT_LOCKED)" begin
    FT = FT_LOCKED
    mp = CMP.Microphysics2MParams(
        FT; with_ice = true, is_limited = true, aerosol = CMP.PrescribedAerosol(FT),
    )
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    states = generate_comprehensive_states(FT; n_lhs = 80, seed = 1)
    @test length(states) > 0

    counts = Dict{Symbol, Int}()
    for s in states
        (s.q_ice > 0 && s.n_ice > 0) || continue
        state0, logλ = context_at(mp, s, FT)
        pp0 = pp_at(mp, tps, s, FT, logλ)
        ppS = pp_at(mp, tps, s, FT, logλ; n_ice_scale = FT(SCALE))

        # ice_aggregation, ice_melting: gated on `ice_population_is_present` at BOTH
        # the base and the scaled n_ice (the predicate is not monotone in n_ice)
        state_scaled = P3.state_from_prognostic(
            mp.ice.scheme, FT(s.ρ) * FT(s.q_ice), FT(s.ρ) * FT(s.n_ice) * FT(SCALE),
            FT(s.ρ) * FT(s.q_rim), FT(s.ρ) * FT(s.b_rim),
        )
        if P3.ice_population_is_present(state0) && P3.ice_population_is_present(state_scaled)
            check_degree!(
                counts, :ice_aggregation_n_ice,
                pp0.ice_aggregation.n_ice, ppS.ice_aggregation.n_ice, SCALE, 2; rtol = 1.0e-6,
            )
            for (slot, d) in ((:q_rai, 1), (:n_rai, 1), (:q_ice, 1), (:n_ice, 1))
                check_degree!(
                    counts, Symbol(:ice_melting_, slot),
                    getfield(pp0.ice_melting, slot), getfield(ppS.ice_melting, slot), SCALE, d;
                    rtol = 1.0e-6,
                )
            end
            for slot in (:q_rim, :b_rim)  # degree 0: `melt_frac` is n_ice-invariant
                check_degree!(
                    counts, Symbol(:ice_melting_, slot),
                    getfield(pp0.ice_melting, slot), getfield(ppS.ice_melting, slot), SCALE, 0;
                    rtol = 1.0e-6,
                )
            end
        end

        # ice_depsub: gated on the deposition timescale not sitting at its degeneracy
        # cap at the (unscaled) base state
        τ_i0 = P3.ice_deposition_timescale(
            mp.ice.terminal_velocity, mp.warm_rain.air_properties, tps,
            FT(s.T), FT(s.ρ), state0, logλ; quad = mp.ice.quad,
        )
        if !P3.ice_deposition_is_degenerate(τ_i0)
            check_degree!(
                counts, :ice_depsub_q_ice,
                pp0.ice_depsub.q_ice, ppS.ice_depsub.q_ice, SCALE, 1; rtol = 1.0e-6,
            )
            check_degree!(
                counts, :ice_depsub_n_ice,
                pp0.ice_depsub.n_ice, ppS.ice_depsub.n_ice, SCALE, 2; rtol = 1.0e-6,
            )
        end
    end

    # every declared (process, slot) pair must actually have been exercised by at least
    # one battery state, or the loop above's gating silently made the whole table vacuous
    for key in (
        :ice_aggregation_n_ice, :ice_depsub_q_ice, :ice_depsub_n_ice,
        :ice_melting_q_rai, :ice_melting_n_rai, :ice_melting_q_ice, :ice_melting_n_ice,
        :ice_melting_q_rim, :ice_melting_b_rim,
    )
        @test get(counts, key, 0) > 0
    end
end

@testset "A2 finiteness under single-species scaling ($FT2)" for FT2 in (Float32, Float64)
    mp = CMP.Microphysics2MParams(
        FT2; with_ice = true, is_limited = true, aerosol = CMP.PrescribedAerosol(FT2),
    )
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT2)
    states = generate_comprehensive_states(FT2; n_lhs = 40, seed = 2)
    species = (:q_lcl, :n_lcl, :q_rai, :n_rai, :q_ice, :n_ice, :q_rim, :b_rim)
    scales = (FT2(0), FT2(0.1), FT2(2), FT2(10))

    for s in states
        state, logλ = context_at(mp, s, FT2)
        base = (;
            q_tot = FT2(s.q_tot), q_lcl = FT2(s.q_lcl), n_lcl = FT2(s.n_lcl),
            q_rai = FT2(s.q_rai), n_rai = FT2(s.n_rai),
            q_ice = FT2(s.q_ice), n_ice = FT2(s.n_ice),
            q_rim = FT2(s.q_rim), b_rim = FT2(s.b_rim),
        )
        thermo = (; ρ = FT2(s.ρ), T = FT2(s.T), w = zero(FT2), p = zero(FT2), logλ)
        for sp in species, scale in scales
            micro = merge(base, (; sp => getfield(base, sp) * scale))
            pp, rs = BMT.p3_2m_process_rates(mp, tps, micro, thermo)
            @test all(isfinite, sum(values(pp)))
            @test isfinite(rs.f_shd)
        end
    end
end
