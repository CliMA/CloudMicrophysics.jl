using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import StaticArrays: SVector

# `Verbose(mode)` augments the `RosenbrockAverage` averaged tendency with a post-solve
# per-process attribution. These tests check that (a) the per-process instantaneous
# parts sum to the instantaneous total, (b) the realized per-process tendencies plus the
# correction reconstruct the verbose net to the per-substep linear-solve roundoff, (c) the
# verbose net equals the same mode's own non-verbose net, and (d) the correction is at
# roundoff, not merely small, on a state constructed so the floor cannot have engaged.
#
# PORTED FROM the campaign's `test/rosenbrock_verbose_tests.jl` (`campaign/he/p3-density-
# floors` in this project's CM clone). Ported here as this design's own tests, not a
# mechanical carry-over: the campaign's `Verbose` was a SEPARATE substep implementation
# ("previously the verbose path ran its own substep loop that omitted the limiter, the
# water bound, the projection and the shape refresh, and refused the manual Jacobian
# outright" - `src/BMT_diagnostics.jl`'s own account of what it replaced), so the
# campaign's "verbose net equals non-verbose net" claim only held for the one artificially
# unlimited, donor-free mode it was restricted to, and needed its own bespoke test to say
# so. The rebuilt design runs ONE substep body for every mode, parameterized by a record
# sink; production's sink is a no-op that compiles away. `Verbose`/`Trace`/production are
# therefore bit-identical BY CONSTRUCTION, for every mode, not by two implementations
# happening to agree - so the net-equality identity below is asserted across several
# modes, not the one the campaign was restricted to, and would be a genuinely stronger
# regression guard than the campaign's own version once it can run (see below).
#
# THE KNOWN-FALSE ASSERTION the team lead named is not carried over: the campaign version
# restricted its net-equality check to one hand-picked unlimited configuration because
# that was the only case where its two SEPARATE implementations happened to agree; a
# general "verbose net equals non-verbose net for any mode" claim did not hold for that
# design and is not what is asserted here. What holds NOW, by construction, is exactly
# that general claim, so it is what is asserted below, across the modes this file already
# exercises rather than one special-cased one.
#
# STALE REFERENCES REPAIRED:
#   - `v.clamp_correction` renamed to `v.correction` (`_verbose_extras`,
#     `src/BMT_diagnostics.jl`: `(; processes, correction, substeps)`).
#   - `BMT.rosenbrock_donor()` (deleted, D4 ruling) replaced with `BMT.RosenbrockAverage()`
#     (bare defaults reproduce the same donor-based configuration) in the throws-test.
#   - The whole 1M section (`test_rosenbrock_verbose_1m`, `BMT.Raw1MTendency`,
#     `BMT.Verbose1MTendency`, `net_vec_1m`) is DROPPED, not repaired, for the same reason
#     as in `rosenbrock_framework_tests.jl`: the `RosenbrockAverage`-on-`Microphysics1Moment`
#     framework it needs does not exist on this branch (`he/p3-cmLUT-tables`) - it lives on
#     the sibling `he/p3-cm10b-1m`, which has not merged onto this trunk. The 1M JET
#     inference check at the bottom is dropped with it.
#
# NOT REPAIRABLE WITHOUT A SOURCE FIX - this is the headline fact about this file, not a
# footnote: `Verbose(mode)`/`Trace(mode)` for `Microphysics2Moment` are NOT WIRED to the
# production substep on this branch. Traced in full (see the report to the team lead):
#   - `bulk_microphysics_tendencies(v.mode, cm, mp, tps, ..., Δt, nsub, w, p; sink)`
#     (`src/BMT_diagnostics.jl`) passes `sink` as a keyword to a call that resolves to
#     `bulk_microphysics_tendencies(mode::RosenbrockAverage{...}, ::Microphysics2Moment, ...)`
#     (`src/BMT_2mp3_march.jl:593-599`, `:626-632`) - NEITHER method accepts a `sink`
#     keyword, so this throws `MethodError` before reaching a single substep.
#   - Even granting the keyword got through: `_rosenbrock_average_2mp3`
#     (`src/BMT_2mp3_march.jl:542`) calls `_march_2mp3` with no sink argument at all, so
#     `_march_2mp3` always uses its own unrelated `sink = nothing` default and its
#     `record!` call passes `(; x, x_prev, T_pre, T_post, logλ, h, diag)` - not the
#     `(; pp, W, h, α_w, α_l, δ_acc)` shape `record!(::VerboseSink, ctx)` destructures.
#     Two independent mismatches, not one; `git show --stat` on the commit that added
#     `BMT_diagnostics.jl` touches only that file plus a one-line include - `_march_2mp3`
#     was never updated to build the richer context or accept the sink.
# So EVERY test below that calls `bulk_microphysics_tendencies(Verbose(mode), ...)` throws
# on its first call, for any input, until `_march_2mp3`/`_rosenbrock_substep_diag` are
# changed to build and thread `_record_context` through. The two tests that do NOT go
# through that path - "2M verbose instantaneous parts sum to total" (calls
# `Verbose2MP3Tendency`/`Instantaneous2MP3Tendency` directly, a primal-level pair
# unrelated to the record-sink machinery) and "the number-without-mass regime
# discriminates" (calls no microphysics tendency function at all) - are NOT affected and
# should pass now. Every other testset in this file is marked individually below.
#
# The unrelated `_jacobian_2mp3_manual` `τ_act` field-access bug that used to additionally
# block `rosenbrock_manual()` (see `rosenbrock_framework_tests.jl`) is FIXED as of the tip
# this file was last checked against - `rosenbrock_manual()`'s inclusion in the modes below
# now depends only on the wiring gap above, not on a second bug. The branch is being
# restacked continuously; re-check before relying on either fact.
#
# AN EXPOSURE GAP, NOW AUTHORIZED TO CLOSE (team lead, on review): `VerboseSink`'s
# per-substep record is `(; pieces, correction, total)` (`record!(sink::VerboseSink, ctx)`,
# `src/BMT_diagnostics.jl`) - it does NOT carry `α_w`/`α_l` (the water bound's and
# limiter's rescale factors) or any other signal of whether the positivity floor /
# rime-pair projection actually engaged on that substep, even though `TraceSink`'s own
# per-substep record (`(; δ_acc, α_w, α_l)`) carries `α_w`/`α_l` for the SAME substep. The
# contract's own docstring claims "`correction` vanishes only where no floor engaged" - a
# testable, meaningful claim - but nothing exposed in a `Verbose` return lets a caller
# confirm INDEPENDENTLY, from the returned data alone, whether a given substep's floor was
# the reason `correction` is nonzero. This was judged a defect in the design, not a
# limitation to test around, and `VerboseSink`'s records are being extended to carry the
# same scale factors `TraceSink` already does.
#
# The testset below ("correction is roundoff on a floor-inactive state") tests the CLAIM
# by construction (a state built far from every threshold) and stays as it stands - it
# is honest about testing the claim rather than the mechanism, and the team lead confirmed
# to keep it. TODO once `VerboseSink`'s records carry the scale factors: add a companion
# assertion that CORRELATES `correction`'s magnitude with the mechanism the contract
# attributes it to directly - e.g. at a state where the floor is deliberately forced to
# engage (a species driven toward zero within the substep), assert the exposed
# floor-engaged signal is set and that `correction` is where the size of it lives, rather
# than inferring engagement from the state's construction alone as the testset below still
# has to.

# Reduce a 2M+P3 net-tendency NamedTuple to the eight prognostic-species vector.
net_vec_2m(t) = SVector(
    t.dq_lcl_dt, t.dn_lcl_dt, t.dq_rai_dt, t.dn_rai_dt,
    t.dq_ice_dt, t.dn_ice_dt, t.dq_rim_dt, t.db_rim_dt,
)

function test_rosenbrock_verbose_2m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme

    consistent_logλ(ρ, x) =
        P3.get_distribution_logλ(P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))

    # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
    regimes = (
        (; ρ = FT(1.05), T = FT(288), q_tot = FT(0.015),  # warm rain
            x = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0], logλ = FT(-Inf)),
        (; ρ = FT(0.78), T = FT(273.5), q_tot = FT(0.009),  # mixed phase
            x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8], logλ = nothing),
        (; ρ = FT(0.45), T = FT(253), q_tot = FT(4e-4),  # ice sublimation
            x = FT[0, 0, 0, 0, 8e-4, 5e5, 5e-4, 9e-7], logλ = nothing),
        # ORPHAN MASS in all three species: mass present, number exactly zero. Every species in
        # the three regimes above is either fully populated or fully empty, so no orphan
        # predicate can fire in any of them, and this assertion passed without ever exercising
        # the drains - presence is not activity, in a test. It therefore could not see that the
        # drains lived on the per-process path and not on the instantaneous one, which is the
        # divergence this regime exists to catch.
        #
        # Subsaturated on purpose: the drains are identically zero at or above saturation, so a
        # saturated orphan state would reproduce the same vacuity one level down.
        (; ρ = FT(0.9), T = FT(265), q_tot = FT(1.0e-3),  # orphan mass, subsaturated
            x = FT[1e-4, 0, 5e-5, 0, 1e-4, 0, 3e-5, 5e-8], logλ = nothing),
        # NUMBER WITHOUT MASS in SUPERSATURATED air - an orphan's mirror image, and the one corner
        # none of the five regimes above reaches. `number_bounded_by_mass_limits` consults its
        # `sat_excess` argument in exactly one branch, mass absent, where it keeps the number if
        # the air is supersaturated and drops it to zero if it is not. So this is the only kind of
        # state at which passing that argument and omitting it disagree, and the disagreement
        # propagates into every warm-rain rate built from the bounded droplet number.
        #
        # Both halves of the condition are load-bearing: q_lcl = 0 selects the branch, and the
        # SUPERSATURATION is what makes the two answers differ within it. A subsaturated version
        # of this state agrees on both paths and would be one more vacuous regime.
    )

    # KEPT OUT OF `regimes`: mass exactly zero with number present is precisely where the
    # substep drivers' documented difference used to bite under the campaign's separate-
    # implementation design. The sum-to-total invariant below is INSTANTANEOUS - one
    # evaluation, no substepping - so it is unaffected either way, and it is the invariant
    # this regime exists to make discriminating.
    numberless = (; ρ = FT(1.0), T = FT(290), q_tot = FT(0.015),
        x = FT[0, 1e6, 0, 0, 0, 0, 0, 0], logλ = FT(-Inf))
    sum_regimes = (regimes..., numberless)

    # THE SIXTH REGIME VERIFIES ITSELF, because picking numbers that happen to land in a branch is
    # how a regime silently stops discriminating. Both halves are asserted: the air must be
    # SUPERSATURATED (a subsaturated twin agrees on both paths and would be vacuous), and the two
    # bounding calls must actually DISAGREE there (the property the regime exists to exercise).
    #
    # NOT affected by the record-sink wiring gap: calls no `bulk_microphysics_tendencies`
    # at all.
    @testset "the number-without-mass regime discriminates ($FT)" begin
        r = numberless
        sx = BMT._liquid_sat_excess(tps, r.ρ, r.T, r.q_tot, r.x[1], r.x[3], r.x[5])
        @test sx > 0                              # supersaturated, or the regime is vacuous
        @test r.x[1] == 0 && r.x[2] > 0           # mass absent, number present: the branch
        xc = (; x_min = mp.warm_rain.seifert_beheng.pdf_c.xc_min,
            x_max = mp.warm_rain.seifert_beheng.pdf_c.xc_max)
        with = CM.Microphysics2M.number_bounded_by_mass_limits(xc, r.x[1], r.x[2], sx; invent_from_zero = false)
        without = CM.Microphysics2M.number_bounded_by_mass_limits(xc, r.x[1], r.x[2]; invent_from_zero = false)
        @test with != without                     # the two bounding calls genuinely differ here
    end

    # NOT affected by the record-sink wiring gap: `Verbose2MP3Tendency` and
    # `Instantaneous2MP3Tendency` are both primal-level callables (one evaluation of
    # `p3_2m_process_rates`, no substep, no sink) unrelated to `Verbose(mode)`'s
    # substep-averaged attribution below.
    @testset "2M verbose instantaneous parts sum to total ($FT)" begin
        for r in sum_regimes
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            x = SVector{8, FT}(r.x...)
            g = BMT.Instantaneous2MP3Tendency(mp, tps, r.ρ, r.T, r.q_tot, logλ)
            gv = BMT.Verbose2MP3Tendency(mp, tps, r.ρ, r.T, r.q_tot, logλ)
            full = SVector(g(x)...)
            psum = SVector(sum(values(gv(x)))...)
            @test all(isfinite, psum)
            @test psum == full
        end
    end

    # KNOWN TO CURRENTLY THROW (record-sink wiring gap, see file header): `Verbose(mode)`
    # never reaches a substep, so `v.processes`/`v.correction` never get built.
    @testset "2M verbose attribution reconstructs net ($FT)" begin
        # Σ_p (per-process realized tendency) + correction == net realized tendency, to
        # the per-substep linear-solve roundoff (relative to the net scale; F32 carries
        # the larger number-species roundoff)
        Δt = FT(60)
        rtol = FT == Float64 ? FT(1e-10) : FT(1e-3)
        for mode in (BMT.rosenbrock_exact(), BMT.rosenbrock_manual()), r in regimes, nsub in (1, 4, 16)
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            v = BMT.bulk_microphysics_tendencies(
                BMT.Verbose(mode), BMT.Microphysics2Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., logλ, Δt, nsub,
            )
            net = net_vec_2m(v)
            recon = SVector((sum(values(v.processes)) + v.correction)...)
            @test all(isfinite, recon)
            scale = maximum(abs.(net)) + eps(FT)
            @test maximum(abs.(recon - net)) ≤ rtol * scale
        end
    end

    # KNOWN TO CURRENTLY THROW (record-sink wiring gap, see file header). The contract's own
    # claim ("`correction` vanishes only where no floor engaged" - `VerboseSink`'s
    # docstring): at a state comfortably far from every threshold the floor and rime-pair
    # projection can act on (every species well clear of zero, admissible rime density
    # already, no near-saturation limiter edge), `correction` should be at roundoff, not
    # merely small relative to the net. See the file header's account of why this is tested
    # by CONSTRUCTING such a state rather than by reading a signal that confirms the floor
    # was inactive - the attribution contract exposes no such signal.
    @testset "correction is roundoff on a floor-inactive state ($FT)" begin
        r = (; ρ = FT(1.0), T = FT(280), q_tot = FT(0.02),
            x = FT[2e-3, 5e7, 1e-3, 2e4, 2e-3, 2e5, 4e-4, 1e-6], logλ = nothing)
        logλ = consistent_logλ(r.ρ, r.x)
        Δt = FT(2)  # a short step, so the increment is small relative to every species above
        # RELATIVE to the sum being formed, not absolute. `correction` is
        # `total - sum(values(pieces))`, so its roundoff scale is the magnitude of those
        # summands, and this state spans eleven orders of magnitude across species, from
        # `1e-6` to `5e7`. A single absolute tolerance over all of them asks the number
        # concentrations to cancel to a part in 1e17, far below their own representation
        # error, which is what the two failures here were: exactly 2^-32 and 2^-30. The
        # testset's own name and the comment above it say the claim is that the correction is
        # roundoff, and roundoff is relative by definition.
        rtol = FT == Float64 ? FT(1e-10) : FT(1e-4)
        for mode in (BMT.rosenbrock_exact(), BMT.rosenbrock_manual())
            v = BMT.bulk_microphysics_tendencies(
                BMT.Verbose(mode), BMT.Microphysics2Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., logλ, Δt, 1,
            )
            scale = sum(abs.(p) for p in values(v.processes))
            @test maximum(abs.(v.correction) ./ max.(scale, floatmin(FT))) ≤ rtol
        end
    end

    # KNOWN TO CURRENTLY THROW (record-sink wiring gap, see file header). This is the
    # identity the team lead named as what should replace the campaign's narrower,
    # one-mode-only claim: `Verbose`, `Trace` and production run the same substep body
    # with only the sink differing, so the verbose net equals the SAME mode's own
    # non-verbose net BY CONSTRUCTION - not just for one hand-picked unlimited
    # configuration, which is why it is checked here across several modes, including
    # ones WITH their own limiter (`rosenbrock_exact()`) and the manual Jacobian
    # (`rosenbrock_manual()`), neither of which the campaign's version dared include.
    @testset "2M verbose net equals non-verbose net ($FT)" begin
        modes = (
            BMT.rosenbrock_exact(),
            BMT.rosenbrock_manual(),
            BMT.RosenbrockAverage(BMT.ExactJacobian(), BMT.ExplicitGrowthDiagonal(), BMT.NoLimiter()),
        )
        Δt = FT(60)
        for mode in modes, r in regimes, nsub in (1, 4)
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            v = BMT.bulk_microphysics_tendencies(
                BMT.Verbose(mode), BMT.Microphysics2Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., logλ, Δt, nsub,
            )
            nv = BMT.bulk_microphysics_tendencies(
                mode, BMT.Microphysics2Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., logλ, Δt, nsub,
            )
            @test net_vec_2m(v) == net_vec_2m(nv)
        end
    end

    # KNOWN TO CURRENTLY THROW (record-sink wiring gap, see file header) - the throw this
    # test expects is `ArgumentError` from the mode-mismatch catch-all
    # (`src/BMT_2mp3_march.jl:554`), reached BEFORE the sink keyword would matter, so this
    # one test's own assertion is unaffected by the wiring gap; it is marked anyway
    # because it is easy to misread as passing "for the right reason" once the gap is
    # fixed elsewhere, when actually the ArgumentError path never depends on the sink.
    @testset "2M verbose on a non-Exact Jacobian throws ($FT)" begin
        r = regimes[2]
        logλ = consistent_logλ(r.ρ, r.x)
        @test_throws ArgumentError BMT.bulk_microphysics_tendencies(
            BMT.Verbose(BMT.RosenbrockAverage()), BMT.Microphysics2Moment(), mp, tps,
            r.ρ, r.T, r.q_tot, r.x..., logλ, FT(60), 4,
        )
    end
end

test_rosenbrock_verbose_2m(Float64)
test_rosenbrock_verbose_2m(Float32)

# JET report-freedom on the verbose entry (compiler-version sensitive, like the other
# perf/inference assertions; see rosenbrock_mode_tests.jl). KNOWN TO CURRENTLY THROW
# before JET even gets to analyze it (record-sink wiring gap, see file header) - kept as
# the check this entry should pass once wired, not weakened.
if VERSION >= v"1.12"
    import JET
    @testset "verbose entry inference" begin
        for FT in (Float64, Float32)
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
            p3 = mp.ice.scheme
            st = P3.state_from_prognostic(
                p3, FT(0.78) * FT(1e-4), FT(0.78) * FT(2e5),
                FT(0.78) * FT(4e-5), FT(0.78) * FT(6e-8),
            )
            logλ = P3.get_distribution_logλ(st)
            rep2 = JET.report_call(
                BMT.bulk_microphysics_tendencies,
                typeof.((
                    BMT.Verbose(BMT.rosenbrock_exact()), BMT.Microphysics2Moment(), mp, tps,
                    FT(0.78), FT(273.5), FT(0.009),
                    FT(2e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8),
                    logλ, FT(60), 4,
                )),
            )
            @test isempty(JET.get_reports(rep2))
        end
    end
end
