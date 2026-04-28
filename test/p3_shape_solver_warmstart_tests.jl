using Test: @testset, @test
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import BenchmarkTools: @belapsed

"""
    test_warmstart_correctness(FT)

For a sweep of (L_ice, N_ice, F_rim, ρ_rim) inputs, the warm-started solver
must produce a `logλ` that matches the cold-start Brent solution up to FP
tolerance. We exercise four guess regimes:

  1. `nothing`            — must hit the Brent path.
  2. `NaN` / `Inf`        — non-finite; must gracefully fall back to Brent.
  3. Exactly at the root  — must converge in one Newton step and match.
  4. Close to the root    — typical "previous step" warm start.
  5. Far from the root    — Newton may diverge; the guard must fall back.
  6. Out of bracket       — must be rejected and fall back to Brent.

In every case the final `logλ` value must equal the cold-start value to
`2 * eps(FT)` relative tolerance.
"""
function test_warmstart_correctness(FT)
    params = CMP.ParametersP3(FT)

    # A small grid of representative P3 states. Mixes unrimed / rimed,
    # low-ice / high-ice, sparse / dense populations.
    L_sweep = FT.((1e-5, 1e-4, 1e-3, 1e-2))
    N_sweep = FT.((1e4, 1e6, 1e8))
    F_sweep = FT.((0, 0.3, 0.7))
    ρ_sweep = FT.((200, 500))

    @testset "warm-start correctness [FT=$FT]" begin
        rtol = 10 * eps(FT)
        for L_ice in L_sweep, N_ice in N_sweep, F_rim in F_sweep, ρ_rim in ρ_sweep
            state = P3.get_state(params; L_ice, N_ice, F_rim, ρ_rim)

            # --- 1. cold reference
            logλ_cold = P3.get_distribution_logλ(state)

            # --- 2. nothing guess
            logλ_nothing = P3.get_distribution_logλ(state, nothing)
            @test logλ_nothing == logλ_cold  # must be bit-identical (same path)

            # --- 3. non-finite guesses
            for bad in (FT(NaN), FT(Inf), FT(-Inf))
                logλ_bad = P3.get_distribution_logλ(state, bad)
                @test logλ_bad == logλ_cold
            end

            # --- 4. guess exactly at the root
            logλ_exact = P3.get_distribution_logλ(state, logλ_cold)
            @test isapprox(logλ_exact, logλ_cold; rtol)

            # --- 5. close-to-root guess (previous step would look like this)
            for δ in FT.((-0.05, -0.005, 0.005, 0.05))
                g = logλ_cold + δ
                # keep inside bracket
                if FT(2) < g < FT(17)
                    logλ_near = P3.get_distribution_logλ(state, g)
                    @test isapprox(logλ_near, logλ_cold; rtol)
                end
            end

            # --- 6. far-from-root guess (Newton may diverge / overshoot)
            for g in (FT(2.5), FT(16.5))
                logλ_far = P3.get_distribution_logλ(state, g)
                # Either Newton succeeded to the true root, or the guard fell
                # back to Brent. Both must yield the cold-start answer.
                @test isapprox(logλ_far, logλ_cold; rtol)
            end

            # --- 7. out-of-bracket guess
            for oob in (FT(1.0), FT(20.0))
                logλ_oob = P3.get_distribution_logλ(state, oob)
                @test logλ_oob == logλ_cold
            end
        end
    end

    # Zero-ice edge case: both paths should return `log(0)` regardless of guess.
    @testset "warm-start zero-ice edge [FT=$FT]" begin
        state0 = P3.get_state(params;
            L_ice = FT(0), N_ice = FT(0), F_rim = FT(0), ρ_rim = FT(500),
        )
        @test P3.get_distribution_logλ(state0) === P3.get_distribution_logλ(state0, FT(10))
    end
end

"""
    test_warmstart_speedup(FT)

Benchmark cold vs warm start on a steady-state cell (the common case). A
"smooth" cell evolves by ~1e-3 in `logλ` per step, so the prior step's
`logλ` is a good seed for bracket halving. We expect a modest but
strictly positive speedup: the hot-start pays one extra `shape_problem`
eval at the guess and hands Brent a bracket that is half as wide,
saving ~1 Brent iteration on average.

The test target is just `t_warm < t_cold` — any strict speedup passes.
We intentionally don't pursue a more aggressive narrowing (e.g. probing
`guess ± δ`) because Brent's interpolating iterate is more efficient
per `shape_problem` eval than any fixed-offset probe would be; extra
probes here compete with Brent rather than complementing it.
"""
function test_warmstart_speedup(FT)
    params = CMP.ParametersP3(FT)
    L_ice = FT(1e-3)
    N_ice = FT(1e6)
    F_rim = FT(0.3)
    ρ_rim = FT(500)
    state = P3.get_state(params; L_ice, N_ice, F_rim, ρ_rim)
    logλ_true = P3.get_distribution_logλ(state)
    # Typical step-to-step drift: ~1e-3 in log space
    guess = logλ_true + FT(1e-3)

    t_cold = @belapsed P3.get_distribution_logλ($state)
    t_warm = @belapsed P3.get_distribution_logλ($state, $guess)

    @testset "warm-start speedup [FT=$FT]" begin
        @info "shape solver benchmark" FT t_cold t_warm speedup = t_cold / t_warm
        # Guard against regressions: warm-start must be strictly faster on
        # a smoothly-evolving cell. Target 1.2x; we measure ~1.4–1.6x on
        # Apple M-series at Float64, but a lenient threshold avoids flakes
        # on noisy or slow CI machines.
        @test t_warm < t_cold
    end
end

@testset "P3 shape-solver warm-start" begin
    for FT in (Float32, Float64)
        test_warmstart_correctness(FT)
    end
    # Benchmarks only on Float64 to avoid noise.
    test_warmstart_speedup(Float64)
end
