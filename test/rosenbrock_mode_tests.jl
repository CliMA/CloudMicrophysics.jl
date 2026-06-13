using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import StaticArrays: SVector

# `RosenbrockAverage` substeps the raw instantaneous pointwise 2M+P3 tendency
# with a linearized-implicit (Rosenbrock-Euler) update. The reference for
# accuracy tests is a finely-resolved forward-Euler integration of the same
# raw tendency (identical T update and frozen logλ/q_tot semantics).

function explicit_reference(mp, tps, ρ, T, q_tot, x0, logλ, Δt, nsub)
    FT = typeof(q_tot)
    h = Δt / FT(nsub)
    x = SVector{8, FT}(x0...)
    Tsub = T
    cp_d = TDI.TD.Parameters.cp_d(tps)
    for _ in 1:nsub
        g = BMT.Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        f = g(x)
        xp = x
        x = max.(x .+ h .* f, zero(FT))
        Ts = max(FT(150), Tsub)
        Tsub +=
            (
                TDI.Lᵥ(tps, Ts) * ((x[1] - xp[1]) + (x[3] - xp[3])) +
                TDI.Lₛ(tps, Ts) * (x[5] - xp[5])
            ) / cp_d
    end
    return x
end

function test_rosenbrock_mode(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme

    function consistent_logλ(ρ, x)
        st = P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8])
        return P3.get_distribution_logλ(st)
    end
    function step(x0, ρ, T, q_tot, logλ, Δt, nsub)
        t = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics2Moment(), mp, tps,
            ρ, T, q_tot, x0..., logλ, Δt, nsub,
        )
        return SVector{8, FT}(x0...) .+ Δt .* SVector(values(t)...)
    end

    # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
    regimes = (
        (; name = "warm rain", ρ = FT(1.05), T = FT(288), q_tot = FT(0.015),
            x = FT[4e-4, 8e7, 2.1e-3, 5e4, 0, 0, 0, 0], logλ = FT(-Inf), tol = 0.002),
        (; name = "mixed phase", ρ = FT(0.78), T = FT(273.5), q_tot = FT(0.009),
            x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8], logλ = nothing, tol = 0.07),
        (; name = "ice sublimation", ρ = FT(0.45), T = FT(253), q_tot = FT(4e-4),
            x = FT[0, 0, 0, 0, 8e-4, 5e5, 5e-4, 9e-7], logλ = nothing, tol = 0.012),
    )

    # Per-species scales with physical floors: collapsing species (reference
    # depleted to ~0) would otherwise dominate a plain relative metric. The
    # rime channels (q_rim, b_rim) carry genesis-scale floors: once ice has
    # melted away, single-precision leaves an O(1e-7) coupled volume/mass
    # residual that is physically negligible but unbounded against an
    # eps-scale floor.
    floors = SVector{8, FT}(1e-9, 1e-1, 1e-9, 1e-1, 1e-9, 1e-1, 1e-8, 1e-6)
    err_metric(x, x_ref, x0) = maximum(abs.(x .- x_ref) ./ (abs.(x0) .+ abs.(x_ref) .+ floors))

    @testset "RosenbrockAverage vs fine explicit reference ($FT)" begin
        Δt = FT(10)
        for r in regimes
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            x_ref = explicit_reference(mp, tps, r.ρ, r.T, r.q_tot, r.x, logλ, Δt, 2048)
            x0 = SVector{8, FT}(r.x...)
            err(x) = err_metric(x, x_ref, x0)
            errs = [err(step(r.x, r.ρ, r.T, r.q_tot, logλ, Δt, n)) for n in (1, 4, 16)]
            @test all(isfinite, errs)
            # accuracy improves under substep refinement...
            @test errs[3] ≤ max(errs[1], FT(1e-3)) * (1 + sqrt(eps(FT)))
            # ...to within a regime-calibrated distance of the reference
            @test errs[3] < (FT == Float64 ? r.tol : 2 * r.tol)
        end
    end

    @testset "degenerate and trivial states ($FT)" begin
        # all-zero state: degenerate gate -> explicit substeps -> exactly zero
        t0 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics2Moment(), mp, tps,
            FT(1), FT(273), FT(0),
            FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0), FT(0),
            FT(-Inf), FT(60), 4,
        )
        @test all(iszero, values(t0))
        # nsub defaults to 1 and accepts the trailing-argument form
        r = regimes[2]
        logλ = consistent_logλ(r.ρ, r.x)
        t1 = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics2Moment(), mp, tps,
            r.ρ, r.T, r.q_tot, r.x..., logλ, FT(60),
        )
        @test all(isfinite, values(t1))
    end

    @testset "substeps stay finite and non-negative ($FT)" begin
        # hostile stress state: strongly supersaturated over liquid at 233 K
        # (fast condensation-freezing cascade), near-empty rain alongside ice
        x_stress = FT[1e-6, 1e6, 1e-12, 1e-2, 8e-4, 5e5, 5e-4, 9e-7]
        logλ = consistent_logλ(FT(0.45), x_stress)
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            t = BMT.bulk_microphysics_tendencies(
                BMT.RosenbrockAverage(), BMT.Microphysics2Moment(), mp, tps,
                FT(0.45), FT(233), FT(0.003), x_stress..., logλ, Δt, nsub,
            )
            x1 = SVector{8, FT}(x_stress...) .+ Δt .* SVector(values(t)...)
            @test all(isfinite, x1)
            # non-negative up to the roundoff of the host-side x + Δt * t
            # reconstruction (internally the state is floored at zero)
            tol = eps(FT) .* (abs.(SVector{8, FT}(x_stress...)) .+ Δt .* abs.(SVector(values(t)...)))
            @test all(x1 .>= -tol)
        end
    end

    @testset "near-empty channels take the explicit path ($FT)" begin
        # condensed masses in (eps, 1e-10) produce finite but enormous
        # Jacobian rows from the steep dependence of the raw process rates
        # on a near-zero channel; the linearized update fabricates phantom
        # droplet number that substep refinement cannot remove. The channel
        # mask must route such channels to forward Euler: the result has to
        # track the explicit reference. spin-up from a band state is a fast
        # transient, so same-nsub trajectories legitimately differ at first
        # order; the bug signature is droplet number orders of magnitude
        # beyond the fine reference (~2.6e5 here; the unmasked Jacobian gave
        # ~7e9)
        x_band = FT[1e-13, 1e2, 0, 0, 0, 0, 0, 0]
        Δt = FT(60)
        x_ref = explicit_reference(mp, tps, FT(1), FT(288), FT(0.02), x_band, FT(-Inf), Δt, 2048)
        x16 = step(x_band, FT(1), FT(288), FT(0.02), FT(-Inf), Δt, 16)
        @test all(isfinite, x16)
        @test x16[2] < 10 * x_ref[2]
        # a band-mass channel alongside a healthy ice channel must not
        # poison the implicit update of the healthy channels either
        x_mixb = FT[1e-13, 1e2, 0, 0, 8e-4, 5e5, 5e-4, 9e-7]
        logλ_m = consistent_logλ(FT(0.45), x_mixb)
        xm = step(x_mixb, FT(0.45), FT(253), FT(4e-4), logλ_m, FT(10), 16)
        @test all(isfinite, xm)
        @test xm[2] < FT(1e6)  # no phantom droplet number orders of magnitude above the reference
    end
end

test_rosenbrock_mode(Float64)
test_rosenbrock_mode(Float32)

# Allocation check on the hot call (compiler-version sensitive, like the other
# perf assertions; see performance_tests.jl)
if VERSION >= v"1.12"
    @testset "RosenbrockAverage allocations" begin
        FT = Float64
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
        p3 = mp.ice.scheme
        st = P3.state_from_prognostic(
            p3,
            FT(0.78) * FT(1e-4),
            FT(0.78) * FT(2e5),
            FT(0.78) * FT(4e-5),
            FT(0.78) * FT(6e-8),
        )
        logλ = P3.get_distribution_logλ(st)
        call() = BMT.bulk_microphysics_tendencies(
            BMT.RosenbrockAverage(), BMT.Microphysics2Moment(), mp, tps,
            FT(0.78), FT(273.5), FT(0.009),
            FT(2e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8),
            logλ, FT(60), 4,
        )
        call()
        # On <= 1.11 the differentiated path boxes wholesale (the known
        # inference-depth limit behind the other >= 1.12 perf gates), hence
        # the version gate on this testset.
        @test (@allocated call()) == 0
    end
end
