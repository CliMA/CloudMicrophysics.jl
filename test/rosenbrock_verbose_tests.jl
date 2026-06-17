using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import StaticArrays: SVector

# `Verbose(mode)` augments the `RosenbrockAverage` averaged tendency with a
# post-solve per-process attribution. These tests check that (a) the per-process
# instantaneous parts sum to the instantaneous total and (b) the realized
# per-process tendencies plus the clamp correction reconstruct the verbose net to
# the per-substep linear-solve roundoff.

# Reduce a 2M+P3 net-tendency NamedTuple to the eight prognostic-species vector.
net_vec_2m(t) = SVector(
    t.dq_lcl_dt, t.dn_lcl_dt, t.dq_rai_dt, t.dn_rai_dt,
    t.dq_ice_dt, t.dn_ice_dt, t.dq_rim_dt, t.db_rim_dt,
)
# Reduce a 1M net-tendency NamedTuple to the four prognostic-species vector.
net_vec_1m(t) = SVector(t.dq_lcl_dt, t.dq_icl_dt, t.dq_rai_dt, t.dq_sno_dt)

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
    )

    @testset "2M verbose instantaneous parts sum to total ($FT)" begin
        # the per-process decomposition uses the same physics as the summed
        # entry, so the sum over processes equals the full raw tendency
        for r in regimes
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

    @testset "2M verbose attribution reconstructs net ($FT)" begin
        # Σ_p (per-process realized tendency) + clamp-correction == net realized
        # tendency, to the per-substep linear-solve roundoff (relative to the
        # net scale; F32 carries the larger number-species roundoff)
        Δt = FT(60)
        rtol = FT == Float64 ? FT(1e-10) : FT(1e-3)
        for r in regimes, nsub in (1, 4, 16)
            logλ = isnothing(r.logλ) ? consistent_logλ(r.ρ, r.x) : r.logλ
            v = BMT.bulk_microphysics_tendencies(
                BMT.Verbose(BMT.rosenbrock_exact()), BMT.Microphysics2Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., logλ, Δt, nsub,
            )
            net = net_vec_2m(v)
            recon = SVector((sum(values(v.processes)) + v.clamp_correction)...)
            @test all(isfinite, recon)
            scale = maximum(abs.(net)) + eps(FT)
            @test maximum(abs.(recon - net)) ≤ rtol * scale
        end
    end
end

function test_rosenbrock_verbose_1m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)
    T_frz = TDI.T_freeze(tps)

    # x = [q_lcl, q_icl, q_rai, q_sno]
    regimes = (
        (; ρ = FT(1.0), T = T_frz + FT(17), q_tot = FT(0.018),
            x = FT[2e-3, 0, 5e-4, 0]),  # warm rain
        (; ρ = FT(1.2), T = T_frz + FT(5), q_tot = FT(0.012),
            x = FT[5e-4, 2e-4, 3e-4, 3e-4]),  # mixed warm
        (; ρ = FT(1.2), T = T_frz - FT(10), q_tot = FT(0.012),
            x = FT[3e-4, 5e-4, 2e-4, 4e-4]),  # mixed cold
        (; ρ = FT(1.2), T = T_frz - FT(15), q_tot = FT(0.008),
            x = FT[0, 0, 0, 1e-4]),  # snow cold deposition
    )

    @testset "1M verbose instantaneous parts sum to total ($FT)" begin
        for r in regimes
            x = SVector{4, FT}(r.x...)
            g = BMT.Raw1MTendency(mp, tps, r.ρ, r.T, r.q_tot)
            gv = BMT.Verbose1MTendency(mp, tps, r.ρ, r.T, r.q_tot)
            full = SVector(g(x)...)
            psum = SVector(sum(values(gv(x)))...)
            @test all(isfinite, psum)
            @test psum == full
        end
    end

    @testset "1M verbose net equals non-verbose net ($FT)" begin
        # rosenbrock_donor() carries no limiter, so its verbose net is the same
        # unlimited solve net as the non-verbose mode
        Δt = FT(20)
        for r in regimes, nsub in (1, 4, 16)
            v = BMT.bulk_microphysics_tendencies(
                BMT.Verbose(BMT.rosenbrock_donor()), BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, nsub,
            )
            nv = BMT.bulk_microphysics_tendencies(
                BMT.rosenbrock_donor(), BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, nsub,
            )
            @test net_vec_1m(v) == net_vec_1m(nv)
        end
    end

    @testset "1M verbose attribution reconstructs net ($FT)" begin
        Δt = FT(20)
        rtol = FT == Float64 ? FT(1e-10) : FT(1e-4)
        for r in regimes, nsub in (1, 4, 16)
            v = BMT.bulk_microphysics_tendencies(
                BMT.Verbose(BMT.rosenbrock_donor()), BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, nsub,
            )
            net = net_vec_1m(v)
            recon = SVector((sum(values(v.processes)) + v.clamp_correction)...)
            @test all(isfinite, recon)
            scale = maximum(abs.(net)) + eps(FT)
            @test maximum(abs.(recon - net)) ≤ rtol * scale
        end
    end
end

test_rosenbrock_verbose_2m(Float64)
test_rosenbrock_verbose_2m(Float32)
test_rosenbrock_verbose_1m(Float64)
test_rosenbrock_verbose_1m(Float32)

# JET report-freedom on the verbose entries (compiler-version sensitive, like the
# other perf/inference assertions; see rosenbrock_mode_tests.jl).
if VERSION >= v"1.12"
    import JET
    @testset "verbose entries inference" begin
        for FT in (Float64, Float32)
            tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
            # 2M
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
            # 1M
            mp1 = CMP.Microphysics1MParams(FT)
            rep1 = JET.report_call(
                BMT.bulk_microphysics_tendencies,
                typeof.((
                    BMT.Verbose(BMT.rosenbrock_donor()), BMT.Microphysics1Moment(), mp1, tps,
                    FT(1.2), FT(278), FT(0.012),
                    FT(5e-4), FT(2e-4), FT(3e-4), FT(3e-4), FT(20), 4,
                )),
            )
            @test isempty(JET.get_reports(rep1))
        end
    end
end
