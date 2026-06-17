using Test

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI
import StaticArrays: SVector

# The unified `RosenbrockAverage{Jacobian, GrowthTreatment, TendencyLimiter}`
# framework: presets (`rosenbrock_donor`, `rosenbrock_coupled`,
# `rosenbrock_exact`), the keyword constructor, the `LinearizedAverage` ≡ donor
# equivalence on the 1M model, and the 2M+P3 `ExactJacobian`-only contract.

net_vec_1m(t) = SVector(t.dq_lcl_dt, t.dq_icl_dt, t.dq_rai_dt, t.dq_sno_dt)

function test_framework_1m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)
    T_frz = TDI.T_freeze(tps)

    # x = [q_lcl, q_icl, q_rai, q_sno]
    regimes = (
        (; ρ = FT(1.0), T = T_frz + FT(17), q_tot = FT(0.018), x = FT[2e-3, 0, 5e-4, 0]),  # warm rain
        (; ρ = FT(1.2), T = T_frz + FT(5), q_tot = FT(0.012), x = FT[5e-4, 2e-4, 3e-4, 3e-4]),  # mixed warm
        (; ρ = FT(1.2), T = T_frz - FT(10), q_tot = FT(0.012), x = FT[3e-4, 5e-4, 2e-4, 4e-4]),  # mixed cold
    )

    @testset "rosenbrock_donor() ≡ LinearizedAverage() ($FT)" begin
        for r in regimes, nsub in (1, 4, 16)
            Δt = FT(20)
            donor = BMT.bulk_microphysics_tendencies(
                BMT.rosenbrock_donor(), BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, nsub,
            )
            lin = BMT.bulk_microphysics_tendencies(
                BMT.LinearizedAverage(), BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, nsub,
            )
            @test net_vec_1m(donor) == net_vec_1m(lin)
        end
    end

    @testset "keyword constructor matches rosenbrock_donor() ($FT)" begin
        kw = BMT.RosenbrockAverage(
            jacobian = BMT.DonorJacobian(),
            growth = BMT.ImplicitGrowth(),
            limiter = BMT.NoLimiter(),
        )
        @test kw == BMT.rosenbrock_donor()
        for r in regimes
            Δt = FT(20)
            t_kw = BMT.bulk_microphysics_tendencies(
                kw, BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, 4,
            )
            t_preset = BMT.bulk_microphysics_tendencies(
                BMT.rosenbrock_donor(), BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, 4,
            )
            @test net_vec_1m(t_kw) == net_vec_1m(t_preset)
        end
    end

    @testset "1M presets give finite tendencies ($FT)" begin
        presets = (BMT.rosenbrock_donor(), BMT.rosenbrock_coupled(), BMT.rosenbrock_exact())
        for mode in presets, r in regimes, nsub in (1, 2, 8), Δt in (FT(20), FT(120))
            t = BMT.bulk_microphysics_tendencies(
                mode, BMT.Microphysics1Moment(), mp, tps,
                r.ρ, r.T, r.q_tot, r.x..., Δt, nsub,
            )
            @test all(isfinite, net_vec_1m(t))
        end
    end
end

# The unified `RosenbrockAverage` framework on the two-moment + P3 model:
# `rosenbrock_exact()` gives finite tendencies, and a non-Exact
# `RosenbrockAverage` throws (only `ExactJacobian` is supported there).

function test_framework_2m(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme

    consistent_logλ(ρ, x) =
        P3.get_distribution_logλ(P3.state_from_prognostic(p3, ρ * x[5], ρ * x[6], ρ * x[7], ρ * x[8]))

    # x = [q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim]
    ρ = FT(0.78)
    T = FT(273.5)
    q_tot = FT(0.009)
    x = FT[2e-4, 5e7, 1e-4, 4e4, 1e-4, 2e5, 4e-5, 6e-8]
    logλ = consistent_logλ(ρ, x)

    @testset "rosenbrock_exact() on Microphysics2Moment works ($FT)" begin
        for nsub in (1, 2, 8), Δt in (FT(60), FT(300))
            t = BMT.bulk_microphysics_tendencies(
                BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps,
                ρ, T, q_tot, x..., logλ, Δt, nsub,
            )
            @test all(isfinite, values(t))
        end
    end

    @testset "non-Exact RosenbrockAverage on Microphysics2Moment throws ($FT)" begin
        non_exact = (
            BMT.RosenbrockAverage(BMT.DonorJacobian(), BMT.ImplicitGrowth(), BMT.NoLimiter()),
            BMT.RosenbrockAverage(BMT.CoupledDonorJacobian(), BMT.ImplicitGrowth(), BMT.NoLimiter()),
        )
        for mode in non_exact
            @test_throws ArgumentError BMT.bulk_microphysics_tendencies(
                mode, BMT.Microphysics2Moment(), mp, tps,
                ρ, T, q_tot, x..., logλ, FT(60), 4,
            )
        end
    end
end

test_framework_1m(Float64)
test_framework_1m(Float32)
test_framework_2m(Float64)
test_framework_2m(Float32)
