using Test

import JET
import BenchmarkTools as BT

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.ThermodynamicsInterface as TDI

# The unified `RosenbrockAverage{Jacobian, GrowthTreatment, TendencyLimiter}`
# framework on the two-moment + P3 model: `rosenbrock_exact()` gives finite
# tendencies, and a non-Exact `RosenbrockAverage` throws (only `ExactJacobian`
# is supported there).

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

    @testset "RosenbrockAverage on warm-rain-only parameters throws ($FT)" begin
        mp_warm = CMP.Microphysics2MParams(FT; with_ice = false)
        @test_throws "requires P3 ice parameters" BMT.bulk_microphysics_tendencies(
            BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp_warm, tps,
            ρ, T, q_tot, x..., logλ, FT(60), 4,
        )
    end
end

function _framework_exact_args(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics2MParams(FT; with_ice = true, is_limited = true)
    p3 = mp.ice.scheme
    st = P3.state_from_prognostic(
        p3, FT(0.78) * FT(1e-4), FT(0.78) * FT(2e5), FT(0.78) * FT(4e-5), FT(0.78) * FT(6e-8),
    )
    logλ = P3.get_distribution_logλ(st)
    return (
        BMT.rosenbrock_exact(), BMT.Microphysics2Moment(), mp, tps,
        FT(0.78), FT(273.5), FT(0.009),
        FT(2e-4), FT(5e7), FT(1e-4), FT(4e4), FT(1e-4), FT(2e5), FT(4e-5), FT(6e-8),
        logλ, FT(60), 4,
    )
end

function test_framework_exact_inference(FT)
    args = _framework_exact_args(FT)
    @testset "rosenbrock_exact() inference and allocations ($FT)" begin
        @test (@inferred BMT.bulk_microphysics_tendencies(args...)) isa NamedTuple
        JET.@test_opt BMT.bulk_microphysics_tendencies(args...)
        trail = BT.@benchmark $(splat(BMT.bulk_microphysics_tendencies))($args) samples = 100 evals = 1
        @test trail.memory == 0
    end
end

test_framework_2m(Float64)
test_framework_2m(Float32)

test_framework_exact_inference(Float64)
test_framework_exact_inference(Float32)
