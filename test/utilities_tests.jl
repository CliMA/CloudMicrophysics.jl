using Test
import ForwardDiff as FD
import CloudMicrophysics.Utilities as UT

# `rime_density` and `rime_mass_fraction` must agree between `Float32` and `Float64`
# to ordinary single-precision rounding across the physical parameter range, not to a
# precision-dependent regularization artifact.
@testset "rime_density / rime_mass_fraction: Float32-Float64 consistency" begin
    @testset "a state with a small, proportionally consistent (q_rim, b_rim) pair" begin
        q_rim, b_rim = 1.2176073047297676e-4, 1.2132270859560276e-6
        ρ_rim64 = UT.rime_density(q_rim, b_rim)
        ρ_rim32 = UT.rime_density(Float32(q_rim), Float32(b_rim))
        @test isapprox(Float64(ρ_rim32), ρ_rim64; rtol = 1e-5)
    end

    @testset "sweep over a physically plausible (F_rim, ρ_rim) grid" begin
        for q_ice in (1e-6, 1e-4, 1e-2), F_rim in (0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99),
            ρ_rim_true in (50.0, 100.0, 400.0, 700.0, 900.0)

            q_rim = F_rim * q_ice
            b_rim = q_rim / ρ_rim_true
            ρ_rim64 = UT.rime_density(q_rim, b_rim)
            ρ_rim32 = UT.rime_density(Float32(q_rim), Float32(b_rim))
            relerr = abs(Float64(ρ_rim32) - ρ_rim64) / ρ_rim64
            @test relerr < 1e-4
        end
    end

    @testset "sweep over a physically plausible F_rim grid (rime_mass_fraction)" begin
        for q_ice in (1e-6, 1e-4, 1e-2), F_rim_true in (0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99)
            q_rim = F_rim_true * q_ice
            F_rim64 = UT.rime_mass_fraction(q_rim, q_ice)
            F_rim32 = UT.rime_mass_fraction(Float32(q_rim), Float32(q_ice))
            relerr = abs(Float64(F_rim32) - F_rim64) / F_rim64
            @test relerr < 1e-4
        end
    end

    @testset "a genuinely negligible rime/ice pair still regularizes to zero" begin
        for FT in (Float64, Float32)
            @test isfinite(UT.rime_density(FT(0), FT(0)))
            @test UT.rime_density(FT(0), FT(0)) == 0
            @test isfinite(UT.rime_mass_fraction(FT(0), FT(0)))
            @test UT.rime_mass_fraction(FT(0), FT(0)) == 0
        end
    end
end

# A physically admissible (q_rim, b_rim) pair must not be suppressed to zero merely
# because `b_rim` is small, as long as the pair sits above `rime_density`'s cutoff.
# Only a pair genuinely below that cutoff regularizes to zero.
@testset "rime_density: admissible pairs are not suppressed at small b_rim" begin
    @testset "a specific admissible state at the rho_min boundary" begin
        q_rim64, b_rim64 = 3.1287404e-19, 1.961572e-21
        @test isapprox(UT.rime_density(q_rim64, b_rim64), 159.5; rtol = 1e-4)
        q_rim32, b_rim32 = Float32(q_rim64), Float32(b_rim64)
        @test isapprox(UT.rime_density(q_rim32, b_rim32), 159.5f0; rtol = 1.0f-4)
    end

    @testset "sweep: b_rim spans small scales, q_rim proportionally small" begin
        for FT in (Float64, Float32)
            for ρ_rim_true in FT[159.5, 300.0, 916.7], scale in FT[1e-10, 1e-15, 1e-20]
                b_rim = scale
                q_rim = ρ_rim_true * scale
                ρ_rim = UT.rime_density(q_rim, b_rim)
                @test isfinite(ρ_rim)
                @test isapprox(ρ_rim, ρ_rim_true; rtol = FT === Float64 ? 1e-8 : 1.0f-4)
            end
        end
    end

    @testset "a pair genuinely below the presence cutoff still regularizes to zero" begin
        for FT in (Float64, Float32)
            b_below = FT(1e-30) / 2
            @test UT.rime_density(FT(1e-10), b_below) == 0
        end
    end
end

@testset "guarded_quotient" begin
    for FT in (Float64, Float32)
        # A `Dual` whose value is zero and whose partials are not: the value decides
        # presence, while the `Dual` comparison reads the same number as present.
        rate = FD.Dual{Nothing}(FT(0), FT(1), FT(0))
        q = FD.Dual{Nothing}(FT(0), FT(0), FT(1))
        @test !(FD.value(q) > zero(FD.value(q)))
        @test q > zero(q)

        f = UT.guarded_quotient(rate, q)
        @test iszero(FD.value(f))
        @test all(iszero, FD.partials(f))

        # The exact quotient, with finite partials, where the denominator is present.
        qp = FD.Dual{Nothing}(FT(2), FT(0), FT(1))
        rp = FD.Dual{Nothing}(FT(3), FT(1), FT(0))
        @test FD.value(UT.guarded_quotient(rp, qp)) == FT(1.5)
        @test all(isfinite, FD.partials(UT.guarded_quotient(rp, qp)))

        # `absent` is returned unchanged, and defaults to zero.
        @test UT.guarded_quotient(rate, q, oftype(FD.value(q), Inf)) == FT(Inf)
        @test UT.guarded_quotient(FT(1), FT(0)) == 0
    end
end

@testset "nearest_admissible_b" begin
    for FT in (Float64, Float32)
        ρ_min, ρ_max, q = FT(50), FT(900), FT(1e-4)

        # A ratio already in range returns the trial value itself.
        b_ok = q / FT(400)
        @test UT.nearest_admissible_b(q, b_ok, ρ_min, ρ_max) === b_ok

        # A ratio out of range returns the nearer endpoint.
        @test UT.nearest_admissible_b(q, q / FT(2000), ρ_min, ρ_max) ≈ q / ρ_max
        @test UT.nearest_admissible_b(q, q / FT(10), ρ_min, ρ_max) ≈ q / ρ_min

        # No mass pairs with no volume.
        @test UT.nearest_admissible_b(zero(FT), b_ok, ρ_min, ρ_max) == 0

        for ρ_trial in FT[1, 49, 50, 400, 900, 2000]
            b = UT.nearest_admissible_b(q, q / ρ_trial, ρ_min, ρ_max)
            @test q / b >= ρ_min * (1 - 4 * eps(FT))
            @test q / b <= ρ_max * (1 + 4 * eps(FT))
        end
    end
end
nothing
