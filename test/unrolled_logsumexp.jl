

using Test
using Random
import LogExpFunctions: logsumexp
import CloudMicrophysics.Utilities: unrolled_logsumexp

Random.seed!(42)

@testset "logsumexp NTuple{4}" begin

    @testset "random Float64 values" begin
        for _ = 1:100
            for N in 4:16:84
                t = Tuple(randn(N))
                @test @inferred(unrolled_logsumexp(t)) ≈ logsumexp(t)
            end
        end
    end

    @testset "random Float32 values" begin
        for _ = 1:100
            for N in 4:16:84
                t = Tuple(randn(Float32, N))
                @test @inferred(unrolled_logsumexp(t)) ≈ logsumexp(t)
            end
        end
    end

    @testset "large values (overflow regime)" begin
        t = (1000.0, 1001.0, 999.0, 998.0)
        @test unrolled_logsumexp(t) ≈ logsumexp(t)

        t2 = (1e3, 1e3, 1e3, 1e3)
        @test unrolled_logsumexp(t2) ≈ logsumexp(t2)
    end

    @testset "small values (underflow regime)" begin
        t = (-1000.0, -1001.0, -999.0, -998.0)
        @test unrolled_logsumexp(t) ≈ logsumexp(t)

        t2 = (-1e3, -1e3, -1e3, -1e3)
        @test unrolled_logsumexp(t2) ≈ logsumexp(t2)
    end

    @testset "all equal values" begin
        for v in [-5.0, 0.0, 1.0, 100.0]
            t = (v, v, v, v)
            @test unrolled_logsumexp(t) ≈ logsumexp(t)
            @test unrolled_logsumexp(t) ≈ v + log(4)
        end
    end

    @testset "mixed sign values" begin
        t = (-10.0, 5.0, -3.0, 8.0)
        @test unrolled_logsumexp(t) ≈ logsumexp(t)
    end

    @testset "NaN propagation" begin
        @test isnan(unrolled_logsumexp((NaN, 1.0, 2.0, 3.0)))
        @test isnan(unrolled_logsumexp((1.0, NaN, 2.0, 3.0)))
        @test isnan(unrolled_logsumexp((1.0, 2.0, NaN, 3.0)))
        @test isnan(unrolled_logsumexp((1.0, 2.0, 3.0, NaN)))
        @test isnan(unrolled_logsumexp((NaN, NaN, NaN, NaN)))
    end

    @testset "Inf handling" begin
        @test unrolled_logsumexp((Inf, 1.0, 2.0, 3.0)) == Inf
        @test unrolled_logsumexp((1.0, Inf, 2.0, 3.0)) == Inf
        @test unrolled_logsumexp((Inf, Inf, 1.0, 2.0)) == Inf
    end

    @testset "-Inf handling" begin
        @test unrolled_logsumexp((-Inf, -Inf, -Inf, -Inf)) == -Inf
        t = (-Inf, 0.0, -Inf, -Inf)
        @test unrolled_logsumexp(t) ≈ logsumexp(t)
    end

    @testset "+Inf -Inf handling" begin
        @test unrolled_logsumexp((Inf, -Inf)) == Inf
        @test unrolled_logsumexp((-Inf, Inf)) == Inf
        @test unrolled_logsumexp((Inf, -Inf, 1.0, 2.0)) == Inf
        @test unrolled_logsumexp((-Inf, Inf, 1.0, 2.0)) == Inf
    end

    @testset "NaN with Inf" begin
        @test isnan(unrolled_logsumexp((NaN, Inf, 1.0, 2.0)))
        @test isnan(unrolled_logsumexp((NaN, -Inf, 1.0, 2.0)))
    end

    @testset "type promotion (Int input)" begin
        t = (1, 2, 3, 4)
        result = @inferred unrolled_logsumexp(t)
        @test result ≈ logsumexp(t)
        @test result isa Float64
    end

    @testset "type stability" begin
        @test @inferred(unrolled_logsumexp((1.0f0, 2.0f0, 3.0f0, 4.0f0))) isa Float32
        @test @inferred(unrolled_logsumexp((1.0, 2.0, 3.0, 4.0))) isa Float64
    end
end
