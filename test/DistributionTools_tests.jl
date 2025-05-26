using Test
import CloudMicrophysics.DistributionTools as DT

@testset "DistributionTools" begin
    @testset "Generalized Gamma Distribution" begin
        # Test parameters
        ν = 2.0
        μ = 3.0
        B = 1.0

        # Test CDF and quantile correspondence
        for Y in [0.1, 0.25, 0.5, 0.75, 0.9]
            x = DT.generalized_gamma_quantile(ν, μ, B, Y)
            p = DT.generalized_gamma_cdf(ν, μ, B, x)
            @test p ≈ Y rtol = 1e-10
        end

        # Test edge cases
        @test DT.generalized_gamma_cdf(ν, μ, B, 0) == 0
        @test DT.generalized_gamma_cdf(ν, μ, B, -1) == 0

        # Test error handling
        @test_throws DomainError DT.generalized_gamma_cdf(ν, -1, B, 1)  # μ < 0
        @test_throws DomainError DT.generalized_gamma_cdf(ν, μ, -1, 1)  # B < 0
    end

    @testset "Exponential Distribution" begin
        # Test parameters
        D_mean = 2.0

        # Test CDF and quantile correspondence
        for Y in [0.1, 0.25, 0.5, 0.75, 0.9]
            D = DT.exponential_quantile(D_mean, Y)
            p = DT.exponential_cdf(D_mean, D)
            @test p ≈ Y rtol = 1e-10
        end

        # Test edge cases
        @test DT.exponential_cdf(D_mean, 0) == 0
        @test DT.exponential_cdf(D_mean, -1) == 0
        @test DT.exponential_cdf(D_mean, Inf) ≈ 1 rtol = 1e-10

        # Test error handling
        @test_throws DomainError DT.exponential_cdf(-1, 1)  # D_mean < 0
        @test_throws DomainError DT.exponential_quantile(D_mean, -0.1)  # Y < 0
        @test_throws DomainError DT.exponential_quantile(D_mean, 1.1)   # Y > 1
        @test_throws DomainError DT.exponential_quantile(-1, 0.5)     # D_mean < 0
    end
end
