import Test as TT
import CloudMicrophysics.DistributionTools as DT

TT.@testset "DistributionTools" begin
    TT.@testset "Generalized Gamma Distribution" begin
        # Test parameters
        ν = 2.0
        μ = 3.0
        B = 2.0

        # Test CDF and quantile correspondence
        for Y in [0.1, 0.25, 0.5, 0.75, 0.9]
            x = DT.generalized_gamma_quantile(ν, μ, B, Y)
            p = DT.generalized_gamma_cdf(ν, μ, B, x)
            TT.@test p ≈ Y rtol = 1e-10
        end

        # Test edge cases
        TT.@test DT.generalized_gamma_cdf(ν, μ, B, 0) == 0
        TT.@test DT.generalized_gamma_cdf(ν, μ, B, -1) == 0

        # Invalid parameters return NaN (device-kernel-safe; no DomainError)
        TT.@test isnan(DT.generalized_gamma_cdf(ν, -1, B, 1))  # μ < 0
        TT.@test isnan(DT.generalized_gamma_cdf(ν, μ, -1, 1))  # B < 0
        # A non-positive B (degenerate distribution) is the fractional-power site
        # that aborted GPU kernels; it must return NaN, not raise.
        TT.@test isnan(DT.generalized_gamma_quantile(ν, μ, -1.0, 0.5))  # B < 0
        TT.@test isnan(DT.generalized_gamma_Mⁿ(ν, μ, -1.0, 1.0, 2))     # B < 0
    end

    TT.@testset "Exponential Distribution" begin
        # Test parameters
        D_mean = 2.0

        # Test CDF and quantile correspondence
        for Y in [0.1, 0.25, 0.5, 0.75, 0.9]
            D = DT.exponential_quantile(D_mean, Y)
            p = DT.exponential_cdf(D_mean, D)
            TT.@test p ≈ Y rtol = 1e-10
        end

        # Test edge cases
        TT.@test DT.exponential_cdf(D_mean, 0) == 0
        TT.@test DT.exponential_cdf(D_mean, -1) == 0
        TT.@test DT.exponential_cdf(D_mean, Inf) ≈ 1 rtol = 1e-10

        # Invalid parameters return NaN (device-kernel-safe; no DomainError)
        TT.@test isnan(DT.exponential_cdf(-1, 1))  # D_mean < 0
        TT.@test isnan(DT.exponential_quantile(D_mean, -0.1))  # Y < 0
        TT.@test isnan(DT.exponential_quantile(D_mean, 1.1))   # Y > 1
        TT.@test isnan(DT.exponential_quantile(-1, 0.5))     # D_mean < 0
    end
end
nothing
