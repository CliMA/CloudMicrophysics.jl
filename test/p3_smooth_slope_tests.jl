import Test: @testset, @test
import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import ForwardDiff as FD

# P3 parameters with the smoothed slope law; κ is read from the repo-local override TOML.
function smooth_slope_params(::Type{FT}) where {FT}
    override = joinpath(pkgdir(CM), "src", "parameters", "toml", "P3_smooth_slope.toml")
    toml_dict = CP.create_toml_dict(FT; override_file = override)
    return CMP.ParametersP3(toml_dict; slope_law = :smooth_powerlaw)
end

function test_smooth_slope_law(FT)
    params = smooth_slope_params(FT)
    slope = params.slope
    # Name the hard law explicitly. Reading the unqualified default here was a latent bug: it
    # silently became the SMOOTH law when :smooth_powerlaw became the default, so the sharp-limit
    # testset compared smooth(κ=50) against smooth(κ=2.68) while calling the latter "hard".
    hard = CMP.ParametersP3(FT; slope_law = :powerlaw).slope
    (; a, b, c, μ_max, κ) = slope
    logλs = FT.(range(2, 17; length = 400))

    @testset "bounds and interior limit [FT=$FT]" begin
        for logλ in logλs
            μ = P3.get_μ(slope, logλ)
            @test isfinite(μ)
            @test μ ≤ μ_max + eps(FT)
            @test μ ≥ -sqrt(eps(FT))
            μh = P3.get_μ(hard, logλ)
            if FT(0.5) < μh < μ_max - FT(0.5)
                @test isapprox(μ, μh; atol = 3 / κ)
            end
        end
    end

    @testset "monotone in logλ [FT=$FT]" begin
        μprev = P3.get_μ(slope, logλs[1])
        for logλ in logλs[2:end]
            μ = P3.get_μ(slope, logλ)
            @test μ ≥ μprev - sqrt(eps(FT))
            μprev = μ
        end
    end

    @testset "sharp limit converges to hard clamp [FT=$FT]" begin
        for κ_big in FT.((50, 200, 1000))
            sharp = CMP.SmoothSlopePowerLaw(; a, b, c, μ_max, κ = κ_big)
            for logλ in logλs
                @test isapprox(
                    P3.get_μ(sharp, logλ), P3.get_μ(hard, logλ); atol = 5 / κ_big,
                )
            end
        end
    end

    @testset "ForwardDiff derivative finite and nonnegative [FT=$FT]" begin
        for logλ in logλs
            d = FD.derivative(x -> P3.get_μ(slope, x), logλ)
            @test isfinite(d)
            @test d ≥ -sqrt(eps(FT))
            h = sqrt(eps(FT))
            fd = (P3.get_μ(slope, logλ + h) - P3.get_μ(slope, logλ - h)) / (2h)
            @test isapprox(d, fd; atol = FT(1e-2), rtol = FT(1e-2))
        end
    end

    @testset "shape solver is single-valued [FT=$FT]" begin
        N_ice = FT(1e8)
        F_rims = FT.((0, 0.2, 0.5, 0.8, 0.95))
        ρ_rims = FT.((100, 400, 800))
        logλ_trues = FT.(range(5, 12; length = 40))
        for F_rim in F_rims, ρ_rim in ρ_rims, logλ_true in logλ_trues
            dummy = P3.P3State(params, FT(1), FT(1), F_rim, ρ_rim)
            x_ice = exp(P3.logLdivN(dummy, logλ_true))
            state = P3.P3State(params, x_ice * N_ice, N_ice, F_rim, ρ_rim)
            # Skip degenerate states with mass below the representable minimum.
            isfinite(P3.get_distribution_logλ(state)) || continue
            sols = P3.get_distribution_logλ_all_solutions(state)
            @test length(sols) == 1
            @test isapprox(only(sols), logλ_true; atol = FT(0.1))
        end
    end
end

@testset "P3 SmoothSlopePowerLaw" begin
    for FT in (Float32, Float64)
        test_smooth_slope_law(FT)
    end
end
