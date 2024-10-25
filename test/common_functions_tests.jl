import Test as TT

import ClimaParams
import Thermodynamics as TD
import CloudMicrophysics as CM
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP

@info "Common Functions Tests"

TT.@testset "logistic_function unit tests" begin

    TT.@test CO.logistic_function(-1.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function(0.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function(1.0, 1.0, 2.0) == 0.5
    TT.@test CO.logistic_function(2.0, 1.0, 2.0) ≈ 0.9525 atol = 1e-4

    TT.@test CO.logistic_function(1.0, 0.0, 2.0) == 1.0
    TT.@test CO.logistic_function(0.0, 0.0, 2.0) == 0.0
    TT.@test_throws AssertionError("x_0 >= 0") CO.logistic_function(
        1.0,
        -1.0,
        2.0,
    )
    TT.@test_throws AssertionError("k > 0") CO.logistic_function(1.0, 1.0, 0.0)

end

TT.@testset "logistic_function_integral unit tests" begin

    TT.@test CO.logistic_function_integral(-1.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function_integral(0.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function_integral(1.0, 1.0, 2.0) ≈ 0.3115 atol = 1e-4
    TT.@test CO.logistic_function_integral(3.0, 1.0, 2.0) ≈ 2.0 atol = 1e-2

    TT.@test CO.logistic_function_integral(1.0, 0.0, 2.0) == 1.0
    TT.@test CO.logistic_function_integral(0.0, 0.0, 2.0) == 0.0
    TT.@test_throws AssertionError("x_0 >= 0") CO.logistic_function_integral(
        1.0,
        -1.0,
        2.0,
    )
    TT.@test_throws AssertionError("k > 0") CO.logistic_function_integral(
        1.0,
        1.0,
        0.0,
    )

end

function test_H2SO4_soln_saturation_vapor_pressure(FT)

    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

    TT.@testset "H2SO4 solution saturated vapor pressure" begin

        T_warm = FT(225.0)
        T_cold = FT(200.0)
        T_too_warm = FT(240)
        T_too_cold = FT(180)
        x_sulph = FT(0.1)

        # If T out of range
        TT.@test_throws AssertionError("T < T_max") CO.H2SO4_soln_saturation_vapor_pressure(
            H2SO4_prs,
            x_sulph,
            T_too_warm,
        )
        TT.@test_throws AssertionError("T > T_min") CO.H2SO4_soln_saturation_vapor_pressure(
            H2SO4_prs,
            x_sulph,
            T_too_cold,
        )

        # p_sol should be higher at warmer temperatures
        TT.@test CO.H2SO4_soln_saturation_vapor_pressure(
            H2SO4_prs,
            x_sulph,
            T_warm,
        ) > CO.H2SO4_soln_saturation_vapor_pressure(
            H2SO4_prs,
            x_sulph,
            T_cold,
        )
    end
end

function test_a_w_xT(FT)

    tps = TD.Parameters.ThermodynamicsParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

    TT.@testset "a_w_xT" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph_low = FT(0.06)
        x_sulph_high = FT(0.1)

        # a_w greater at warmer temperatures
        for x_sulph in [x_sulph_high, x_sulph_low]
            TT.@test CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_cold) <
                     CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_warm)
        end

        # a_w greater at lower sulphuric acid concentration
        for T in [T_warm, T_cold]
            TT.@test CO.a_w_xT(H2SO4_prs, tps, x_sulph_high, T) <
                     CO.a_w_xT(H2SO4_prs, tps, x_sulph_low, T)
        end
    end
end

function test_a_w_eT(FT)

    tps = TD.Parameters.ThermodynamicsParameters(FT)

    TT.@testset "a_w_eT" begin

        T_warm = FT(285)
        T_cold = FT(251)
        e_high = FT(1088)
        e_low = FT(544)

        # a_w greater at higher altitudes
        TT.@test CO.a_w_eT(tps, e_low, T_cold) > CO.a_w_eT(tps, e_high, T_warm)

        # a_w greater at greater partial pressures
        for T in [T_warm, T_cold]
            TT.@test CO.a_w_eT(tps, e_low, T) < CO.a_w_eT(tps, e_high, T)
        end
    end
end

function test_a_w_ice(FT)

    tps = TD.Parameters.ThermodynamicsParameters(FT)

    TT.@testset "a_w_ice" begin

        T_warm = FT(240)
        T_cold = FT(230)

        # a_w greater at warmer temperatures
        TT.@test CO.a_w_ice(tps, T_cold) < CO.a_w_ice(tps, T_warm)

    end
end

function test_Chen_coefficients(FT)
    ρ = FT(1.2)
    tol = 10 * eps(FT)
    Ch2022 = CMP.Chen2022VelType(FT)
    snow = CMP.Snow(FT)
    ice = CMP.CloudIce(FT)

    TT.@testset "Chen terminal velocity rain (B1)" begin
        aiu, bi, ciu = CO.Chen2022_vel_coeffs_B1(Ch2022.rain, ρ)

        TT.@test all(
            isapprox.(
                aiu,
                [
                    FT(286768.02047954104),
                    FT(-1.6916433443360287e6),
                    FT(9843.240767655458),
                ],
                rtol = tol,
            ),
        )
        TT.@test all(
            isapprox.(
                bi,
                [FT(2.249342), FT(2.249342), FT(1.098942)],
                rtol = tol,
            ),
        )
        TT.@test all(
            isapprox.(ciu, [FT(0), FT(184.325), FT(184.325)], rtol = tol),
        )
    end

    TT.@testset "Chen terminal velocity small ice (B2)" begin
        aiu, bi, ciu = CO.Chen2022_vel_coeffs_B2(Ch2022.small_ice, ρ, ice.ρᵢ)

        TT.@test all(
            isapprox.(
                aiu,
                [380.9523818928577, -378.94964461870484],
                rtol = tol,
            ),
        )
        TT.@test all(
            isapprox.(bi, [0.7065511618279822, 0.7065511618279822], rtol = tol),
        )
        TT.@test all(isapprox.(ciu, [0.0, 5194.870484289714], rtol = tol))
    end

    TT.@testset "Chen terminal velocity large ice (B4)" begin
        aiu, bi, ciu = CO.Chen2022_vel_coeffs_B4(Ch2022.large_ice, ρ, snow.ρᵢ)

        TT.@test all(
            isapprox.(
                aiu,
                [144.32735209844674, -0.9732523067775996],
                rtol = tol,
            ),
        )
        TT.@test all(
            isapprox.(
                bi,
                [0.5375612365666844, 0.020917783512726773],
                rtol = tol,
            ),
        )
        TT.@test all(isapprox.(ciu, [0.0, 86.47212246687042], rtol = tol))
    end
end

println("Testing Float64")
test_H2SO4_soln_saturation_vapor_pressure(Float64)
test_a_w_xT(Float64)
test_a_w_eT(Float64)
test_a_w_ice(Float64)
test_Chen_coefficients(Float64)

println("Testing Float32")
test_H2SO4_soln_saturation_vapor_pressure(Float32)
test_a_w_xT(Float32)
test_a_w_eT(Float32)
test_a_w_ice(Float32)
test_Chen_coefficients(Float32)
