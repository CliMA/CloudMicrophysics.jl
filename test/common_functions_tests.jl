import Test as TT

import CloudMicrophysics
import CloudMicrophysics.Common as CO

import Thermodynamics as TD
import CLIMAParameters as CP

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

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
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "H2SO4 solution saturated vapor pressure" begin

        T_warm = FT(225.0)
        T_cold = FT(200.0)
        T_too_warm = FT(240)
        T_too_cold = FT(180)
        x_sulph = FT(0.1)

        # If T out of range
        TT.@test_throws AssertionError("T < T_max") CO.H2SO4_soln_saturation_vapor_pressure(
            prs,
            x_sulph,
            T_too_warm,
        )
        TT.@test_throws AssertionError("T > T_min") CO.H2SO4_soln_saturation_vapor_pressure(
            prs,
            x_sulph,
            T_too_cold,
        )

        # p_sol should be higher at warmer temperatures
        TT.@test CO.H2SO4_soln_saturation_vapor_pressure(prs, x_sulph, T_warm) >
                 CO.H2SO4_soln_saturation_vapor_pressure(prs, x_sulph, T_cold)
    end
end

function test_a_w_xT(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "a_w_xT" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph_low = FT(0.06)
        x_sulph_high = FT(0.1)

        # a_w greater at warmer temperatures
        for x_sulph in [x_sulph_high, x_sulph_low]
            TT.@test CO.a_w_xT(prs, x_sulph, T_cold) <
                     CO.a_w_xT(prs, x_sulph, T_warm)
        end

        # a_w greater at lower sulphuric acid concentration
        for T in [T_warm, T_cold]
            TT.@test CO.a_w_xT(prs, x_sulph_high, T) <
                     CO.a_w_xT(prs, x_sulph_low, T)
        end
    end
end

function test_a_w_eT(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "a_w_eT" begin

        T_warm = FT(285)
        T_cold = FT(251)
        e_high = FT(1088)
        e_low = FT(544)

        # a_w greater at higher altitudes
        TT.@test CO.a_w_eT(prs, e_low, T_cold) > CO.a_w_eT(prs, e_high, T_warm)

        # a_w greater at greater partial pressures
        for T in [T_warm, T_cold]
            TT.@test CO.a_w_eT(prs, e_low, T) < CO.a_w_eT(prs, e_high, T)
        end
    end
end

function test_a_w_ice(FT)
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)

    TT.@testset "a_w_ice" begin

        T_warm = FT(240)
        T_cold = FT(230)

        # a_w greater at warmer temperatures
        TT.@test CO.a_w_ice(prs, T_cold) < CO.a_w_ice(prs, T_warm)

    end
end

println("Testing Float64")
test_H2SO4_soln_saturation_vapor_pressure(Float64)
test_a_w_xT(Float64)
test_a_w_eT(Float64)
test_a_w_ice(Float64)

println("Testing Float32")
test_H2SO4_soln_saturation_vapor_pressure(Float32)
test_a_w_xT(Float32)
test_a_w_eT(Float32)
test_a_w_ice(Float32)
