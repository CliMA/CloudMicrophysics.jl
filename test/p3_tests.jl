import Test as TT
import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP

@info "P3 Scheme Tests"

function test_p3_thresholds(FT)

    p3 = CMP.ParametersP3(FT)

    TT.@testset "thresholds (nonlinear solver function)" begin

        # initialize test values:
        ρ_r = FT(400)
        F_r = FT(0.8)
        ρ_r_good = (FT(200), FT(400), FT(800)) # representative ρ_r values
        F_r_good = (FT(0.5), FT(0.8), FT(0.95)) # representative F_r values

        # test asserts
        for _ρ_r in (FT(0), FT(-1))
            TT.@test_throws AssertionError("ρ_r > FT(0)") P3.thresholds(
                p3,
                _ρ_r,
                F_r,
            )
        end
        for _ρ_r in (FT(0), FT(-1))
            TT.@test_throws AssertionError("ρ_r <= p3.ρ_l") P3.thresholds(
                p3,
                FT(1200),
                F_r,
            )
        end

        for _F_r in (FT(-1 * eps(FT)), FT(-1))
            TT.@test_throws AssertionError("F_r >= FT(0)") P3.thresholds(
                p3,
                ρ_r,
                _F_r,
            )
        end
        for _F_r in (FT(1), FT(1.5))
            TT.@test_throws AssertionError("F_r < FT(1)") P3.thresholds(
                p3,
                ρ_r,
                _F_r,
            )
        end

        # Test if the P3 scheme solution satisifies the conditions
        # from eqs. 14-17 in Morrison and Milbrandt 2015
        for F_r in F_r_good
            for ρ_r in ρ_r_good
                sol = P3.thresholds(p3, ρ_r, F_r)
                atol = 5e3 * eps(FT)
                TT.@test sol.D_cr ≈ P3.D_cr_helper(p3, F_r, sol[3]) atol = atol
                TT.@test sol.D_gr ≈ P3.D_gr_helper(p3, sol[3]) atol = atol
                TT.@test sol.ρ_g ≈ P3.ρ_g_helper(ρ_r, F_r, sol[4]) atol = atol
                TT.@test sol.ρ_d ≈ P3.ρ_d_helper(p3, sol[1], sol[2]) atol = atol
            end
        end

        # Check that the P3 scheme solution matches the published values
        function diff(ρ_r, F_r, el, gold, rtol = 2e-2)
            TT.@test P3.thresholds(p3, ρ_r, F_r)[el] * 1e3 ≈ gold rtol = rtol
        end
        # D_cr and D_gr vs Fig. 1a Morrison and Milbrandt 2015
        D_cr_fig_1a_ref = [FT(0.4946323381999426), FT(1.0170979628696817)]
        D_gr_fig_1a_ref = [FT(0.26151186272014415), FT(0.23392868352755775)]
        for val in [1, 2]
            diff(ρ_r_good[2], F_r_good[val], 1, D_cr_fig_1a_ref[val])
            diff(ρ_r_good[2], F_r_good[val], 2, D_gr_fig_1a_ref[val])
        end
        # D_cr and D_gr vs Fig. 1b Morrison and Milbrandt 2015
        #! format: off
        D_cr_fig_1b_ref = [FT(6.152144691917768), FT(3.2718818175768405), FT(1.7400778369620664)]
        D_gr_fig_1b_ref = [FT(0.39875043123651077), FT(0.2147085163169669), FT(0.11516682512848)]
        #! format: on
        for val in [1, 2, 3]
            diff(ρ_r_good[val], F_r_good[3], 1, D_cr_fig_1b_ref[val])
            diff(ρ_r_good[val], F_r_good[3], 2, D_gr_fig_1b_ref[val])
        end
    end
end

function test_p3_shape_solver(FT)

    p3 = CMP.ParametersP3(FT)

    TT.@testset "shape parameters (nonlinear solver function)" begin

        # initialize test values:
        ep = 1 #1e4 * eps(FT)
        N_test = (FT(1e7), FT(1e8), FT(1e9), FT(1e10))                             # N values
        λ_test = (FT(1e1), FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))                # test λ values in range also do 15000, 20000
        ρ_r_test = (FT(200), FT(400), FT(600), FT(800))    # representative ρ_r values
        F_r_test = (FT(0), FT(0.5), FT(0.8), FT(0.95))        # representative F_r values

        # check that the shape solution solves to give correct values
        for N in N_test
            for λ_ex in λ_test
                for ρ_r in ρ_r_test
                    for F_r in F_r_test
                        # Compute the shape parameters that correspond to the
                        # input test values
                        μ_ex = P3.DSD_μ(p3, λ_ex)
                        N₀_ex = P3.DSD_N₀(p3, N, λ_ex)
                        # Find the P3 scheme  thresholds
                        th = P3.thresholds(p3, ρ_r, F_r)
                        # Convert λ to ensure it remains positive
                        x = log(λ_ex)
                        # Compute mass density based on input shape parameters
                        q_calc = N * P3.q_over_N_gamma(p3, F_r, x, μ_ex, th)

                        if q_calc < FT(1)
                            # Solve for shape parameters
                            (λ, N₀) = P3.distribution_parameter_solver(
                                p3,
                                q_calc,
                                N,
                                ρ_r,
                                F_r,
                            )

                            # Compare solved values with the input expected values
                            TT.@test λ ≈ λ_ex rtol = ep
                            TT.@test N₀ ≈ N₀_ex rtol = ep
                        end
                    end
                end
            end
        end
    end
end

function test_velocities(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    p3 = CMP.ParametersP3(FT)
    q = FT(0.22)
    N = FT(1e6)
    ρ_a = FT(1.2)
    ρ_rs = [FT(200), FT(400), FT(600), FT(800)]
    F_rs = [FT(0), FT(0.2), FT(0.4), FT(0.6), FT(0.8)]

    TT.@testset "Mass and number weighted terminal velocities" begin
        paper_vals = [
            [1.5, 1.5, 1.5, 1.5, 1.5],
            [1.5, 1.5, 2.5, 2.5, 2.5],
            [1.5, 2.5, 2.5, 2.5, 2.5],
            [1.5, 2.5, 3.5, 3.5, 3.5],
        ]
        expected_vals = [
            [1.52, 1.46, 1.41, 1.36, 1.24],
            [1.52, 1.47, 1.44, 1.42, 1.35],
            [1.52, 1.47, 1.45, 1.44, 1.42],
            [1.52, 1.47, 1.45, 1.45, 1.45],
        ]
        for i in 1:length(ρ_rs)
            for j in 1:length(F_rs)
                ρ_r = ρ_rs[i]
                F_r = F_rs[j]

                calculated_vel = P3.terminal_velocity(
                    p3,
                    Chen2022.snow_ice,
                    q,
                    N,
                    ρ_r,
                    F_r,
                    ρ_a,
                )

                # number weighted
                TT.@test calculated_vel[2] > 0
                TT.@test expected_vals[i][j] ≈ calculated_vel[1] atol = 0.1

                # mass weighted
                TT.@test calculated_vel[1] > 0
                TT.@test paper_vals[i][j] ≈ calculated_vel[2] atol = 3.14
            end
        end
    end

    TT.@testset "Mass-weighted mean diameters" begin
        paper_vals = [
            [5, 5, 5, 5, 5],
            [4.5, 4.5, 4.5, 4.5, 4.5],
            [3.5, 3.5, 3.5, 3.5, 3.5],
            [3.5, 3.5, 2.5, 2.5, 2.5],
        ]
        for i in 1:length(ρ_rs)
            for j in 1:length(F_rs)
                ρ_r = ρ_rs[i]
                F_r = F_rs[j]

                calculated_dm = P3.D_m(p3, q, N, ρ_r, F_r) * 1e3

                TT.@test calculated_dm > 0
                TT.@test paper_vals[i][j] ≈ calculated_dm atol = 3.14

            end
        end
    end
end

function test_tendencies(FT)

    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    TT.@testset "Collision Tendencies Smoke Test" begin
        N = FT(1e8)
        ρ_a = FT(1.2)
        ρ_r = FT(500)
        F_r = FT(0.5)
        T_warm = FT(300)
        T_cold = FT(200)

        qs = range(0.001, stop = 0.005, length = 5)
        q_const = FT(0.05)

        cloud_expected_warm =
            [5.78e-5, 0.00019256, 0.00039239, 0.00065184, 0.00096698]
        cloud_expected_cold =
            [0.0016687, 0.0026921, 0.0035912, 0.0044255, 0.00522]
        rain_expected_warm = [0.0003392, 0.000713, 0.001103, 0.001506, 0.00192]
        rain_expected_cold = [0.2905, 0.2982, 0.3033, 0.3072, 0.3104]

        for i in axes(qs, 1)
            cloud_warm = P3.ice_collisions(
                "cloud",
                p3,
                Chen2022,
                qs[i],
                N,
                q_const,
                N,
                ρ_a,
                F_r,
                ρ_r,
                T_warm,
            )
            cloud_cold = P3.ice_collisions(
                "cloud",
                p3,
                Chen2022,
                qs[i],
                N,
                q_const,
                N,
                ρ_a,
                F_r,
                ρ_r,
                T_cold,
            )

            TT.@test cloud_warm >= 0
            TT.@test cloud_warm ≈ cloud_expected_warm[i] rtol = 1e-3
            TT.@test cloud_cold >= 0
            TT.@test cloud_cold ≈ cloud_expected_cold[i] rtol = 1e-3

            rain_warm = P3.ice_collisions(
                "rain",
                p3,
                Chen2022,
                qs[i],
                N,
                q_const,
                N,
                ρ_a,
                F_r,
                ρ_r,
                T_warm,
            )
            rain_cold = P3.ice_collisions(
                "rain",
                p3,
                Chen2022,
                qs[i],
                N,
                q_const,
                N,
                ρ_a,
                F_r,
                ρ_r,
                T_cold,
            )

            TT.@test rain_warm >= 0
            TT.@test rain_warm ≈ rain_expected_warm[i] rtol = 1e-3
            TT.@test rain_cold >= 0
            TT.@test rain_cold ≈ rain_expected_cold[i] rtol = 1e-3
        end
    end

end

println("Testing Float32")
test_p3_thresholds(Float32)
#TODO - only works for Float64 now. We should switch the units inside the solver
# from SI base to something more managable
#test_p3_shape_solver(Float32)

println("Testing Float64")
test_p3_thresholds(Float64)
test_p3_shape_solver(Float64)
test_velocities(Float64)
test_tendencies(Float64)
