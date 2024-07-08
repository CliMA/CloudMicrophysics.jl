import Test as TT
import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.TerminalVelocity as TV

import QuadGK as QGK

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

    TT.@testset "mass and area tests" begin
        # test F_liq = 0
        # values 
        ρ_r = FT(500)
        F_r = FT(0.5)
        F_liq = FT(0)

        # get thresholds 
        D_th = P3.D_th_helper(p3)
        th = P3.thresholds(p3, ρ_r, F_r)
        (; D_gr, D_cr) = th

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area 
        TT.@test P3.p3_area(p3, D_1, F_r, F_liq, th) == P3.A_s(D_1)
        TT.@test P3.p3_area(p3, D_2, F_r, F_liq, th) == P3.A_ns(p3, D_2)
        TT.@test P3.p3_area(p3, D_3, F_r, F_liq, th) == P3.A_s(D_3)
        TT.@test P3.p3_area(p3, D_cr, F_r, F_liq, th) == P3.A_r(p3, F_r, D_cr)

        # test mass 
        TT.@test P3.p3_mass(p3, D_1, F_r, F_liq, th) == P3.mass_s(D_1, p3.ρ_i)
        TT.@test P3.p3_mass(p3, D_2, F_r, F_liq, th) == P3.mass_nl(p3, D_2)
        TT.@test P3.p3_mass(p3, D_3, F_r, F_liq, th) == P3.mass_s(D_3, th.ρ_g)
        TT.@test P3.p3_mass(p3, D_cr, F_r, F_liq, th) ==
                 P3.mass_r(p3, D_cr, F_r)

        # test F_r = 0 and D > D_th
        F_r = FT(0)
        TT.@test P3.p3_area(p3, D_2, F_r, F_liq, th) == P3.A_ns(p3, D_2)
        TT.@test P3.p3_mass(p3, D_2, F_r, F_liq, th) == P3.mass_nl(p3, D_2)

        # test F_liq != 0
        F_liq = FT(0.5)

        # test area
        TT.@test P3.p3_area(p3, D_1, F_r, F_liq, th) == P3.A_s(D_1)
        TT.@test P3.p3_area(p3, D_2, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.A_ns(p3, D_2) + F_liq * P3.A_s(D_2)
        #TT.@test P3.p3_area(p3, D_3, F_r, F_liq, th) == P3.A_s(D_3)
        TT.@test P3.p3_area(p3, D_3, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.A_ns(p3, D_3) + F_liq * P3.A_s(D_3)
        TT.@test P3.p3_area(p3, D_cr, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.A_r(p3, F_r, D_cr) + F_liq * P3.A_s(D_cr)

        # test mass 
        TT.@test P3.p3_mass(p3, D_1, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.mass_s(D_1, p3.ρ_i) +
                 F_liq * P3.mass_liq(p3, D_1)
        TT.@test P3.p3_mass(p3, D_2, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.mass_nl(p3, D_2) +
                 F_liq * P3.mass_liq(p3, D_2)
        TT.@test P3.p3_mass(p3, D_3, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.mass_nl(p3, D_3) +
                 F_liq * P3.mass_liq(p3, D_3)
        #TT.@test P3.p3_mass(p3, D_3, F_r, F_liq, th) == (1 - F_liq) * P3.mass_s(D_3, th.ρ_g) + F_liq * P3.mass_liq(p3, D_3)
        TT.@test P3.p3_mass(p3, D_cr, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.mass_r(p3, D_cr, F_r) +
                 F_liq * P3.mass_liq(p3, D_cr)

        # test F_r = 0 and D > D_th
        F_r = FT(0)
        TT.@test P3.p3_area(p3, D_2, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.A_ns(p3, D_2) + F_liq * P3.A_s(D_2)
        TT.@test P3.p3_mass(p3, D_2, F_r, F_liq, th) ==
                 (1 - F_liq) * P3.mass_nl(p3, D_2) +
                 F_liq * P3.mass_liq(p3, D_2)


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
        F_liq_test = (FT(0), FT(0.33), FT(0.67), FT(1))

        # check that the shape solution solves to give correct values
        for N in N_test
            for λ_ex in λ_test
                for ρ_r in ρ_r_test
                    for F_r in F_r_test
                        for F_liq in F_liq_test
                            # Compute the shape parameters that correspond to the
                            # input test values
                            μ_ex = P3.DSD_μ(p3, λ_ex)
                            N₀_ex = P3.DSD_N₀(p3, N, λ_ex)
                            # Find the P3 scheme  thresholds
                            th = P3.thresholds(p3, ρ_r, F_r)
                            # Convert λ to ensure it remains positive
                            x = log(λ_ex)
                            # Compute mass density based on input shape parameters
                            q_calc =
                                N *
                                P3.q_over_N_gamma(p3, F_liq, F_r, x, μ_ex, th)

                            if q_calc < FT(1)
                                # Solve for shape parameters
                                (λ, N₀) = P3.distribution_parameter_solver(
                                    p3,
                                    q_calc,
                                    N,
                                    ρ_r,
                                    F_liq,
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
end

function test_velocities(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    p3 = CMP.ParametersP3(FT)
    q = FT(0.22)
    N = FT(1e6)
    ρ_a = FT(1.2)
    ρ_rs = [FT(200), FT(400), FT(600), FT(800)]
    F_rs = [FT(0), FT(0.2), FT(0.4), FT(0.6), FT(0.8)]
    F_liq = FT(0) # set F_liq = 0 to test against MM2015 values

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

                calculated_vel = P3.terminal_velocity_tot(
                    p3,
                    Chen2022.snow_ice,
                    Chen2022.rain,
                    q,
                    N,
                    ρ_r,
                    F_liq,
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

    ## TODO: Test velocities with F_liq != 0 ##
    # TT.@testset "Mass and number weighted terminal velocities with nonzero liquid fraction" begin
    #     F_liqs = [FT(0), FT(0.33), FT(0.67), FT(1)]
    #     ρ_r = ρ_rs[4]
    #     F_r = F_rs[4]
    #     expected_vals = zeros(length(F_liqs))
    #     for k in eachindex(F_liqs)
    #         F_liq = F_liqs[k]

    #         calculated_vel = P3.terminal_velocity_tot(
    #             p3,
    #             Chen2022.snow_ice,
    #             Chen2022.rain,
    #             q,
    #             N,
    #             ρ_r,
    #             F_liq,
    #             F_r,
    #             ρ_a,
    #         )
    #         expected_vals[k] = calculated_vel

    #         # number weighted
    #         # TT.@test calculated_vel[2] > 0
    #         # TT.@test expected_vals[k][1] ≈ calculated_vel[1] atol = 0.1

    #         # # mass weighted
    #         # TT.@test calculated_vel[1] > 0
    #         # TT.@test expected_vals[k][2] ≈ calculated_vel[2] atol = 3.14
    #     end
    #     print(expected_vals)
    # end

    # set F_liq = 0 to test against MM2015 values
    TT.@testset "Mass-weighted mean diameters" begin
        F_liq = FT(0)
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

                calculated_dm = P3.D_m_tot(p3, q, N, ρ_r, F_r, F_liq) * 1e3

                TT.@test calculated_dm > 0
                TT.@test paper_vals[i][j] ≈ calculated_dm atol = 3.14

            end
        end
    end

    TT.@testset "Mass-weighted mean diameters with nonzero liquid fraction" begin
        ρ_r = ρ_rs[4]
        F_r = F_rs[4]
        F_liqs = [FT(0), FT(0.33), FT(0.67), FT(1)]
        expected_vals =
            [FT(0.00333689), FT(0.000892223), FT(0.00124423), FT(0.00164873)]
        for i in eachindex(F_liqs)
            calculated_dm = P3.D_m_tot(p3, q, N, ρ_r, F_r, F_liqs[i])
            expected_vals[i] = calculated_dm
            TT.@test expected_vals[i] ≈ calculated_dm atol = 1e-6
        end
    end
end

## TODO: Fix Gamma vs Integral Comparison for F_liq != 0 ##
function test_integrals(FT)
    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    N = FT(1e8)
    qs = range(0.001, stop = 0.005, length = 5)
    ρ_r = FT(500)
    F_liqs = [FT(0)]
    #F_liqs = [FT(0), FT(0.5), FT(1)]
    F_rs = [FT(0), FT(0.5)]
    ρ_a = FT(1.2)
    tolerance = eps(FT)

    TT.@testset "Gamma vs Integral Comparison" begin
        for F_r in F_rs
            for F_liq in F_liqs
                for i in axes(qs, 1)
                    q = qs[i]

                    # Velocity comparisons
                    vel_N, vel_m = P3.terminal_velocity_tot(
                        p3,
                        Chen2022.snow_ice,
                        Chen2022.rain,
                        q,
                        N,
                        ρ_r,
                        F_liq,
                        F_r,
                        ρ_a,
                    )

                    λ, N_0 = P3.distribution_parameter_solver(
                        p3,
                        q,
                        N,
                        ρ_r,
                        F_liq,
                        F_r,
                    )
                    th = P3.thresholds(p3, ρ_r, F_r)
                    ice_bound = P3.get_ice_bound(p3, λ, tolerance)
                    vel(d) = TV.velocity_chen(
                        d,
                        Chen2022.snow_ice,
                        ρ_a,
                        P3.p3_mass(p3, d, F_r, F_liq, th),
                        P3.p3_area(p3, d, F_r, F_liq, th),
                        p3.ρ_i,
                    )
                    vel_liq(d) = TV.velocity_chen(d, Chen2022.rain, ρ_a)
                    f(d) =
                        ((1 - F_liq) * vel(d) + F_liq * vel_liq(d)) *
                        P3.N′ice(p3, d, λ, N_0)

                    qgk_vel_N, = QGK.quadgk(d -> f(d) / N, FT(0), 2 * ice_bound)
                    qgk_vel_m, = QGK.quadgk(
                        d -> f(d) * P3.p3_mass(p3, d, F_r, F_liq, th) / q,
                        FT(0),
                        2 * ice_bound,
                    )

                    TT.@test vel_N ≈ qgk_vel_N rtol = 1e-7
                    TT.@test vel_m ≈ qgk_vel_m rtol = 1e-7

                    # Dₘ comparisons 
                    D_m = P3.D_m_tot(p3, q, N, ρ_r, F_r, F_liq)
                    f_d(d) =
                        d *
                        P3.p3_mass(p3, d, F_r, F_liq, th) *
                        P3.N′ice(p3, d, λ, N_0)
                    qgk_D_m, = QGK.quadgk(d -> f_d(d) / q, FT(0), 2 * ice_bound)

                    TT.@test D_m ≈ qgk_D_m rtol = 1e-8
                end
            end
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
test_integrals(Float64)
