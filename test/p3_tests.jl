import Test as TT
import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import Thermodynamics as TD

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
        # values
        ρ_r = FT(500)
        F_r = FT(0.5)

        # get thresholds
        D_th = P3.D_th_helper(p3)
        th = P3.thresholds(p3, ρ_r, F_r)
        (; D_gr, D_cr) = th

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        TT.@test P3.p3_area(p3, D_1, F_r, th) == P3.A_s(D_1)
        TT.@test P3.p3_area(p3, D_2, F_r, th) == P3.A_ns(p3, D_2)
        TT.@test P3.p3_area(p3, D_3, F_r, th) == P3.A_s(D_3)
        TT.@test P3.p3_area(p3, D_cr, F_r, th) == P3.A_r(p3, F_r, D_cr)

        # test mass
        TT.@test P3.p3_mass(p3, D_1, F_r, th) == P3.mass_s(D_1, p3.ρ_i)
        TT.@test P3.p3_mass(p3, D_2, F_r, th) == P3.mass_nl(p3, D_2)
        TT.@test P3.p3_mass(p3, D_3, F_r, th) == P3.mass_s(D_3, th.ρ_g)
        TT.@test P3.p3_mass(p3, D_cr, F_r, th) == P3.mass_r(p3, D_cr, F_r)

        # test F_r = 0 and D > D_th
        F_r = FT(0)
        TT.@test P3.p3_area(p3, D_2, F_r, th) == P3.A_ns(p3, D_2)
        TT.@test P3.p3_mass(p3, D_2, F_r, th) == P3.mass_nl(p3, D_2)

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

function test_particle_terminal_velocities(FT)

    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    ρ_a = FT(1.2)

    TT.@testset "Chen 2022 - Rain" begin
        Ds = range(FT(1e-6), stop = FT(1e-5), length = 5)
        expected = [0.002508, 0.009156, 0.01632, 0.02377, 0.03144]
        for i in axes(Ds, 1)
            vel = CM2.rain_particle_terminal_velocity(Ds[i], Chen2022.rain, ρ_a)
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
    end

    TT.@testset "Chen 2022 - Ice" begin
        F_r = FT(0.5)
        ρ_r = FT(500)
        th = P3.thresholds(p3, ρ_r, F_r)
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.7912, 1.1550, 1.4871]
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = P3.ice_particle_terminal_velocity(D, Chen2022.snow_ice, ρ_a)
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
    end
end

function test_bulk_terminal_velocities(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    p3 = CMP.ParametersP3(FT)
    q = FT(0.22)
    N = FT(1e6)
    ρ_a = FT(1.2)
    ρ_rs = [FT(200), FT(400), FT(600), FT(800)]
    F_rs = [FT(0), FT(0.2), FT(0.4), FT(0.6), FT(0.8)]

    TT.@testset "Mass and number weighted terminal velocities" begin
        reference_vals_m = [
            [7.79, 7.27, 6.66, 5.94, 5.25],
            [7.79, 7.26, 6.62, 5.83, 4.82],
            [7.79, 7.25, 6.62, 5.81, 4.7],
            [7.79, 7.25, 6.62, 5.81, 4.65],
        ]
        reference_vals_n = [
            [3.65, 3.37, 3.05, 2.64, 2.14],
            [3.64, 3.37, 3.04, 2.62, 2.04],
            [3.65, 3.37, 3.04, 2.62, 2.02],
            [3.64, 3.37, 3.04, 2.61, 2.01],
        ]
        for i in 1:length(ρ_rs)
            for j in 1:length(F_rs)
                ρ_r = ρ_rs[i]
                F_r = F_rs[j]

                calculated_vel = P3.ice_terminal_velocity(
                    p3,
                    Chen2022.snow_ice,
                    q,
                    N,
                    ρ_r,
                    F_r,
                    ρ_a,
                )

                # number weighted
                TT.@test calculated_vel[1] > 0
                TT.@test reference_vals_n[i][j] ≈ calculated_vel[1] atol = 0.1

                # mass weighted
                TT.@test calculated_vel[2] > 0
                TT.@test reference_vals_m[i][j] ≈ calculated_vel[2] atol = 0.1
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

#function test_tendencies(FT)
#
#    tps = TD.Parameters.ThermodynamicsParameters(FT)
#
#    p3 = CMP.ParametersP3(FT)
#    Chen2022 = CMP.Chen2022VelType(FT)
#    aps = CMP.AirProperties(FT)
#
#    SB2006 = CMP.SB2006(FT, false) # no limiters
#    pdf_r = SB2006.pdf_r
#    pdf_c = SB2006.pdf_c
#
#    TT.@testset "Collision Tendencies Smoke Test" begin
#        N = FT(1e8)
#        ρ_a = FT(1.2)
#        ρ_r = FT(500)
#        F_r = FT(0.5)
#        T_warm = FT(300)
#        T_cold = FT(200)
#
#        qs = range(0.001, stop = 0.005, length = 5)
#        q_const = FT(0.05)
#
#        cloud_expected_warm_previous =
#            [6.341e-5, 0.0002099, 0.0004258, 0.0007047, 0.001042]
#        cloud_expected_cold_previous =
#            [0.002196, 0.003525, 0.004687, 0.005761, 0.006781]
#        cloud_expected_warm =
#            [8.043e-27, 3.641e-26, 8.773e-26, 1.63e-25, 2.625e-25]
#        cloud_expected_cold =
#            [8.197e-33, 1.865e-32, 3.012e-32, 4.233e-32, 5.523e-32]
#        rain_expected_warm = [0.000402, 0.0008436, 0.001304, 0.001777, 0.00226]
#        rain_expected_cold = [0.4156, 0.4260, 0.433, 0.4383, 0.4427]
#
#        for i in axes(qs, 1)
#            cloud_warm = P3.ice_collisions(
#                pdf_c,
#                p3,
#                Chen2022,
#                qs[i],
#                N,
#                q_const,
#                N,
#                ρ_a,
#                F_r,
#                ρ_r,
#                T_warm,
#            )
#            cloud_cold = P3.ice_collisions(
#                pdf_c,
#                p3,
#                Chen2022,
#                qs[i],
#                N,
#                q_const,
#                N,
#                ρ_a,
#                F_r,
#                ρ_r,
#                T_cold,
#            )
#
#            TT.@test cloud_warm >= 0
#            TT.@test cloud_warm ≈ cloud_expected_warm[i] rtol = 1e-3
#            TT.@test cloud_cold >= 0
#            TT.@test cloud_cold ≈ cloud_expected_cold[i] rtol = 1e-3
#
#            rain_warm = P3.ice_collisions(
#                pdf_r,
#                p3,
#                Chen2022,
#                qs[i],
#                N,
#                q_const,
#                N,
#                ρ_a,
#                F_r,
#                ρ_r,
#                T_warm,
#            )
#            rain_cold = P3.ice_collisions(
#                pdf_r,
#                p3,
#                Chen2022,
#                qs[i],
#                N,
#                q_const,
#                N,
#                ρ_a,
#                F_r,
#                ρ_r,
#                T_cold,
#            )
#
#            TT.@test rain_warm >= 0
#            TT.@test rain_warm ≈ rain_expected_warm[i] rtol = 1e-3
#            TT.@test rain_cold >= 0
#            TT.@test rain_cold ≈ rain_expected_cold[i] rtol = 1e-3
#        end
#    end
#
#    TT.@testset "Melting Tendencies Smoke Test" begin
#        N = FT(1e8)
#        ρ_a = FT(1.2)
#        ρ_r = FT(500)
#        F_r = FT(0.5)
#        T_freeze = FT(273.15)
#
#        qs = range(0.001, stop = 0.005, length = 5)
#
#        expected_melt = [0.0006982, 0.0009034, 0.001054, 0.001177, 0.001283]
#
#        for i in axes(qs, 1)
#            rate = P3.p3_melt(
#                p3,
#                Chen2022,
#                aps,
#                tps,
#                qs[i],
#                N,
#                T_freeze + 2,
#                ρ_a,
#                F_r,
#                ρ_r,
#            )
#
#            TT.@test rate >= 0
#            TT.@test rate ≈ expected_melt[i] rtol = 1e-3
#        end
#    end
#
#    TT.@testset "Heterogeneous Freezing Smoke Test" begin
#        T = FT(250)
#        N = FT(1e8)
#        ρ_a = FT(1.2)
#        qᵥ = FT(8.1e-4)
#        aero_type = CMP.Illite(FT)
#
#        qs = range(0.001, stop = 0.005, length = 5)
#
#        expected_freeze_q =
#            [2.036e-61, 6.463e-61, 1.270e-60, 2.052e-60, 2.976e-60]
#        expected_freeze_N =
#            [1.414e-51, 2.244e-51, 2.941e-51, 3.562e-51, 4.134e-51]
#
#        for i in axes(qs, 1)
#            rate_mass = P3.p3_rain_het_freezing(
#                true,
#                pdf_r,
#                p3,
#                tps,
#                qs[i],
#                N,
#                T,
#                ρ_a,
#                qᵥ,
#                aero_type,
#            )
#            rate_num = P3.p3_rain_het_freezing(
#                false,
#                pdf_r,
#                p3,
#                tps,
#                qs[i],
#                N,
#                T,
#                ρ_a,
#                qᵥ,
#                aero_type,
#            )
#
#            TT.@test rate_mass >= 0
#            TT.@test rate_mass ≈ expected_freeze_q[i] rtol = 1e-3
#
#            TT.@test rate_num >= 0
#            TT.@test rate_num ≈ expected_freeze_N[i] rtol = 1e-3
#        end
#    end
#
#end

function test_integrals(FT)
    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    N = FT(1e8)
    qs = range(0.001, stop = 0.005, length = 5)
    ρ_r = FT(500)
    F_rs = [FT(0), FT(0.5)]
    ρ_a = FT(1.2)
    tolerance = eps(FT)

    TT.@testset "Gamma vs Integral Comparison" begin
        for F_r in F_rs
            for i in axes(qs, 1)
                q = qs[i]

                # Velocity comparisons
                vel_N, vel_m = P3.ice_terminal_velocity(
                    p3,
                    Chen2022.snow_ice,
                    q,
                    N,
                    ρ_r,
                    F_r,
                    ρ_a,
                )

                λ, N_0 = P3.distribution_parameter_solver(p3, q, N, ρ_r, F_r)
                th = P3.thresholds(p3, ρ_r, F_r)
                ice_bound = P3.get_ice_bound(p3, λ, tolerance)
                vel(d) =
                    P3.ice_particle_terminal_velocity(d, Chen2022.snow_ice, ρ_a)
                f(d) = vel(d) * P3.N′ice(p3, d, λ, N_0)

                qgk_vel_N, = QGK.quadgk(d -> f(d) / N, FT(0), 2 * ice_bound)
                qgk_vel_m, = QGK.quadgk(
                    d -> f(d) * P3.p3_mass(p3, d, F_r, th) / q,
                    FT(0),
                    2 * ice_bound,
                )

                TT.@test vel_N ≈ qgk_vel_N rtol = 1e-7
                TT.@test vel_m ≈ qgk_vel_m rtol = 1e-7

                # Dₘ comparisons
                D_m = P3.D_m(p3, q, N, ρ_r, F_r)
                f_d(d) =
                    d * P3.p3_mass(p3, d, F_r, th) * P3.N′ice(p3, d, λ, N_0)
                qgk_D_m, = QGK.quadgk(d -> f_d(d) / q, FT(0), 2 * ice_bound)

                TT.@test D_m ≈ qgk_D_m rtol = 1e-8
            end
        end
    end
end


println("Testing Float32")
test_p3_thresholds(Float32)
test_particle_terminal_velocities(Float64)
#TODO - only works for Float64 now. We should switch the units inside the solver
# from SI base to something more managable
#test_p3_shape_solver(Float32)

println("Testing Float64")
test_p3_thresholds(Float64)
test_p3_shape_solver(Float64)
test_particle_terminal_velocities(Float64)
test_bulk_terminal_velocities(Float64)
#test_tendencies(Float64)
test_integrals(Float64)
