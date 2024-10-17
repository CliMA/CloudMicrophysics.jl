import Test as TT
import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import Thermodynamics as TD

import QuadGK as QGK

@info("P3 Scheme Tests")

function test_incomplete_gamma_function_approximation(FT)
    # A smoke tests to see if the reference values are not changed
    TT.@testset "Incomplete Γ function approximation" begin

        Γ_lower_reference = [
            [0.6171173028812539, 0.9999647375921944, 0.9999999999999561],
            [66.43238525128415, 192195.93434392574, 362553.66f0],
            [9.801186821098075e11, 8.309791f14, 1.1845235f17],
        ]

        a = [FT(1), FT(10), FT(20)]
        z = [FT(1), FT(10), FT(30)]

        for it in [1, 2, 3]
            for jt in [1, 2, 3]
                TT.@test P3.Γ_lower(a[it], z[jt]) ≈ Γ_lower_reference[it][jt] rtol =
                    1e-3
            end
        end
    end
end

function test_thresholds_solver(FT)

    p3 = CMP.ParametersP3(FT)

    TT.@testset "Thresholds - nonlinear solver" begin

        # initialize test values:
        ρ_r = FT(400)
        F_rim = FT(0.8)
        ρ_r_good = (FT(200), FT(400), FT(800)) # representative ρ_r values
        F_rim_good = (FT(0.5), FT(0.8), FT(0.95)) # representative F_rim values

        # test asserts
        for _ρ_r in (FT(0), FT(-1))
            TT.@test_throws AssertionError("ρ_r > FT(0)") P3.thresholds(
                p3,
                _ρ_r,
                F_rim,
            )
        end
        for _ρ_r in (FT(0), FT(-1))
            TT.@test_throws AssertionError("ρ_r <= p3.ρ_l") P3.thresholds(
                p3,
                FT(1200),
                F_rim,
            )
        end

        for _F_rim in (FT(-1 * eps(FT)), FT(-1))
            TT.@test_throws AssertionError("F_rim >= FT(0)") P3.thresholds(
                p3,
                ρ_r,
                _F_rim,
            )
        end
        for _F_rim in (FT(1), FT(1.5))
            TT.@test_throws AssertionError("F_rim < FT(1)") P3.thresholds(
                p3,
                ρ_r,
                _F_rim,
            )
        end

        # Test if the P3 scheme solution satisifies the conditions
        # from eqs. 14-17 in Morrison and Milbrandt 2015
        for F_rim in F_rim_good
            for ρ_r in ρ_r_good
                sol = P3.thresholds(p3, ρ_r, F_rim)
                atol = 5e3 * eps(FT)
                TT.@test sol.D_cr ≈ P3.D_cr_helper(p3, F_rim, sol[3]) atol =
                    atol
                TT.@test sol.D_gr ≈ P3.D_gr_helper(p3, sol[3]) atol = atol
                TT.@test sol.ρ_g ≈ P3.ρ_g_helper(ρ_r, F_rim, sol[4]) atol = atol
                TT.@test sol.ρ_d ≈ P3.ρ_d_helper(p3, sol[1], sol[2]) atol = atol
            end
        end

        # Check that the P3 scheme solution matches the published values
        function diff(ρ_r, F_rim, el, gold, rtol = 2e-2)
            TT.@test P3.thresholds(p3, ρ_r, F_rim)[el] * 1e3 ≈ gold rtol = rtol
        end
        # D_cr and D_gr vs Fig. 1a Morrison and Milbrandt 2015
        D_cr_fig_1a_ref = [FT(0.4946323381999426), FT(1.0170979628696817)]
        D_gr_fig_1a_ref = [FT(0.26151186272014415), FT(0.23392868352755775)]
        for val in [1, 2]
            diff(ρ_r_good[2], F_rim_good[val], 1, D_cr_fig_1a_ref[val])
            diff(ρ_r_good[2], F_rim_good[val], 2, D_gr_fig_1a_ref[val])
        end
        # D_cr and D_gr vs Fig. 1b Morrison and Milbrandt 2015
        #! format: off
        D_cr_fig_1b_ref = [FT(6.152144691917768), FT(3.2718818175768405), FT(1.7400778369620664)]
        D_gr_fig_1b_ref = [FT(0.39875043123651077), FT(0.2147085163169669), FT(0.11516682512848)]
        #! format: on
        for val in [1, 2, 3]
            diff(ρ_r_good[val], F_rim_good[3], 1, D_cr_fig_1b_ref[val])
            diff(ρ_r_good[val], F_rim_good[3], 2, D_gr_fig_1b_ref[val])
        end
    end

    TT.@testset "Thresholds - mass, area, density, aspect ratio" begin
        # values
        ρ_r = FT(500)
        F_rim = FT(0.5)
        F_liq = FT(0) # testing F_liq = 0

        # get thresholds
        D_th = P3.D_th_helper(p3)
        th = P3.thresholds(p3, ρ_r, F_rim)
        (; D_gr, D_cr) = th

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        TT.@test P3.p3_area(p3, D_1, F_rim, F_liq, th) == P3.A_s(D_1)
        TT.@test P3.p3_area(p3, D_2, F_rim, F_liq, th) == P3.A_ns(p3, D_2)
        TT.@test P3.p3_area(p3, D_3, F_rim, F_liq, th) == P3.A_s(D_3)
        TT.@test P3.p3_area(p3, D_cr, F_rim, F_liq, th) ==
                 P3.A_r(p3, F_rim, D_cr)

        # test mass
        TT.@test P3.p3_mass(p3, D_1, F_rim, F_liq, th) == P3.mass_s(D_1, p3.ρ_i)
        TT.@test P3.p3_mass(p3, D_2, F_rim, F_liq, th) == P3.mass_nl(p3, D_2)
        TT.@test P3.p3_mass(p3, D_3, F_rim, F_liq, th) == P3.mass_s(D_3, th.ρ_g)
        TT.@test P3.p3_mass(p3, D_cr, F_rim, F_liq, th) ==
                 P3.mass_r(p3, D_cr, F_rim)

        # test density
        TT.@test P3.p3_density(p3, D_1, F_rim, th) ≈ p3.ρ_i
        TT.@test P3.p3_density(p3, D_2, F_rim, th) ≈ 544.916989830
        TT.@test P3.p3_density(p3, D_3, F_rim, th) ≈ th.ρ_g
        TT.@test P3.p3_density(p3, D_cr, F_rim, th) ≈ 383.33480937

        # test aspect ratio
        TT.@test P3.ϕᵢ(p3, D_1, F_rim, th) ≈ 1
        TT.@test P3.ϕᵢ(p3, D_2, F_rim, th) ≈ 1
        TT.@test P3.ϕᵢ(p3, D_3, F_rim, th) ≈ 1
        TT.@test P3.ϕᵢ(p3, D_cr, F_rim, th) ≈ 1

        # test F_rim = 0 and D > D_th
        F_rim = FT(0)
        TT.@test P3.p3_area(p3, D_2, F_rim, F_liq, th) == P3.A_ns(p3, D_2)
        TT.@test P3.p3_mass(p3, D_2, F_rim, F_liq, th) == P3.mass_nl(p3, D_2)

        # test F_liq != 0
        F_liq = FT(0.5)

        # get D_i values again
        # values
        ρ_r = FT(500)
        F_rim = FT(0.5)

        # get thresholds
        D_th = P3.D_th_helper(p3)
        th = P3.thresholds(p3, ρ_r, F_rim)
        (; D_gr, D_cr) = th

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        TT.@test P3.p3_area(p3, D_1, F_rim, F_liq, th) == P3.A_s(D_1)
        TT.@test P3.p3_area(p3, D_2, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.A_ns(p3, D_2) + F_liq * P3.A_s(D_2)
        TT.@test P3.p3_area(p3, D_3, F_rim, F_liq, th) == P3.A_s(D_3)
        TT.@test P3.p3_area(p3, D_cr, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.A_r(p3, F_rim, D_cr) + F_liq * P3.A_s(D_cr)

        # test mass
        TT.@test P3.p3_mass(p3, D_1, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.mass_s(D_1, p3.ρ_i) +
                 F_liq * P3.mass_s(D_1, p3.ρ_l)
        TT.@test P3.p3_mass(p3, D_2, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.mass_nl(p3, D_2) +
                 F_liq * P3.mass_s(D_2, p3.ρ_l)
        TT.@test P3.p3_mass(p3, D_3, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.mass_s(D_3, th.ρ_g) +
                 F_liq * P3.mass_s(D_3, p3.ρ_l)
        TT.@test P3.p3_mass(p3, D_cr, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.mass_r(p3, D_cr, F_rim) +
                 F_liq * P3.mass_s(D_cr, p3.ρ_l)

        # test F_rim = 0 and D > D_th
        F_rim = FT(0)
        TT.@test P3.p3_area(p3, D_2, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.A_ns(p3, D_2) + F_liq * P3.A_s(D_2)
        TT.@test P3.p3_mass(p3, D_2, F_rim, F_liq, th) ==
                 (1 - F_liq) * P3.mass_nl(p3, D_2) +
                 F_liq * P3.mass_s(D_2, p3.ρ_l)
    end
end

function test_shape_solver(FT)

    p3 = CMP.ParametersP3(FT)

    TT.@testset "Shape parameters - nonlinear solver" begin

        # initialize test values:
        ep = 1 #1e4 * eps(FT)
        N_test = (FT(1e7), FT(1e8), FT(1e9), FT(1e10))                         # N values
        λ_test = (FT(1e1), FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))        # test λ values in range also do 15000, 20000
        ρ_r_test = (FT(200), FT(400), FT(600), FT(800))                        # representative ρ_r values
        F_rim_test = (FT(0), FT(0.5), FT(0.8), FT(0.95))                       # representative F_rim values
        F_liq_test = (FT(0), FT(0.33), FT(0.67), FT(1))                        # representative F_rim values

        # check that the shape solution solves to give correct values
        for N in N_test
            for λ_ex in λ_test
                for ρ_r in ρ_r_test
                    for F_rim in F_rim_test
                        for F_liq in F_liq_test
                            # Compute the shape parameters that correspond to the
                            # input test values
                            μ_ex = P3.DSD_μ(p3, λ_ex)
                            N₀_ex = P3.DSD_N₀(p3, N, λ_ex)
                            # Find the P3 scheme  thresholds
                            th = P3.thresholds(p3, ρ_r, F_rim)
                            # Convert λ to ensure it remains positive
                            x = log(λ_ex)
                            # Compute mass density based on input shape parameters
                            L_calc =
                                N *
                                P3.L_over_N_gamma(p3, F_liq, F_rim, x, μ_ex, th)

                            if L_calc < FT(1)
                                # Solve for shape parameters
                                (λ, N₀) = P3.distribution_parameter_solver(
                                    p3,
                                    L_calc,
                                    N,
                                    ρ_r,
                                    F_liq,
                                    F_rim,
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

function test_particle_terminal_velocities(FT)

    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    ρ_a = FT(1.2)

    TT.@testset "Smoke tests for rain particle terminal vel from Chen 2022" begin
        Ds = range(FT(1e-6), stop = FT(1e-5), length = 5)
        expected = [0.002508, 0.009156, 0.01632, 0.02377, 0.03144]
        for i in axes(Ds, 1)
            vel = CM2.rain_particle_terminal_velocity(Ds[i], Chen2022.rain, ρ_a)
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
    end

    TT.@testset "Smoke tests for ice particle terminal vel from Chen 2022" begin
        F_rim = FT(0.5)
        ρ_r = FT(500)
        F_liq = FT(0)
        th = P3.thresholds(p3, ρ_r, F_rim)
        use_aspect_ratio = false
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.7912, 1.1550, 1.4871]
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = P3.ice_particle_terminal_velocity(
                p3,
                D,
                Chen2022.snow_ice,
                ρ_a,
                F_rim,
                th,
                use_aspect_ratio,
            )
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end

        use_aspect_ratio = true
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.79121, 1.155, 1.487]
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = P3.ice_particle_terminal_velocity(
                p3,
                D,
                Chen2022.snow_ice,
                ρ_a,
                F_rim,
                th,
                use_aspect_ratio,
            )
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
    end

    TT.@testset "Smoke tests for mixed phase particle terminal velocity" begin
        F_rim = FT(0.5)
        F_liq = FT(0.5)
        ρ_r = FT(500)
        th = P3.thresholds(p3, ρ_r, F_rim)
        use_aspect_ratio = true
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13192, 0.50457, 0.90753, 1.3015, 1.6757]
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = P3.p3_particle_terminal_velocity(
                p3,
                D,
                Chen2022,
                ρ_a,
                F_rim,
                F_liq,
                th,
                use_aspect_ratio,
            )
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
        use_aspect_ratio = false
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13191, 0.50457, 0.90753, 1.301499, 1.67569]
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = P3.p3_particle_terminal_velocity(
                p3,
                D,
                Chen2022,
                ρ_a,
                F_rim,
                F_liq,
                th,
                use_aspect_ratio,
            )
            TT.@test vel >= 0
            TT.@test vel ≈ expected[i] rtol = 1e-3
        end
    end
end

function test_bulk_terminal_velocities(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    p3 = CMP.ParametersP3(FT)
    L = FT(0.22)
    N = FT(1e6)
    ρ_a = FT(1.2)
    ρ_r = FT(800)
    F_rims = [FT(0), FT(0.6)]
    F_liqs = [FT(0.5), FT(1)]

    TT.@testset "Mass and number weighted terminal velocities" begin

        # Liquid fraction = 0

        ref_v_n = [3.637320375996216, 2.6130243256593313]
        ref_v_n_ϕ = [3.637320375996216, 2.6130243256593313]
        #ref_v_n_ϕ = [2.2058692424028923, 2.24224820764301]
        ref_v_m = [7.7208191246659705, 5.751352472270078]
        ref_v_m_ϕ = [7.7208191246659705, 5.751352472270078]
        #ref_v_m_ϕ = [4.193799087495265, 4.773592179767702]

        for k in 1:length(F_rims)
            F_rim = F_rims[k]
            F_liq = FT(0)
            vel = P3.ice_terminal_velocity(
                p3,
                Chen2022,
                L,
                N,
                ρ_r,
                F_rim,
                F_liq,
                ρ_a,
                false,
            )
            vel_ϕ = P3.ice_terminal_velocity(
                p3,
                Chen2022,
                L,
                N,
                ρ_r,
                F_rim,
                F_liq,
                ρ_a,
                true,
            )

            # number weighted
            TT.@test vel[1] > 0
            TT.@test vel_ϕ[1] > 0
            TT.@test vel[1] ≈ ref_v_n[k] rtol = 1e-6
            TT.@test vel_ϕ[1] ≈ ref_v_n_ϕ[k] rtol = 1e-6

            # mass weighted
            TT.@test vel[2] > 0
            TT.@test vel_ϕ[2] > 0
            TT.@test vel[2] ≈ ref_v_m[k] rtol = 1e-6
            TT.@test vel_ϕ[2] ≈ ref_v_m_ϕ[k] rtol = 1e-6

            # slower with aspect ratio
            TT.@test vel_ϕ[1] <= vel[1]
            TT.@test vel_ϕ[2] <= vel[2]
        end

        # Liquid fraction != 0
        ref_v_n = [1.674591925057434, 1.6180970319460353]
        ref_v_n_ϕ = [1.674591925057434, 1.6180970319460353]
        #ref_v_n_ϕ = [1.549777478756061, 1.6180970319460353]
        ref_v_m = [5.126648558302173, 5.416679316254198]
        ref_v_m_ϕ = [5.126648558302173, 5.416679316254198]
        #ref_v_m_ϕ = [4.6358422594886495, 5.416679316254198]

        for k in 1:length(F_liqs)
            F_liq = F_liqs[k]
            F_rim = FT(0.4)
            vel = P3.ice_terminal_velocity(
                p3,
                Chen2022,
                L,
                N,
                ρ_r,
                F_rim,
                F_liq,
                ρ_a,
                false,
            )
            vel_ϕ = P3.ice_terminal_velocity(
                p3,
                Chen2022,
                L,
                N,
                ρ_r,
                F_rim,
                F_liq,
                ρ_a,
                true,
            )
            # number weighted
            TT.@test vel[1] > 0
            TT.@test vel_ϕ[1] > 0
            TT.@test vel[1] ≈ ref_v_n[k] rtol = 1e-6
            TT.@test vel_ϕ[1] ≈ ref_v_n_ϕ[k] rtol = 1e-6

            # mass weighted
            TT.@test vel[2] > 0
            TT.@test vel_ϕ[2] > 0
            TT.@test vel[2] ≈ ref_v_m[k] rtol = 1e-6
            TT.@test vel_ϕ[2] ≈ ref_v_m_ϕ[k] rtol = 1e-6

            # slower with aspect ratio
            TT.@test vel_ϕ[1] <= vel[1]
            TT.@test vel_ϕ[2] <= vel[2]
        end
    end
    TT.@testset "Mass-weighted mean diameters" begin
        F_liq = FT(0)
        ref_vals = [0.00536934536778844, 0.003321590762128599]
        for j in 1:length(F_rims)
            Dₘ = P3.D_m(p3, L, N, ρ_r, F_rims[j], F_liq)
            TT.@test Dₘ > 0
            TT.@test Dₘ ≈ ref_vals[j]
        end

        # nonzero F_liq
        F_rim = F_rims[2]
        F_liqs = [FT(0.33), FT(1)]
        ref_vals = [FT(0.0021371920600012184), FT(0.0016487352655895715)]
        for i in eachindex(F_liqs)
            Dₘ = P3.D_m(p3, L, N, ρ_r, F_rim, F_liqs[i])
            TT.@test Dₘ ≈ ref_vals[i]
        end
    end
end

function test_numerical_integrals(FT)
    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    N = FT(1e8)
    Ls = range(FT(0.001), stop = FT(0.005), length = 5)
    ρ_r = FT(500)
    F_rims = [FT(0), FT(0.5)]
    ρ_a = FT(1.2)
    use_aspect_ratio = false
    tolerance = eps(FT)
    F_liq = FT(0) # test F_liq = 0

    TT.@testset "Numerical integrals sanity checks for N, velocity and diameter" begin
        for F_rim in F_rims
            for i in axes(Ls, 1)
                L = Ls[i]

                # Get shape parameters, thresholds and intergal bounds
                λ, N_0 = P3.distribution_parameter_solver(
                    p3,
                    L,
                    N,
                    ρ_r,
                    F_rim,
                    F_liq,
                )
                th = P3.thresholds(p3, ρ_r, F_rim)
                ice_bound = P3.get_ice_bound(p3, λ, N, tolerance)

                # Number concentration comparison
                f(d) = P3.N′ice(p3, d, λ, N_0)
                qgk_N, = QGK.quadgk(d -> f(d), FT(0), ice_bound)
                TT.@test N ≈ qgk_N rtol = sqrt(eps(FT))

                # Bulk velocity comparison
                vel_N, vel_m = P3.ice_terminal_velocity(
                    p3,
                    Chen2022,
                    L,
                    N,
                    ρ_r,
                    F_rim,
                    F_liq,
                    ρ_a,
                    use_aspect_ratio,
                )

                particle_vel(d) = P3.ice_particle_terminal_velocity(
                    p3,
                    d,
                    Chen2022.snow_ice,
                    ρ_a,
                    F_rim,
                    th,
                    use_aspect_ratio,
                )
                g(d) = particle_vel(d) * P3.N′ice(p3, d, λ, N_0)

                qgk_vel_N, = QGK.quadgk(d -> g(d) / N, FT(0), ice_bound)
                qgk_vel_m, = QGK.quadgk(
                    d -> g(d) * P3.p3_mass(p3, d, F_rim, F_liq, th) / L,
                    FT(0),
                    ice_bound,
                )
                TT.@test vel_N ≈ qgk_vel_N rtol = 1e-6
                TT.@test vel_m ≈ qgk_vel_m rtol = 1e-6

                # Dₘ comparisons
                D_m = P3.D_m(p3, L, N, ρ_r, F_rim, F_liq)
                h(d) =
                    d *
                    P3.p3_mass(p3, d, F_rim, F_liq, th) *
                    P3.N′ice(p3, d, λ, N_0)
                qgk_D_m, = QGK.quadgk(d -> h(d) / L, FT(0), ice_bound)
                TT.@test D_m ≈ qgk_D_m rtol = 1e-2
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
#end

function test_p3_het_freezing(FT)

    TT.@testset "Heterogeneous Freezing Smoke Test" begin
        tps = TD.Parameters.ThermodynamicsParameters(FT)
        p3 = CMP.ParametersP3(FT)
        SB2006 = CMP.SB2006(FT, false) # no limiters
        pdf_c = SB2006.pdf_c
        aerosol = CMP.Illite(FT)

        N_liq = FT(1e8)
        T = FT(244)
        p = FT(500 * 1e2)

        dt = FT(1)

        expected_freeze_L =
            [1.544e-22, 1.068e-6, 0.0001428, 0.0001428, 0.0001428, 0.0001428]
        expected_freeze_N = [1.082e-10, 747647.5, N_liq, N_liq, N_liq, N_liq]
        qᵥ_range = range(FT(0.5e-3), stop = FT(1.5e-3), length = 6)

        for it in range(1, 6)
            qₚ = TD.PhasePartition(FT(qᵥ_range[it]), FT(2e-4), FT(0))
            eᵥ_sat = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
            ϵ = tps.molmass_water / tps.molmass_dryair
            eᵥ = p * qᵥ_range[it] / (ϵ + qᵥ_range[it] * (1 - ϵ))
            RH = eᵥ / eᵥ_sat
            ρₐ = TD.air_density(tps, T, p, qₚ)
            rate = P3.het_ice_nucleation(aerosol, tps, qₚ, N_liq, RH, T, ρₐ, dt)

            TT.@test rate.dNdt >= 0
            TT.@test rate.dLdt >= 0

            TT.@test rate.dNdt ≈ expected_freeze_N[it] rtol = 1e-2
            TT.@test rate.dLdt ≈ expected_freeze_L[it] rtol = 1e-2
        end
    end
end

function test_p3_melting(FT)

    TT.@testset "Melting Smoke Test" begin

        p3 = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        tps = TD.Parameters.ThermodynamicsParameters(FT)

        ρₐ = FT(1.2)
        qᵢ = FT(1e-4)
        Lᵢ = qᵢ * ρₐ
        Nᵢ = FT(2e5) * ρₐ
        F_rim = FT(0.8)
        ρ_rim = FT(800)
        dt = FT(1)
        F_liq = FT(0) # process not dependent on F_liq

        T_cold = FT(273.15 - 0.01)
        rate = P3.ice_melt(
            p3,
            vel.snow_ice,
            aps,
            tps,
            Lᵢ,
            Nᵢ,
            T_cold,
            ρₐ,
            F_rim,
            ρ_rim,
            dt,
        )

        TT.@test rate.dNdt == 0
        TT.@test rate.dLdt == 0

        T_warm = FT(273.15 + 0.01)
        rate = P3.ice_melt(
            p3,
            vel.snow_ice,
            aps,
            tps,
            Lᵢ,
            Nᵢ,
            T_warm,
            ρₐ,
            F_rim,
            ρ_rim,
            dt,
        )

        TT.@test rate.dNdt >= 0
        TT.@test rate.dLdt >= 0

        TT.@test rate.dNdt ≈ 172032.67844162227 rtol = 1e-6
        TT.@test rate.dLdt ≈ 8.601633922081113e-5 rtol = 1e-6

        T_vwarm = FT(273.15 + 0.1)
        rate = P3.ice_melt(
            p3,
            vel.snow_ice,
            aps,
            tps,
            Lᵢ,
            Nᵢ,
            T_vwarm,
            ρₐ,
            F_rim,
            ρ_rim,
            dt,
        )

        TT.@test rate.dNdt == Nᵢ
        TT.@test rate.dLdt == Lᵢ
    end
end

for FT in [Float32, Float64]
    @info("Testing " * string(FT))

    # numerics
    test_incomplete_gamma_function_approximation(FT)
    test_thresholds_solver(FT)

    # velocity
    test_particle_terminal_velocities(FT)

    # processes
    test_p3_het_freezing(FT)
end

# TODO - make work for Float32
for FT in [Float64]
    @info("Tests that only work for Float64")

    # numerics
    test_shape_solver(FT)
    test_numerical_integrals(FT)

    # velocity
    test_bulk_terminal_velocities(FT)

    # processes
    test_p3_melting(FT)
end
