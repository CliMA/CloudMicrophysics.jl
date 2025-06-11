using Test: @testset, @test, @test_throws, @test_broken
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Common as CO
import Thermodynamics as TD
import ClimaParams as CP

@info("P3 Scheme Tests")

function test_p3_state_creation(FT)
    @testset "P3State Creation and Properties" begin
        # Test creating a state with valid parameters
        params = CMP.ParametersP3(FT)

        # Test unrimed state
        state_unrimed = P3.get_state(params; F_rim = FT(0.0), ρ_r = FT(400))
        @test P3.isunrimed(state_unrimed)

        # Test rimed state
        state_rimed = P3.get_state(params; F_rim = FT(0.5), ρ_r = FT(400))
        @test !P3.isunrimed(state_rimed)

        # Test thresholds for unrimed state
        @test length(P3.threshold_tuple(state_unrimed)) == 1
        @test P3.threshold_tuple(state_unrimed)[1] == state_unrimed.D_th

        # Test thresholds for rimed state
        @test length(P3.threshold_tuple(state_rimed)) == 3
        @test P3.threshold_tuple(state_rimed) ==
              (state_rimed.D_th, state_rimed.D_gr, state_rimed.D_cr)

        # Test parameter boundary validation
        @test_throws AssertionError P3.get_state(
            params;
            F_rim = FT(-0.1),
            ρ_r = FT(400),
        )
        @test_throws AssertionError P3.get_state(
            params;
            F_rim = FT(1),
            ρ_r = FT(400),
        )
        @test_throws AssertionError P3.get_state(
            params;
            F_rim = FT(0.5),
            ρ_r = FT(-400),
        )
    end
end

function test_thresholds_solver(FT)

    params = CMP.ParametersP3(FT)

    @testset "Thresholds - nonlinear solver" begin

        # initialize test values:
        ρ_r = FT(400)
        F_rim = FT(0.8)
        ρ_r_good = (FT(200), FT(400), FT(800)) # representative ρ_r values
        F_rim_good = (FT(0.5), FT(0.8), FT(0.95)) # representative F_rim values

        # test asserts
        for _ρ_r in (FT(0), FT(-1), params.ρ_l + 1)
            @test_throws AssertionError("0 < ρ_r <= ρ_l") P3.get_state(
                params;
                F_rim,
                ρ_r = _ρ_r,
            )
        end

        for _F_rim in (FT(-eps(FT)), FT(-1), FT(1), FT(1.5))
            @test_throws AssertionError(
                "Rime mass fraction, `F_rim`, must be between 0 and 1",
            ) P3.get_state(params; F_rim = _F_rim, ρ_r)
        end

        # Test if the P3 scheme solution satisifies the conditions
        # from eqs. 14-17 in Morrison and Milbrandt 2015
        function get_ρ_d_paper((; α_va, β_va)::CMP.MassPowerLaw; D_cr, D_gr)
            # This is Eq. 17 in Morrison and Milbrandt 2015
            βm2 = β_va - 2
            num = 6 * α_va * (D_cr^βm2 - D_gr^βm2)
            den = π * βm2 * (D_cr - D_gr)
            return num / den
        end

        for F_rim in F_rim_good
            for ρ_r in ρ_r_good
                state = P3.get_state(params; F_rim, ρ_r)

                ρ_d = P3.get_ρ_d(params.mass, F_rim, ρ_r)
                @test get_ρ_d_paper(
                    params.mass;
                    D_cr = state.D_cr,
                    D_gr = state.D_gr,
                ) ≈ ρ_d
                @test P3.get_ρ_g(ρ_r, F_rim, ρ_d) ≈ state.ρ_g
                @test P3.get_D_cr(params.mass, state.ρ_g, F_rim) ≈ state.D_cr
                @test P3.get_D_gr(params.mass, state.ρ_g) ≈ state.D_gr
            end
        end

        # Check that the P3 scheme solution matches the published values
        function diff(ρ_r, F_rim, el, gold, rtol = 2e-2)
            state = P3.get_state(params; F_rim, ρ_r)
            @test getfield(state, el) * 1e3 ≈ gold rtol = rtol
        end
        # D_cr and D_gr vs Fig. 1a Morrison and Milbrandt 2015
        D_cr_fig_1a_ref = [FT(0.4946323381999426), FT(1.0170979628696817)]
        D_gr_fig_1a_ref = [FT(0.26151186272014415), FT(0.23392868352755775)]
        for val in [1, 2]
            diff(ρ_r_good[2], F_rim_good[val], :D_cr, D_cr_fig_1a_ref[val])
            diff(ρ_r_good[2], F_rim_good[val], :D_gr, D_gr_fig_1a_ref[val])
        end
        # D_cr and D_gr vs Fig. 1b Morrison and Milbrandt 2015
        #! format: off
        D_cr_fig_1b_ref = [FT(6.152144691917768), FT(3.2718818175768405), FT(1.7400778369620664)]
        D_gr_fig_1b_ref = [FT(0.39875043123651077), FT(0.2147085163169669), FT(0.11516682512848)]
        #! format: on
        for val in [1, 2, 3]
            # TODO: fix this. Where do the reference values come from? They are close to one digit only.
            # diff(ρ_r_good[val], F_rim_good[3], :D_cr, D_cr_fig_1b_ref[val])
            # diff(ρ_r_good[val], F_rim_good[3], :D_gr, D_gr_fig_1b_ref[val])
        end
    end

    @testset "Thresholds - mass, area, density, aspect ratio" begin
        # values
        ρ_r = FT(500)
        F_rim = FT(0.5)

        # get thresholds
        state = P3.get_state(params; F_rim, ρ_r)
        (; D_th, D_gr, D_cr) = state

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        @test P3.ice_area(state, D_1) == P3.area_spherical(D_1)
        @test P3.ice_area(state, D_2) == P3.area_nonspherical(params.area, D_2)
        @test P3.ice_area(state, D_3) == P3.area_spherical(D_3)
        @test P3.ice_area(state, D_cr) ==
              P3.area_rimed(params.area, F_rim, D_cr)

        # test mass
        @test P3.ice_mass(state, D_1) == P3.mass_spherical(params.ρ_i, D_1)
        @test P3.ice_mass(state, D_2) == P3.mass_nonspherical(params.mass, D_2)
        @test P3.ice_mass(state, D_3) == P3.mass_spherical(state.ρ_g, D_3)
        @test P3.ice_mass(state, D_cr) ==
              P3.mass_rimed(params.mass, D_cr, F_rim)

        # test density
        @test P3.ice_density(state, D_1) ≈ params.ρ_i
        @test P3.ice_density(state, D_2) ≈ 544.916989830
        @test P3.ice_density(state, D_3) ≈ state.ρ_g
        @test P3.ice_density(state, D_cr) ≈ 383.33480937

        # test aspect ratio
        @test P3.ϕᵢ(state, D_1) ≈ 1
        @test P3.ϕᵢ(state, D_2) ≈ 1
        @test P3.ϕᵢ(state, D_3) ≈ 1
        @test P3.ϕᵢ(state, D_cr) ≈ 1

        # test F_rim = 0 and D > D_th
        state = P3.get_state(params; F_rim = FT(0), ρ_r)
        @test P3.ice_area(state, D_2) == P3.area_nonspherical(params.area, D_2)
        @test P3.ice_mass(state, D_2) == P3.mass_nonspherical(params.mass, D_2)

        # TODO: Add tests for F_liq != 0
    end
end

function test_shape_solver(FT)

    slope_laws = (:constant, :powerlaw)
    for slope_law in slope_laws
        params = CMP.ParametersP3(FT; slope_law)

        @testset "Shape parameters - nonlinear solver" begin

            # initialize test values:
            ep = 1 #1e4 * eps(FT)
            N_test = (FT(1e7), FT(1e8), FT(1e9), FT(1e10))                         # N values
            λ_test = (FT(1e1), FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))        # test λ values in range also do 15000, 20000
            ρ_r_test = (FT(200), FT(400), FT(600), FT(800))                        # representative ρ_r values
            F_rim_test = (FT(0), FT(0.5), FT(0.8), FT(0.95))                       # representative F_rim values

            # TODO: Add tests for F_liq != 0
            # F_liq_test = (FT(0), FT(0.33), FT(0.67), FT(1))                        # representative F_rim values

            # check that the shape solution solves to give correct values
            for N in N_test
                for λ_ex in λ_test
                    for ρ_r in ρ_r_test
                        for F_rim in F_rim_test
                            # for F_liq in F_liq_test

                            state = P3.get_state(params; F_rim, ρ_r)
                            # Compute the shape parameters that correspond to the input test values
                            logλ_ex = log(λ_ex)
                            logN₀_ex = P3.get_log_N₀(state; N, log_λ = logλ_ex)
                            # Compute mass density based on input shape parameters
                            L_calc = exp(log(N) + P3.log_L_div_N(state, logλ_ex))

                            if L_calc < FT(1)
                                # Solve for shape parameters
                                (; log_λ, log_N₀) =
                                    P3.get_distribution_parameters(
                                        state;
                                        L = L_calc,
                                        N,
                                    )

                                # Compare solved values with the input expected values
                                @test log_λ ≈ logλ_ex rtol = ep
                                @test log_N₀ ≈ logN₀_ex rtol = ep
                            end
                        end
                    end
                end
            end
        end
    end
end

function test_particle_terminal_velocities(FT)

    params = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    ρ_a = FT(1.2)

    @testset "Smoke tests for cloud/rain particle terminal vel from Chen 2022" begin
        Ds = range(FT(1e-6), stop = FT(1e-5), length = 5)  # TODO: Add tests for larger sizes
        expected = [0.002508, 0.009156, 0.01632, 0.02377, 0.03144]
        v_term = CO.particle_terminal_velocity(Chen2022.rain, ρ_a)
        for i in axes(Ds, 1)
            vel = v_term(Ds[i])
            @test vel >= 0
            @test vel ≈ expected[i] rtol = 1e-3
        end
    end

    @testset "Smoke tests for ice particle terminal vel from Chen 2022" begin
        F_rim = FT(0.5)
        ρ_r = FT(500)
        state = P3.get_state(params; F_rim, ρ_r)
        use_aspect_ratio = false
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.7912, 1.1550, 1.4871]
        v_term = P3.ice_particle_terminal_velocity(state, Chen2022, ρ_a; use_aspect_ratio)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test vel ≈ expected[i] rtol = 1e-3
        end

        use_aspect_ratio = true
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.79121, 1.155, 1.487]
        v_term = P3.ice_particle_terminal_velocity(state, Chen2022, ρ_a; use_aspect_ratio)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test vel ≈ expected[i] rtol = 1e-3
        end
    end

    @testset "Smoke tests for mixed phase particle terminal velocity" begin
        F_rim = FT(0.5)
        F_liq = FT(0.5)  # TODO: Broken test since it assumes `F_liq != 0`
        ρ_r = FT(500)
        state = P3.get_state(params; F_rim, ρ_r)
        use_aspect_ratio = true
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13192, 0.50457, 0.90753, 1.3015, 1.6757]
        v_term = P3.ice_particle_terminal_velocity(state, Chen2022, ρ_a; use_aspect_ratio)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test_broken vel ≈ expected[i] rtol = 1e-3  # TODO: Implement `F_liq != 0`
        end
        use_aspect_ratio = false
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13191, 0.50457, 0.90753, 1.301499, 1.67569]
        v_term = P3.ice_particle_terminal_velocity(state, Chen2022, ρ_a; use_aspect_ratio)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test_broken vel ≈ expected[i] rtol = 1e-3  # TODO: Implement `F_liq != 0`
        end
    end
end

function test_bulk_terminal_velocities(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    params = CMP.ParametersP3(FT)
    L = FT(0.22)
    N = FT(1e6)
    ρ_a = FT(1.2)
    ρ_r = FT(800)
    F_rims = FT[0, 0.6]

    # TODO: Implement `F_liq != 0`. The tests break below since they expect `F_liq != 0`
    # F_liqs = [FT(0.5), FT(1)]

    @testset "Mass and number weighted terminal velocities" begin

        state₀ = P3.get_state(params; F_rim = FT(0.5), ρ_r)
        dist₀ = P3.get_distribution_parameters(state₀; L = FT(0), N)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(dist₀, Chen2022, ρ_a)
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(dist₀, Chen2022, ρ_a)
        @test iszero(vel_n₀)
        @test iszero(vel_m₀)

        state₀ = P3.get_state(params; F_rim = FT(0.5), ρ_r)
        dist₀ = P3.get_distribution_parameters(state₀; L, N = FT(0))
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(dist₀, Chen2022, ρ_a)
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(dist₀, Chen2022, ρ_a)
        @test iszero(vel_n₀)
        @test iszero(vel_m₀)

        # NOTE: All reference values are output from the code.
        # A failing test indicates that the code has changed.
        # But if the changes are intentional, the reference values can be updated.

        # Liquid fraction = 0

        ref_v_n = [3.6459501724701835, 2.619088950728837]
        ref_v_n_ϕ = [3.6459501724701835, 2.619088950728837]
        ref_v_m = [7.7881075425985085, 5.797674122909204]
        ref_v_m_ϕ = [7.7881075425985085, 5.797674122909204]

        for (k, F_rim) in enumerate(F_rims)
            state = P3.get_state(params; F_rim, ρ_r)
            dist = P3.get_distribution_parameters(state; L, N)
            args = (dist, Chen2022, ρ_a)
            accurate = true
            vel_n = P3.ice_terminal_velocity_number_weighted(args...; use_aspect_ratio = false, accurate)
            vel_m = P3.ice_terminal_velocity_mass_weighted(args...; use_aspect_ratio = false, accurate)
            vel_n_ϕ = P3.ice_terminal_velocity_number_weighted(args...; use_aspect_ratio = true, accurate)
            vel_m_ϕ = P3.ice_terminal_velocity_mass_weighted(args...; use_aspect_ratio = true, accurate)

            # number weighted
            @test vel_n > 0
            @test vel_n_ϕ > 0
            @test vel_n ≈ ref_v_n[k] rtol = 1e-6
            @test vel_n_ϕ ≈ ref_v_n_ϕ[k] rtol = 1e-6

            # mass weighted
            @test vel_m > 0
            @test vel_m_ϕ > 0
            @test vel_m ≈ ref_v_m[k] rtol = 5e-5
            @test vel_m_ϕ ≈ ref_v_m_ϕ[k] rtol = 5e-5

            # slower with aspect ratio (within machine precision)
            @test vel_n_ϕ <= vel_n + eps(vel_n)
            @test vel_m_ϕ <= vel_m + eps(vel_m)
        end

        # Liquid fraction != 0
        ref_v_n = [1.674591925057434, 1.6180970319460353]
        ref_v_n_ϕ = [1.674591925057434, 1.6180970319460353]
        #ref_v_n_ϕ = [1.549777478756061, 1.6180970319460353]
        ref_v_m = [5.126648558302173, 5.416679316254198]
        ref_v_m_ϕ = [5.126648558302173, 5.416679316254198]
        #ref_v_m_ϕ = [4.6358422594886495, 5.416679316254198]

        # TODO: Add tests for F_liq != 0
        # for k = 1:length(F_liqs)
        #     F_liq = F_liqs[k]
        #     F_rim = FT(0.4)
        #     vel =
        #         P3.ice_terminal_velocity(p3, Chen2022, L, N, ρ_r, F_rim, F_liq, ρ_a, false)
        #     vel_ϕ =
        #         P3.ice_terminal_velocity(p3, Chen2022, L, N, ρ_r, F_rim, F_liq, ρ_a, true)
        #     # number weighted
        #     @test vel[1] > 0
        #     @test vel_ϕ[1] > 0
        #     @test vel[1] ≈ ref_v_n[k] rtol = 1e-6
        #     @test vel_ϕ[1] ≈ ref_v_n_ϕ[k] rtol = 1e-6

        #     # mass weighted
        #     @test vel[2] > 0
        #     @test vel_ϕ[2] > 0
        #     @test vel[2] ≈ ref_v_m[k] rtol = 1e-6
        #     @test vel_ϕ[2] ≈ ref_v_m_ϕ[k] rtol = 1e-6

        #     # slower with aspect ratio
        #     @test vel_ϕ[1] <= vel[1]
        #     @test vel_ϕ[2] <= vel[2]
        # end
    end
    @testset "Mass-weighted mean diameters" begin
        ref_vals = [0.005397144197921535, 0.0033368960364578005]
        for (F_rim, ref_val) in zip(F_rims, ref_vals)
            state = P3.get_state(params; F_rim, ρ_r)
            dist = P3.get_distribution_parameters(state; L, N)
            Dₘ = P3.D_m(dist)
            @test Dₘ > 0
            @test Dₘ ≈ ref_val
        end

        # TODO: Add tests for F_liq != 0
        # nonzero F_liq
        # F_rim = F_rims[2]
        # F_liqs = [FT(0.33), FT(1)]
        # ref_vals = [FT(0.0021371920600012184), FT(0.0016487352655895715)]
        # for i in eachindex(F_liqs)
        #     Dₘ = P3.D_m(p3, L, N, ρ_r, F_rim, F_liqs[i])
        #     @test Dₘ ≈ ref_vals[i]
        # end
    end
end

function test_numerical_integrals(FT)
    params = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    N = FT(1e8)
    Ls = range(FT(0.001), stop = FT(0.005), length = 5)
    ρ_r = FT(500)
    F_rims = FT[0, 0.5]
    ρ_a = FT(1.2)
    use_aspect_ratio = false
    ps = [1e-3, 1e-6]

    @testset "Numerical integrals sanity checks for N, velocity and diameter" begin
        for (F_rim, L, p) in Iterators.product(F_rims, Ls, ps)

            # Get shape parameters, thresholds and intergal bounds
            state = P3.get_state(params; F_rim, ρ_r)
            dist = P3.get_distribution_parameters(state; L, N)

            # Number concentration comparison
            # Note: To achieve sufficient accuracy, we need to substantially
            # increase the `order` of the quadrature rule, and set `rtol=0`.
            # The `rtol` settings essentially forces max evaluations of the method.
            # Note 2: For F_rim=0, L=0.002, even higher order quadrature rules are needed.
            N_estim = P3.∫fdD(dist) do D
                P3.N′ice(dist, D)
            end
            @test N ≈ N_estim rtol = 1e-5

            # Bulk velocity comparison
            vel_N = P3.ice_terminal_velocity_number_weighted(dist, Chen2022, ρ_a; use_aspect_ratio, p)
            vel_m = P3.ice_terminal_velocity_mass_weighted(dist, Chen2022, ρ_a; use_aspect_ratio, p)

            v_term = P3.ice_particle_terminal_velocity(state, Chen2022, ρ_a; use_aspect_ratio)
            g(D) = v_term(D) * exp(P3.log_N′ice(dist, D))
            vel_N_estim = P3.∫fdD(g, dist; p) / N
            vel_m_estim = P3.∫fdD(dist; p) do D
                g(D) * P3.ice_mass(state, D) / L
            end

            @test vel_N ≈ vel_N_estim rtol = 1e-6
            @test vel_m ≈ vel_m_estim rtol = 1e-5

            # Dₘ comparisons
            D_m = P3.D_m(dist)
            D_m_estim = P3.∫fdD(dist; moment_order = 1 + 3) do D
                D * P3.ice_mass(state, D) * P3.N′ice(dist, D) / L
            end
            @test D_m ≈ D_m_estim rtol = 5e-4
        end
    end
end

function test_p3_het_freezing(FT)

    @testset "Heterogeneous Freezing Smoke Test" begin
        tps = TD.Parameters.ThermodynamicsParameters(FT)
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

            @test rate.dNdt >= 0
            @test rate.dLdt >= 0

            @test rate.dNdt ≈ expected_freeze_N[it] rtol = 1e-2
            @test rate.dLdt ≈ expected_freeze_L[it] rtol = 1e-2
        end
    end
end

function test_p3_melting(FT)

    @testset "Melting Smoke Test" begin

        params = CMP.ParametersP3(FT)
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

        state = P3.get_state(params; F_rim, ρ_r = ρ_rim)
        dist = P3.get_distribution_parameters(state; L = Lᵢ, N = Nᵢ)

        T_cold = FT(273.15 - 0.01)

        rate = P3.ice_melt(dist, vel, aps, tps, T_cold, ρₐ, dt)

        @test rate.dNdt == 0
        @test rate.dLdt == 0

        T_warm = FT(273.15 + 0.01)
        rate = P3.ice_melt(dist, vel, aps, tps, T_warm, ρₐ, dt)

        @test rate.dNdt >= 0
        @test rate.dLdt >= 0

        # NOTE: All reference values are output from the code.
        # A failing test indicates that the code has changed.
        # But if the changes are intentional, the reference values can be updated.
        if FT == Float64
            ref_dNdt = FT(171951.91644755047)
            ref_dLdt = FT(8.597595822377523e-5)
        else
            ref_dNdt = FT(172119.73f0)
            ref_dLdt = FT(8.605987f-5)
        end
        @test rate.dNdt ≈ ref_dNdt
        @test rate.dLdt ≈ ref_dLdt

        T_vwarm = FT(273.15 + 0.1)
        rate = P3.ice_melt(dist, vel, aps, tps, T_vwarm, ρₐ, dt)

        @test rate.dNdt == Nᵢ
        @test rate.dLdt == Lᵢ
    end
end


for FT in [Float32, Float64]
    @info("Testing " * string(FT))

    # state creation
    test_p3_state_creation(FT)

    # numerics
    test_thresholds_solver(FT)
    test_shape_solver(FT)
    test_numerical_integrals(FT)

    # velocity
    test_particle_terminal_velocities(FT)
    test_bulk_terminal_velocities(FT)

    # processes
    test_p3_het_freezing(FT)
    test_p3_melting(FT)
end
