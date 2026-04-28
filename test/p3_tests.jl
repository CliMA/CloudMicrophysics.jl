using Test: @testset, @test, @test_throws, @test_broken
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Common as CO
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams as CP
import SpecialFunctions as SF
import QuadGK as QGK

function test_p3_state_creation(FT)
    @testset "P3State Creation and Properties" begin
        # Test creating a state with valid parameters
        params = CMP.ParametersP3(FT)
        L_ice = FT(0.22)
        N_ice = FT(1e6)
        F_rim = FT(0.5)
        ρ_rim = FT(400)

        # Test unrimed state
        state_unrimed = P3.get_state(params; F_rim = FT(0), ρ_rim, L_ice, N_ice)
        @test P3.isunrimed(state_unrimed)

        # Test rimed state
        state_rimed = P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
        @test !P3.isunrimed(state_rimed)

        # Test thresholds for unrimed state. Per the docstring on
        # `get_thresholds_ρ_g`, unrimed ice has no graupel → `D_gr = D_cr = Inf`
        # (the "always before graupel regime" sentinel) and `ρ_g = NaN`
        # (should not be used by callers).
        (; D_th, D_gr, D_cr) = P3.get_thresholds_ρ_g(state_unrimed)
        @test isfinite(D_th)
        @test D_gr == Inf
        @test D_cr == Inf

        # Test thresholds for rimed state
        (; D_th, D_gr, D_cr) = P3.get_thresholds_ρ_g(state_rimed)
        @test D_th < D_gr < D_cr

        # Test parameter boundary validation
        @test_throws DomainError P3.get_state(params; F_rim = FT(-0.1), ρ_rim, L_ice, N_ice)
        @test_throws DomainError P3.get_state(params; F_rim = FT(1), ρ_rim, L_ice, N_ice)
        @test_throws DomainError P3.get_state(params; F_rim, ρ_rim = FT(-400), L_ice, N_ice)
    end
end

function test_thresholds_solver(FT)

    params = CMP.ParametersP3(FT)

    @testset "Thresholds - exact solution" begin

        # initialize test values:
        ρ_rim = FT(400)
        F_rim = FT(0.8)
        L_ice = FT(0.22)
        N_ice = FT(1e6)
        ρ_rim_good = (FT(200), FT(400), FT(800)) # representative ρ_rim values
        F_rim_good = (FT(0.5), FT(0.8), FT(0.95)) # representative F_rim values

        # test asserts
        for ρ_rim in (FT(-1), params.ρ_l + 1)
            @test_throws DomainError(ρ_rim,
                "Rime density, ρ_rim, must be 0 ≤ ρ_rim ≤ ρ_l",
            ) P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
        end

        for F_rim in (FT(-eps(FT)), FT(-1), FT(1), FT(1.5))
            @test_throws DomainError(F_rim,
                "Rime mass fraction, F_rim, must be 0 ≤ F_rim < 1",
            ) P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
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

        (; mass, ρ_i) = params
        for F_rim in F_rim_good
            for ρ_rim in ρ_rim_good
                (; D_th, D_gr, D_cr) = P3.get_thresholds_ρ_g(params, F_rim, ρ_rim)
                ρ_d = P3.get_ρ_d(mass, F_rim, ρ_rim)
                @test D_th < D_gr < D_cr
                @test get_ρ_d_paper(mass; D_cr, D_gr) ≈ ρ_d
            end
        end

        # For very high rimed density, the thresholds are ill-defined. TODO: Investigate this
        (; mass, ρ_i) = params
        F_rim_bad = FT(0.93)
        ρ_rim_bad = FT(975)
        (; D_th, D_gr, D_cr) = P3.get_thresholds_ρ_g(params, F_rim_bad, ρ_rim_bad)
        @test_broken D_th_bad < D_gr_bad

        # Check that the P3 scheme solution matches the published values
        # D_cr and D_gr vs Fig. 1a Morrison and Milbrandt 2015
        D_cr_fig_1a_ref = FT[0.4946323381999426, 1.0170979628696817]
        D_gr_fig_1a_ref = FT[0.26151186272014415, 0.23392868352755775]
        for i in 1:2
            (; D_cr, D_gr) = P3.get_thresholds_ρ_g(params, F_rim_good[i], ρ_rim_good[2])
            @test 1000 * D_cr ≈ D_cr_fig_1a_ref[i] rtol = 2e-2
            @test 1000 * D_gr ≈ D_gr_fig_1a_ref[i] rtol = 2e-2
        end
        # D_cr and D_gr vs Fig. 1b Morrison and Milbrandt 2015
        # D_cr_fig_1b_ref = FT[6.152144691917768, 3.2718818175768405, 1.7400778369620664]
        # D_gr_fig_1b_ref = FT[0.39875043123651077, 0.2147085163169669, 0.11516682512848]
        # for val in 1:3
        #     # TODO: fix this. Where do the reference values come from? They are close to one digit only.
        #     D_cr = P3.get_D_cr(mass, F_rim_good[3], ρ_rim_good[val])
        #     D_gr = P3.get_D_gr(mass, ρ_g)
        #     @test 1000 * D_cr ≈ D_cr_fig_1b_ref[val] rtol = 2e-2
        #     @test 1000 * D_gr ≈ D_gr_fig_1b_ref[val] rtol = 2e-2
        # end
    end

    @testset "Thresholds - mass, area, density, aspect ratio" begin
        # values
        ρ_rim = FT(500)
        F_rim = FT(0.5)
        L_ice = FT(0.22)
        N_ice = FT(1e6)

        (; area, mass, ρ_i) = params

        # get thresholds
        (; D_th, D_gr, D_cr, ρ_g) = P3.get_thresholds_ρ_g(params, F_rim, ρ_rim)

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        spherical_area(D) = D^2 * π / 4
        nonspherical_area(D) = area.γ * D^area.σ
        @test P3.ice_area(params, F_rim, ρ_rim, D_1) == spherical_area(D_1)
        @test P3.ice_area(params, F_rim, ρ_rim, D_2) == nonspherical_area(D_2)
        @test P3.ice_area(params, F_rim, ρ_rim, D_3) == spherical_area(D_3)
        @test P3.ice_area(params, F_rim, ρ_rim, D_cr) ==
              F_rim * spherical_area(D_cr) + (1 - F_rim) * nonspherical_area(D_cr)

        # test mass
        spherical_mass(ρ, D) = ρ * π / 6 * D^3
        nonspherical_mass(D) = mass.α_va * D^mass.β_va
        @test P3.ice_mass(params, F_rim, ρ_rim, D_1) == spherical_mass(ρ_i, D_1)
        @test P3.ice_mass(params, F_rim, ρ_rim, D_2) == nonspherical_mass(D_2)
        @test P3.ice_mass(params, F_rim, ρ_rim, D_3) == spherical_mass(ρ_g, D_3)
        @test P3.ice_mass(params, F_rim, ρ_rim, D_cr) == nonspherical_mass(D_cr) / (1 - F_rim)

        # test density
        @test P3.ice_density(params, F_rim, ρ_rim, D_1) ≈ ρ_i
        @test P3.ice_density(params, F_rim, ρ_rim, D_2) ≈ 544.916989830
        @test P3.ice_density(params, F_rim, ρ_rim, D_3) ≈ ρ_g
        @test P3.ice_density(params, F_rim, ρ_rim, D_cr) ≈ 383.33480937

        # test aspect ratio
        @test P3.ϕᵢ(params, F_rim, ρ_rim, D_1) ≈ 1
        @test P3.ϕᵢ(params, F_rim, ρ_rim, D_2) ≈ 1
        @test P3.ϕᵢ(params, F_rim, ρ_rim, D_3) ≈ 1
        @test P3.ϕᵢ(params, F_rim, ρ_rim, D_cr) ≈ 1

        # test F_rim = 0 and D > D_th
        @test P3.ice_area(params, FT(0), ρ_rim, D_2) == nonspherical_area(D_2)
        @test P3.ice_mass(params, FT(0), ρ_rim, D_2) == nonspherical_mass(D_2)

        # TODO: Add tests for F_liq != 0
    end
end

function test_shape_solver(FT)

    slope_laws = (:constant, :powerlaw)
    for slope_law in slope_laws
        params = CMP.ParametersP3(FT; slope_law)

        @testset "Shape parameters - nonlinear solver" begin
            # -- First, test limiting behavior: `N_ice = L_ice = 0` --
            state = P3.get_state(params; F_rim = FT(0.5), ρ_rim = FT(500), L_ice = FT(0), N_ice = FT(0))
            logλ = P3.get_distribution_logλ(state)
            @test logλ == -Inf
            # --

            # initialize test values:
            ep = 1 #1e4 * eps(FT)
            N_test = (FT(1e7), FT(1e8), FT(1e9), FT(1e10))                         # N values
            λ_test = (FT(1e1), FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))        # test λ values in range also do 15000, 20000
            ρ_rim_test = (FT(200), FT(400), FT(600), FT(800))                        # representative ρ_rim values
            F_rim_test = (FT(0), FT(0.5), FT(0.8), FT(0.95))                       # representative F_rim values

            # TODO: Add tests for F_liq != 0
            # F_liq_test = (FT(0), FT(0.33), FT(0.67), FT(1))                        # representative F_rim values

            # check that the shape solution solves to give correct values
            for N_ice in N_test
                for λ_ex in λ_test
                    for ρ_rim in ρ_rim_test
                        for F_rim in F_rim_test
                            # for F_liq in F_liq_test

                            state = P3.get_state(params; F_rim, ρ_rim, L_ice = FT(0), N_ice = FT(0)) # L_ice, N_ice not used in this test
                            # Compute the shape parameters that correspond to the input test values
                            logλ_ex = log(λ_ex)
                            μ = P3.get_μ(params.slope, logλ_ex)
                            logN₀_ex = P3.get_logN₀(N_ice, μ, logλ_ex)
                            # Compute mass density based on input shape parameters
                            L_calc = exp(log(N_ice) + P3.logLdivN(state, logλ_ex))

                            if L_calc < FT(1)
                                # Solve for shape parameters
                                logλ = P3.get_distribution_logλ(params, L_calc, N_ice, F_rim, ρ_rim)
                                log_N₀ = P3.get_logN₀(N_ice, μ, logλ)

                                # Compare solved values with the input expected values
                                @test logλ ≈ logλ_ex rtol = ep
                                @test log_N₀ ≈ logN₀_ex rtol = ep
                            end
                        end
                    end
                end
            end
        end

        @testset "Shape solver - robustness across physical inputs" begin
            params = CMP.ParametersP3(FT)

            # Regression test: this specific `(L_ice, N_ice, F_rim, ρ_rim)`
            # triggered a NaN return under the previous `SecantMethod`-based
            # solver because a secant step extrapolated outside the search
            # interval into a region where `logLdivN` is not finite. The
            # bracketing `BrentsMethod` must return a finite, positive
            # `logλ` strictly inside the search bounds.
            logλ = P3.get_distribution_logλ(
                params, FT(2.366e-5), FT(16461.6), FT(0.2), FT(800),
            )
            @test isfinite(logλ)
            @test FT(2) < logλ < FT(17)

            # Broader sweep covering typical P3 microphysics inputs.
            # All entries must give a finite `logλ` within the search bounds.
            for L_ice in (FT(1e-6), FT(1e-5), FT(2.366e-5), FT(1e-4), FT(1e-3))
                for N_ice in (FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))
                    for F_rim in (FT(0), FT(0.2), FT(0.5), FT(0.8), FT(0.95))
                        for ρ_rim in (FT(200), FT(400), FT(600), FT(800))
                            logλ = P3.get_distribution_logλ(
                                params, L_ice, N_ice, F_rim, ρ_rim,
                            )
                            @test isfinite(logλ)
                            @test FT(2) ≤ logλ ≤ FT(17)
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
        ρ_rim = FT(500)
        state = P3.get_state(params; F_rim, ρ_rim, L_ice = FT(0), N_ice = FT(0))
        use_aspect_ratio = false
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.7912, 1.1550, 1.4871]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state; use_aspect_ratio)
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
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state; use_aspect_ratio)
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
        ρ_rim = FT(500)
        state = P3.get_state(params; F_rim, ρ_rim, L_ice = FT(0), N_ice = FT(0))
        use_aspect_ratio = true
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13192, 0.50457, 0.90753, 1.3015, 1.6757]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state; use_aspect_ratio)
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
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state; use_aspect_ratio)
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
    L_ice = FT(0.22)
    N_ice = FT(1e6)
    ρ_a = FT(1.2)
    ρ_rim = FT(800)
    F_rims = FT[0, 0.6]

    # TODO: Implement `F_liq != 0`. The tests break below since they expect `F_liq != 0`
    # F_liqs = [FT(0.5), FT(1)]

    @testset "Mass and number weighted terminal velocities" begin

        state₀ = P3.get_state(params; F_rim = FT(0.5), ρ_rim, L_ice = FT(0), N_ice)
        logλ = P3.get_distribution_logλ(state₀)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state₀, logλ)
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state₀, logλ)
        @test iszero(vel_n₀)
        @test iszero(vel_m₀)

        state₀ = P3.get_state(params; F_rim = FT(0.5), ρ_rim, L_ice, N_ice = FT(0))
        logλ = P3.get_distribution_logλ(state₀)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state₀, logλ)
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state₀, logλ)
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
            state = P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
            logλ = P3.get_distribution_logλ(state)
            args = (Chen2022, ρ_a, state, logλ)
            vel_n = P3.ice_terminal_velocity_number_weighted(args...; use_aspect_ratio = false)
            vel_m = P3.ice_terminal_velocity_mass_weighted(args...; use_aspect_ratio = false)
            vel_n_ϕ = P3.ice_terminal_velocity_number_weighted(args...; use_aspect_ratio = true)
            vel_m_ϕ = P3.ice_terminal_velocity_mass_weighted(args...; use_aspect_ratio = true)

            # number weighted
            @test vel_n > 0
            @test vel_n_ϕ > 0
            @test vel_n ≈ ref_v_n[k] rtol = 1e-4
            @test vel_n_ϕ ≈ ref_v_n_ϕ[k] rtol = 1e-4

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
        #         P3.ice_terminal_velocity(p3, Chen2022, L, N, ρ_rim, F_rim, F_liq, ρ_a, false)
        #     vel_ϕ =
        #         P3.ice_terminal_velocity(p3, Chen2022, L, N, ρ_rim, F_rim, F_liq, ρ_a, true)
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
            state = P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
            logλ = P3.get_distribution_logλ(state)
            Dₘ = P3.D_m(state, logλ)
            @test Dₘ > 0
            @test Dₘ ≈ ref_val
        end

        # TODO: Add tests for F_liq != 0
        # nonzero F_liq
        # F_rim = F_rims[2]
        # F_liqs = [FT(0.33), FT(1)]
        # ref_vals = [FT(0.0021371920600012184), FT(0.0016487352655895715)]
        # for i in eachindex(F_liqs)
        #     Dₘ = P3.D_m(p3, L, N, ρ_rim, F_rim, F_liqs[i])
        #     @test Dₘ ≈ ref_vals[i]
        # end
    end
end

function test_numerical_integrals(FT)
    params = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)

    N_ice = FT(1e8)
    L_ices = range(FT(0.001), stop = FT(0.005), length = 5)
    ρ_rim = FT(500)
    F_rims = FT[0, 0.5]
    ρ_a = FT(1.2)
    use_aspect_ratio = false
    ps = [1e-3, 1e-6]

    @testset "Chebyshev-Gauss quadrature" begin
        quad = P3.ChebyshevGauss(10)
        f(x) = x^4
        # test that integration gives the correct result
        num_int = P3.integrate(f, 0, 1, quad)
        @test num_int ≈ 0.2 rtol = 0.1
        # test that increasing the number of points improves the accuracy
        num_int2 = P3.integrate(f, 0, 1, P3.ChebyshevGauss(100))
        @test abs(num_int2 - 0.2) < abs(num_int - 0.2)
    end

    @testset "Numerical integrals sanity checks for N, velocity and diameter" begin
        for (F_rim, L_ice, p) in Iterators.product(F_rims, L_ices, ps)

            # Get shape parameters, thresholds and intergal bounds
            state = P3.get_state(params; F_rim, ρ_rim, L_ice, N_ice)
            logλ = P3.get_distribution_logλ(state)

            # Number concentration comparison
            # Note: To achieve sufficient accuracy, we need to substantially
            # increase the `order` of the quadrature rule, and set `rtol=0`.
            # The `rtol` settings essentially forces max evaluations of the method.
            # Note 2: For F_rim=0, L=0.002, even higher order quadrature rules are needed.
            N′ = P3.size_distribution(state, logλ)
            bnds = P3.integral_bounds(state, logλ; p = 1e-6, moment_order = 0)
            N_estim_cheb = P3.integrate(N′, bnds)
            @test N_ice ≈ N_estim_cheb rtol = 1e-5

            # Compare with quadgk
            N_estim_qgk = QGK.quadgk(N′, bnds...)[1]
            @test N_estim_cheb ≈ N_estim_qgk rtol = 1e-5


            # Bulk velocity comparison
            vel_N = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state, logλ; use_aspect_ratio, p)
            vel_m = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state, logλ; use_aspect_ratio, p)

            v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state; use_aspect_ratio)
            g(D) = v_term(D) * N′(D)
            gm(D) = g(D) * P3.ice_mass(state, D)
            vel_N_estim_cheb = P3.integrate(g, bnds, P3.ChebyshevGauss(10)) / N_ice
            vel_m_estim_cheb = P3.integrate(gm, bnds, P3.ChebyshevGauss(10)) / L_ice
            @test vel_N ≈ vel_N_estim_cheb rtol = 0.005
            @test vel_m ≈ vel_m_estim_cheb rtol = 0.05

            # Compare with quadgk
            vel_N_estim_qgk = QGK.quadgk(g, bnds...)[1] / N_ice
            vel_m_estim_qgk = QGK.quadgk(gm, bnds...)[1] / L_ice

            @test vel_N_estim_cheb ≈ vel_N_estim_qgk rtol = 0.005
            @test vel_m_estim_cheb ≈ vel_m_estim_qgk rtol = 0.05


            # Dₘ comparisons
            D_m = P3.D_m(state, logλ)
            D_m_func(D) = D * P3.ice_mass(state, D) * N′(D) / L_ice
            D_m_estim_cheb = P3.integrate(D_m_func, bnds, P3.ChebyshevGauss(100))
            @test D_m ≈ D_m_estim_cheb rtol = 5e-4

            # Compare with quadgk
            D_m_estim_qgk = QGK.quadgk(D_m_func, bnds...)[1]
            @test D_m_estim_cheb ≈ D_m_estim_qgk rtol = 5e-4
        end
    end
end

function test_p3_het_freezing(FT)

    @testset "Heterogeneous Freezing Smoke Test" begin
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        aerosol = CMP.Illite(FT)

        N_lcl = FT(1e8)
        T = FT(244)
        p = FT(500 * 1e2)

        expected_freeze_L =
            [1.4953923796668346e-22, 1.0365387091217913e-6, 0.0001428, 0.0001428, 0.0001428, 0.0001428]
        expected_freeze_N = [1.0473022910416716e-10, 726031.0899744622, N_lcl, N_lcl, N_lcl, N_lcl]
        qᵥ_range = range(FT(0.5e-3), stop = FT(1.5e-3), length = 6)

        for it in range(1, 6)
            q_lcl = FT(2e-4)
            eᵥ_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)
            ϵ = TDI.Rd_over_Rv(tps)
            eᵥ = p * qᵥ_range[it] / (ϵ + qᵥ_range[it] * (1 - ϵ))
            RH = eᵥ / eᵥ_sat
            ρₐ = TDI.air_density(tps, T, p, qᵥ_range[it] + q_lcl, q_lcl, FT(0))
            rate = P3.het_ice_nucleation(aerosol, tps, q_lcl, N_lcl, RH, T, ρₐ)

            @test rate.dNdt >= 0
            @test rate.dLdt >= 0

            @test rate.dNdt ≈ expected_freeze_N[it] rtol = 2e-2
            @test rate.dLdt ≈ expected_freeze_L[it] rtol = 2e-2
        end
    end
end

function test_p3_melting(FT)

    @testset "Melting Smoke Test" begin

        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

        ρₐ = FT(1.2)
        qᵢ = FT(1e-4)
        Lᵢ = qᵢ * ρₐ
        Nᵢ = FT(2e5) * ρₐ
        F_rim = FT(0.8)
        ρ_rim = FT(800)

        state = P3.get_state(params; F_rim, ρ_rim, L_ice = Lᵢ, N_ice = Nᵢ)
        logλ = P3.get_distribution_logλ(state)

        T_cold = FT(273.15 - 0.01)

        rate = P3.ice_melt(vel, aps, tps, T_cold, ρₐ, state, logλ)

        @test rate.dNdt == 0
        @test rate.dLdt == 0

        T_warm = FT(273.15 + 0.01)
        rate = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state, logλ)

        @test rate.dNdt >= 0
        @test rate.dLdt >= 0

        # NOTE: All reference values are output from the code.
        # A failing test indicates that the code has changed.
        # But if the changes are intentional, the reference values can be updated.
        if FT == Float64
            ref_dNdt = FT(171964.6981857982)
            ref_dLdt = FT(8.59823490928991e-5)
        else
            ref_dNdt = FT(172119.73f0)
            ref_dLdt = FT(8.605987f-5)
        end
        @test rate.dNdt ≈ ref_dNdt
        @test rate.dLdt ≈ ref_dLdt

        T_vwarm = FT(273.15 + 0.1)
        rate = P3.ice_melt(vel, aps, tps, T_vwarm, ρₐ, state, logλ)

        @test rate.dNdt == Nᵢ
        @test rate.dLdt == Lᵢ
    end
end

function test_p3_bulk_liquid_ice_collisions(FT)
    params = CMP.ParametersP3(FT)
    vel_params = CMP.Chen2022VelType(FT)
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    (; T_freeze) = params

    ρₐ = FT(1.2)
    qᵢ = FT(1e-4)
    Lᵢ = qᵢ * ρₐ
    Nᵢ = FT(2e5) * ρₐ
    F_rim = FT(0.8)
    ρ_rim = FT(800)

    state = P3.get_state(params; F_rim, ρ_rim, L_ice = Lᵢ, N_ice = Nᵢ)
    logλ = P3.get_distribution_logλ(state)
    D̄ = exp(-logλ)

    @testset "maximum dry freezing rate" begin
        # Below freezing, max freeze rate is non-zero (check against reference value)
        Tₐ = T_freeze - 1 // 10
        max_rate = P3.compute_max_freeze_rate(aps, tps, vel_params, ρₐ, Tₐ, state)
        @test max_rate(D̄) ≈ FT(9.35962884896919e-13) rtol = 2e-4

        # At freezing, max freeze rate is zero
        Tₐ = T_freeze
        max_rate = P3.compute_max_freeze_rate(aps, tps, vel_params, ρₐ, Tₐ, state)
        @test iszero(max_rate(D̄))

        # Above freezing, max freeze rate is zero
        Tₐ = T_freeze + 1 // 10
        max_rate = P3.compute_max_freeze_rate(aps, tps, vel_params, ρₐ, Tₐ, state)
        @test iszero(max_rate(D̄))
    end

    @testset "local rime density" begin
        Tₐ = T_freeze - 1 // 10
        ρ′_rim_func = P3.compute_local_rime_density(vel_params, ρₐ, Tₐ, state)
        @test ρ′_rim_func(D̄, D̄) ≈ FT(159.5) rtol = 1e-6

        a, b, c = 51, 114, -11 // 2 # coeffs for Eq. 17 in Cober and List (1993), converted to [kg / m³]
        ρ′_rim_CL93(Rᵢ) = a + b * Rᵢ + c * Rᵢ^2  # Eq. 17 in Cober and List (1993), in [kg / m³], valid for 1 ≤ Rᵢ ≤ 8
        ρ_ice = FT(916.7)  # density of solid bulk ice

        ρ_rim_local = params.ρ_rim_local

        @test ρ_rim_local(1) == ρ′_rim_CL93(1)
        @test ρ_rim_local(8) == ρ′_rim_CL93(8)
        @test ρ_rim_local(12) == ρ_ice
    end

    @testset "∫liquid_ice_collisions" begin
        # Test liquid_integrals function in isolation
        # Mock simple functions for analytical comparison
        ∂ₜV(Dᵢ, D) = Dᵢ * D  # Simple collision rate
        n(D) = exp(-D)       # Simple size distribution
        n_c = n_r = n_i = n  # Mock cloud, rain and ice size distributions
        m_l(D) = D^3         # Simple mass function
        ρ′_rim(Dᵢ, D) = 500     # Constant rime density
        liq_bounds = ice_bounds = (FT(0), FT(1))
        Dᵢ = FT(2.5)
        cloud_integrals = P3.get_liquid_integrals(n_c, ∂ₜV, m_l, ρ′_rim, liq_bounds; quad = P3.ChebyshevGauss(100))
        rain_integrals = P3.get_liquid_integrals(n_r, ∂ₜV, m_l, ρ′_rim, liq_bounds; quad = P3.ChebyshevGauss(100))

        # Test with known analytical result, noting e.g. that:
        # ∫₀¹ Dᵢ * D * exp(-D) * D³ dD = ∫₀¹ D⁴ * exp(-D) dD = γ(5, 1) [lower incomplete gamma function]
        γ(a, x) = SF.gamma_inc(a, x)[1] * SF.gamma(a)

        (∫∂ₜVn, ∫∂ₜVnm, ∫∂ₜVnm_ρ′) = cloud_integrals(Dᵢ)
        @test all(x -> x isa FT, (∫∂ₜVn, ∫∂ₜVnm, ∫∂ₜVnm_ρ′))  # check type stability
        @test ∫∂ₜVn ≈ γ(2, 1) * Dᵢ rtol = 5e-5
        @test ∫∂ₜVnm ≈ γ(5, 1) * Dᵢ rtol = 1e-4
        @test ∫∂ₜVnm_ρ′ ≈ γ(5, 1) * Dᵢ / 500 rtol = 1e-4

        # Test edge cases for liquid_integrals

        # Zero ice diameter
        result = cloud_integrals(FT(0))
        @test all(iszero, result)

        # Zero liquid content (n(D) = 0)
        n_zero(D) = FT(0)
        integrals_n0 = P3.get_liquid_integrals(n_zero, ∂ₜV, m_l, ρ′_rim, liq_bounds)
        result = integrals_n0(Dᵢ)
        @test all(iszero, result)

        # Mass rate should be related to number rate through mass
        # ∂ₜM_col should be approximately ∫ n(D) * m_l(D) * ∂ₜV dD
        # This is a simplified check - in reality it's more complex

        # Test the full ∫liquid_ice_collisions function
        # Mock functions
        ∂ₜM_max(Dᵢ) = FT(0.04)  # Small max freeze rate
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds,
        )
        @test all(x -> x isa FT, rates)  # check type stability
        @test all(>(0), rates)  # check positivity

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col = rates

        # Mass conservation: QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col

        # Wet growth indicator should be ≤ total collision rate
        @test ∫𝟙_wet_M_col <= ∫M_col

        # Since we specified identical size distributions, we expect:
        @test QCFRZ == QRFRZ
        @test QCSHD == QRSHD
        @test NCCOL == NRCOL
        @test BCCOL == BRCOL

        # Test edge cases for full collision integration

        # Zero ice content
        n_i_zero(Dᵢ) = FT(0)
        rates = P3.∫liquid_ice_collisions(
            n_i_zero, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds,
        )
        @test all(iszero, rates)

        # Zero liquid content
        n_zero = Returns(FT(0))
        zero_liq_integrals = P3.get_liquid_integrals(n_zero, ∂ₜV, m_l, ρ′_rim, liq_bounds)
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max, zero_liq_integrals, zero_liq_integrals, ice_bounds,
        )
        @test all(iszero, rates)

        # No freezing (above freezing temperature)
        ∂ₜM_max_zero(Dᵢ) = FT(0)
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max_zero, cloud_integrals, rain_integrals, ice_bounds,
        )
        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col = rates
        @test QCFRZ == 0  # No cloud freezing
        @test QRFRZ == 0  # No rain freezing
        @test QCSHD > 0  # All cloud particles should shed
        @test QRSHD > 0  # All rain particles should freeze
        @test QCSHD + QRSHD == ∫M_col  # All collisions should result shedding
        @test ∫𝟙_wet_M_col == ∫M_col  # Above freezing, collisions at all sizes are wet
    end

    @testset "Bulk liquid-ice collisions" begin
        # Test the high-level interface with real P3 parameters
        state = P3.get_state(params; F_rim, ρ_rim, L_ice = Lᵢ, N_ice = Nᵢ)
        logλ = P3.get_distribution_logλ(state)

        # Create mock particle size distributions
        toml_dict = CP.create_toml_dict(FT)
        psd_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        psd_r = CMP.RainParticlePDF_SB2006_limited(toml_dict)

        # Test parameters
        L_c = FT(1e-3)  # 1 g/m³ cloud water
        N_c = FT(1e8)   # 100 million cloud droplets per m³
        L_r = FT(1e-4)  # 0.1 g/m³ rain water
        N_r = FT(1e6)   # 1 million raindrops per m³
        T = T_freeze - FT(5)  # 5K below freezing

        # Liquid particle mass function
        ρw = psd_c.ρw
        m_l(Dₗ) = ρw * CO.volume_sphere_D(Dₗ)

        # Test the high-level interface
        rates = P3.∫liquid_ice_collisions(
            state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T, m_l;
            quad = P3.ChebyshevGauss(50),
        )
        @test eltype(rates) == FT  # check type stability

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col = rates

        # Basic sanity checks
        @test all(rates .>= 0)
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col
        @test ∫𝟙_wet_M_col <= ∫M_col

        # Smoke tests, aka: Check that rates don't change with new commits.
        @test QCFRZ ≈ 5.89686152717295e-7
        @test QCSHD ≈ 2.075334534409237e-9
        @test NCCOL ≈ 60226.258f0
        @test QRFRZ ≈ 6.646808312782509e-5
        @test QRSHD ≈ 3.656428833353944e-6
        @test NRCOL ≈ 172.79896499770385
        @test ∫M_col ≈ 7.071627344843075e-5
        @test BCCOL ≈ 3.6970918665661123e-9
        @test BRCOL ≈ 4.1672779390485956e-7
        @test ∫𝟙_wet_M_col ≈ 1.561091379329206e-5

        ### Test the bulk source function
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        rates = P3.bulk_liquid_ice_collision_sources(
            state, logλ,
            psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T,
        )
        @test eltype(rates) == FT  # check type stability
    end
end

function test_p3_ice_self_collection(FT)
    params = CMP.ParametersP3(FT)
    vel_params = CMP.Chen2022VelType(FT)
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    (; T_freeze) = params

    ρₐ = FT(1.2)
    qᵢ = FT(1e-4)
    Lᵢ = qᵢ * ρₐ
    Nᵢ = FT(2e5) * ρₐ
    F_rim = FT(0.8)
    ρ_rim = FT(800)

    state = P3.get_state(params; F_rim, ρ_rim, L_ice = Lᵢ, N_ice = Nᵢ)
    logλ = P3.get_distribution_logλ(state)
    T = T_freeze - FT(5)  # 5K below freezing

    @testset "ice self-collection rate" begin
        # Call the new ice self-collection parameterization
        rates = P3.ice_self_collection(state, logλ, aps, tps, vel_params, ρₐ, T; quad = P3.ChebyshevGauss(50))
        @test eltype(rates) == FT  # check type stability

        # Self-collection should represent a positive loss rate
        @test rates.dNdt > 0

        # Test edge case with virtually zero L_ice and N_ice
        state_zero = P3.get_state(params; F_rim, ρ_rim, L_ice = FT(0), N_ice = FT(0))
        logλ_zero = P3.get_distribution_logλ(state_zero)
        rates_zero =
            P3.ice_self_collection(state_zero, logλ_zero, aps, tps, vel_params, ρₐ, T; quad = P3.ChebyshevGauss(50))
        @test rates_zero.dNdt == 0

        # TODO: compare against an analytically derived reference
        # For a simple size distribution and uniform velocity difference, one could compute analytical dNdt.
    end
end

@testset "P3 tests ($FT)" for FT in (Float64, Float32)
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

    # bulk liquid-ice collisions and related processes
    test_p3_bulk_liquid_ice_collisions(FT)
    test_p3_ice_self_collection(FT)
end
nothing
