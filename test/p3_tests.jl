using Test: @testset, @test, @test_throws, @test_broken, @inferred
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams as CP
import SpecialFunctions as SF
import QuadGK as QGK
import ForwardDiff as FD

function test_p3_state_creation(FT)
    @testset "P3State Creation and Properties" begin
        # Test creating a state with valid parameters
        params = CMP.ParametersP3(FT)
        L_ice = FT(0.22)
        N_ice = FT(1e6)
        F_rim = FT(0.5)
        ρ_rim = FT(400)

        # Test unrimed state
        state_unrimed = P3.P3State(params, L_ice, N_ice, FT(0), ρ_rim)
        @test P3.isunrimed(state_unrimed)

        # Test rimed state
        state_rimed = P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)
        @test !P3.isunrimed(state_rimed)

        # Test thresholds for unrimed state. Per the `P3State` constructor,
        # unrimed ice has no graupel → `D_gr = D_cr = Inf` (the "always before
        # graupel regime" sentinel) and `ρ_g = NaN` (should not be used).
        (; D_th, D_gr, D_cr) = state_unrimed
        @test isfinite(D_th)
        @test D_gr == Inf
        @test D_cr == Inf

        # Test thresholds for rimed state
        (; D_th, D_gr, D_cr) = state_rimed
        @test D_th < D_gr < D_cr
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
        D_th = P3.get_D_th(mass, ρ_i)
        for F_rim in F_rim_good
            for ρ_rim in ρ_rim_good
                ρ_d = P3.get_ρ_d(mass, F_rim, ρ_rim)
                ρ_g = P3.get_ρ_g(F_rim, ρ_rim, ρ_d)
                D_gr = P3.get_D_gr(mass, ρ_g)
                D_cr = P3.get_D_cr(mass, F_rim, ρ_g)
                @test D_th < D_gr < D_cr
                @test get_ρ_d_paper(mass; D_cr, D_gr) ≈ ρ_d
            end
        end

        # For very high rimed density, the thresholds are ill-defined. TODO: Investigate this
        F_rim_bad = FT(0.93)
        ρ_rim_bad = FT(975)
        ρ_g_bad = P3.get_ρ_g(mass, F_rim_bad, ρ_rim_bad)
        D_gr_bad = P3.get_D_gr(mass, ρ_g_bad)
        @test_broken D_th < D_gr_bad

        # Check that the P3 scheme solution matches the published values
        # D_cr and D_gr vs Fig. 1a Morrison and Milbrandt 2015
        D_cr_fig_1a_ref = FT[0.4946323381999426, 1.0170979628696817]
        D_gr_fig_1a_ref = FT[0.26151186272014415, 0.23392868352755775]
        for i in 1:2
            ρ_g = P3.get_ρ_g(mass, F_rim_good[i], ρ_rim_good[2])
            D_gr = P3.get_D_gr(mass, ρ_g)
            D_cr = P3.get_D_cr(mass, F_rim_good[i], ρ_g)
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
        ρ_g = P3.get_ρ_g(mass, F_rim, ρ_rim)
        D_th = P3.get_D_th(mass, ρ_i)
        D_gr = P3.get_D_gr(mass, ρ_g)
        D_cr = P3.get_D_cr(mass, F_rim, ρ_g)
        state = P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        spherical_area(D) = D^2 * π / 4
        nonspherical_area(D) = area.γ * D^area.σ
        @test P3.ice_area(state, D_1) == spherical_area(D_1)
        @test P3.ice_area(state, D_2) == nonspherical_area(D_2)
        @test P3.ice_area(state, D_3) == spherical_area(D_3)
        @test P3.ice_area(state, D_cr) == F_rim * spherical_area(D_cr) + (1 - F_rim) * nonspherical_area(D_cr)

        # test mass
        spherical_mass(ρ, D) = ρ * π / 6 * D^3
        nonspherical_mass(D) = mass.α_va * D^mass.β_va
        @test P3.ice_mass(state, D_1) == spherical_mass(ρ_i, D_1)
        @test P3.ice_mass(state, D_2) == nonspherical_mass(D_2)
        @test P3.ice_mass(state, D_3) == spherical_mass(ρ_g, D_3)
        @test P3.ice_mass(state, D_cr) == nonspherical_mass(D_cr) / (1 - F_rim)

        # test density
        @test P3.ice_density(state, D_1) ≈ ρ_i
        @test P3.ice_density(state, D_2) ≈ 544.916989830
        @test P3.ice_density(state, D_3) ≈ ρ_g
        @test P3.ice_density(state, D_cr) ≈ 383.33480937

        # test aspect ratio (oblate ϕ = 3√π m / (4 ρ A^{3/2}), ρ the per-regime
        # material density, see `P3.ϕᵢ`)
        aspect_ratio_closed(ρ, D) = 3 * sqrt(FT(π)) * P3.ice_mass(state, D) /
                                    (4 * ρ * P3.ice_area(state, D)^FT(1.5))
        @test P3.ϕᵢ(state, D_1) ≈ 1                              # D < D_th: spherical
        @test P3.ϕᵢ(state, D_2) ≈ aspect_ratio_closed(ρ_i, D_2)  # dense nonspherical
        @test P3.ϕᵢ(state, D_2) < 1                              # oblate
        @test P3.ϕᵢ(state, D_3) ≈ 1                              # graupel: spherical (ρ_g)
        @test P3.ϕᵢ(state, D_cr) ≈ aspect_ratio_closed(ρ_i, D_cr)  # partially rimed
        @test P3.ϕᵢ(state, D_cr) < 1                             # oblate
        # residual ϕ > 1 band just above D_th (area discontinuity, see `P3.ϕᵢ`)
        @test 1 < P3.ϕᵢ(state, D_th * FT(1.001)) < FT(1.3)

        # test F_rim = 0 and D > D_th
        state′ = P3.P3State(params, L_ice, N_ice, FT(0), ρ_rim)
        @test P3.ice_area(state′, D_2) == nonspherical_area(D_2)
        @test P3.ice_mass(state′, D_2) == nonspherical_mass(D_2)

        # TODO: Add tests for F_liq != 0
    end
end

function test_shape_solver(FT)

    slope_laws = (:constant, :powerlaw)
    for slope_law in slope_laws
        params = CMP.ParametersP3(FT; slope_law)

        @testset "Shape parameters - nonlinear solver" begin
            # -- First, test limiting behavior: `N_ice = L_ice = 0` --
            state = P3.P3State(params, FT(0), FT(0), FT(0.5), FT(500))
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

                            state = P3.P3State(params, FT(0), FT(0), F_rim, ρ_rim) # L_ice, N_ice not used in this test
                            # Compute the shape parameters that correspond to the input test values
                            logλ_ex = log(λ_ex)
                            μ = P3.get_μ(params.slope, logλ_ex)
                            logN₀_ex = P3.get_logN₀(N_ice, μ, logλ_ex)
                            # Compute mass density based on input shape parameters
                            L_calc = exp(log(N_ice) + P3.logLdivN(state, logλ_ex))

                            if L_calc < FT(1)
                                # Solve for shape parameters
                                state′ = P3.P3State(params, L_calc, N_ice, F_rim, ρ_rim)
                                logλ = P3.get_distribution_logλ(state′)
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
                P3.P3State(params, FT(2.366e-5), FT(16461.6), FT(0.2), FT(800)),
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
                                P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim),
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
        params_noar = CMP.ParametersP3(FT; aspect_ratio = CMP.NoAspectRatio())
        state = P3.P3State(params_noar, FT(0), FT(0), F_rim, ρ_rim)
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.7912, 1.1550, 1.4871]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test vel ≈ expected[i] rtol = 1e-3
        end

        state = P3.P3State(params, FT(0), FT(0), F_rim, ρ_rim)  # aspect_ratio = Oblate() default
        # one D per P3 regime; `cbrt(ϕ) ≤ 1` slows the nonspherical sizes
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.38381, 0.79121, 1.155, 1.1477]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state)
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
        state = P3.P3State(params, FT(0), FT(0), F_rim, ρ_rim)  # aspect_ratio = Oblate() default
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13192, 0.50457, 0.90753, 1.3015, 1.6757]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test_broken vel ≈ expected[i] rtol = 1e-3  # TODO: Implement `F_liq != 0`
        end
        state = P3.P3State(CMP.ParametersP3(FT; aspect_ratio = CMP.NoAspectRatio()), FT(0), FT(0), F_rim, ρ_rim)
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13191, 0.50457, 0.90753, 1.301499, 1.67569]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state)
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

        state₀ = P3.P3State(params, FT(0), N_ice, FT(0.5), ρ_rim)
        logλ = P3.get_distribution_logλ(state₀)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state₀, logλ; quad = P3.GaussLegendre(FT, 12))
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state₀, logλ; quad = P3.GaussLegendre(FT, 12))
        @test iszero(vel_n₀)
        @test iszero(vel_m₀)

        state₀ = P3.P3State(params, L_ice, FT(0), FT(0.5), ρ_rim)
        logλ = P3.get_distribution_logλ(state₀)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state₀, logλ; quad = P3.GaussLegendre(FT, 12))
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state₀, logλ; quad = P3.GaussLegendre(FT, 12))
        @test iszero(vel_n₀)
        @test iszero(vel_m₀)

        # NOTE: All reference values are output from the code.
        # A failing test indicates that the code has changed.
        # But if the changes are intentional, the reference values can be updated.

        # Liquid fraction = 0. The `_ϕ` (aspect-ratio-on) references are below
        # their aspect-off counterparts (`cbrt(ϕ) < 1`).
        ref_v_n = [3.646059575504377, 2.6191026241691695]
        ref_v_n_ϕ = [1.5223915218714987, 1.4656564581919258]
        ref_v_m = [7.788114224053879, 5.797675366222473]
        ref_v_m_ϕ = [2.427666066669716, 2.3683439025452544]

        params_noar = CMP.ParametersP3(FT; aspect_ratio = CMP.NoAspectRatio())
        for (k, F_rim) in enumerate(F_rims)
            state = P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)
            state_noar = P3.P3State(params_noar, L_ice, N_ice, F_rim, ρ_rim)
            logλ = P3.get_distribution_logλ(state)
            quad = P3.GaussLegendre(FT, 12)
            vel_n = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state_noar, logλ; quad)
            vel_m = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state_noar, logλ; quad)
            vel_n_ϕ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state, logλ; quad)
            vel_m_ϕ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state, logλ; quad)

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
            state = P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)
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
    params = CMP.ParametersP3(FT; aspect_ratio = CMP.NoAspectRatio())
    Chen2022 = CMP.Chen2022VelType(FT)

    N_ice = FT(1e8)
    L_ices = range(FT(0.001), stop = FT(0.005), length = 5)
    ρ_rim = FT(500)
    F_rims = FT[0, 0.5]
    ρ_a = FT(1.2)
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

    @testset "Gauss-Legendre quadrature" begin
        quad = P3.GaussLegendre(16)
        # exact for polynomials up to degree 2n-1 (here deg 4 ≤ 31)
        @test P3.integrate(x -> x^4, 0, 1, quad) ≈ 0.2 rtol = 1e-12
        @test P3.integrate(x -> x^7, 0.0, 1.0, P3.GaussLegendre(16)) ≈ 0.125 rtol = 1e-12
        # higher order remains spectrally accurate on a non-polynomial integrand
        ref = exp(1) - 1                          # ∫₀¹ eˣ dx
        e_lo = abs(P3.integrate(exp, 0.0, 1.0, P3.GaussLegendre(16)) - ref)
        e_hi = abs(P3.integrate(exp, 0.0, 1.0, P3.GaussLegendre(40)) - ref)
        @test e_lo < 1e-12 && e_hi < 1e-12
        # nested NTuple form (used by the P3 collision integrals)
        @test P3.integrate(x -> x^2, (0.0, 1.0, 2.0), P3.GaussLegendre(32)) ≈ 8 / 3 rtol = 1e-12
        # GPU-safety invariants: the constructed rule is isbits / concrete, so
        # it ships to GPU kernels with no per-call construction.
        @test isbits(P3.GaussLegendre(40))
        @test isconcretetype(typeof(P3.GaussLegendre(40)))
        @test eltype(P3.GaussLegendre(Float32, 32).nodes) == Float32
        @test eltype(P3.GaussLegendre(Float64, 32).nodes) == Float64
        # Arbitrary orders are supported now (no baked tables) — including
        # orders the old table-based design rejected, e.g. 37.
        @test P3.integrate(x -> x^4, 0.0, 1.0, P3.GaussLegendre(37)) ≈ 0.2 rtol = 1e-12
        @test sum(P3.GaussLegendre(37).weights) ≈ 2 rtol = 1e-12
        # The reused-once design: reading nodes/weights in the hot loop is a
        # static SVector lookup, type-stable for a concretely-typed rule.
        let q = P3.GaussLegendre(40)
            @test (@inferred P3.node(q, 1.0, q.n)) isa Float64
            @test (@inferred P3.weight(q, 1.0, q.n)) isa Float64
        end
    end

    @testset "Numerical integrals sanity checks for N, velocity and diameter" begin
        for (F_rim, L_ice, p) in Iterators.product(F_rims, L_ices, ps)

            # Get shape parameters, thresholds and intergal bounds
            state = P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)
            logλ = P3.get_distribution_logλ(state)

            # Number concentration comparison
            N′ = P3.size_distribution(state, logλ)
            bnds = P3.integral_bounds(state, logλ; p = 1e-6, moment_order = 0)
            N_estim_gl = P3.integrate(N′, bnds, P3.GaussLegendre(FT, 32))
            N_tol = FT == Float32 ? 2e-5 : 1e-5  # native-FT gamma_inc slightly less precise than Float64-backed SF
            @test N_ice ≈ N_estim_gl rtol = N_tol

            # Compare with quadgk
            N_estim_qgk = QGK.quadgk(N′, bnds...)[1]
            @test N_estim_gl ≈ N_estim_qgk rtol = 1e-5


            # Bulk velocity comparison
            vel_N = P3.ice_terminal_velocity_number_weighted(
                Chen2022, ρ_a, state, logλ;
                p, quad = P3.GaussLegendre(FT, 12),
            )
            vel_m = P3.ice_terminal_velocity_mass_weighted(
                Chen2022, ρ_a, state, logλ;
                p, quad = P3.GaussLegendre(FT, 12),
            )

            v_term = P3.ice_particle_terminal_velocity(Chen2022, ρ_a, state)
            g(D) = v_term(D) * N′(D)
            gm(D) = g(D) * P3.ice_mass(state, D)
            vel_N_estim_gl = P3.integrate(g, bnds, P3.GaussLegendre(FT, 32)) / N_ice
            vel_m_estim_gl = P3.integrate(gm, bnds, P3.GaussLegendre(FT, 32)) / L_ice
            @test vel_N ≈ vel_N_estim_gl rtol = 0.005
            @test vel_m ≈ vel_m_estim_gl rtol = 0.05

            # Compare with quadgk
            vel_N_estim_qgk = QGK.quadgk(g, bnds...)[1] / N_ice
            vel_m_estim_qgk = QGK.quadgk(gm, bnds...)[1] / L_ice

            @test vel_N_estim_gl ≈ vel_N_estim_qgk rtol = 0.005
            @test vel_m_estim_gl ≈ vel_m_estim_qgk rtol = 0.05


            # Dₘ comparisons
            D_m = P3.D_m(state, logλ)
            D_m_func(D) = D * P3.ice_mass(state, D) * N′(D) / L_ice
            D_m_estim_gl = P3.integrate(D_m_func, bnds, P3.GaussLegendre(FT, 32))
            @test D_m ≈ D_m_estim_gl rtol = 5e-4

            # Compare with quadgk
            D_m_estim_qgk = QGK.quadgk(D_m_func, bnds...)[1]
            @test D_m_estim_gl ≈ D_m_estim_qgk rtol = 5e-4
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

        # Reference values are output from the code. These are now the *uncapped*
        # instantaneous rates (the availability/dt cap was removed). The `qᵥ` sweep
        # is held below ~RH 1.16 so the raw ABIFM rate stays finite and consistent
        # across Float32/Float64; higher RH overflows `J` in Float32 (→ 0 via the
        # `isfinite` guard) while Float64 explodes to ~1e69. Update if the rate changes.
        expected_freeze_N = [
            1.0473022910416842e-10, 5.925723559806242e-6, 0.33501487392087853,
            18925.187757721098, 1.0682422440661902e9, 6.0249407658238766e13,
        ]
        expected_freeze_L = [
            1.4953923796668527e-22, 8.460745965684499e-18, 4.783166694522096e-13,
            2.701940516414268e-8, 0.0015250690076232267, 86.01153323961839,
        ]
        qᵥ_range = range(FT(0.5e-3), stop = FT(0.8e-3), length = 6)

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

        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        logλ = P3.get_distribution_logλ(state)

        T_cold = FT(273.15 - 0.01)

        rate = P3.ice_melt(vel, aps, tps, T_cold, ρₐ, state, logλ; quad = P3.GaussLegendre(FT, 12))

        @test rate.dNdt == 0
        @test rate.dLdt == 0

        T_warm = FT(273.15 + 0.01)
        rate = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state, logλ; quad = P3.GaussLegendre(FT, 12))

        @test rate.dNdt >= 0
        @test rate.dLdt >= 0

        # NOTE: All reference values are output from the code.
        # A failing test indicates that the code has changed.
        # But if the changes are intentional, the reference values can be updated.
        if FT == Float64
            ref_dNdt = FT(172084.75278912345)
            ref_dLdt = FT(8.604237639456172e-5)
        else
            ref_dNdt = FT(172265.67f0)
            ref_dLdt = FT(8.613284f-5)
        end
        @test rate.dNdt ≈ ref_dNdt
        @test rate.dLdt ≈ ref_dLdt

        T_vwarm = FT(273.15 + 0.1)
        rate = P3.ice_melt(vel, aps, tps, T_vwarm, ρₐ, state, logλ; quad = P3.GaussLegendre(FT, 12))
        if FT == Float64
            ref_vwarm_dNdt = FT(1.7198680382990765e6)
            ref_vwarm_dLdt = FT(8.599340191495382e-4)
        else
            ref_vwarm_dNdt = FT(1.7201018f6)
            ref_vwarm_dLdt = FT(8.6005084f-4)
        end
        @test rate.dNdt ≈ ref_vwarm_dNdt
        @test rate.dLdt ≈ ref_vwarm_dLdt
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

    state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
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
        cloud_integrals = P3.get_liquid_integrals(n_c, ∂ₜV, m_l, ρ′_rim, liq_bounds; quad = P3.GaussLegendre(FT, 12))
        rain_integrals = P3.get_liquid_integrals(n_r, ∂ₜV, m_l, ρ′_rim, liq_bounds; quad = P3.GaussLegendre(FT, 12))

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
        integrals_n0 = P3.get_liquid_integrals(n_zero, ∂ₜV, m_l, ρ′_rim, liq_bounds; quad = P3.GaussLegendre(FT, 12))
        result = integrals_n0(Dᵢ)
        @test all(iszero, result)

        # Mass rate should be related to number rate through mass
        # ∂ₜM_col should be approximately ∫ n(D) * m_l(D) * ∂ₜV dD
        # This is a simplified check - in reality it's more complex

        # Test the full ∫liquid_ice_collisions function
        # Mock functions
        ∂ₜM_max(Dᵢ) = FT(0.04)  # Small max freeze rate
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad = P3.GaussLegendre(FT, 12),
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
            n_i_zero, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad = P3.GaussLegendre(FT, 12),
        )
        @test all(iszero, rates)

        # Zero liquid content
        n_zero = Returns(FT(0))
        zero_liq_integrals =
            P3.get_liquid_integrals(n_zero, ∂ₜV, m_l, ρ′_rim, liq_bounds; quad = P3.GaussLegendre(FT, 12))
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max, zero_liq_integrals, zero_liq_integrals, ice_bounds; quad = P3.GaussLegendre(FT, 12),
        )
        @test all(iszero, rates)

        # No freezing (above freezing temperature)
        ∂ₜM_max_zero(Dᵢ) = FT(0)
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max_zero, cloud_integrals, rain_integrals, ice_bounds; quad = P3.GaussLegendre(FT, 12),
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
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
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
            quad = P3.GaussLegendre(FT, 12),
        )
        @test eltype(rates) == FT  # check type stability

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col = rates

        # Basic sanity checks
        @test all(rates .>= 0)
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col
        @test ∫𝟙_wet_M_col <= ∫M_col

        # Smoke tests, aka: Check that rates don't change with new commits.
        # `rtol = 5e-4` admits both Float32 and Float64 against these (Float64)
        # reference values.
        @test QCFRZ ≈ 5.942471550989089e-7 rtol = 5e-4
        @test QCSHD ≈ 2.0728862241368704e-9 rtol = 5e-4
        @test NCCOL ≈ 60651.35670910096 rtol = 5e-4
        @test QRFRZ ≈ 6.642674674038379e-5 rtol = 5e-4
        @test QRSHD ≈ 3.6526001759370415e-6 rtol = 5e-4
        @test NRCOL ≈ 172.61819652435105 rtol = 5e-4
        @test ∫M_col ≈ 7.067566695764388e-5 rtol = 5e-4
        @test BCCOL ≈ 3.725687492783128e-9 rtol = 5e-4
        @test BRCOL ≈ 4.164686317018988e-7 rtol = 5e-4
        @test ∫𝟙_wet_M_col ≈ 1.5520362321253953e-5 rtol = 5e-4

        ### Test the bulk source function
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        rates = P3.bulk_liquid_ice_collision_sources(
            state, logλ,
            psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T;
            quad = P3.GaussLegendre(FT, 12),
        )
        @test eltype(rates) == FT  # check type stability
    end
end

function test_p3_ice_self_collection(FT)
    params = CMP.ParametersP3(FT)
    vel_params = CMP.Chen2022VelType(FT)

    ρₐ = FT(1.2)
    qᵢ = FT(1e-4)
    Lᵢ = qᵢ * ρₐ
    Nᵢ = FT(2e5) * ρₐ
    F_rim = FT(0.8)
    ρ_rim = FT(800)

    state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
    logλ = P3.get_distribution_logλ(state)

    @testset "ice self-collection rate" begin
        # Call the new ice self-collection parameterization
        rates = P3.ice_self_collection(state, logλ, vel_params, ρₐ; quad = P3.GaussLegendre(FT, 12))
        @test eltype(rates) == FT  # check type stability

        # Self-collection should represent a positive loss rate
        @test rates.dNdt > 0

        # Test edge case with virtually zero L_ice and N_ice
        state_zero = P3.P3State(params, FT(0), FT(0), F_rim, ρ_rim)
        logλ_zero = P3.get_distribution_logλ(state_zero)
        rates_zero =
            P3.ice_self_collection(state_zero, logλ_zero, vel_params, ρₐ; quad = P3.GaussLegendre(FT, 12))
        @test rates_zero.dNdt == 0

        # Cross-check the triangular domain against the full-square double
        # integral, where the ½ factor counts each unordered pair once
        quad32 = P3.GaussLegendre(FT, 32)
        rates32 = P3.ice_self_collection(state, logλ, vel_params, ρₐ; quad = quad32)
        n_i = DT.size_distribution(state, logλ)
        v_i = P3.ice_particle_terminal_velocity(vel_params, ρₐ, state)
        bnds = P3.velocity_integral_bounds(state, logλ, vel_params.small_ice.cutoff; p = eps(one(ρₐ)))
        square = P3.integrate(
            D₁ -> begin
                v₁ = v_i(D₁)
                collision_rate =
                    D₂ -> P3.collision_cross_section_ice_ice(state, D₁, D₂) * abs(v₁ - v_i(D₂)) * n_i(D₂)
                # Split the inner integral at D₂ = D₁, where |v₁ - v(D₂)| is not smooth
                inner =
                    P3.integrate(collision_rate, (first(bnds), D₁), quad32) +
                    P3.integrate(collision_rate, (D₁, last(bnds)), quad32)
                n_i(D₁) * inner
            end,
            bnds, quad32,
        )
        @test rates32.dNdt ≈ square / 2 rtol = 0.05
    end
end

function test_p3_closed_form_rain_inner(FT)
    # The closed form is the exact analytic reduction of the rain inner integral.
    # Testset checks:
    #  - N and M matches an adaptive (QuadGK) reference numerical quadrature
    #  - correctness of the rime volume quadrature
    #  - smoke test for collision cross section
    @testset "P3 closed-form rain inner (N, M, B) + cross-section coeffs" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        psd_r = CMP.SB2006(FT).pdf_r
        ρₐ = FT(1)
        ρ_w = psd_r.ρw
        p = eltype(params)(1e-5)
        m_liq(D) = ρ_w * FT(π) / 6 * D^3
        rtol = FT == Float64 ? FT(1e-10) : FT(1e-3)
        qrtol = FT == Float64 ? FT(1e-12) : FT(1e-7)
        v_l = CO.particle_terminal_velocity(vel.rain, ρₐ)
        for (L_ice, N_ice, F_rim, ρ_rim) in (
                (1e-3, 1e6, 0.5, 500),
                (1e-2, 1e8, 0.95, 800),
                (1e-5, 1e4, 0.0, 200),
            ),
            (L_r, N_r) in ((1e-6, 1e4), (1e-4, 1e3), (2e-3, 5e2))

            state = P3.P3State(params, FT(L_ice), FT(N_ice), FT(F_rim), FT(ρ_rim))
            n_r = DT.size_distribution(psd_r, FT(L_r) / ρₐ, ρₐ, FT(N_r))
            ∂ₜV = P3.volumetric_collision_rate_integrand(vel, ρₐ, state)
            ρ′_rim = P3.compute_local_rime_density(vel, ρₐ, FT(270), state)
            D_min, D_max = bnds = CM2.get_size_distribution_bounds(psd_r, FT(L_r) / ρₐ, ρₐ, FT(N_r), p)
            D_max > D_min || continue
            v_i = P3.ice_particle_terminal_velocity(vel, ρₐ, state)
            rc = P3.get_liquid_integrals_rain_closed(
                psd_r, vel, n_r, ρₐ, FT(L_r), FT(N_r), state, ∂ₜV,
                m_liq, ρ′_rim, bnds; quad = P3.ChebyshevGauss(40), v_i, v_l,
            )
            rn = P3.get_liquid_integrals(  # numerical fallback
                n_r, ∂ₜV, m_liq, ρ′_rim, bnds;
                quad = P3.ChebyshevGauss(40), v_i, v_l,
            )
            for Dᵢ in FT.(10 .^ range(-5, -2; length = 5))
                vi = v_i(Dᵢ)
                Dstar = P3.crossover_diameter(vi, v_l, D_min, D_max)

                # N, M: closed form vs adaptive reference
                Nc, Mc, Bc = rc(Dᵢ)
                σ(D) = P3.collision_cross_section_ice_liquid(state, Dᵢ, D)
                gN(D) = σ(D) * abs(vi - v_l(D)) * n_r(D)
                Nref = QGK.quadgk(gN, D_min, Dstar, D_max; rtol = qrtol)[1]
                Mref = QGK.quadgk(D -> gN(D) * m_liq(D), D_min, Dstar, D_max; rtol = qrtol)[1]
                @test isapprox(Nc, Nref; rtol)
                @test isapprox(Mc, Mref; rtol)

                # B: closed quadrature vs numerical fallback at the same order
                @test isapprox(Bc, rn(Dᵢ)[3]; rtol = sqrt(eps(FT)))

                # smoke test: collision cross section has form: `π(rᵢ + Dₗ/2)²`, with rᵢ derived from `P3.ice_area`
                rᵢ = sqrt(P3.ice_area(state, Dᵢ) / FT(π))
                K = P3.collision_cross_section_ice_liquid_coeffs(state, Dᵢ)
                @test K == P3.collision_cross_section_ice_liquid_coeffs(rᵢ)
                @test evalpoly(Dstar, K) ≈ FT(π) * (rᵢ + Dstar / 2)^2
                @test evalpoly(D_max, K) ≈ FT(π) * (rᵢ + D_max / 2)^2
            end
        end
    end

    # AD smoke test for the closed form.
    # find D where ice/liquid sedimentation velocities are equal,
    # separate the integrals, then check that closed form is correct.
    @testset "P3 closed-form ForwardDiff AD smoke (v_i,r_i,Dr,N₀r)" begin
        vel = CMP.Chen2022VelType(FT)
        psd_r = CMP.SB2006(FT).pdf_r
        ρₐ = FT(1)
        ρ_w = psd_r.ρw
        ai, bi, ci = CO.Chen2022_vel_coeffs(vel.rain, ρₐ)
        v_l = CO.particle_terminal_velocity(vel.rain, ρₐ)
        L_r, N_r = FT(1e-4), FT(1e3)
        (; N₀r, Dr_mean) =
            CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
        rᵢ0 = FT(1e-3)
        vi0 = FT(3.0)
        D_min, D_max = FT(1e-5), FT(5e-3)
        rtolAD = FT(1e-4)
        Dstar0 = P3.crossover_diameter(vi0, v_l, D_min, D_max)
        cases = (
            ("v_i", vi0,
                x -> P3.closed_rain_inner_NM(x, Dstar0, rᵢ0, ρ_w,
                    ai, bi, ci, D_min, D_max, N₀r, Dr_mean)),
            ("r_i", rᵢ0,
                x -> P3.closed_rain_inner_NM(vi0, Dstar0, x, ρ_w,
                    ai, bi, ci, D_min, D_max, N₀r, Dr_mean)),
            ("Dr", Dr_mean,
                x -> P3.closed_rain_inner_NM(vi0, Dstar0, rᵢ0, ρ_w,
                    ai, bi, ci, D_min, D_max, N₀r, x)),
            ("N₀r", N₀r,
                x -> P3.closed_rain_inner_NM(vi0, Dstar0, rᵢ0, ρ_w,
                    ai, bi, ci, D_min, D_max, x, Dr_mean)),
        )
        for (_, x0, g) in cases, idx in (1, 2)
            f(x) = g(x)[idx]
            d_ad = FD.derivative(f, x0)
            @test isfinite(d_ad)  # check that AD works
            if FT == Float64
                h = max(abs(x0) * FT(1e-5), FT(1e-12))
                d_fd = (f(x0 + h) - f(x0 - h)) / (2h)
                @test isapprox(d_ad, d_fd;
                    rtol = rtolAD,
                    atol = 100 * eps(FT) * max(abs(d_ad), abs(d_fd), one(FT)),
                )
            end
        end
        # Check that the crossover diameter is differentiable
        v_tgt = (v_l(D_min) + v_l(D_max)) / 2
        dDstar_bisect = FD.derivative(
            vt -> P3.crossover_diameter(vt, v_l, D_min, D_max), v_tgt,
        )
        @test isfinite(dDstar_bisect)            # frozen-root: == 0
        Dstar = P3.crossover_diameter(v_tgt, v_l, D_min, D_max)
        vlp = FD.derivative(v_l, Dstar)          # v_l′ elementary, AD-clean
        @test isfinite(vlp) && vlp > 0
        # Check that it matches a numerical derivative
        if FT == Float64
            hv = abs(v_tgt) * FT(1e-6)
            dDstar_fd =
                (
                    P3.crossover_diameter(v_tgt + hv, v_l, D_min, D_max) -
                    P3.crossover_diameter(v_tgt - hv, v_l, D_min, D_max)
                ) /
                (2hv)
            @test isapprox(dDstar_fd, inv(vlp); rtol = FT(1e-3))
        end
    end

    # Check edge cases (e.g. ice velocity > liquid velocity for all sizes)
    # still returns the correct result
    @testset "P3 closed-form D*-at-bracket-end (v_i outside band)" begin
        vel = CMP.Chen2022VelType(FT)
        v_l = CO.particle_terminal_velocity(vel.rain, FT(1))
        D_min, D_max = FT(1e-5), FT(5e-3)
        # v_i far below v_l(D_min) ⇒ D* = D_min ; far above ⇒ D* = D_max.
        @test P3.crossover_diameter(FT(-1), v_l, D_min, D_max) == D_min
        @test P3.crossover_diameter(FT(1e6), v_l, D_min, D_max) == D_max
        # An interior target lands strictly inside the bracket.
        v_mid = (v_l(D_min) + v_l(D_max)) / 2
        Dstar = P3.crossover_diameter(v_mid, v_l, D_min, D_max)
        @test D_min <= Dstar <= D_max
        @test isapprox(v_l(Dstar), v_mid; rtol = FT(1e-4))
        # Full closed path with v_i outside the band stays finite and the
        # absent-crossover side collapses (one piece zero).
        r_i = FT(2e-4)
        ai, bi, ci = CO.Chen2022_vel_coeffs(vel.rain, FT(1))
        for v_i in (FT(-1), FT(1e6))
            Dstar_i = P3.crossover_diameter(v_i, v_l, D_min, D_max)
            N, M = P3.closed_rain_inner_NM(
                v_i, Dstar_i, r_i, FT(1000), ai, bi, ci,
                D_min, D_max, FT(1e7), FT(5e-4),
            )
            @test isfinite(N) && isfinite(M)
        end
    end

    @testset "P3 closed-form dispatch fallback (non-Chen / non-SB2006)" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        psd_r = CMP.SB2006(FT).pdf_r
        state = P3.P3State(params, FT(1e-3), FT(1e6), FT(0.5), FT(500))
        rest = (
            identity, (a, b) -> a, identity, identity, identity,
            FT(1), FT(1e-4), FT(1e3), state,
        )
        m_closed = which(P3._rain_inner_integrals, typeof((psd_r, vel, rest...)))
        m_fallback = which(P3._rain_inner_integrals, typeof((1.0, 2.0, rest...)))
        # closed-form method: first two params are the typed bundle
        @test m_closed.sig.parameters[2] <: CMP.RainParticlePDF_SB2006
        @test m_closed.sig.parameters[3] <: CMP.Chen2022VelType
        # fallback method: first two params are ::Any (Any === Any)
        @test m_fallback.sig.parameters[2] === Any
        @test m_fallback.sig.parameters[3] === Any
        # and the two methods are distinct (no accidental ambiguity merge)
        @test m_closed !== m_fallback
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
    test_p3_closed_form_rain_inner(FT)
end
nothing
