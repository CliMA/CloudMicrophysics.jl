using Test: @testset, @test, @test_throws, @test_broken, @inferred
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Common as CO
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Utilities as UT
import ClimaParams as CP
import SpecialFunctions as SF
import QuadGK as QGK
import ForwardDiff as FD

# P3 state from the rime mass fraction and the rime density
rimed_state(params, ρq_ice, ρn_ice, F_rim, ρ_rim; quad = P3.GaussLegendre(typeof(params.ρ_i), 12)) =
    P3.P3State(params, quad, ρq_ice, ρn_ice, F_rim * ρq_ice, F_rim * ρq_ice / ρ_rim)

# Number- and mass-weighted means of `f(D)` over the size distribution of `state` on `bnds`, by
# adaptive quadrature
function adaptive_weighted_means(f, state, bnds)
    n = P3.size_distribution(state.dist)
    m(D) = P3.ice_mass(state.particles, D)
    ∫(g) = QGK.quadgk(g, bnds...)[1]
    return (∫(D -> n(D) * f(D)) / ∫(n), ∫(D -> n(D) * m(D) * f(D)) / ∫(D -> n(D) * m(D)))
end

function test_p3_particles(FT)
    @testset "P3Particles from the rime mass fraction and the rime density" begin
        params = CMP.ParametersP3(FT)
        ρ_rim = FT(400)
        # Unrimed ice has no graupel regime
        unrimed = P3.P3Particles(params, FT(0), ρ_rim)
        @test P3.isunrimed(unrimed)
        @test isfinite(unrimed.D_th) && unrimed.D_gr == unrimed.D_cr == Inf
        rimed = P3.P3Particles(params, FT(0.5), ρ_rim)
        @test !P3.isunrimed(rimed) && isbits(rimed)
        @test rimed.D_th < rimed.D_gr < rimed.D_cr
        # The state builds the same particles from its bulk rime moments
        particles = rimed_state(params, FT(0.22), FT(1e6), FT(0.5), ρ_rim).particles
        @test all(f -> getfield(particles, f) ≈ getfield(rimed, f), (:F_rim, :ρ_rim, :ρ_g, :D_th, :D_gr, :D_cr))
    end
end

function test_state_construction(FT)
    @testset "P3State projects the rime pair onto the admissible densities" begin
        params = CMP.ParametersP3(FT)
        ρ_min, ρ_max = params.ρ_rim_local(one(FT)), params.ρ_i
        ρq_ice, ρn_ice = FT(1e-3), FT(1e5)
        # (ρq_rim, ρb_rim): admissible, too dense, too soft, mass without volume, negative
        # volume, volume without mass, negative mass, no rime
        rime_pairs = (
            (FT(3e-4), FT(1e-6)), (FT(3e-4), FT(1e-8)), (FT(3e-4), FT(1e-5)), (FT(3e-4), FT(0)),
            (FT(3e-4), FT(-1e-6)), (FT(0), FT(1e-6)), (FT(-3e-4), FT(1e-6)), (FT(0), FT(0)),
        )
        for (ρq_rim, ρb_rim) in rime_pairs
            state = @inferred P3.P3State(params, P3.GaussLegendre(FT, 12), ρq_ice, ρn_ice, ρq_rim, ρb_rim)
            if ρq_rim > 0
                # The rime mass is kept, and the volume moves to the nearest admissible value
                @test state.particles.F_rim ≈ ρq_rim / ρq_ice
                @test state.particles.ρ_rim ≈ (ρb_rim > 0 ? clamp(ρq_rim / ρb_rim, ρ_min, ρ_max) : ρ_max)
            else
                @test state.particles.F_rim == 0 && state.particles.ρ_rim == 0
            end
            # The state holds the prognostic rime volume, clamped to be nonnegative
            @test state.ρb_rim == max(ρb_rim, 0)
            f_q(q) = P3.P3State(params, P3.GaussLegendre(FT, 12), ρq_ice, ρn_ice, q, ρb_rim).particles.ρ_rim
            f_b(b) = P3.P3State(params, P3.GaussLegendre(FT, 12), ρq_ice, ρn_ice, ρq_rim, b).particles.ρ_rim
            @test isfinite(FD.derivative(f_q, ρq_rim))
            @test isfinite(FD.derivative(f_b, ρb_rim))
        end

        # Negative ice mass and number are clamped to zero
        state = P3.P3State(params, P3.GaussLegendre(FT, 12), -ρq_ice, -ρn_ice, FT(0), FT(0))
        @test state.ρq_ice == 0 && state.ρn_ice == 0
    end

    @testset "the distribution bounds the mean ice mass and solves a consistent logλ" begin
        params = CMP.ParametersP3(FT)
        quad = P3.GaussLegendre(FT, 64)
        x_min = P3.ice_mean_particle_mass_min(params)
        x_max = P3.ice_mean_particle_mass_max(params)
        ρq_ice = FT(1e-4)
        # (ρn_ice, bounded ρn_ice): too few particles, too many particles, within bounds
        cases = ((FT(1e-2), ρq_ice / x_max), (FT(1e14), ρq_ice / x_min), (FT(1e5), FT(1e5)))
        for (ρn_ice, ρn_bounded) in cases, (ρq_rim, ρb_rim) in ((FT(0), FT(0)), (FT(5e-5), FT(1e-7)))
            st = P3.P3State(params, quad, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
            dist = @inferred P3.P3Distribution(st.particles, quad, st.ρq_ice, st.ρn_ice)
            @test st.ρn_ice == ρn_ice
            @test dist.ρn ≈ ρn_bounded rtol = 4 * eps(FT)
            # The size distribution integrates back to the bounded moments
            N′ = P3.size_distribution(dist)
            N = P3.integrate(N′, P3.integral_bounds(st.particles, st.dist; p = 1e-6), quad)
            bnds_L = P3.integral_bounds(st.particles, st.dist; p = 1e-6, moment_order = 2)
            L = P3.integrate(D -> N′(D) * P3.ice_mass(st.particles, D), bnds_L, quad)
            rtol = FT === Float64 ? 1e-5 : 1e-3
            @test N ≈ dist.ρn rtol = rtol
            @test L ≈ ρq_ice rtol = rtol
        end
        # Without ice mass or without ice number, the number is zero
        @test P3.P3State(params, quad, FT(0), FT(1e5), FT(0), FT(0)).dist.ρn == 0
        @test P3.P3State(params, quad, ρq_ice, FT(0), FT(0), FT(0)).dist.ρn == 0
        state = P3.P3State(params, quad, ρq_ice, FT(1e5), FT(0), FT(0))
        @test isbits(state) && isbits(state.particles) && isbits(state.dist)
        (; particles, ρq_ice, ρn_ice) = state
        @test P3.P3Distribution(particles, quad, ρq_ice, ρn_ice) === state.dist
    end

    @testset "admissibility of the ice moments and the rime density bounds" begin
        params = CMP.ParametersP3(FT)
        x_min = P3.ice_mean_particle_mass_min(params)
        x_max = P3.ice_mean_particle_mass_max(params)
        @test x_max == params.m_max
        ρq_ice = FT(1e-4)
        @test P3.ice_moments_are_admissible(params, ρq_ice, FT(1e5))
        @test P3.ice_moments_are_admissible(params, ρq_ice, ρq_ice / x_max)
        @test P3.ice_moments_are_admissible(params, ρq_ice, ρq_ice / x_min)
        @test !P3.ice_moments_are_admissible(params, ρq_ice, ρq_ice / x_max / 2)
        @test !P3.ice_moments_are_admissible(params, ρq_ice, 2 * ρq_ice / x_min)
        @test !P3.ice_moments_are_admissible(params, FT(0), FT(1e5))
        @test !P3.ice_moments_are_admissible(params, ρq_ice, FT(0))
        # Admissible moments are the moments the distribution leaves unchanged
        for ρn_ice in (FT(1e-2), FT(1e5), FT(1e14))
            dist = P3.P3State(params, P3.GaussLegendre(FT, 12), ρq_ice, ρn_ice, FT(0), FT(0)).dist
            @test P3.ice_moments_are_admissible(params, ρq_ice, ρn_ice) == (dist.ρn == ρn_ice)
        end
        (ρ_rim_min, ρ_rim_max) = P3.rime_density_bounds(params)
        @test ρ_rim_min == params.ρ_rim_local(1)
        @test ρ_rim_max == params.ρ_i
        # State construction moves the rime density into the bounds
        for ρb_rim in (FT(1e-9), FT(1e-6), FT(1e-3))
            state = P3.P3State(params, P3.GaussLegendre(FT, 12), ρq_ice, FT(1e5), FT(5e-5), ρb_rim)
            @test ρ_rim_min * (1 - 4eps(FT)) <= state.particles.ρ_rim <= ρ_rim_max * (1 + 4eps(FT))
        end
    end
end

function test_velocities_from_nonphysical_prognostic(FT)
    @testset "Mean fall speeds from nonphysical prognostic input" begin
        params = CMP.ParametersP3(FT)
        Chen = CMP.Chen2022VelType(FT)
        # (ρq_ice, ρn_ice, ρq_rim, ρb_rim) with negative moments, a rime mass above the ice
        # mass, and a rime density above solid ice
        nonphysical = (
            (FT(1e-3), FT(1e5), FT(-3e-4), FT(1e-6)),
            (FT(-1e-3), FT(1e5), FT(3e-4), FT(1e-6)),
            (FT(1e-3), FT(-2e4), FT(3e-4), FT(1e-6)),
            (FT(1e-10), FT(1e5), FT(1e-2), FT(1e-6)),
            (FT(1e-3), FT(1e5), FT(1e-2), FT(1e-8)),
            (FT(-1e-3), FT(-2e4), FT(-3e-4), FT(-1e-6)),
        )
        for prognostic in nonphysical
            state = P3.P3State(params, P3.GaussLegendre(FT, 12), prognostic...)
            @test all(v -> isfinite(v) && v ≥ 0, P3.ice_terminal_velocities(Chen, state, FT(1)))
        end
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
                ρ_g = P3.get_ρ_g(mass, F_rim, ρ_rim)
                D_gr = P3.get_D_gr(mass, ρ_g)
                D_cr = P3.get_D_cr(mass, F_rim, ρ_g)
                @test D_th < D_gr < D_cr
                @test get_ρ_d_paper(mass; D_cr, D_gr) ≈ ρ_d
            end
        end

        # A rime density above solid ice inverts the thresholds, D_gr < D_th. The `P3State`
        # constructor bounds ρ_rim by ρ_i, which keeps them ordered.
        F_rim_bad = FT(0.93)
        @test P3.get_D_gr(mass, P3.get_ρ_g(mass, F_rim_bad, FT(975))) < D_th
        for ρ_rim in (FT(50), FT(400), ρ_i, FT(975), 2 * ρ_i)
            state = rimed_state(params, FT(1e-4), FT(1e6), F_rim_bad, ρ_rim)
            @test state.particles.ρ_rim ≈ clamp(ρ_rim, P3.rime_density_bounds(params)...)
            @test state.particles.D_th ≤ state.particles.D_gr ≤ state.particles.D_cr
        end

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
        state = rimed_state(params, L_ice, N_ice, F_rim, ρ_rim)

        # define in between values
        D_1 = D_th / 2
        D_2 = (D_th + D_gr) / 2
        D_3 = (D_gr + D_cr) / 2

        # test area
        spherical_area(D) = D^2 * π / 4
        nonspherical_area(D) = area.γ * D^area.σ
        @test P3.ice_area(state.particles, D_1) == spherical_area(D_1)
        @test P3.ice_area(state.particles, D_2) == nonspherical_area(D_2)
        @test P3.ice_area(state.particles, D_3) == spherical_area(D_3)
        @test P3.ice_area(state.particles, D_cr) == F_rim * spherical_area(D_cr) + (1 - F_rim) * nonspherical_area(D_cr)

        # test mass
        spherical_mass(ρ, D) = ρ * π / 6 * D^3
        nonspherical_mass(D) = mass.α_va * D^mass.β_va
        @test P3.ice_mass(state.particles, D_1) == spherical_mass(ρ_i, D_1)
        @test P3.ice_mass(state.particles, D_2) == nonspherical_mass(D_2)
        @test P3.ice_mass(state.particles, D_3) == spherical_mass(ρ_g, D_3)
        @test P3.ice_mass(state.particles, D_cr) == nonspherical_mass(D_cr) / (1 - F_rim)

        # test density
        @test P3.ice_density(state.particles, D_1) ≈ ρ_i
        @test P3.ice_density(state.particles, D_2) ≈ 544.916989830
        @test P3.ice_density(state.particles, D_3) ≈ ρ_g
        @test P3.ice_density(state.particles, D_cr) ≈ 383.33480937

        # test aspect ratio (oblate ϕ = 3√π m / (4 ρ A^{3/2}), ρ the per-regime
        # material density, see `P3.ϕᵢ`)
        aspect_ratio_closed(ρ, D) =
            3 * sqrt(FT(π)) * P3.ice_mass(state.particles, D) /
            (4 * ρ * P3.ice_area(state.particles, D)^FT(1.5))
        @test P3.ϕᵢ(state.particles, D_1) ≈ 1                              # D < D_th: spherical
        @test P3.ϕᵢ(state.particles, D_2) ≈ aspect_ratio_closed(ρ_i, D_2)  # dense nonspherical
        @test P3.ϕᵢ(state.particles, D_2) < 1                              # oblate
        @test P3.ϕᵢ(state.particles, D_3) ≈ 1                              # graupel: spherical (ρ_g)
        @test P3.ϕᵢ(state.particles, D_cr) ≈ aspect_ratio_closed(ρ_i, D_cr)  # partially rimed
        @test P3.ϕᵢ(state.particles, D_cr) < 1                             # oblate
        # residual ϕ > 1 band just above D_th (area discontinuity, see `P3.ϕᵢ`)
        @test 1 < P3.ϕᵢ(state.particles, D_th * FT(1.001)) < FT(1.3)

        # test F_rim = 0 and D > D_th
        state′ = rimed_state(params, L_ice, N_ice, FT(0), ρ_rim)
        @test P3.ice_area(state′.particles, D_2) == nonspherical_area(D_2)
        @test P3.ice_mass(state′.particles, D_2) == nonspherical_mass(D_2)

        # TODO: Add tests for F_liq != 0
    end
end

function test_shape_solver(FT)

    slope_laws = (CMP.SlopeConstant(FT), CMP.SlopePowerLaw(FT), CMP.SmoothSlopePowerLaw(FT))
    for slope in slope_laws
        params = CMP.ParametersP3(FT; slope)

        @testset "Shape parameters - nonlinear solver" begin
            # -- First, test limiting behavior: `N_ice = L_ice = 0` gives the upper bound --
            state = rimed_state(params, FT(0), FT(0), FT(0.5), FT(500))
            logλ = state.dist.logλ
            @test logλ == P3._logλ_max(params)
            # --

            # initialize test values:
            ep = 1 #1e4 * eps(FT)
            x_min = P3.ice_mean_particle_mass_min(params)
            x_max = P3.ice_mean_particle_mass_max(params)
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

                            state = rimed_state(params, FT(1e-4), FT(1e5), F_rim, ρ_rim) # L_ice, N_ice not used in this test
                            # Compute the shape parameters that correspond to the input test values
                            logλ_ex = log(λ_ex)
                            μ = P3.get_μ(params.slope, logλ_ex)
                            logN₀_ex = P3.get_logN₀(N_ice, μ, logλ_ex)
                            # Compute mass density based on input shape parameters
                            L_calc = exp(log(N_ice) + P3.logLdivN(state.particles, logλ_ex))

                            if x_min ≤ L_calc / N_ice ≤ x_max
                                # Solve for shape parameters
                                state′ = rimed_state(params, L_calc, N_ice, F_rim, ρ_rim)
                                logλ = state′.dist.logλ
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
            state_regr = rimed_state(params, FT(2.366e-5), FT(16461.6), FT(0.2), FT(800))
            logλ = state_regr.dist.logλ
            @test FT(2) < logλ < P3._logλ_max(params)

            # Broader sweep covering typical P3 microphysics inputs.
            # All entries must give a finite `logλ` within the search bounds.
            for L_ice in (FT(1e-6), FT(1e-5), FT(2.366e-5), FT(1e-4), FT(1e-3))
                for N_ice in (FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))
                    for F_rim in (FT(0), FT(0.2), FT(0.5), FT(0.8), FT(0.95))
                        for ρ_rim in (FT(200), FT(400), FT(600), FT(800))
                            state_sweep = rimed_state(params, L_ice, N_ice, F_rim, ρ_rim)
                            logλ = state_sweep.dist.logλ
                            @test FT(2) ≤ logλ ≤ P3._logλ_max(params)
                        end
                    end
                end
            end
        end
    end

    @testset "Upper bound of the shape solve" begin
        # At the bound, the mean particle mass is that of a newly nucleated crystal
        for slope in (CMP.SlopeConstant(FT), CMP.SlopePowerLaw(FT), CMP.SmoothSlopePowerLaw(FT))
            params = CMP.ParametersP3(FT; slope)
            state = rimed_state(params, FT(1), FT(1), FT(0), FT(500))
            m_nuc = params.ρ_i * FT(π) / 6 * params.D_nuc^3
            @test P3.logLdivN(state.particles, P3._logλ_max(params)) ≈ log(m_nuc) rtol = 8 * eps(FT)
        end
    end

    @testset "Shape solver and size distribution without ice mass or number" begin
        params = CMP.ParametersP3(FT)
        (lo, hi) = (FT(2), P3._logλ_max(params))
        state(q, n) = rimed_state(params, q, n, FT(0.5), FT(500))
        logλ(q, n) = state(q, n).dist.logλ
        N′(q, n) = P3.size_distribution(state(q, n).dist)(FT(1e-4))

        # No mass gives the smallest particles, and mass without number the largest
        @test logλ(FT(0), FT(0)) == hi
        @test logλ(FT(0), FT(1e6)) == hi
        @test logλ(FT(1e-4), FT(0)) == lo
        # A subnormal mass counts as zero
        @test logλ(nextfloat(FT(0)), FT(1e3)) == hi
        # A mean mass out of bounds is moved to the nearer bound
        @test logλ(FT(1e-20), FT(1e6)) ≈ hi rtol = 16 * eps(FT)
        @test logλ(FT(1e-4), FT(1e-20)) == logλ(FT(1e-4), FT(1e-4) / P3.ice_mean_particle_mass_max(params))
        # Without number, the distribution is zero
        @test iszero(N′(FT(1e-4), FT(0)))

        # A dual number with a zero value is also empty, although ForwardDiff compares it as
        # greater than zero when its partials are positive
        dual0 = FD.Dual(FT(0), FT(1))
        @test FD.value(logλ(dual0, FT(1e6))) == hi
        @test FD.value(logλ(FT(1e-4), dual0)) == lo
        @test isfinite(FD.derivative(x -> logλ(x, FT(1e6)), FT(0)))
        @test isfinite(FD.derivative(x -> logλ(FT(1e-4), x), FT(0)))
        @test isfinite(FD.derivative(x -> N′(FT(1e-4), x), FT(0)))
    end

    # The `N ≈ ∫N′ dD` checks hold for any root, since `logN₀` follows from `logλ`, so the
    # residual is tested directly. `logLdivN` is the log of the mean particle mass, so `expm1` of
    # the residual is its relative error. The states are heavily rimed small ice at low rime density.
    @testset "Shape solver residual at the hard states" begin
        L_ice = FT(1e-4)
        for slope in (CMP.SlopeConstant(FT), CMP.SlopePowerLaw(FT), CMP.SmoothSlopePowerLaw(FT))
            params = CMP.ParametersP3(FT; slope)
            for m̄ in FT.((1e-9, 1e-10, 1e-8)),
                F_rim in FT.((0.9, 0.99)),
                ρ_rim in (first(P3.rime_density_bounds(params)), FT(200), FT(900))

                state = rimed_state(params, L_ice, L_ice / m̄, F_rim, ρ_rim)
                dist = state.dist
                logλ = dist.logλ
                @test FT(2) < logλ < P3._logλ_max(params)
                target = log(state.ρq_ice) - log(dist.ρn)
                residual = P3.logLdivN(state.particles, logλ) - target
                @test abs(expm1(residual)) < 256 * eps(FT)
            end
        end
    end

    # Unrimed to lightly rimed ice across the `μ` ramp, where `logLdivN` flattens and then steepens.
    # The hard-clamp `SlopePowerLaw` is left out because `logLdivN` is not monotone at its corner.
    @testset "Shape solver residual across the μ ramp" begin
        L_ice = FT(1e-4)
        # The Float32 rounding floor of `logLdivN` at the knee is close to 256 eps
        rtol = FT === Float32 ? 1024 * eps(FT) : 256 * eps(FT)
        for slope in (CMP.SlopeConstant(FT), CMP.SmoothSlopePowerLaw(FT))
            params = CMP.ParametersP3(FT; slope)
            for F_rim in FT.((0, 0.1, 0.3, 0.5)),
                ρ_rim in FT.((200, 700)),
                logλ_true in FT.(range(9, 11; length = 81))

                state₀ = rimed_state(params, L_ice, FT(1e5), F_rim, ρ_rim)
                N_ice = L_ice / exp(P3.logLdivN(state₀.particles, logλ_true))
                state = rimed_state(params, L_ice, N_ice, F_rim, ρ_rim)
                dist = state.dist
                logλ = dist.logλ
                target = log(state.ρq_ice) - log(dist.ρn)
                residual = P3.logLdivN(state.particles, logλ) - target
                @test abs(expm1(residual)) < rtol
            end
        end
    end

    @testset "loggamma_inc_moment cancellation term" begin
        # For close diameters, `Δq` can round to zero or below
        D₁ = FT(1e-4)
        D₂ = nextfloat(D₁)
        for μ in FT.((0, 2, 6)), logλ in FT.((5, 10, 15)), k in (0, 2)
            val = P3.loggamma_inc_moment(D₁, D₂, μ, logλ, k)
            @test !isnan(val)
            d = FD.derivative(x -> P3.loggamma_inc_moment(D₁, x, μ, logλ, k), D₂)
            @test !isnan(d)
        end

        # An empty segment, such as [D_gr, D_cr] = [Inf, Inf] for unrimed ice, contributes
        # nothing, with zero derivatives
        @test P3.loggamma_inc_moment(D₁, D₁, FT(2), FT(10)) == -Inf
        @test iszero(FD.derivative(x -> P3.loggamma_inc_moment(D₁, D₁, FT(2), x), FT(10)))
    end
end

function test_smooth_slope_law(FT)
    smooth = CMP.SmoothSlopePowerLaw(FT)
    hard = CMP.SlopePowerLaw(FT)
    logλs = FT.(range(2, 17; length = 3001))

    @testset "SmoothSlopePowerLaw: a smooth, nondecreasing approximation of the hard limits" begin
        for logλ in logλs
            μ = P3.get_μ(smooth, logλ)
            @test 0 ≤ μ ≤ smooth.μ_max
            # A softplus corner deviates from the hard limit by at most log(2)/κ
            @test abs(μ - P3.get_μ(hard, logλ)) ≤ log(FT(2)) / smooth.κ
            @test FD.derivative(x -> P3.get_μ(smooth, x), logλ) ≥ 0
        end
    end

    @testset "SmoothSlopePowerLaw: log(L/N) decreases strictly in logλ, so the shape solve has one root" begin
        params = CMP.ParametersP3(FT; slope = smooth)
        for F_rim in FT.((0, 0.2, 0.5, 0.8, 0.95)), ρ_rim in (first(P3.rime_density_bounds(params)), FT(400), FT(800))
            state = rimed_state(params, FT(1), FT(1), F_rim, ρ_rim)
            @test all(<(0), diff(P3.logLdivN.(state.particles, logλs)))
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
        state = rimed_state(params_noar, FT(1e-4), FT(1e5), F_rim, ρ_rim)
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.4115, 0.7912, 1.1550, 1.4871]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, state.particles, ρ_a)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test vel ≈ expected[i] rtol = 1e-3
        end

        state = rimed_state(params, FT(1e-4), FT(1e5), F_rim, ρ_rim)  # aspect_ratio = Oblate() default
        # one D per P3 regime; `cbrt(ϕ) ≤ 1` slows the nonspherical sizes
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.08109, 0.38381, 0.79121, 1.155, 1.1477]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, state.particles, ρ_a)
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
        state = rimed_state(params, FT(1e-4), FT(1e5), F_rim, ρ_rim)  # aspect_ratio = Oblate() default
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13192, 0.50457, 0.90753, 1.3015, 1.6757]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, state.particles, ρ_a)
        for i in axes(Ds, 1)
            D = Ds[i]
            vel = v_term(D)
            @test vel >= 0
            @test_broken vel ≈ expected[i] rtol = 1e-3  # TODO: Implement `F_liq != 0`
        end
        state = rimed_state(
            CMP.ParametersP3(FT; aspect_ratio = CMP.NoAspectRatio()),
            FT(1e-4),
            FT(1e5),
            F_rim,
            ρ_rim,
        )
        # Allow for a D falling into every regime of the P3 Scheme
        Ds = range(FT(0.5e-4), stop = FT(4.5e-4), length = 5)
        expected = [0.13191, 0.50457, 0.90753, 1.301499, 1.67569]
        v_term = P3.ice_particle_terminal_velocity(Chen2022, state.particles, ρ_a)
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

    @testset "Mass and number weighted terminal velocities" begin

        # Without ice mass or without ice number, both mean velocities vanish
        for (ρq_ice, ρn_ice) in ((FT(0), N_ice), (L_ice, FT(0)))
            state₀ = rimed_state(params, ρq_ice, ρn_ice, FT(0.5), ρ_rim)
            (; v_n, v_m) = P3.ice_terminal_velocities(Chen2022, state₀, ρ_a)
            @test iszero(v_n)
            @test iszero(v_m)
        end

        # Finite ForwardDiff derivatives with respect to an absent quantity
        function velocities(ρq_ice, ρn_ice)
            state = rimed_state(params, ρq_ice, ρn_ice, FT(0.5), ρ_rim)
            (; v_n, v_m) = P3.ice_terminal_velocities(Chen2022, state, ρ_a)
            return (v_n, v_m)
        end
        for ρn_ice in (FT(0), N_ice)
            d = FD.derivative(x -> velocities(L_ice, x)[1], ρn_ice)
            @test !isnan(d)
        end
        for ρn_ice in (FT(0), N_ice)
            d = FD.derivative(x -> velocities(x, ρn_ice)[2], FT(0))
            @test !isnan(d)
        end

        # The aspect ratio slows the mean fall speeds
        params_noar = CMP.ParametersP3(FT; aspect_ratio = CMP.NoAspectRatio())
        for F_rim in F_rims
            v_ϕ = P3.ice_terminal_velocities(Chen2022, rimed_state(params, L_ice, N_ice, F_rim, ρ_rim), ρ_a)
            v = P3.ice_terminal_velocities(Chen2022, rimed_state(params_noar, L_ice, N_ice, F_rim, ρ_rim), ρ_a)
            @test v_ϕ.v_n <= v.v_n + eps(v.v_n) && v_ϕ.v_m <= v.v_m + eps(v.v_m)
        end

        # Output of the code at F_rim = 0.6, without and with the aspect ratio, to detect any change
        for (prm, ref) in (
            (params_noar, [2.6232178141323046, 5.7902652329103335]),
            (params, [1.4679656741018663, 2.3672147498976797]),
        )
            (; v_n, v_m) = P3.ice_terminal_velocities(Chen2022, rimed_state(prm, L_ice, N_ice, F_rims[2], ρ_rim), ρ_a)
            @test [v_n, v_m] ≈ ref rtol = 5e-6
        end
    end

    @testset "Mean fall speeds lie between the fall speeds at the integration bounds" begin
        for aspect_ratio in (CMP.Oblate(), CMP.NoAspectRatio())
            params_ar = CMP.ParametersP3(FT; aspect_ratio)
            # (L_ice, N_ice, F_rim, ρ_rim): small and large unrimed ice, rimed ice, near-graupel,
            # and two states with the mean mass out of bounds
            for (L, N, F_rim, ρ_r) in (
                (FT(1e-5), FT(1e6), FT(0), FT(400)),
                (FT(1e-3), FT(1e3), FT(0), FT(400)),
                (FT(1e-4), FT(1e5), FT(0.5), FT(800)),
                (FT(1e-3), FT(1e4), FT(0.9), FT(900)),
                (FT(1e-20), FT(1e6), FT(0.5), FT(400)),
                (FT(1e-4), FT(1e-20), FT(0.5), FT(400)),
            )
                state = rimed_state(params_ar, L, N, F_rim, ρ_r)
                v_term = P3.ice_particle_terminal_velocity(Chen2022, state.particles, ρ_a)
                bnds = P3.velocity_integral_bounds(state.particles, state.dist, v_term; p = 1e-6)
                (; v_n, v_m) = P3.ice_terminal_velocities(Chen2022, state, ρ_a)
                @test v_term(first(bnds)) ≤ v_n ≤ v_m ≤ v_term(last(bnds))
            end
        end
    end

    @testset "Mass-weighted mean diameters" begin
        ref_vals = [0.005388435466357483, 0.0033291124145735426]
        for (F_rim, ref_val) in zip(F_rims, ref_vals)
            state = rimed_state(params, L_ice, N_ice, F_rim, ρ_rim)
            @test P3.D_m(state.particles, state.dist) ≈ ref_val
        end
    end

    @testset "D_m presence gate at an absent number" begin
        state0 = rimed_state(params, L_ice, FT(0), F_rims[1], ρ_rim)
        @test iszero(P3.D_m(state0.particles, state0.dist))

        # Finite derivatives for unrimed ice, whose segment [D_gr, D_cr] = [Inf, Inf] is empty
        function Dm(ρn_ice)
            state = rimed_state(params, L_ice, ρn_ice, F_rims[1], ρ_rim)
            return P3.D_m(state.particles, state.dist)
        end
        @test isfinite(FD.derivative(Dm, FT(0)))
        @test isfinite(FD.derivative(Dm, N_ice))
    end
end

function test_numerical_integrals(FT)
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
            @test (@inferred P3.node(q, 1)) isa Float64
            @test (@inferred P3.weight(q, 1)) isa Float64
        end
    end

    @testset "Number, mean fall speeds and D_m match adaptive quadrature" begin
        N_tol = FT == Float32 ? 2e-5 : 1e-5  # native-FT gamma_inc slightly less precise than Float64-backed SF
        for aspect_ratio in (CMP.Oblate(), CMP.NoAspectRatio()),
            (F_rim, L_ice, p) in Iterators.product(F_rims, L_ices, ps)

            state = rimed_state(CMP.ParametersP3(FT; aspect_ratio), L_ice, N_ice, F_rim, ρ_rim)
            (; particles, dist) = state
            bnds = P3.integral_bounds(particles, dist; p = 1e-6)
            @test QGK.quadgk(P3.size_distribution(dist), bnds...)[1] ≈ N_ice rtol = N_tol

            v_term = P3.ice_particle_terminal_velocity(Chen2022, particles, ρ_a)
            (; v_n, v_m) = P3.ice_terminal_velocities(Chen2022, state, ρ_a; p)
            v_bnds = P3.velocity_integral_bounds(particles, dist, v_term; p)
            @test [v_n, v_m] ≈ collect(adaptive_weighted_means(v_term, state, v_bnds)) rtol = 2e-4

            @test P3.D_m(particles, dist) ≈ last(adaptive_weighted_means(identity, state, bnds)) rtol = 1e-4
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
    params = CMP.ParametersP3(FT)
    vel = CMP.Chen2022VelType(FT)
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    quad = P3.GaussLegendre(FT, 12)
    (; T_freeze) = params

    ρₐ = FT(1.2)
    Lᵢ = FT(1e-4) * ρₐ
    Nᵢ = FT(2e5) * ρₐ
    F_rim, ρ_rim = FT(0.8), FT(800)
    state = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
    melt(T, st = state) = P3.ice_melt(vel, aps, tps, st, ρₐ, T)

    @testset "Melting by heat conduction" begin
        # No melting at and below freezing
        @test melt(T_freeze - FT(0.01)) == (; dNdt = FT(0), dLdt = FT(0))
        @test melt(T_freeze) == (; dNdt = FT(0), dLdt = FT(0))

        # The mass rate is 2π K_therm ΔT / L_f ∫ D F_v N′ dD, checked with an independent quadrature
        T = T_freeze + FT(1)
        v_term = P3.ice_particle_terminal_velocity(vel, state.particles, ρₐ)
        F_v = CO.ventilation_factor(params.vent, aps, v_term)
        dist = state.dist
        N′ = P3.size_distribution(dist)
        ∫DFvN = QGK.quadgk(D -> D * F_v(D) * N′(D), FT(0), 100 / exp(dist.logλ))[1]
        @test melt(T).dLdt ≈ 2 * FT(π) * aps.K_therm * (T - T_freeze) / TDI.Lf(tps, T) * ∫DFvN rtol = 1e-3

        # The rate is linear in ΔT / L_f(T)
        T₂ = T_freeze + FT(2)
        @test melt(T₂).dLdt / melt(T).dLdt ≈ 2 * TDI.Lf(tps, T) / TDI.Lf(tps, T₂) rtol = sqrt(eps(FT))
    end

    @testset "Mass and number melt at the fractional rate of the size distribution" begin
        T = T_freeze + FT(1)
        rate = melt(T)
        @test rate.dNdt / Nᵢ ≈ rate.dLdt / Lᵢ

        # At a fixed mean particle mass, both rates scale with the moments
        rate_2 = melt(T, rimed_state(params, 2 * Lᵢ, 2 * Nᵢ, F_rim, ρ_rim))
        @test rate_2.dLdt ≈ 2 * rate.dLdt && rate_2.dNdt ≈ 2 * rate.dNdt

        # With a mean particle mass below the nucleation mass, the mass and the prognostic number
        # melt at the fractional rate of the bounded distribution
        L_small = FT(1e-12)
        state_small = rimed_state(params, L_small, Nᵢ, F_rim, ρ_rim)
        rate = melt(T, state_small)
        @test rate.dLdt / L_small ≈ rate.dNdt / Nᵢ

        # Without ice mass or without ice number, nothing melts
        for (ρq_ice, ρn_ice) in ((FT(0), Nᵢ), (Lᵢ, FT(0)))
            state_absent = rimed_state(params, ρq_ice, ρn_ice, F_rim, ρ_rim)
            @test melt(T, state_absent) == (; dNdt = FT(0), dLdt = FT(0))
        end
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

    state = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
    D̄ = exp(-state.dist.logλ)

    @testset "maximum dry freezing rate" begin
        # Below freezing, max freeze rate is non-zero (check against reference value)
        Tₐ = T_freeze - 1 // 10
        max_rate = P3.compute_max_freeze_rate(aps, tps, vel_params, state.particles, ρₐ, Tₐ)
        @test max_rate(D̄) ≈ FT(9.35962884896919e-13) rtol = 2e-4

        # At freezing, max freeze rate is zero
        Tₐ = T_freeze
        max_rate = P3.compute_max_freeze_rate(aps, tps, vel_params, state.particles, ρₐ, Tₐ)
        @test iszero(max_rate(D̄))

        # Above freezing, max freeze rate is zero
        Tₐ = T_freeze + 1 // 10
        max_rate = P3.compute_max_freeze_rate(aps, tps, vel_params, state.particles, ρₐ, Tₐ)
        @test iszero(max_rate(D̄))
    end

    @testset "local rime density" begin
        Tₐ = T_freeze - 1 // 10
        ρ′_rim_func = P3.compute_local_rime_density(vel_params, state.particles, ρₐ, Tₐ)
        # In Float32, `Tₐ - T_freeze` = -0.1 K has a relative rounding error of about 1e-4
        @test ρ′_rim_func(D̄, D̄) ≈ FT(282.8765520969483) rtol = 2e-4

        # Rᵢ > 0 for T < T_freeze, so ρ′_rim densifies toward ρ_ice as T → T_freeze.
        Dₗ = FT(200e-6)
        ρ′_rim(T) = P3.compute_local_rime_density(vel_params, state.particles, ρₐ, FT(T))(D̄, Dₗ)
        @test issorted(ρ′_rim.((240, 250, 260, 265, 270)))
        # At and above T_freeze, the solid-ice (wet-growth) limit
        @test ρ′_rim(T_freeze) ≈ FT(916.7)
        @test ρ′_rim(T_freeze + 5) ≈ FT(916.7)

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
        cloud_integrals = P3.get_liquid_integrals(n_c, ∂ₜV, m_l, ρ′_rim, liq_bounds, P3.GaussLegendre(FT, 12))
        rain_integrals = P3.get_liquid_integrals(n_r, ∂ₜV, m_l, ρ′_rim, liq_bounds, P3.GaussLegendre(FT, 12))

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
        n_zero(D) = zero(D)
        integrals_n0 = P3.get_liquid_integrals(n_zero, ∂ₜV, m_l, ρ′_rim, liq_bounds, P3.GaussLegendre(FT, 12))
        result = integrals_n0(Dᵢ)
        @test all(iszero, result)

        # Mass rate should be related to number rate through mass
        # ∂ₜM_col should be approximately ∫ n(D) * m_l(D) * ∂ₜV dD
        # This is a simplified check - in reality it's more complex

        # Test the full ∫liquid_ice_collisions function
        # Mock functions
        ∂ₜM_max(Dᵢ) = oftype(Dᵢ, 0.04)  # Small max freeze rate
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, P3.GaussLegendre(FT, 12),
        )
        @test all(x -> x isa FT, rates)  # check type stability
        @test all(>(0), rates)  # check positivity

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL = rates

        # Mass conservation: QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col


        # Since we specified identical size distributions, we expect:
        @test QCFRZ == QRFRZ
        @test QCSHD == QRSHD
        @test NCCOL == NRCOL
        @test BCCOL == BRCOL

        # Test edge cases for full collision integration

        # Zero ice content
        n_i_zero(Dᵢ) = zero(Dᵢ)
        rates = P3.∫liquid_ice_collisions(
            n_i_zero, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, P3.GaussLegendre(FT, 12),
        )
        @test all(iszero, rates)

        # Zero liquid content
        n_zero = Returns(FT(0))
        zero_liq_integrals =
            P3.get_liquid_integrals(n_zero, ∂ₜV, m_l, ρ′_rim, liq_bounds, P3.GaussLegendre(FT, 12))
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max, zero_liq_integrals, zero_liq_integrals, ice_bounds, P3.GaussLegendre(FT, 12),
        )
        @test all(iszero, rates)

        # No freezing (above freezing temperature)
        ∂ₜM_max_zero(Dᵢ) = zero(Dᵢ)
        rates = P3.∫liquid_ice_collisions(
            n_i, ∂ₜM_max_zero, cloud_integrals, rain_integrals, ice_bounds, P3.GaussLegendre(FT, 12),
        )
        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL = rates
        @test QCFRZ == 0  # No cloud freezing
        @test QRFRZ == 0  # No rain freezing
        @test QCSHD > 0  # All cloud particles should shed
        @test QRSHD > 0  # All rain particles should freeze
        @test QCSHD + QRSHD == ∫M_col  # All collisions should result shedding
    end

    @testset "Bulk liquid-ice collisions" begin
        # Test the high-level interface with real P3 parameters
        state = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim)

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
            psd_c, psd_r, aps, tps, vel_params, m_l, state, L_c, N_c, L_r, N_r, ρₐ, T,
        )
        @test eltype(rates) == FT  # check type stability

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL = rates

        # Basic sanity checks
        @test all(rates .>= 0)
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col

        # Smoke tests, aka: Check that rates don't change with new commits.
        # `rtol = 5e-4` admits both Float32 and Float64 against these (Float64)
        # reference values.
        @test QCFRZ ≈ 5.942470654139085e-7 rtol = 5e-4
        @test QCSHD ≈ 2.0729761363894677e-9 rtol = 5e-4
        @test NCCOL ≈ 60651.3567152377 rtol = 5e-4
        @test QRFRZ ≈ 6.642660875225302e-5 rtol = 5e-4
        @test QRSHD ≈ 3.6527381638257998e-6 rtol = 5e-4
        @test NRCOL ≈ 172.61819652336334 rtol = 5e-4
        @test ∫M_col ≈ 7.067566695762912e-5 rtol = 5e-4
        @test BCCOL ≈ 3.5089259323110104e-9 rtol = 5e-4
        @test BRCOL ≈ 7.24718229605571e-8 rtol = 5e-4

        ### Test the bulk source function
        rates = P3.bulk_liquid_ice_collision_sources(
            psd_c, psd_r, aps, tps, vel_params, state, L_c, N_c, L_r, N_r, ρₐ, T,
        )
        @test eltype(rates) == FT  # check type stability
    end

    # The reference integrates the outer integral adaptively over the same subintervals
    @testset "freeze/shed partition matches an adaptive reference" begin
        toml_dict = CP.create_toml_dict(Float64)
        psd_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        psd_r = CMP.RainParticlePDF_SB2006_limited(toml_dict)
        L_c, N_c, L_r, N_r = 1e-3, 1e8, 1e-4, 1e6
        params64 = CMP.ParametersP3(Float64)
        vel64 = CMP.Chen2022VelType(Float64)
        aps64 = CMP.AirProperties(Float64)
        tps64 = TDI.TD.Parameters.ThermodynamicsParameters(Float64)
        ρₐ64 = 1.2
        T = params64.T_freeze - 5
        st = rimed_state(params64, 1e-4 * ρₐ64, 2e5 * ρₐ64, 0.8, 800.0)
        (; particles, dist) = st
        m_l(D) = psd_c.ρw * CO.volume_sphere_D(D)

        ∂ₜV = P3.volumetric_collision_rate_integrand(vel64, particles, ρₐ64)
        ρ′_rim = P3.compute_local_rime_density(vel64, particles, ρₐ64, T)
        ∂ₜM_max = P3.compute_max_freeze_rate(aps64, tps64, vel64, particles, ρₐ64, T)
        n_i = DT.size_distribution(dist)
        p = 1e-5
        bounds_c = CM2.get_size_distribution_bounds(psd_c, L_c / ρₐ64, ρₐ64, N_c, p)
        bounds_r = CM2.get_size_distribution_bounds(psd_r, L_r / ρₐ64, ρₐ64, N_r, p)
        quad64 = P3.GaussLegendre(Float64, 32)
        n_c = DT.size_distribution(psd_c, L_c / ρₐ64, ρₐ64, N_c)
        n_r = DT.size_distribution(psd_r, L_r / ρₐ64, ρₐ64, N_r)
        cloud = P3.get_liquid_integrals(n_c, ∂ₜV, m_l, ρ′_rim, bounds_c, quad64)
        rain = P3.get_liquid_integrals(n_r, ∂ₜV, m_l, ρ′_rim, bounds_r, quad64)
        function partition(Dᵢ)
            (N_c_col, M_c_col, B_c_col) = cloud(Dᵢ)
            (N_r_col, M_r_col, B_r_col) = rain(Dᵢ)
            M_col = M_c_col + M_r_col
            f_frz = min(M_col, ∂ₜM_max(Dᵢ)) / M_col
            n = n_i(Dᵢ)
            return P3.SA.SVector(
                n * M_c_col * f_frz, n * M_c_col * (1 - f_frz), n * N_c_col,
                n * M_r_col * f_frz, n * M_r_col * (1 - f_frz), n * N_r_col,
                n * M_col, n * B_c_col * f_frz, n * B_r_col * f_frz,
            )
        end
        ice_bounds = P3.velocity_integral_bounds(particles, dist, ∂ₜV.v_i; p)
        D_lo, D_hi = first(ice_bounds), last(ice_bounds)
        D_wet = P3.wet_growth_onset_diameters(cloud, rain, ∂ₜM_max, D_lo, D_hi)
        @test all(D -> D_lo ≤ D ≤ D_hi, D_wet)
        breakpoints = sort(unique([ice_bounds..., D_wet...]))
        reference, _ = QGK.quadgk(partition, breakpoints...; rtol = 1e-8)
        @test reference[2] + reference[5] > 0  # wet growth sheds liquid

        state_FT =
            rimed_state(params, FT(1e-4 * ρₐ64), FT(2e5 * ρₐ64), FT(0.8), FT(800); quad = P3.GaussLegendre(FT, 16))
        ρw = FT(psd_c.ρw)
        rates = P3.∫liquid_ice_collisions(
            CMP.CloudParticlePDF_SB2006(CP.create_toml_dict(FT)),
            CMP.RainParticlePDF_SB2006_limited(CP.create_toml_dict(FT)),
            aps, tps, vel_params, D -> ρw * CO.volume_sphere_D(D),
            state_FT, FT(L_c), FT(N_c), FT(L_r), FT(N_r), FT(ρₐ64), FT(T),
        )
        for i in 1:9
            @test rates[i] ≈ reference[i] rtol = 1e-3
        end
    end

    # The state clamps ρ_rim to ρ_i, so only the prognostic rime volume shows the excess volume
    @testset "wet-growth densification uses the prognostic rime volume" begin
        toml_dict = CP.create_toml_dict(FT)
        psd_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        psd_r = CMP.RainParticlePDF_SB2006_limited(toml_dict)
        L_c, N_c, L_r, N_r = FT(1e-3), FT(1e8), FT(1e-4), FT(1e6)
        T = T_freeze - FT(5)
        ρ_i = params.ρ_i
        L_rim = Lᵢ * F_rim

        # The state `st` with the prognostic rime volume `B`
        with_ρb_rim(st, B) = P3.P3State(st.particles, st.integration, st.ρq_ice, st.ρn_ice, B, st.dist)
        collide(st, B) = P3.bulk_liquid_ice_collision_sources(
            psd_c, psd_r, aps, tps, vel_params, with_ρb_rim(st, B), L_c, N_c, L_r, N_r, ρₐ, T,
        )
        # The rime volume implied by the state's clamped rime density
        recon_B((; ρq_ice, particles)) = particles.ρ_rim > 0 ? ρq_ice * particles.F_rim / particles.ρ_rim : zero(FT)
        # The collision integrals, to separate QIWET and BIWET from the totals
        integrals(st) = P3.∫liquid_ice_collisions(
            psd_c, psd_r, aps, tps, vel_params, D -> psd_c.ρw * CO.volume_sphere_D(D), st, L_c, N_c, L_r, N_r, ρₐ, T,
        )

        st = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        r = integrals(st)
        f_shd = (r[2] + r[5]) / r[7]
        @test 0 < f_shd <= 1

        # On a consistent state the two volumes agree
        @test st.particles.ρ_rim ≈ ρ_rim rtol = 8 * eps(FT)
        @test all(Tuple(collide(st, L_rim / ρ_rim)) .=== Tuple(collide(st, recon_B(st))))

        # Above solid ice, the prognostic volume is smaller, so the increment's density is below ρ_i
        recon_B_rims = FT[]
        prog_BIWETs = FT[]
        for δ in FT[1e-3, 2e-2]
            B_hard = L_rim / (ρ_i * (1 + δ))
            st_h = P3.P3State(params, P3.GaussLegendre(FT, 12), Lᵢ, Nᵢ, L_rim, B_hard)
            @test st_h.ρb_rim == B_hard
            @test st_h.particles.ρ_rim ≈ ρ_i rtol = 8 * eps(FT)
            r_h = integrals(st_h)
            QCFRZ, QRFRZ, BCCOL, BRCOL = r_h[1], r_h[4], r_h[8], r_h[9]

            prog = collide(st_h, B_hard)
            recon = collide(st_h, recon_B(st_h))
            @test prog.∂ₜB_rim > recon.∂ₜB_rim
            @test prog.∂ₜL_rim == recon.∂ₜL_rim

            QIWET = prog.∂ₜL_rim - (QCFRZ + QRFRZ)
            BIWET_p = prog.∂ₜB_rim - (BCCOL + BRCOL)
            BIWET_r = recon.∂ₜB_rim - (BCCOL + BRCOL)
            @test QIWET > 0 && BIWET_p > 0 && BIWET_r > 0
            @test QIWET / BIWET_r ≈ ρ_i rtol = 1e-3
            @test QIWET / BIWET_p ≈ ρ_i * (1 - F_rim) / (1 - F_rim / (1 + δ)) rtol = 1e-3
            push!(recon_B_rims, recon.∂ₜB_rim)
            push!(prog_BIWETs, prog.∂ₜB_rim)
        end
        # Densification needs freezing, so it vanishes at and above the freezing temperature
        warm = P3.bulk_liquid_ice_collision_sources(
            psd_c, psd_r, aps, tps, vel_params, st, L_c, N_c, L_r, N_r, ρₐ, T_freeze + FT(1),
        )
        @test warm.∂ₜL_rim == 0 && warm.∂ₜB_rim == 0

        # The same clamped density with different prognostic volumes densifies differently
        @test recon_B_rims[1] ≈ recon_B_rims[2] rtol = 1e-3
        @test prog_BIWETs[1] != prog_BIWETs[2]
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

    state = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim)

    @testset "ice self-collection rate" begin
        # Call the new ice self-collection parameterization
        rates = P3.ice_self_collection(vel_params, state, ρₐ)
        @test eltype(rates) == FT  # check type stability

        # Self-collection should represent a positive loss rate
        @test rates.dNdt > 0

        # Test edge case with virtually zero L_ice and N_ice
        state_zero = rimed_state(params, FT(0), FT(0), F_rim, ρ_rim)
        rates_zero = P3.ice_self_collection(vel_params, state_zero, ρₐ)
        @test rates_zero.dNdt == 0

        # The triangle equals half the full square, which counts each unordered pair twice
        quad32 = P3.GaussLegendre(FT, 32)
        (; particles, dist) = state32 = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim; quad = quad32)
        rates32 = P3.ice_self_collection(vel_params, state32, ρₐ)
        n_i = DT.size_distribution(dist)
        v_i = P3.ice_particle_terminal_velocity(vel_params, particles, ρₐ)
        bnds = P3.velocity_integral_bounds(particles, dist, v_i; p = eps(one(ρₐ)))
        r(D) = sqrt(P3.ice_area(particles, D) / π)  # effective radius
        square = P3.integrate(
            D₁ -> begin
                v₁ = v_i(D₁)
                collision_rate = D₂ -> π * (r(D₁) + r(D₂))^2 * abs(v₁ - v_i(D₂)) * n_i(D₂)
                # Split the inner integral at D₂ = D₁, where |v₁ - v(D₂)| is not smooth
                inner =
                    P3.integrate(collision_rate, (first(bnds), D₁), quad32) +
                    P3.integrate(collision_rate, (D₁, last(bnds)), quad32)
                n_i(D₁) * inner
            end,
            bnds, quad32,
        )
        @test rates32.dNdt ≈ square / 2 rtol = 2e-3
    end
end

function test_p3_number_tendencies_with_bounded_number(FT)
    params = CMP.ParametersP3(FT)
    vel = CMP.Chen2022VelType(FT)
    aps = CMP.AirProperties(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    quad = P3.GaussLegendre(FT, 12)
    ρₐ = FT(1.2)
    T = params.T_freeze + FT(2)
    Lᵢ = FT(1e-4) * ρₐ
    F_rim, ρ_rim = FT(0.5), FT(500)
    x_min = P3.ice_mean_particle_mass_min(params)
    x_max = P3.ice_mean_particle_mass_max(params)

    @testset "Number tendencies scale with the prognostic number ($label)" for (label, Nᵢ) in (
        ("too few particles", Lᵢ / x_max / 10),
        ("too many particles", Lᵢ / x_min * 10),
    )
        state = rimed_state(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        dist = state.dist
        state_bounded = rimed_state(params, Lᵢ, dist.ρn, F_rim, ρ_rim)
        @test dist.ρn != Nᵢ

        # Melting removes mass and prognostic number at one fractional rate
        melt = P3.ice_melt(vel, aps, tps, state, ρₐ, T)
        f_melt = melt.dLdt / Lᵢ
        @test melt.dNdt / Nᵢ ≈ f_melt rtol = 16 * eps(FT)
        Δt = FT(0.99) / f_melt
        @test Nᵢ - Δt * melt.dNdt ≥ 0

        # Self-collection removes the prognostic number at the fractional rate of the distribution
        agg = P3.ice_self_collection(vel, state, ρₐ)
        agg_bounded = P3.ice_self_collection(vel, state_bounded, ρₐ)
        f_agg = agg_bounded.dNdt / dist.ρn
        @test agg.dNdt / Nᵢ ≈ f_agg rtol = 16 * eps(FT)
        Δt = FT(0.99) / f_agg
        @test Nᵢ - Δt * agg.dNdt ≥ 0
        # A rate of the bounded number would remove more number than exists
        if Nᵢ < dist.ρn
            @test Nᵢ - Δt * agg_bounded.dNdt < 0
        end
    end
end

function test_p3_rain_inner(FT)
    @testset "P3 rain inner (N, M, B) + cross-section coeffs" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        psd_r = CMP.SB2006(FT).pdf_r
        ρₐ = FT(1)
        ρ_w = psd_r.ρw
        p = eltype(params)(1e-5)
        m_liq(D) = ρ_w * π / 6 * D^3
        rtol = FT(5e-3)
        qrtol = FT == Float64 ? FT(1e-12) : FT(1e-7)
        Brtol = FT(1e-2)
        v_l = CO.particle_terminal_velocity(vel.rain, ρₐ)
        for (L_ice, N_ice, F_rim, ρ_rim) in (
                (1e-3, 1e6, 0.5, 500),
                (1e-2, 1e8, 0.95, 800),
                (1e-5, 1e4, 0.0, 200),
            ),
            (L_r, N_r) in ((1e-6, 1e4), (1e-4, 1e3), (2e-3, 5e2))

            state = rimed_state(params, FT(L_ice), FT(N_ice), FT(F_rim), FT(ρ_rim))
            n_r = DT.size_distribution(psd_r, FT(L_r) / ρₐ, ρₐ, FT(N_r))
            ∂ₜV = P3.volumetric_collision_rate_integrand(vel, state.particles, ρₐ)
            ρ′_rim = P3.compute_local_rime_density(vel, state.particles, ρₐ, FT(270))
            D_min, D_max = bnds = CM2.get_size_distribution_bounds(psd_r, FT(L_r) / ρₐ, ρₐ, FT(N_r), p)
            D_max > D_min || continue
            rn = P3.get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bnds, P3.ChebyshevGauss(FT, 40))
            v_i = P3.ice_particle_terminal_velocity(vel, state.particles, ρₐ)
            for Dᵢ in FT.(10 .^ range(-5, -2; length = 5))
                vi = v_i(Dᵢ)
                Dstar = P3.crossover_diameter(vi, v_l, D_min, D_max)

                # N, M, B: quadrature split at the crossing vs adaptive reference
                Nc, Mc, Bc = rn(Dᵢ)
                σ(D) = P3.collision_cross_section_ice_liquid(state.particles, Dᵢ, D)
                gN(D) = σ(D) * abs(vi - v_l(D)) * n_r(D)
                Nref = QGK.quadgk(gN, D_min, Dstar, D_max; rtol = qrtol)[1]
                Mref = QGK.quadgk(D -> gN(D) * m_liq(D), D_min, Dstar, D_max; rtol = qrtol)[1]
                @test isapprox(Nc, Nref; rtol)
                @test isapprox(Mc, Mref; rtol)
                Bref = QGK.quadgk(D -> gN(D) * m_liq(D) / ρ′_rim(Dᵢ, D), D_min, Dstar, D_max; rtol = qrtol)[1]
                @test isapprox(Bc, Bref; rtol = Brtol)

                # smoke test: collision cross section has form: `π(rᵢ + Dₗ/2)²`, with rᵢ derived from `P3.ice_area`
                rᵢ = sqrt(P3.ice_area(state.particles, Dᵢ) / FT(π))
                K = P3.collision_cross_section_ice_liquid_coeffs(state.particles, Dᵢ)
                @test K == P3.collision_cross_section_ice_liquid_coeffs(rᵢ)
                @test evalpoly(Dstar, K) ≈ FT(π) * (rᵢ + Dstar / 2)^2
                @test evalpoly(D_max, K) ≈ FT(π) * (rᵢ + D_max / 2)^2
            end
        end
    end

    @testset "P3 crossover diameter ForwardDiff" begin
        vel = CMP.Chen2022VelType(FT)
        v_l = CO.particle_terminal_velocity(vel.rain, FT(1))
        D_min, D_max = FT(1e-5), FT(5e-3)
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

    @testset "P3 crossover diameter at the bracket ends (v_i outside band)" begin
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
    end
end

function test_local_rime_density_wet_growth_limit(FT)
    params = CMP.ParametersP3(FT)
    ρ_rim_local = params.ρ_rim_local
    T_freeze = FT(params.T_freeze)
    ρ_ice = ρ_rim_local(FT(12))
    ρ_lowest = ρ_rim_local(FT(1))
    # The local rime density at a given velocity difference
    free(T°C, v, Dₗ) =
        P3.rime_density_at(P3.CollisionRimeDensity(ρ_rim_local, T°C, identity, identity), v, Dₗ)

    Dₗs = FT[0, 1e-6, 20e-6, 200e-6, 3e-3]
    vs = FT[0, 1e-2, 1e-1, 1, 5, 20]
    T°Cs = FT[-40, -10, -4.2, -1, -1e-3, -1e-4, 0, 1e-4, 1, 5]

    @testset "local rime density reaches the wet-growth limit" begin
        @testset "finite and inside the parameterization's range everywhere" begin
            for T°C in T°Cs, v in vs, Dₗ in Dₗs
                ρ′ = free(T°C, v, Dₗ)
                @test isfinite(ρ′)
                @test ρ_lowest ≤ ρ′ ≤ ρ_ice
            end
        end

        @testset "the melting point and above return the wet-growth limit" begin
            for T°C in FT[0, 1e-4, 1, 5], v in vs, Dₗ in Dₗs
                @test free(T°C, v, Dₗ) === ρ_ice
            end
        end

        @testset "densifies toward solid ice as the supercooling vanishes" begin
            for v in FT[1e-2, 1, 5], Dₗ in FT[1e-6, 20e-6, 200e-6]
                @test issorted([free(T°C, v, Dₗ) for T°C in FT[-40, -20, -10, -1, -1e-3, 0]])
            end
        end

        # The discarded branch divides by zero at the melting point
        @testset "no non-finite value or partial through ForwardDiff" begin
            for v in vs, Dₗ in Dₗs, Δ in FT[-10, -1, -1e-3, 0, 1]
                f(t) = free(t - T_freeze, v, Dₗ)
                @test isfinite(f(T_freeze + Δ))
                @test isfinite(FD.derivative(f, T_freeze + Δ))
            end
        end
    end
end

@testset "P3 tests ($FT)" for FT in (Float64, Float32)
    # state creation
    test_p3_particles(FT)
    test_state_construction(FT)

    # numerics
    test_thresholds_solver(FT)
    test_shape_solver(FT)
    test_smooth_slope_law(FT)
    test_numerical_integrals(FT)

    # velocity
    test_particle_terminal_velocities(FT)
    test_bulk_terminal_velocities(FT)
    test_velocities_from_nonphysical_prognostic(FT)

    # processes
    test_p3_het_freezing(FT)
    test_p3_melting(FT)

    # bulk liquid-ice collisions and related processes
    test_p3_bulk_liquid_ice_collisions(FT)
    test_p3_ice_self_collection(FT)
    test_p3_number_tendencies_with_bounded_number(FT)
    test_p3_rain_inner(FT)
    test_local_rime_density_wet_growth_limit(FT)
end
nothing
