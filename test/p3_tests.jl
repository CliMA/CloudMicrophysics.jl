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

"`ρ_rim` as `state_from_prognostic` computes it upstream of `P3State`'s own `clamp(ρ_rim, 0, ρ_i)` -
the value that clamp receives. Equal to `state.ρ_rim` iff that clamp is a no-op."
function preclamp_ρ_rim(params, ρq_rim, ρb_rim)
    ρ_min, ρ_max = P3.rime_density_bounds(params)
    ρq_rim_c = UT.clamp_to_nonneg(ρq_rim)
    ρb_rim_c = UT.nearest_admissible_b(ρq_rim_c, ρb_rim, ρ_min, ρ_max)
    UT.rime_density(ρq_rim_c, ρb_rim_c)
end

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

function test_p3_nonphysical_state_bounds(FT)
    @testset "P3State bounds on non-physical prognostic input" begin
        params = CMP.ParametersP3(FT)
        Chen = CMP.Chen2022VelType(FT)
        ρₐ = FT(1)
        quad = P3.GaussLegendre(FT, 12)
        Ds = FT.((1e-6, 1e-4, 1e-3, 1e-2))

        # (label, ρq_ice, ρn_ice, ρq_rim, ρb_rim) — states drawn from the
        # non-physical prognostic ice moments observed at ice onset (negative
        # rime/number moments, rime mass fraction ≫ 1, rime density ≫ ρ_l)
        nonphysical = (
            (FT(1e-3), FT(1e5), FT(-3e-4), FT(1e-6)),   # neg ρq_rim, pos ρb_rim → neg ρ_rim
            (FT(1e-3), FT(1e5), FT(3e-4), FT(-1e-6)),   # neg ρb_rim
            (FT(-1e-3), FT(1e5), FT(3e-4), FT(1e-6)),   # neg ρq_ice
            (FT(1e-3), FT(-2e4), FT(3e-4), FT(1e-6)),   # neg ρn_ice
            (FT(1e-10), FT(1e5), FT(1e-2), FT(1e-6)),   # q_rim/q_ice ≈ 1e8
            (FT(1e-3), FT(1e5), FT(1e-2), FT(1e-8)),    # rime density ≈ 1e6
            (FT(1e-3), FT(1e5), FT(-1e-2), FT(1e-6)),   # large neg ρ_rim
            (FT(-1e-3), FT(-2e4), FT(-3e-4), FT(-1e-6)),  # all four negative
        )

        for (ρq_ice, ρn_ice, ρq_rim, ρb_rim) in nonphysical
            state = @inferred P3.state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
            @test FT(0) <= state.F_rim <= FT(1)
            @test FT(0) <= state.ρ_rim <= params.ρ_i
            # `P3State`'s own `clamp(ρ_rim, 0, ρ_i)` is a no-op given the upstream density-cone
            # projection: it never sees a value outside its own range to clamp.
            @test state.ρ_rim == preclamp_ρ_rim(params, ρq_rim, ρb_rim)
            @test state.ρq_ice >= FT(0)
            @test state.ρn_ice >= FT(0)
            for D in Ds
                @test isfinite(P3.ice_mass(state, D))
                @test isfinite(P3.ice_area(state, D))
                @test isfinite(P3.ice_particle_terminal_velocity(Chen, ρₐ, state)(D))
            end
            logλ = @inferred P3.get_distribution_logλ_from_prognostic(
                params, ρq_ice, ρn_ice, ρq_rim, ρb_rim,
            )
            @test isfinite(logλ)
            vn = P3.ice_terminal_velocity_number_weighted_from_prognostic(
                Chen, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; quad,
            )
            vm = P3.ice_terminal_velocity_mass_weighted_from_prognostic(
                Chen, ρₐ, params, ρq_ice, ρn_ice, ρq_rim, ρb_rim, logλ; quad,
            )
            @test isfinite(vn) && vn >= FT(0)
            @test isfinite(vm) && vm >= FT(0)
        end

        # The clamps are inert on physical states: a directly-constructed
        # physical state is unchanged.
        phys = P3.P3State(params, FT(1e-3), FT(1e5), FT(0.3), FT(400))
        @test phys.F_rim == FT(0.3)
        @test phys.ρ_rim == FT(400)
    end
end

function test_p3_state_from_prognostic_rime_pair_projection(FT)
    @testset "state_from_prognostic projects (ρq_rim, ρb_rim) jointly" begin
        params = CMP.ParametersP3(FT)
        ρ_min, ρ_max = P3.rime_density_bounds(params)
        ρq_ice, ρn_ice = FT(1e-3), FT(1e6)
        # A genuinely in-cone pair: ρq_rim / ρb_rim = 500 kg/m³ ∈ [ρ_min, ρ_max].
        ρq_rim, ρb_rim = FT(5e-4), FT(1e-6)
        Ds = FT[1e-5, 5e-5, 1e-4, 5e-4, 1e-3]

        unrimed = P3.state_from_prognostic(params, ρq_ice, ρn_ice, FT(0), FT(0))
        @test unrimed.F_rim == FT(0)
        @test unrimed.ρ_rim == FT(0)
        @test P3.isunrimed(unrimed)

        admissible = P3.state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
        @test admissible.ρ_rim == ρq_rim / ρb_rim
        @test ρ_min <= admissible.ρ_rim <= ρ_max

        # An orphan: positive mass, negative volume. Independent clamping strands
        # ρb_rim at zero and ρ_rim reads 0 despite F_rim > 0; the joint projection
        # instead lands ρ_rim inside the admissible interval.
        orphan = P3.state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, -ρb_rim)
        @test orphan.F_rim == admissible.F_rim
        @test orphan.ρ_rim > FT(0)
        @test ρ_min <= orphan.ρ_rim <= ρ_max
        @test orphan.D_gr < FT(1e-2)
        # The mass-size relation is rime-blind below D_gr by construction:
        # D < D_th is the solid-ice sphere, and the dense-rimed band
        # [D_th, D_gr) shares α_va with the unrimed branch (`ice_mass_coeffs`
        # passes the same coefficient to both slots). Rime first enters the
        # relation at D ≥ D_gr, through ρ_g and F_rim.
        for D in Ds
            if D < orphan.D_gr
                @test P3.ice_mass(orphan, D) == P3.ice_mass(unrimed, D)
            else
                @test P3.ice_mass(orphan, D) != P3.ice_mass(unrimed, D)
            end
        end

        # The mirror orphan: negative mass, positive volume. ρq_rim floors to zero
        # first, which must force the paired volume to zero too.
        mirror = P3.state_from_prognostic(params, ρq_ice, ρn_ice, -ρq_rim, ρb_rim)
        @test mirror.F_rim == FT(0)
        @test mirror.ρ_rim == FT(0)
        @test P3.isunrimed(mirror)

        @testset "P3State's own ρ_rim clamp never binds given the upstream projection" begin
            for (q, b) in ((ρq_rim, ρb_rim), (ρq_rim, -ρb_rim), (-ρq_rim, ρb_rim), (FT(0), FT(0)))
                s = P3.state_from_prognostic(params, ρq_ice, ρn_ice, q, b)
                @test s.ρ_rim == preclamp_ρ_rim(params, q, b)
            end
        end

        @testset "no non-finite partial through ForwardDiff at the projection" begin
            for (q0, b0) in ((ρq_rim, -ρb_rim), (-ρq_rim, ρb_rim), (ρq_rim, ρb_rim), (FT(0), FT(0)))
                f_q(q) = P3.state_from_prognostic(params, ρq_ice, ρn_ice, q, b0).ρ_rim
                f_b(b) = P3.state_from_prognostic(params, ρq_ice, ρn_ice, q0, b).ρ_rim
                @test isfinite(FD.derivative(f_q, q0))
                @test isfinite(FD.derivative(f_b, b0))
            end
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
                ρ_g = P3.get_ρ_g(F_rim, ρ_rim, ρ_d)
                D_gr = P3.get_D_gr(mass, ρ_g)
                D_cr = P3.get_D_cr(mass, F_rim, ρ_g)
                @test D_th < D_gr < D_cr
                @test get_ρ_d_paper(mass; D_cr, D_gr) ≈ ρ_d
            end
        end

        # The raw thresholds invert for a rime density above solid ice: `get_D_gr`
        # decreases with ρ_g, so ρ_g > ρ_i gives D_gr < D_th. The `P3State`
        # constructor bounds ρ_rim ≤ ρ_i (hence ρ_g ≤ ρ_i) to prevent this.
        F_rim_bad = FT(0.93)
        ρ_rim_bad = FT(975)  # unphysical: denser than solid ice ρ_i
        D_gr_bad = P3.get_D_gr(mass, P3.get_ρ_g(mass, F_rim_bad, ρ_rim_bad))
        @test D_gr_bad < D_th  # raw inversion for the unphysical input

        # Constructor clamp: binding above ρ_i, inert below, ordering preserved.
        mk(ρ_rim) = P3.P3State(params, FT(1e-4), FT(1e6), F_rim_bad, FT(ρ_rim))
        @test mk(975).ρ_rim == ρ_i        # binding: clamped to ρ_i
        @test mk(2 * ρ_i).ρ_rim == ρ_i    # binding: far above
        @test mk(400).ρ_rim == FT(400)    # inert: physical value unchanged
        @test mk(ρ_i).ρ_rim == ρ_i        # marginal: exactly at the bound
        for ρ_rim in (FT(50), FT(400), FT(800), ρ_i, FT(975), 2 * ρ_i)
            st = mk(ρ_rim)
            @test st.D_th ≤ st.D_gr ≤ st.D_cr  # ordering never inverts
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
            # -- First, test limiting behavior: `N_ice = L_ice = 0`. With no ice
            # the mass floor freezes the target's mass term at its ϵ-limit and
            # the number term goes to `-Inf`, so the solver returns a finite,
            # bounded logλ (the C0-continuous limit across onset) via the
            # bracket fallback below. --
            state = P3.P3State(params, FT(0), FT(0), FT(0.5), FT(500))
            logλ = P3.get_distribution_logλ(state)
            (dlo, dhi) = P3._derived_logλ_bracket(state)
            @test isfinite(logλ)
            @test dlo <= logλ <= dhi
            # `ρq_ice = ρn_ice = 0` returns the derived `logλ_max` exactly: the target is
            # `NaN`, read off directly rather than inferred from `NaN ≤ NaN` (false, which
            # selected the same endpoint before this was made explicit). Card #9's own
            # ruling: 0/0 carries no information, so it is read at the small-particle
            # (`x_min`) end for consistency with the number-keyed presence mask - a stated
            # choice, not a derivation.
            @test logλ == dhi
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
            state_regr = P3.P3State(params, FT(2.366e-5), FT(16461.6), FT(0.2), FT(800))
            logλ = P3.get_distribution_logλ(state_regr)
            (dlo_regr, dhi_regr) = P3._derived_logλ_bracket(state_regr)
            @test isfinite(logλ)
            @test dlo_regr < logλ < dhi_regr

            # Broader sweep covering typical P3 microphysics inputs.
            # All entries must give a finite `logλ` within the search bounds.
            for L_ice in (FT(1e-6), FT(1e-5), FT(2.366e-5), FT(1e-4), FT(1e-3))
                for N_ice in (FT(1e2), FT(1e3), FT(1e4), FT(1e5), FT(1e6))
                    for F_rim in (FT(0), FT(0.2), FT(0.5), FT(0.8), FT(0.95))
                        for ρ_rim in (FT(200), FT(400), FT(600), FT(800))
                            state_sweep = P3.P3State(params, L_ice, N_ice, F_rim, ρ_rim)
                            logλ = P3.get_distribution_logλ(state_sweep)
                            (dlo_sweep, dhi_sweep) = P3._derived_logλ_bracket(state_sweep)
                            @test isfinite(logλ)
                            @test dlo_sweep ≤ logλ ≤ dhi_sweep
                        end
                    end
                end
            end
        end
    end

    @testset "No-bracket fallback direction" begin
        params = CMP.ParametersP3(FT)
        L_ice = FT(0.22)
        N_ice = FT(1e6)
        F_rim = FT(0.5)
        ρ_rim = FT(800)
        # Every state in this testset shares `(F_rim, ρ_rim)`, and `_derived_logλ_bracket`
        # depends on `params`/`F_rim`/`ρ_rim` only (not on `q_ice`/`n_ice`), so the derived
        # bounds are computed once and reused.
        (dlo, dhi) = P3._derived_logλ_bracket(P3.P3State(params, FT(1e-4), FT(1e6), F_rim, ρ_rim))

        # (1) Mass-free, number-carrying: the early return, exact.
        state = P3.P3State(params, FT(0), N_ice, F_rim, ρ_rim)
        @test P3.get_distribution_logλ(state) == dhi

        # (2) The Float32-underflow entrance: `ρq_ice` a positive subnormal such that
        # `ρq_ice / ρn_ice` rounds to exactly zero at this precision, giving the same
        # `-Inf` target as (1) without triggering the early return (`ρq_ice > 0` here).
        # Only reachable at Float32 by construction - `1e-45 / 1e3` is a representable, if
        # tiny, Float64 - so this is scoped to the precision where it occurs.
        if FT == Float32
            q_sub = FT(1e-45)  # subnormal Float32; `q_sub / FT(1e3)` underflows to 0.0
            @test iszero(q_sub / FT(1e3))
            state_sub = P3.P3State(params, q_sub, FT(1e3), F_rim, ρ_rim)
            @test P3.get_distribution_logλ(state_sub) == dhi
        end

        # (3) The opposite direction: populated mass, absent number (`ρn_ice = 0`) gives
        # target `+Inf` and the smallest logλ - the same endpoint the unfixed code already
        # returned here, now explicit rather than an `abs(-Inf) ≤ abs(-Inf)` accident.
        state_noN = P3.P3State(params, L_ice, FT(0), F_rim, ρ_rim)
        @test P3.get_distribution_logλ(state_noN) == dlo

        # (4) A finite but unbracketable target (mean mass far below the nucleation mass)
        # is unaffected by this fix - still the nearest-representable-bound comparison, not
        # the new direction-aware branch. `m̄ = 2e-45 kg` is the campaign's own measured
        # no-valid-shape example.
        m̄ = FT(2e-45)
        state_unbr = P3.P3State(params, N_ice * m̄, N_ice, F_rim, ρ_rim)
        logλ_unbr = P3.get_distribution_logλ(state_unbr)
        @test logλ_unbr == dhi  # nearest bound to a target this far below the bracket's range
    end

    @testset "Derived bracket: asymmetric design (card #9)" begin
        # logλ_max (the small-mass, physical-floor end) is provably state-independent (the
        # small-spherical-ice regime's coefficients do not depend on F_rim/ρ_rim) -
        # regression-pinned across the derivation battery, not merely asserted once.
        # logλ_min stays the literal `2` at every state - a pin census on the gen-2 record
        # found deriving it from `ice_mean_particle_mass_max` (a regularization target, not a
        # physical ceiling) regressed pinning from 0.00% to 0.17%, so this end is retained
        # rather than derived (`notes/logl-derived-bracket-design.md`'s asymmetric revision).
        dhi_vals = FT[]
        dlo_vals = FT[]
        for F_rim in FT.((0, 0.5, 0.9, 0.99)), ρ_rim in FT.((50, 200, 500, 900))
            params = CMP.ParametersP3(FT)
            q_rim = F_rim > 0 ? F_rim / (1 - F_rim) * FT(1e-4) : FT(0)
            b_rim = F_rim > 0 ? q_rim / ρ_rim : FT(0)
            state = P3.state_from_prognostic(params, FT(1e-4), FT(1e6), q_rim, b_rim)
            (dlo, dhi) = P3._derived_logλ_bracket(state)
            push!(dhi_vals, dhi)
            push!(dlo_vals, dlo)
        end
        @test allequal(dhi_vals)
        @test all(==(FT(2)), dlo_vals)
    end

    @testset "C0 consistency at both degenerate corners (card #9)" begin
        # A physical sequence approaching either corner must converge to the SAME value the
        # corner's own degenerate return gives, with no jump - the requirement the derived
        # bracket exists to satisfy at a principled value rather than the bare literal `17`/`2`.
        params = CMP.ParametersP3(FT)
        q_ice, F_rim, ρ_rim = FT(1e-4), FT(0.5), FT(800)

        # Mass-free corner: q_ice -> 0 at fixed n_ice. The degenerate return is the early-return
        # path (`ρq_ice <= 0`); the sequence approaches it from `ρq_ice > 0`.
        n_ice = FT(1e6)
        (dlo_q, dhi_q) = P3._derived_logλ_bracket(P3.P3State(params, q_ice, n_ice, F_rim, ρ_rim))
        for q_ice_test in FT[1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 0]
            st = P3.P3State(params, q_ice_test, n_ice, F_rim, ρ_rim)
            @test P3.get_distribution_logλ(st) == dhi_q
        end

        # Number-free corner: n_ice -> 0 at fixed q_ice. `lo` is the retained literal `2` under
        # the asymmetric design (not state-dependent), but the C0 requirement is unchanged: no
        # early return covers this direction, so the sequence must still converge to `lo` with
        # no jump through the non-finite-target branch (`target_log_LdN = +Inf → lo`).
        (dlo_n, dhi_n) = P3._derived_logλ_bracket(P3.P3State(params, q_ice, FT(1e6), F_rim, ρ_rim))
        for n_ice_test in FT[1e-9, 1e-10, 1e-11, 1e-12, 1e-13, 1e-14, 0]
            st = P3.P3State(params, q_ice, n_ice_test, F_rim, ρ_rim)
            @test P3.get_distribution_logλ(st) == dlo_n
        end
    end

    # The fixed iteration budget in `get_distribution_logλ` is sized so the
    # returned root reaches each precision's own rounding floor. Nothing else in
    # the suite constrains the residual: the `N ≈ ∫N′ dD` checks in
    # `test_numerical_integrals` balance for any root, because `logN₀` is derived
    # from the returned `logλ`. So assert the residual directly, in the units it
    # is consumed in. `logLdivN` is `log(L/N)`, i.e. the log of the mean particle
    # mass, so `expm1` of the residual is the relative error in the mean mass
    # that sets terminal velocity and every size-dependent rate.
    #
    # The states are the measured hard cases: heavily rimed small ice at low rime
    # density, where the worst error was 7.4e-2 (hard law) / 2.3e-2 (smoothed) at
    # the old Float32 budget of 8 iterations. The 1% bound has margin over the
    # measured Float32 floor of 1.28e-3, which both slope laws share. Bounding
    # the residual, not the root, is deliberate: the root is at the rounding
    # floor and its low bits are not a guarantee.
    @testset "Shape solver residual at the hard states" begin
        L_ice = FT(1e-4)
        for slope_law in (:constant, :powerlaw, :smooth_powerlaw)
            params = CMP.ParametersP3(FT; slope_law)
            for m̄ in FT.((1e-9, 1e-10, 1e-8)),
                F_rim in FT.((0.9, 0.99)),
                ρ_rim in FT.((50, 200, 900))

                state = P3.P3State(params, L_ice, L_ice / m̄, F_rim, ρ_rim)
                logλ = P3.get_distribution_logλ(state)
                (dlo_hard, dhi_hard) = P3._derived_logλ_bracket(state)
                # Skip the states the bracket does not contain, where the solver
                # returns the nearer endpoint by design rather than a root.
                (dlo_hard < logλ < dhi_hard) || continue
                # Read the target off the state, as the solver does, rather than
                # from `m̄`, so the assertion does not also depend on the state
                # constructor reproducing the requested moments exactly.
                target = log(state.ρq_ice) - log(state.ρn_ice)
                residual = P3.logLdivN(state, logλ) - target
                @test abs(expm1(residual)) < FT(0.01)
            end
        end
    end

    @testset "loggamma_inc_moment cancellation term" begin
        # `Δq` is a difference of two regularized incomplete gamma
        # evaluations that can round to zero or slightly negative when the
        # two diameters are close together. Value and derivative must stay
        # finite regardless of which branch the clamp takes.
        D₁ = FT(1e-4)
        D₂ = nextfloat(D₁)
        for μ in FT.((0, 2, 6)), logλ in FT.((5, 10, 15)), k in (0, 2)
            val = P3.loggamma_inc_moment(D₁, D₂, μ, logλ, k)
            @test !isnan(val)
            d = FD.derivative(x -> P3.loggamma_inc_moment(D₁, x, μ, logλ, k), D₂)
            @test !isnan(d)
        end

        # The genuinely zero-width segment (`D₁ = D₂`) returns `log(0)`
        # through the early exit above the clamp, unaffected by it.
        @test P3.loggamma_inc_moment(D₁, D₁, FT(2), FT(10)) == log(FT(0))
    end

    @testset "Shape solver - number term at absent and trace population" begin
        params = CMP.ParametersP3(FT)
        ρq_ice = FT(1e-4)
        ρ_rim = FT(500)
        logλ(ρn_ice) = P3.get_distribution_logλ(
            P3.P3State(params, ρq_ice, ρn_ice, FT(0.5), ρ_rim),
        )

        # Absent number: the target diverges and the solver falls back to a
        # bracket endpoint (value and, since the fallback returns a
        # zero-partial constant, derivative both finite).
        (dlo_at, dhi_at) = P3._derived_logλ_bracket(P3.P3State(params, ρq_ice, FT(1), FT(0.5), ρ_rim))
        logλ0 = logλ(FT(0))
        @test isfinite(logλ0) && dlo_at <= logλ0 <= dhi_at
        @test !isnan(FD.derivative(logλ, FT(0)))

        # A trace population far below `eps(FT)`, where an `eps(FT)`-tied
        # floor would substitute a value independent of `ρn_ice` and would
        # differ by tens of orders of magnitude in `log` between precisions.
        # `floatmin(FT)` does not bind here, so the target reads the true
        # value at both precisions.
        ρn_trace = FT(1e-20)
        @test max(ρn_trace, floatmin(FT)) == ρn_trace
        @test isfinite(logλ(ρn_trace))
        @test !isnan(FD.derivative(logλ, ρn_trace))
    end

    @testset "size distribution presence gate at an absent number" begin
        # `get_logN₀` requires a present `N_ice` (its own docstring states the
        # precondition); this tests its caller's gate directly, at the level
        # of the size distribution itself, isolated from the velocity and
        # quadrature machinery built on top of it.
        params = CMP.ParametersP3(FT)
        state0 = P3.P3State(params, FT(1e-4), FT(0), FT(0.5), FT(500))
        logλ0 = FT(10)
        n = P3.size_distribution(state0, logλ0)
        @test iszero(n(FT(1e-4)))

        nD(ρn_ice) = P3.size_distribution(
            P3.P3State(params, FT(1e-4), ρn_ice, FT(0.5), FT(500)), logλ0,
        )(FT(1e-4))
        d = FD.derivative(nD, FT(0))
        @test !isnan(d)
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

        # Zero mass with nonzero number: the mean velocity is the finite
        # smallest-particle limit (C0-continuous across onset).
        quad = P3.GaussLegendre(FT, 12)
        state₀ = P3.P3State(params, FT(0), N_ice, FT(0.5), ρ_rim)
        logλ = P3.get_distribution_logλ(state₀)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state₀, logλ; quad)
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state₀, logλ; quad)
        @test isfinite(vel_n₀) && vel_n₀ >= 0
        @test isfinite(vel_m₀) && vel_m₀ >= 0

        # The smallest-particle limit at `L_ice = 0` must fall within 2x of the same
        # velocity evaluated at decreasing but still-populated `L_ice`, same
        # `N_ice`/`F_rim`/`ρ_rim`, converging to the same limit.
        small_masses = FT[1e-9, 1e-12, 1e-15]
        vel_n_small = map(small_masses) do L
            state = P3.P3State(params, L, N_ice, FT(0.5), ρ_rim)
            logλ_s = P3.get_distribution_logλ(state)
            P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state, logλ_s; quad)
        end
        band = vel_n_small[end]  # L_ice = 1e-15, the closest sampled point to the degenerate limit
        @test vel_n₀ <= 2 * band
        @test vel_n₀ >= band / 2

        # Zero number: no particles, so both mean velocities vanish.
        state₀ = P3.P3State(params, L_ice, FT(0), FT(0.5), ρ_rim)
        logλ = P3.get_distribution_logλ(state₀)
        vel_n₀ = P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state₀, logλ; quad = P3.GaussLegendre(FT, 12))
        vel_m₀ = P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state₀, logλ; quad = P3.GaussLegendre(FT, 12))
        @test iszero(vel_n₀)
        @test iszero(vel_m₀)

        # Value and ForwardDiff derivative both stay finite through
        # `guarded_quotient` at the same absent and trace states, w.r.t. the
        # quantity that is absent. `logλ` is solved once at a populated
        # reference state and held fixed, matching how a host consumes it
        # (a separately cached field, not re-solved per differentiation):
        # re-solving it from a state carrying the differentiated quantity
        # would also differentiate through the shape solve's own `μ`, and
        # `gamma_inc` does not support differentiating its shape parameter.
        logλ0 = P3.get_distribution_logλ(P3.P3State(params, L_ice, N_ice, FT(0.5), ρ_rim))
        function velocities(ρq_ice, ρn_ice)
            state = P3.P3State(params, ρq_ice, ρn_ice, FT(0.5), ρ_rim)
            quad = P3.GaussLegendre(FT, 12)
            return (
                P3.ice_terminal_velocity_number_weighted(Chen2022, ρ_a, state, logλ0; quad),
                P3.ice_terminal_velocity_mass_weighted(Chen2022, ρ_a, state, logλ0; quad),
            )
        end
        for ρn_ice in (FT(0), N_ice)
            d = FD.derivative(x -> velocities(L_ice, x)[1], ρn_ice)
            @test !isnan(d)
        end
        # At `ρn_ice = N_ice` (populated), differentiating the mass-weighted
        # velocity w.r.t. `ρq_ice` calls `logLdivN`, whose `segment_boundaries`
        # (no `D_min`/`D_max` given) always reaches to `D_max = Inf`, exercising
        # `gamma_inc`'s x-derivative rule at `x = Inf`.
        for ρn_ice in (FT(0), N_ice)
            d = FD.derivative(x -> velocities(x, ρn_ice)[2], FT(0))
            @test !isnan(d)
        end

        # NOTE: All reference values are output from the code.
        # A failing test indicates that the code has changed.
        # But if the changes are intentional, the reference values can be updated.

        # Liquid fraction = 0. The `_ϕ` (aspect-ratio-on) references are below
        # their aspect-off counterparts (`cbrt(ϕ) < 1`).
        # Reference values REGENERATED 2026-07-26 when SmoothSlopePowerLaw became the default
        # slope law (CM bc74de71). The smoothed law is not a pure regularization: a softplus differs
        # from a hard clamp everywhere, by log(2)/kappa = 0.259 in mu at the corners, and mu_hard sits
        # AT a corner over 91% of the logl in [2,17] bracket - so every size-distribution moment moves.
        # Measured shift here: 4.3e-4 to 2.3e-3 relative, against rtols of 5e-5 to 1e-4.
        # kappa cannot be raised to shrink it: at kappa = 4 the shape map log(L/N) is already flat at
        # Float32 and monotone only by 1e-3 at Float64, which is the multiple-solution failure the
        # smoothed law exists to remove. Old values, for audit:
        #   ref_v_n   = [3.646059575504377,  2.6191026241691695]
        #   ref_v_n_ϕ = [1.5223915218714987, 1.4656564581919258]
        #   ref_v_m   = [7.788114224053879,  5.797675366222473]
        #   ref_v_m_ϕ = [2.427666066669716,  2.3683439025452544]
        ref_v_n = [3.6498119615119333, 2.623217220648736]
        ref_v_n_ϕ = [1.5237982072911043, 1.467965341985524]
        ref_v_m = [7.780999436865279, 5.789734504936903]
        ref_v_m_ϕ = [2.4266298047731265, 2.3669977741571313]

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
        # Regenerated with the same 2026-07-26 default-slope-law change documented in the
        # terminal-velocity testset above. Old values, for audit:
        #   ref_vals = [0.005397144197921535, 0.0033368960364578005]
        ref_vals = [0.005388435466357483, 0.0033291124145735426]
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

    @testset "D_m presence gate at an absent number" begin
        # Same gate as `logN′ice` (`P3_size_distribution.jl`), at `D_m`'s own
        # call to `get_logN₀` (`P3_integral_properties.jl`).
        state0 = P3.P3State(params, L_ice, FT(0), F_rims[1], ρ_rim)
        logλ0 = P3.get_distribution_logλ(state0)
        @test iszero(P3.D_m(state0, logλ0))

        # Independent of the presence gate above: `D_m` calls
        # `logmass_gamma_moment` unconditionally before the gate runs. At
        # `F_rims[1] = 0`, `segment_boundaries` collapses `D_gr = D_cr =
        # D_max = Inf`, and `loggamma_inc_moment`'s `D₁ < D₂ || return
        # log(FT(0))` early exit sees two coinciding `Dual(Inf, 0)`
        # boundaries there - a separate site from the fixed `gamma_inc`
        # x-derivative rule (this one never reaches it, since the interval is
        # rejected before either endpoint is evaluated), with the identical
        # `log` of an exact-zero `Dual` hazard. Still open; its own future
        # card, not this one's `Utilities.jl` fix.
        Dm(ρn_ice) = P3.D_m(
            P3.P3State(params, L_ice, ρn_ice, F_rims[1], ρ_rim), logλ0,
        )
        d = FD.derivative(Dm, FT(0))
        @test_broken !isnan(d)
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
        quad = P3.GaussLegendre(FT, 12)

        T_cold = FT(273.15 - 0.01)

        rate = P3.ice_melt(vel, aps, tps, T_cold, ρₐ, state, logλ; quad)

        @test rate.dNdt == 0
        @test rate.dLdt == 0

        T_warm = FT(273.15 + 0.01)
        rate = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state, logλ; quad)

        @test rate.dNdt >= 0
        @test rate.dLdt >= 0

        # Smallest-particle fractional bound: 3 K_therm ΔT / (ρᵢ r² L_f) at r = 5e-6 m.
        # Reference values: dLdt = 1.1819847754471943e-7 (F64) / 1.18313906e-7 (F32),
        # dNdt = 236.39695508943885 (F64) / 236.62782 (F32).
        ΔT = T_warm - params.T_freeze
        L_f = TDI.Lf(tps, T_warm)
        frac_bound = 3 * aps.K_therm * ΔT / (params.ρ_i * FT(5e-6)^2 * L_f)
        frac_rate = rate.dLdt / Lᵢ
        @test frac_rate < frac_bound
        @test frac_rate > frac_bound / 10^4

        # The number melting rate is `ρn_ice * melt_frac`, `melt_frac` bounded by
        # `ice_melt_fraction_limit`. At this state the unbounded fraction `dLdt / ρq_ice`
        # sits about two orders below the bound, so the bound is dormant and `melt_frac` is
        # the plain quotient.
        lim = P3.ice_melt_fraction_limit(aps, tps, params, T_warm)
        @test rate.dNdt == state.ρn_ice * rate.melt_frac
        @test rate.melt_frac == rate.dLdt / state.ρq_ice
        @test rate.melt_frac < lim.inv_τ

        # The melt integral is temperature independent, so rates at two temperatures
        # differ exactly by the prefactor ratio ΔT / L_f(T).
        T_vwarm = FT(273.15 + 0.1)
        # Reference values before the melt-rate correction of 2026-08-06:
        #   dLdt = 8.599340191495382e-4 (F64) / 8.6005084e-4 (F32)
        #   dNdt = 1.7198680382990765e6 (F64) / 1.7201018e6 (F32)
        rate_vwarm = P3.ice_melt(vel, aps, tps, T_vwarm, ρₐ, state, logλ; quad)
        ΔT_ratio = (T_vwarm - params.T_freeze) / ΔT
        L_f_ratio = L_f / TDI.Lf(tps, T_vwarm)
        @test rate_vwarm.dLdt / rate.dLdt ≈ ΔT_ratio * L_f_ratio rtol = sqrt(eps(FT))
        @test rate_vwarm.dNdt / rate.dNdt ≈ ΔT_ratio * L_f_ratio rtol = sqrt(eps(FT))

        # The bound at ten times the excess is ten times larger and the corrected `dLdt`
        # keeps the unbounded fraction well below it there too, so the bound stays dormant.
        lim_vwarm = P3.ice_melt_fraction_limit(aps, tps, params, T_vwarm)
        @test rate_vwarm.dNdt == state.ρn_ice * rate_vwarm.melt_frac
        @test rate_vwarm.melt_frac == rate_vwarm.dLdt / state.ρq_ice
        @test rate_vwarm.melt_frac < lim_vwarm.inv_τ
    end

    @testset "Melt conduction kernel matches the deposition capacitance integral" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        quad = P3.GaussLegendre(FT, 12)

        ρₐ = FT(1.2)
        Lᵢ = FT(1e-4) * ρₐ
        Nᵢ = FT(2e5) * ρₐ
        state = P3.P3State(params, Lᵢ, Nᵢ, FT(0.8), FT(800))
        logλ = P3.get_distribution_logλ(state)

        T_warm = FT(273.15 + 0.01)
        ΔT = T_warm - params.T_freeze
        rate = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state, logλ; quad)

        # `ice_deposition_timescale` computes the same capacitance integral
        # `∫ D F_v(D) N′(D) dD` on the same state; invert it and check that `dLdt`
        # applies the Mason prefactor `2π K_therm ΔT / L_f` to the same integral.
        τ_dep = P3.ice_deposition_timescale(vel, aps, tps, T_warm, ρₐ, state, logλ; quad)
        @test !P3.ice_deposition_is_degenerate(τ_dep)
        G = CO.G_func_ice(aps, tps, T_warm)
        qᵥ_sat = TDI.saturation_vapor_specific_content_over_ice(tps, T_warm, ρₐ)
        ∫DFvN = ρₐ * qᵥ_sat / (2 * FT(π) * G * τ_dep)
        L_f = TDI.Lf(tps, T_warm)
        @test rate.dLdt ≈ 2 * FT(π) * aps.K_therm * ΔT / L_f * ∫DFvN rtol = sqrt(eps(FT))

        # The implemented prefactor against hand Mason arithmetic,
        # dm/dt = 2π D K_therm ΔT F_v / L_f at F_v = 1, ΔT = 0.01 K,
        # K_therm = 0.024 W/m/K, L_f(273.16 K) = 333600 J/kg:
        #   2π ⋅ 0.024 / 333600 = 4.5203e-7 kg m⁻¹ K⁻¹ s⁻¹
        #   D =  20 μm: dm/dt = 4.5203e-7 ⋅ 0.01 ⋅ 20e-6  = 9.0406e-14 kg/s
        #   D = 100 μm: dm/dt = 4.5203e-7 ⋅ 0.01 ⋅ 100e-6 = 4.5203e-13 kg/s
        #   D =   1 mm: dm/dt = 4.5203e-7 ⋅ 0.01 ⋅ 1e-3   = 4.5203e-12 kg/s
        c_perK = rate.dLdt / ∫DFvN / ΔT  # the prefactor `2π K_therm / L_f` as implemented
        for (D, dm_dt) in (
            (FT(20e-6), FT(9.0406e-14)),
            (FT(100e-6), FT(4.5203e-13)),
            (FT(1e-3), FT(4.5203e-12)),
        )
            @test c_perK * ΔT * D ≈ dm_dt rtol = 1e-3
        end
    end

    @testset "Melting rate scaling" begin
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        quad = P3.GaussLegendre(FT, 12)

        ρₐ = FT(1.2)
        Lᵢ = FT(1e-4) * ρₐ
        Nᵢ = FT(2e5) * ρₐ
        F_rim = FT(0.8)
        ρ_rim = FT(800)
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        logλ = P3.get_distribution_logλ(state)

        # Doubling the temperature excess doubles the rate, up to the L_f(T) factor
        ΔT = FT(1)
        T₁ = params.T_freeze + ΔT
        T₂ = params.T_freeze + 2 * ΔT
        rate₁ = P3.ice_melt(vel, aps, tps, T₁, ρₐ, state, logλ; quad)
        rate₂ = P3.ice_melt(vel, aps, tps, T₂, ρₐ, state, logλ; quad)
        @test rate₂.dLdt / rate₁.dLdt ≈ 2 * TDI.Lf(tps, T₁) / TDI.Lf(tps, T₂) rtol = sqrt(eps(FT))

        # The rate is linear in the distribution amplitude at fixed shape: doubling
        # both moments at the same `logλ` doubles `dLdt`
        state2x = P3.P3State(params, 2 * Lᵢ, 2 * Nᵢ, F_rim, ρ_rim)
        rate2x = P3.ice_melt(vel, aps, tps, T₁, ρₐ, state2x, logλ; quad)
        @test rate2x.dLdt ≈ 2 * rate₁.dLdt rtol = sqrt(eps(FT))
    end

    @testset "Melting number rate follows the bounded shared fraction" begin
        # `dNdt = ρn_ice * melt_frac` with `melt_frac = min(dLdt / ρq_ice,
        # ice_melt_fraction_limit)`: the zero-mass state melts no number (previously
        # `dLdt / x_min` through the floored mean mass), and populated states satisfy the
        # identity at the true mean mass, bounded.
        params = CMP.ParametersP3(FT)
        vel = CMP.Chen2022VelType(FT)
        aps = CMP.AirProperties(FT)
        tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
        quad = P3.GaussLegendre(FT, 12)

        ρₐ = FT(1.2)
        Nᵢ = FT(2e5) * ρₐ
        F_rim = FT(0.8)
        ρ_rim = FT(800)
        T_warm = FT(273.15 + 0.01)
        lim = P3.ice_melt_fraction_limit(aps, tps, params, T_warm)
        logλ = P3.get_distribution_logλ(P3.P3State(params, FT(1e-4) * ρₐ, Nᵢ, F_rim, ρ_rim))

        # Ice mass underflows to zero while number and the slope state survive: with no mass
        # there is no fraction, so the number rate is zero rather than `dLdt / x_min`.
        state₀ = P3.P3State(params, FT(0), Nᵢ, F_rim, ρ_rim)
        rate = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state₀, logλ; quad)
        @test rate.melt_frac == 0
        @test rate.dNdt == 0
        @test isfinite(rate.dLdt)

        # A populated state melts number in proportion to the fraction of its mass melted,
        # with the fraction bounded by the conduction limit of the nucleation size.
        state₁ = P3.P3State(params, FT(1.2e-3), FT(12), F_rim, ρ_rim)
        logλ₁ = P3.get_distribution_logλ(state₁)
        rate = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state₁, logλ₁; quad)
        @test isfinite(rate.dNdt)
        @test rate.dNdt == state₁.ρn_ice * rate.melt_frac
        @test rate.melt_frac <= lim.inv_τ
        @test rate.melt_frac ==
              min(rate.dLdt / state₁.ρq_ice, lim.inv_τ)

        # At a depleted ρn_ice the same identity holds: the number rate tracks the fraction,
        # not a floored or ceilinged mean mass.
        state_lown = P3.P3State(params, FT(1.2e-3), FT(1e-3), F_rim, ρ_rim)
        logλ_lown = P3.get_distribution_logλ(state_lown)
        rate_lown = P3.ice_melt(vel, aps, tps, T_warm, ρₐ, state_lown, logλ_lown; quad)
        @test isfinite(rate_lown.dNdt)
        @test rate_lown.dNdt == state_lown.ρn_ice * rate_lown.melt_frac
        @test rate_lown.melt_frac <= lim.inv_τ
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
        # Corrected for the Cober-List sign fix (previously pinned at the Rᵢ = 1 floor, 159.5).
        @test ρ′_rim_func(D̄, D̄) ≈ FT(282.8765520969483) rtol = 2e-4

        # Rᵢ > 0 for T < T_freeze, so ρ′_rim densifies toward ρ_ice as T → T_freeze.
        Dₗ = FT(200e-6)
        ρ′_rim(T) = P3.compute_local_rime_density(vel_params, ρₐ, FT(T), state)(D̄, Dₗ)
        @test issorted(ρ′_rim.((240, 250, 260, 265, 270)))

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

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL = rates

        # Mass conservation: QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col

        # The shed fraction, which replaced the wet growth indicator as the densification
        # driver, is bounded by the total collision rate at every state
        @test QCSHD + QRSHD <= ∫M_col

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
        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL = rates
        @test QCFRZ == 0  # No cloud freezing
        @test QRFRZ == 0  # No rain freezing
        @test QCSHD > 0  # All cloud particles should shed
        @test QRSHD > 0  # All rain particles should freeze
        @test QCSHD + QRSHD == ∫M_col  # All collisions should result shedding
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

        QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL = rates

        # Basic sanity checks
        @test all(rates .>= 0)
        @test QCFRZ + QCSHD + QRFRZ + QRSHD ≈ ∫M_col
        @test QCSHD + QRSHD <= ∫M_col

        # Smoke tests, aka: Check that rates don't change with new commits.
        # `rtol = 5e-4` admits both Float32 and Float64 against these (Float64)
        # reference values.
        @test QCFRZ ≈ 5.942471550989089e-7 rtol = 5e-4
        @test QCSHD ≈ 2.07611985935298e-9 rtol = 5e-4
        @test NCCOL ≈ 60651.35670910096 rtol = 5e-4
        @test QRFRZ ≈ 6.642674674038379e-5 rtol = 5e-4
        @test QRSHD ≈ 3.64983632601479e-6 rtol = 5e-4
        @test NRCOL ≈ 172.61819652435105 rtol = 5e-4
        @test ∫M_col ≈ 7.067566695764388e-5 rtol = 5e-4
        # BCCOL, BRCOL updated for the Cober-List sign fix in compute_local_rime_density.
        @test BCCOL ≈ 3.50892649473301e-9 rtol = 5e-4
        @test BRCOL ≈ 7.247197349759124e-8 rtol = 5e-4

        ### Test the bulk source function
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        rates = P3.bulk_liquid_ice_collision_sources(
            state, logλ,
            psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T;
            B_rim = Lᵢ * F_rim / ρ_rim, quad = P3.GaussLegendre(FT, 12),
        )
        @test eltype(rates) == FT  # check type stability
    end

    # The entry states that shed drops re-enter rain at `D_shd`. What makes that the assumption the
    # code applies, rather than one it merely states, is that the drops it adds carry a mean mass of
    # `m(D_shd)`. That is one identity and it is asserted here without a tolerance beyond floating
    # point, because it is the whole content of the assumption.
    @testset "the shed stream re-enters rain at the shedding diameter" begin
        toml_dict = CP.create_toml_dict(FT)
        psd_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        # The SHIPPED rain distribution, which `RainParticlePDF_SB2006` builds in both branches of
        # `is_limited`. The testset above this one names `RainParticlePDF_SB2006_limited`, the
        # retired Eq. 94-97 clamp cascade that nothing constructs by default, and copying that here
        # would characterise the entry on a distribution no run uses.
        psd_r = CMP.RainParticlePDF_SB2006(toml_dict)
        ρw = psd_c.ρw
        m_l(Dₗ) = ρw * CO.volume_sphere_D(Dₗ)
        m_shd = m_l(FT(1e-3))
        quad = P3.GaussLegendre(FT, 12)
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        B = Lᵢ * F_rim / ρ_rim

        rates(L_c, N_c, L_r, N_r, T) = P3.bulk_liquid_ice_collision_sources(
            state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T; B_rim = B, quad,
        )
        chan(L_c, N_c, L_r, N_r, T) = P3.∫liquid_ice_collisions(
            state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T, m_l; quad,
        )

        # Both donors present, above and below freezing, and cloud-heavy as well as rain-heavy, so
        # the identity is checked where the two shares are very different rather than only where
        # they are comparable.
        discriminated = 0
        for (L_c, N_c, L_r, N_r) in (
            (FT(1e-3), FT(1e8), FT(1e-4), FT(1e6)),
            (FT(5e-3), FT(1e8), FT(1e-6), FT(1e3)),
            (FT(1e-5), FT(1e7), FT(3e-3), FT(1e5)),
        )
            for T in (T_freeze - FT(5), T_freeze + FT(2))
                r = rates(L_c, N_c, L_r, N_r, T)
                (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, M_col, BCCOL, BRCOL) =
                    chan(L_c, N_c, L_r, N_r, T)
                # The shedding source is asserted through the tendency it enters rather than
                # recovered from it. `∂ₜN_r = -NRCOL + NRSHD` and `NRCOL` is the larger of the two
                # at most states, so `∂ₜN_r + NRCOL` cancels and reports the identity to about ten
                # units in the last place rather than to one. Comparing the whole tendency against
                # its own definition carries the same cancellation on both sides, and the absolute
                # tolerance is scaled by the larger term because the difference is what cancels.
                shed_mass = QCSHD + QRSHD
                NRSHD = shed_mass / m_shd
                # `∂ₜN_r = -NRCOL + NRSHD`, and `NRCOL` is the larger term at most states, so
                # every assertion here undoes a cancellation whose noise floor is `eps` times that
                # larger term.
                noise = eps(FT) * max(NRCOL, NRSHD, one(FT))
                @test r.∂ₜN_r ≈ -NRCOL + NRSHD atol = 8 * noise
                # The assertion above passes for any NRSHD; this one fixes WHICH mass it is formed
                # from. Dividing only the rain donor's share, which is what the entry did before,
                # gives `-NRCOL + QRSHD/m_shd`, so the two forms differ by exactly `QCSHD/m_shd`.
                # That difference is only resolvable where it clears the cancellation floor, which
                # at Float32 with a large NRCOL it does not always do; the count below requires the
                # loop to contain states where it does, rather than each state to be one.
                rain_only = -NRCOL + QRSHD / m_shd
                signal = QCSHD / m_shd
                if signal > 8 * noise
                    @test !isapprox(r.∂ₜN_r, rain_only; atol = 8 * noise)
                    @test r.∂ₜN_r - rain_only ≈ signal rtol = sqrt(eps(FT)) atol = 8 * noise
                    discriminated += 1
                end
            end
        end

        # The loop must contain at least one state where the two forms are told apart, or the
        # `signal > 8 * noise` guard has silently disabled the only assertion that fixes which mass
        # the rain number is formed from.
        @test discriminated > 0

        # THE WIRING, which no other assertion reaches. `bulk_liquid_ice_collision_sources` gains an
        # `assembly` keyword and `BMT_2mp3` resolves it from `P3IceParams.liqice_partition`; without
        # this, reverting that one line leaves every host model on the per-particle closure and the
        # suite green. The two closures are compared through the parameter set, at a state where
        # they differ.
        let toml2 = CP.create_toml_dict(FT)
            mp_pw = CMP.Microphysics2MParams(toml2; with_ice = true, quadrature_order = 12)
            mp_bk = CMP.Microphysics2MParams(toml2; with_ice = true, quadrature_order = 12,
                liqice_partition = P3.BulkPartition())
            @test mp_pw.ice.liqice_partition === nothing
            @test mp_bk.ice.liqice_partition isa P3.BulkPartition
            src(mp) = P3.bulk_liquid_ice_collision_sources(
                state, logλ, mp.ice.cloud_pdf, mp.ice.rain_pdf,
                FT(8e-3), FT(1e8), FT(4e-3), FT(1e6),
                aps, tps, mp.ice.terminal_velocity, ρₐ, T_freeze - FT(2); B_rim = B,
                quad = mp.ice.quad,
                assembly = P3._liqice_partition(mp.ice.liqice_partition, mp.ice.quad))
            pw = src(mp_pw)
            bk = src(mp_bk)
            # `∫min(a,b) ≤ min(∫a,∫b)`: the bulk form freezes at least as much, so it sheds at most
            # as much. The two must not be the same number at this state, or the field is inert.
            @test bk.∂ₜL_ice >= (1 - sqrt(eps(FT))) * pw.∂ₜL_ice
            @test bk.f_shd <= pw.f_shd
            @test bk.f_shd != pw.f_shd
        end

        # The case the earlier form could not express. With no rain population every rain term
        # vanishes, so the shed mass is entirely the cloud donor's; the entry must still add rain
        # mass and the rain number that goes with it, at the same diameter. Above freezing nothing
        # freezes, so the whole collection is shed and the state is unambiguous.
        let L_c = FT(1e-3), N_c = FT(1e8), T = T_freeze + FT(2)
            r = rates(L_c, N_c, zero(FT), zero(FT), T)
            (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, M_col, BCCOL, BRCOL) =
                chan(L_c, N_c, zero(FT), zero(FT), T)
            @test QRSHD == 0 && NRCOL == 0      # nothing to collect from an absent population
            @test QCSHD > 0                     # and the cloud donor is shedding
            @test r.∂ₜq_r > 0                   # so rain gains mass
            @test r.∂ₜN_r > 0                   # and must gain number with it
            # `NRCOL` is exactly zero here, so `∂ₜN_r` IS the shedding source and the mean mass of
            # the drops the entry adds can be read off it with nothing to cancel. This is the one
            # place the identity is exact, and it is the case the earlier form could not express.
            @test NRCOL == 0
            @test QCSHD / r.∂ₜN_r ≈ m_shd rtol = 8 * eps(FT)
        end
    end

    # The reference P3 code compares the collected mass with the freezing capacity ONCE for the
    # population; this entry compares them at every ice diameter. `BulkPartition` is the first
    # closure written as a selectable option, so what is asserted here is the property that
    # separates the two rather than a reference value: `∫min(a,b) ≤ min(∫a,∫b)` at every state,
    # which is why the bulk form freezes at least as much and sheds at most as much.
    @testset "the bulk partition is the reference closure and freezes at least as much" begin
        toml_dict = CP.create_toml_dict(FT)
        psd_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        psd_r = CMP.RainParticlePDF_SB2006_limited(toml_dict)
        m_l(Dₗ) = psd_c.ρw * CO.volume_sphere_D(Dₗ)
        quad = P3.GaussLegendre(FT, 12)
        state = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)

        chan(L_c, N_c, L_r, N_r, T, asm) = P3.∫liquid_ice_collisions(
            state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T, m_l; quad, assembly = asm,
        )

        # A liquid loading heavy enough that the capacity binds, so the two closures are
        # compared where they differ rather than where both freeze everything.
        for (L_c, N_c, L_r, N_r) in (
            (FT(1e-3), FT(1e8), FT(1e-4), FT(1e6)),
            (FT(8e-3), FT(1e8), FT(4e-3), FT(1e6)),
        )
            for ΔT in (FT(1), FT(5), FT(20))
                T = T_freeze - ΔT
                bulk = chan(L_c, N_c, L_r, N_r, T, P3.BulkPartition())
                pw = chan(L_c, N_c, L_r, N_r, T, P3.PartitionedOuter())

                # Every channel is a rate out of a liquid species, so none is negative, and the
                # two halves of each species sum back to what that species contributed.
                @test all(bulk .>= 0)
                @test bulk[7] ≈ bulk[1] + bulk[2] + bulk[4] + bulk[5] rtol = sqrt(eps(FT))

                # The partition is applied outside the integral, so the collected mass itself is
                # the same integral in both closures and differs only by the quadrature.
                @test bulk[7] ≈ pw[7] rtol = 1e-3

                # `∫min(a,b) ≤ min(∫a,∫b)`: the inequality that makes these two closures
                # different physics rather than two evaluations of one.
                @test bulk[1] + bulk[4] >= (1 - sqrt(eps(FT))) * (pw[1] + pw[4])
                @test bulk[2] + bulk[5] <= (1 + sqrt(eps(FT))) * (pw[2] + pw[5])
            end
        end

        # At and above freezing the capacity is exactly zero, so the whole collection is shed and
        # no rime volume is deposited - the same limit the per-particle closure reaches, and the
        # one the wet-growth densification gate depends on.
        warm = chan(FT(1e-3), FT(1e8), FT(1e-4), FT(1e6), T_freeze + FT(1), P3.BulkPartition())
        @test warm[1] == 0 && warm[4] == 0 && warm[8] == 0 && warm[9] == 0
        @test warm[2] + warm[5] ≈ warm[7] rtol = sqrt(eps(FT))

        # The bulk freezing capacity is the ventilation integral times a scalar, so it is the
        # population integral of the per-particle rate and not a separate parameterization.
        for ΔT in (FT(1), FT(5), FT(20))
            T = T_freeze - ΔT
            W = P3.bulk_max_freeze_rate(aps, tps, vel_params, ρₐ, T, state, logλ; quad)
            ∂ₜM_max = P3.compute_max_freeze_rate(aps, tps, vel_params, ρₐ, T, state)
            n_i = P3.size_distribution(state, logλ)
            bnds = P3.velocity_integral_bounds(
                state, logλ, P3.ice_particle_terminal_velocity(vel_params, ρₐ, state); p = 1e-6)
            @test W ≈ P3.integrate(D -> n_i(D) * ∂ₜM_max(D), bnds, quad) rtol = 1e-5
        end
        @test P3.bulk_max_freeze_rate(
            aps, tps, vel_params, ρₐ, T_freeze + FT(1), state, logλ; quad) == 0
    end

    # Wet-growth densification relaxes the `(L_rim, B_rim)` pair toward the fully-soaked solid
    # endpoint `(ρq_ice, ρq_ice/ρ_i)`. It took the rime volume as `ρq_ice·F_rim/ρ_rim`,
    # reconstructed from the state's CLAMPED and tapered quotient, rather than the prognostic
    # volume it was called with. On a consistent state the two are identical; where the
    # `ρ_rim ≤ ρ_i` clamp binds the reconstruction returns the volume the CLAMP implies, which
    # is the endpoint's own volume, so the increment's implied density is pinned at exactly ρ_i
    # and the term has no excess left to remove. With the prognostic volume the increment's
    # implied density falls strictly BELOW ρ_i there, making wet growth the one process that
    # actively pulls an above-ρ_i quotient back down whenever it fires.
    @testset "wet-growth densification uses the prognostic rime volume" begin
        toml_dict = CP.create_toml_dict(FT)
        psd_c = CMP.CloudParticlePDF_SB2006(toml_dict)
        psd_r = CMP.RainParticlePDF_SB2006_limited(toml_dict)
        L_c, N_c, L_r, N_r = FT(1e-3), FT(1e8), FT(1e-4), FT(1e6)
        T = T_freeze - FT(5)
        ρ_i = params.ρ_i
        τ_wet = params.τ_wet
        quad = P3.GaussLegendre(FT, 12)
        L_rim = Lᵢ * F_rim

        collide(st, lλ, B) = P3.bulk_liquid_ice_collision_sources(
            st, lλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
            aps, tps, vel_params, ρₐ, T; B_rim = B, quad,
        )
        # the replaced expression, verbatim, as the control
        recon_B(st) = st.ρ_rim > 0 ? st.ρq_ice * st.F_rim / st.ρ_rim : zero(FT)
        # the collision integrals a state carries, to split QIWET/BIWET off the returned totals
        integrals(st, lλ) = P3.∫liquid_ice_collisions(
            st, lλ, psd_c, psd_r, L_c, N_c, L_r, N_r, aps, tps, vel_params, ρₐ, T,
            D -> psd_c.ρw * CO.volume_sphere_D(D); quad,
        )

        # the convexity condition the mediant argument needs, at production settings: the
        # one-step update is a convex combination of the pair and the endpoint only for
        # `f_shd·h/τ_wet < 1`, and `f_shd ≤ 1` because the shed mass is a part of the
        # collected mass
        st = P3.P3State(params, Lᵢ, Nᵢ, F_rim, ρ_rim)
        lλ = P3.get_distribution_logλ(st)
        r = integrals(st, lλ)
        f_shd = (r[2] + r[5]) / r[7]
        @test 0 < f_shd <= 1
        h = FT(2)   # the production box step
        @test f_shd * h / τ_wet < 1
        @test h / τ_wet ≈ FT(0.02) rtol = 8 * eps(FT)

        # bit-identical on a consistent state: the reconstruction IS the prognostic volume there
        @test st.ρ_rim ≈ ρ_rim rtol = 8 * eps(FT)   # neither clamp nor taper binds at 800 kg/m³
        @test all(Tuple(collide(st, lλ, L_rim / ρ_rim)) .=== Tuple(collide(st, lλ, recon_B(st))))

        # above solid ice: the constructor clamps ρ_rim to ρ_i, so the reconstruction is blind to
        # the excess. The prognostic volume is smaller, so more volume is added and the
        # increment's implied density is strictly below ρ_i - it pulls the quotient down.
        recon_B_rims = FT[]
        prog_BIWETs = FT[]
        for δ in FT[1e-3, 2e-2]
            B_hard = L_rim / (ρ_i * (1 + δ))
            st_h = P3.state_from_prognostic(params, Lᵢ, Nᵢ, L_rim, B_hard)
            @test st_h.ρ_rim ≈ ρ_i rtol = 8 * eps(FT)     # the clamp binds ...
            @test st_h.ρ_rim < ρ_i * (1 + δ / 2)          # ... and it is the clamp, not the raw
            lλ_h = P3.get_distribution_logλ(st_h)
            r_h = integrals(st_h, lλ_h)
            QCFRZ, QRFRZ, BCCOL, BRCOL = r_h[1], r_h[4], r_h[8], r_h[9]

            prog = collide(st_h, lλ_h, B_hard)
            recon = collide(st_h, lλ_h, recon_B(st_h))
            @test prog.∂ₜB_rim > recon.∂ₜB_rim    # strictly more restoring
            @test prog.∂ₜL_rim == recon.∂ₜL_rim   # and the mass slot does not move

            QIWET = prog.∂ₜL_rim - (QCFRZ + QRFRZ)
            BIWET_p = prog.∂ₜB_rim - (BCCOL + BRCOL)
            BIWET_r = recon.∂ₜB_rim - (BCCOL + BRCOL)
            @test QIWET > 0 && BIWET_p > 0 && BIWET_r > 0
            # the reconstruction pins the increment at exactly the endpoint density ...
            @test QIWET / BIWET_r ≈ ρ_i rtol = 1e-3
            # ... while the prognostic volume puts it strictly below, by a margin that grows
            # with the excess: ρ_i(1 − F_rim)/(1 − F_rim/(1 + δ))
            @test QIWET / BIWET_p ≈ ρ_i * (1 - F_rim) / (1 - F_rim / (1 + δ)) rtol = 1e-3
            @test QIWET / BIWET_p < ρ_i * (1 - δ / 8)
            push!(recon_B_rims, recon.∂ₜB_rim)
            push!(prog_BIWETs, prog.∂ₜB_rim)
        end
        # the sharpest form of the defect: two states with DIFFERENT prognostic rime volumes but
        # the same clamped quotient densify identically under the reconstruction and differently
        # under the prognostic volume. The clamp, not the state, was setting the rate.
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
        bnds = P3.velocity_integral_bounds(state, logλ, v_i; p = eps(one(ρₐ)))
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
            v_i = ∂ₜV.v_i
            rc = P3.get_liquid_integrals_rain_closed(
                psd_r, n_r, ρₐ, FT(L_r), FT(N_r), state, ∂ₜV,
                m_liq, ρ′_rim, bnds; quad = P3.GaussLegendre(FT, 6),
            )
            rn = P3.get_liquid_integrals(  # numerical fallback
                n_r, ∂ₜV, m_liq, ρ′_rim, bnds;
                quad = P3.GaussLegendre(FT, 6),
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
        ∂ₜV = P3.volumetric_collision_rate_integrand(vel, FT(1), state)
        rest = (identity, identity, (FT(0), FT(1)), FT(1), FT(1e-4), FT(1e3), state)
        # closed-form eligibility: SB2006 rain PSD with a Chen velocity curve on the kernel
        m_closed = which(P3._rain_inner_integrals, typeof((psd_r, identity, ∂ₜV, rest...)))
        m_fallback = which(P3._rain_inner_integrals, typeof((1.0, identity, identity, rest...)))
        @test m_closed.sig.parameters[2] <: CMP.RainParticlePDF_SB2006
        @test m_closed.sig.parameters[4] <: P3.VolumetricCollisionRate{<:Any, <:Any, <:CO.Chen2022VelocityCurve}
        @test m_fallback.sig.parameters[2] === Any
        # and the two methods are distinct (no accidental ambiguity merge)
        @test m_closed !== m_fallback
    end
end

# The host-side admissibility repair. Every case is a state the 48 km record actually carries,
# and the assertions are about the CONTRACT (what comes back is admissible, and what was already
# admissible is untouched bit for bit) rather than about the particular numbers.
function test_admissible_ice_moments(FT)
    @testset "admissible_ice_moments" begin
        p3 = CMP.ParametersP3(FT)
        x_min = P3.ice_mean_particle_mass_min(p3)
        x_max = P3.ice_mean_particle_mass_max(FT)
        (ρ_rim_min, ρ_rim_max) = P3.rime_density_bounds(p3)
        A(q, n, qr, br) = P3.admissible_ice_moments(p3, FT(q), FT(n), FT(qr), FT(br))
        is_admissible((q, n, qr, br)) =
            (q == 0 && n == 0 && qr == 0 && br == 0) || (
                q > 0 && n > 0 && n * x_min <= q && q <= n * x_max &&
                0 <= qr <= q && (qr == 0 ? br == 0 : ρ_rim_min <= qr / br <= ρ_rim_max)
            )

        # A healthy in-cone state passes through untouched, bit for bit. This is the assertion
        # that keeps the repair from being a trajectory change in the cells that carry the ice.
        healthy = (8.0352e-4, 4.2395e6, 5.0772e-4, 1.2282e-6)
        @test A(healthy...) === map(FT, healthy)
        @test is_admissible(A(healthy...))

        # ORPHAN MASS: number destroyed by transport, mass real. The whole quartet goes, and the
        # rime goes with it rather than being left behind ice that no longer exists.
        orphan = A(1.112e-5, -2.6966e-1, 1.0891e-5, 1.3275e-8)
        @test all(iszero, orphan)
        @test is_admissible(orphan)

        # ORPHAN NUMBER: the mirror corner, fifty times more common on the record. Nothing is
        # deleted that has mass, because there is none.
        @test all(iszero, A(0, 9.965e4, 0, 0))
        @test all(iszero, A(-1e-12, 1e3, 0, 0))

        # Mean mass ABOVE the regularisation target: evacuated, and the boundary is inclusive so
        # a state sitting exactly on it survives.
        n_ref = FT(1e3)
        @test all(iszero, A(2 * n_ref * x_max, n_ref, 0, 0))
        on_max = A(n_ref * x_max, n_ref, 0, 0)
        @test on_max[1] > 0 && is_admissible(on_max)

        # Mean mass BELOW the nucleation mass: also inadmissible, and a freshly nucleated
        # population sits exactly ON that bound, where the test must be inert. This is the same
        # property the lower bound's own docstring claims, asserted here against the repair.
        @test all(iszero, A(n_ref * x_min / 2, n_ref, 0, 0))
        on_min = A(n_ref * x_min, n_ref, 0, 0)
        @test on_min[1] > 0 && is_admissible(on_min)

        # THE RIME PAIR IS NEVER THE TRIGGER. A negative rime mass on an otherwise healthy state
        # is projected, and the ice that carries it survives; deleting it instead would have cost
        # sixteen times more ice mass on the record than the pair test does.
        q, n = FT(1e-5), FT(1e3)
        (qi, ni, qr, br) = A(q, n, -1e-9, -1e-11)
        @test qi == q && ni == n            # the ice is untouched
        @test qr == 0 && br == 0            # the rime partition is emptied, not the ice
        # a rime mass over the ice mass it partitions is clipped to it, not passed on
        (_, _, qr2, br2) = A(q, n, 10 * q, 1e-8)
        @test qr2 == q
        @test ρ_rim_min <= qr2 / br2 <= ρ_rim_max

        # IDEMPOTENT, which is what makes it safe to call every step: the output of the repair is
        # a fixed point of it.
        for st in ((1.112e-5, -2.6966e-1, 1.0891e-5, 1.3275e-8), healthy,
            (1e-5, 1e3, -1e-9, -1e-11), (0, 9.965e4, 0, 0), (2e-2, 1e3, 0, 0))
            once = A(st...)
            @test P3.admissible_ice_moments(p3, once...) === once
            @test is_admissible(once)
        end

        # THE PREDICATE AND THE REPAIR ARE ONE QUESTION. A host masks its prognostic state with
        # the predicate and the kernel repairs a quartet with the function; if they could
        # disagree, the host would hand the kernel a state the kernel would then change. Asserted
        # over a grid that straddles both bounds and both signs rather than at chosen points.
        for q in FT[-1e-9, 0, 1e-20, 1e-8, 1e-5, 1e-2]
            for n in FT[-1e3, 0, 1e-20, 1e-3, 1e3, 1e6]
                ok = P3.ice_moments_are_admissible(p3, q, n)
                (qa, na, _, _) = P3.admissible_ice_moments(p3, q, n, FT(0), FT(0))
                @test ok == (qa > 0)
                @test ok ? (qa === q && na === n) : (qa == 0 && na == 0)
            end
        end

        # A NON-FINITE MOMENT SURVIVES THE REPAIR. It reports a numerical failure rather than
        # an inadmissible physical state, and the host's whole-state non-finite check is what
        # should see it; a select would write a real zero over it and destroy the only evidence
        # that something upstream broke. This is the assertion that stops the multiply being
        # turned back into a select by a reader who sees only that both give zero on a finite
        # inadmissible state.
        for bad in (FT(NaN), FT(Inf), FT(-Inf))
            (q1, n1, qr1, br1) =
                P3.admissible_ice_moments(p3, bad, FT(1e3), FT(1e-6), FT(1e-8))
            @test !isfinite(q1)
            (q2, n2, qr2, br2) =
                P3.admissible_ice_moments(p3, FT(1e-5), bad, FT(1e-6), FT(1e-8))
            @test !isfinite(n2)
        end
        # And the predicate rejects a non-finite moment, so the mask is zero there: it is the
        # multiply and not the classification that carries the value through.
        @test !P3.ice_moments_are_admissible(p3, FT(NaN), FT(1e3))
        @test !P3.ice_moments_are_admissible(p3, FT(1e-5), FT(NaN))

        # No allocation and no branch divergence: the repair is a multiply by a mask, not a
        # conditional. The mask is also what leaves a non-finite moment non-finite.
        @test 0 == @allocated P3.admissible_ice_moments(
            p3, FT(1e-5), FT(1e3), FT(1e-6), FT(1e-8))
    end
end

# The rime-density interval `[ρ_rim_min, ρ_i]` derives entirely from parameters the scheme already
# carries. Its endpoints are the implied densities of the admitted sources: freezing and
# wet-growth deposit at exactly `ρ_i`, and the only source below solid ice is dry riming through
# the Cober-List law, whose clamped domain `Rᵢ ∈ [1, 12]` puts its infimum at `ρ′_rim(1)`. Nothing
# in the scheme believes softer rime can be created, so with ray-form sinks nothing softer is
# reachable, and the interval needs no new ClimaParams key.
function test_p3_rime_density_bounds(FT)
    @testset "rime density bounds" begin
        params = CMP.ParametersP3(FT)
        (; ρ_rim_local, ρ_i) = params
        ρ_min, ρ_max = P3.rime_density_bounds(params)

        # both endpoints are functions of existing parameters, not literals
        @test ρ_max === ρ_i
        @test ρ_min === ρ_rim_local(one(FT))
        @test ρ_min === ρ_rim_local(zero(FT))     # `Rᵢ` is clamped to [1, 12] by the law itself
        @test ρ_min ≈ ρ_rim_local.a + ρ_rim_local.b + ρ_rim_local.c rtol = 8 * eps(FT)
        # the shipped CL93 defaults: 51 + 114 - 5.5
        @test ρ_min ≈ FT(159.5) rtol = 1e-6
        @test ρ_max ≈ FT(916.7) rtol = 1e-6
        @test ρ_min < ρ_max
        @test eltype((ρ_min, ρ_max)) == FT

        # it IS the infimum of the deposition law over its whole domain, and the law never
        # exceeds solid ice - so no dry-riming source pair can leave the interval
        for Rᵢ in FT[-1, 0, 1, 2, 4, 8, 8.5, 11, 12, 20]
            @test ρ_min <= ρ_rim_local(Rᵢ) <= ρ_max
        end
        @test ρ_rim_local(FT(12)) ≈ ρ_max rtol = 8 * eps(FT)   # the upper clamp reaches solid ice
    end
end

# `compute_local_rime_density` no longer bounds the surface temperature below zero. As the
# supercooling vanishes the riming index diverges, and `LocalRimeDensity` saturates it at
# `RIME_DENSITY_Rᵢ_MAX`, which returns the solid bulk ice density; at and above the melting point
# the same limit is selected directly. `rime_density_at` is exercised with the velocity difference
# supplied explicitly, so the sweep is over (Dₗ, v, T°C) rather than over the velocity laws.
function test_local_rime_density_wet_growth_limit(FT)
    params = CMP.ParametersP3(FT)
    ρ_rim_local = params.ρ_rim_local
    T_freeze = FT(params.T_freeze)
    ρ_ice = ρ_rim_local(FT(CMP.RIME_DENSITY_Rᵢ_MAX))
    ρ_lowest = ρ_rim_local(FT(CMP.RIME_DENSITY_Rᵢ_MIN))
    μm = 1_000_000

    # The retired bound, as the reference the sweep is compared against.
    T°C_bound = -FT(1e-3)
    bounded(T°C, v, Dₗ) = ρ_rim_local(-(Dₗ * μm * v) / (2 * min(T°C, T°C_bound)))
    free(T°C, v, Dₗ) =
        P3.rime_density_at(P3.RimeDensityRate(ρ_rim_local, T°C, identity, identity), v, Dₗ)

    # Below this the bounded form has not yet reached RIME_DENSITY_Rᵢ_MAX, so the two forms can
    # differ; above it the consumer saturates first and the bound decides nothing.
    Dv_shadow = 2 * abs(T°C_bound) * CMP.RIME_DENSITY_Rᵢ_MAX / μm

    Dₗs = FT[0, 1e-6, 20e-6, 200e-6, 3e-3]
    vs = FT[0, 1e-2, 1e-1, 1, 5, 20]
    T°Cs = FT[-40, -10, -4.2, -1, -1e-3, -1e-4, 0, 1e-4, 1, 5]

    @testset "local rime density reaches the wet-growth limit" begin
        @test Dv_shadow ≈ FT(2.4e-8)

        @testset "finite and inside the parameterization's range everywhere" begin
            for T°C in T°Cs, v in vs, Dₗ in Dₗs
                ρ′ = free(T°C, v, Dₗ)
                @test isfinite(ρ′)
                @test ρ_lowest ≤ ρ′ ≤ ρ_ice
            end
        end

        @testset "identical to the bounded form wherever the bound decided nothing" begin
            for T°C in T°Cs, v in vs, Dₗ in Dₗs
                if T°C ≤ T°C_bound || Dₗ * v ≥ Dv_shadow
                    @test free(T°C, v, Dₗ) === bounded(T°C, v, Dₗ)
                end
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

        # A finite returned value does not imply finite partials; the discarded arm of the select
        # divides by zero at the melting point.
        @testset "no non-finite value or partial through ForwardDiff" begin
            for v in vs, Dₗ in Dₗs, Δ in FT[-10, -1, -1e-3, 0, 1]
                f(t) = free(t - T_freeze, v, Dₗ)
                @test isfinite(f(T_freeze + Δ))
                @test isfinite(FD.derivative(f, T_freeze + Δ))
            end
        end
    end
end

function test_gamma_inc_Q_chain_ice_channel(FT)
    @testset "gamma_inc Q-chain: the ice channel's Q(1,x) is its own boundary term" begin
        # The chain's recurrence denominators, `invΓ[k] = 1/Γ(zf0 + k)`.
        invΓ_at(zf0) = ntuple(k -> FT(1) / FT(SF.gamma(Float64(zf0) + k)), 5)
        # α·D over the range the collision channels actually reach.
        xs = FT[1e-3, 0.1, 1, 2.5, 8, 11.5]

        @testset "at zf0 == 1 it returns exp(-x) exactly, not approximately" begin
            zf0 = one(FT)
            for x in xs
                q = P3.gamma_inc_Q_chain(zf0, x, invΓ_at(zf0))
                # Deliberately `===` rather than `≈`: Q(1,x) = e^{-x} is an algebraic identity, and
                # the iterative evaluation it replaces lands a few ulp away, so exact equality is
                # the assertion that distinguishes the two. `≈` would pass on either.
                @test q[1] === exp(-x)
            end
        end

        @testset "the remaining orders keep Q increasing in order" begin
            zf0 = one(FT)
            for x in xs
                q = P3.gamma_inc_Q_chain(zf0, x, invΓ_at(zf0))
                @test all(isfinite, q)
                # Q(a,x) is the upper tail of a Gamma(a,1), which grows with a at fixed x.
                @test issorted(collect(q))
            end
        end

        @testset "a non-integer zf0 takes the general branch untouched" begin
            zf0 = FT(1.35)   # a rain channel: b + 1 with b a Chen velocity exponent
            for x in xs
                q = P3.gamma_inc_Q_chain(zf0, x, invΓ_at(zf0))
                @test q[1] === UT.gamma_inc(zf0, x)[2]
                @test q[1] != exp(-x)
            end
        end
    end
end

function test_ice_sticking_efficiency(FT)
    @testset "the ice sticking efficiency reproduces the reference at its kinks and its ends" begin
        p3 = CMP.ParametersP3(FT)
        e = (T, F) -> P3.ice_sticking_efficiency(p3, FT(T), FT(F))
        (; e_cold, e_warm, T_cold, F_rim_lo, F_rim_hi) = p3.sticking
        T_frz = p3.T_freeze

        @testset "the temperature ramp, both saturated ends exact" begin
            @test e(T_cold - 20, 0) == e_cold        # saturated cold: exact, not merely close
            @test e(T_cold, 0) == e_cold             # the cold kink itself
            @test e(T_frz, 0) == e_warm              # the warm kink
            @test e(T_frz + 20, 0) == e_warm         # saturated warm
            # the reference's own arithmetic, written as it writes it: a linear ramp over the
            # 20 K span. It hard-codes the reciprocal span as 0.05; we derive it, so this asserts
            # the two agree at the values the reference ships.
            for T in (255.0, 260.0, 265.0, 270.0)
                ref = 0.001 + (T - 253.15) * (0.3 - 0.001) * 0.05
                @test e(T, 0) ≈ FT(ref) rtol = sqrt(eps(FT))
            end
        end

        @testset "the rime shutoff, zero above the cutoff exactly" begin
            @test e(T_frz, 0) == e_warm                       # unrimed: no reduction
            @test e(T_frz, F_rim_lo) == e_warm                # the lower kink
            @test e(T_frz, F_rim_hi) == 0                     # the upper kink: EXACTLY zero
            @test e(T_frz, 0.95) == 0                         # above it: exactly zero, not small
            @test e(T_frz, 1) == 0
            mid = (F_rim_lo + F_rim_hi) / 2                   # halfway down the ramp
            @test e(T_frz, mid) ≈ e_warm / 2 rtol = sqrt(eps(FT))
        end

        @testset "the two factors multiply, and the shutoff wins everywhere" begin
            # a cold, heavily rimed particle is shut off by EITHER factor; the product must be zero
            # rather than merely small, because the shutoff is exact.
            @test e(T_cold - 10, 0.95) == 0
            # and away from the ends the product is the two factors' own product
            @test e(263.15, 0.75) ≈ e(263.15, 0) * FT(0.5) rtol = sqrt(eps(FT))
        end

        @testset "it is bounded by construction over a wide sweep" begin
            for T in range(FT(200), FT(300); length = 21), F in range(FT(0), FT(1); length = 11)
                v = e(T, F)
                @test 0 <= v <= e_warm
            end
        end
    end
end

@testset "P3 tests ($FT)" for FT in (Float64, Float32)
    # state creation
    test_p3_state_creation(FT)
    test_admissible_ice_moments(FT)
    test_p3_nonphysical_state_bounds(FT)
    test_p3_state_from_prognostic_rime_pair_projection(FT)
    test_p3_rime_density_bounds(FT)

    # numerics
    test_thresholds_solver(FT)
    test_shape_solver(FT)
    test_numerical_integrals(FT)
    test_gamma_inc_Q_chain_ice_channel(FT)
    test_ice_sticking_efficiency(FT)

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
    test_local_rime_density_wet_growth_limit(FT)
end
nothing
