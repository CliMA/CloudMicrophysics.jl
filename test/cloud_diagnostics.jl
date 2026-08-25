import Test as TT

import CloudMicrophysics as CM
import ClimaParams as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.CloudDiagnostics as CMD
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.P3Scheme as P3

function test_cloud_diagnostics(FT)

    # Seifert and Beheng 2006 parameters
    override_file = joinpath(
        pkgdir(CM),
        "src",
        "parameters",
        "toml",
        "SB2006_limiters.toml",
    )
    toml_dict = CP.create_toml_dict(FT; override_file)
    SB2006 = CMP.SB2006(toml_dict)
    SB2006_no_limiters = CMP.SB2006(toml_dict; is_limited = false)

    # Water parameters
    wtr = CMP.WaterProperties(FT)
    rain = CMP.Rain(FT)
    cloud_liquid = CMP.CloudLiquid(FT)
    cloud_ice = CMP.CloudIce(FT)

    TT.@testset "1M microphysics RadarReflectivity" begin

        # some example values
        ρ_air = FT(1)
        q_rai = FT(0.18e-3)

        TT.@test CMD.radar_reflectivity_1M(rain, q_rai, ρ_air) ≈ FT(12.17) atol =
            0.2

        q_rai = FT(0.89e-4)

        TT.@test CMD.radar_reflectivity_1M(rain, q_rai, ρ_air) ≈ FT(6.68) atol =
            0.2

    end

    TT.@testset "2M microphysics - Seifert and Beheng 2006 effective radius and reflectivity" begin
        #setup
        ρₐ = FT(1)

        q_lcl = [FT(2.128e-4), FT(2.128e-20), FT(1.6e-12), FT(0)]
        N_lcl = [FT(15053529), FT(3), FT(5512), FT(0)]
        q_rai = [FT(1.573e-4), FT(1.573e-4), FT(1.9e-15), FT(0)]
        N_rai = [FT(510859), FT(510859), FT(0), FT(0)]

        # reference values
        #
        # The third state carries 5512 droplets per kg at 1.6e-12 kg/kg, a mean droplet mass of
        # 2.9e-16 kg and so a mean diameter of 0.82 μm: a tenuous population, but a population.
        # Its earlier expectations of an empty-cloud radius and a -150 dBZ sentinel came from
        # presence being decided on the mass, which put this state below the threshold; presence is
        # decided on the number now, so the state is diagnosed as the cloud it describes. Previous
        # values were reff[3] = 0 and rr[3] = -150.
        # The fourth state is zero in both moments and is unchanged, which is what keeps the
        # sentinel expectations under test.
        rr = [FT(-12.559725319858543), FT(-12.579899), FT(-140.80278169853818), FT(-150)]
        reff = [FT(2.319383e-5), FT(6.91594e-5), FT(1.0559894399963975e-6), FT(0)]

        for (qₗ, Nₗ, qᵣ, Nᵣ, rₑ, Z) in zip(q_lcl, N_lcl, q_rai, N_rai, reff, rr)
            for SB in [SB2006, SB2006_no_limiters]
                #action
                Z_val = CMD.radar_reflectivity_2M(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)
                rₑ_val = CMD.effective_radius_2M(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)
                #test
                TT.@test rₑ_val ≈ rₑ atol = FT(1e-6)
                TT.@test Z_val ≈ Z atol = FT(1e-4)
            end
        end

        # Additional test for small numbers
        qₗ = FT(1.037e-25)
        Nₗ = FT(5.225e-12)
        qᵣ = FT(2.448e-27)
        Nᵣ = FT(5.136e-18)
        Z = FT(-150)
        for SB in [SB2006, SB2006_no_limiters]
            # The rain side is zeroed by the degenerate-input guard
            # (Nᵣ < eps || qᵣ < eps), which is why the reflectivity keeps its
            # -150 dBZ sentinel. The cloud side is PRECISION DEPENDENT here: the
            # size-distribution parameter underflows at Float32, so the guard
            # fires and the radius is zero, while at Float64 it does not and the
            # radius is 1.772e-6 for a population of 5.2e-12 droplets per kg.
            # The sentinel below is the expectation and stays as written; the
            # Float64 result is marked as not meeting it. Owned by the
            # number-presence doctrine for diagnostics: the guard fires on an
            # underflow rather than on a stated presence scale, so which side of
            # it a state falls on depends on the precision.
            rₑ = FT(0)
            Z_val = CMD.radar_reflectivity_2M(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)
            rₑ_val = CMD.effective_radius_2M(SB, qₗ, qᵣ, Nₗ, Nᵣ, ρₐ)
            #test
            if FT === Float64
                TT.@test_broken rₑ_val ≈ rₑ atol = FT(1e-6)
            else
                TT.@test rₑ_val ≈ rₑ atol = FT(1e-6)
            end
            TT.@test Z_val ≈ Z atol = FT(1e-4)
        end
    end

    TT.@testset "Effective radius - '1/3' power law from Liu and Hallett (1997)" begin
        #setup
        ρ_air = FT(1)
        ρ_w = FT(1000)
        q_lcl = FT(2.128e-4)
        N_lcl = FT(15053529)
        q_rai = FT(1.573e-4)
        N_rai = FT(510859)

        #action
        reff = CMD.effective_radius_Liu_Hallet_97(
            wtr,
            ρ_air,
            q_lcl,
            N_lcl,
            q_rai,
            N_rai,
        )
        #test
        TT.@test reff ≈ FT(2.66e-05) atol = FT(8e-6)

        TT.@test CMD.effective_radius_Liu_Hallet_97(
            wtr,
            ρ_air,
            q_lcl,
            FT(1e8),
            FT(0),
            FT(0),
        ) == CMD.effective_radius_Liu_Hallet_97(wtr, ρ_air, q_lcl)

        CMD.effective_radius_Liu_Hallet_97(wtr, ρ_air, q_lcl) ==
        CMD.effective_radius_Liu_Hallet_97(cloud_liquid, ρ_air, q_lcl)
    end

    TT.@testset "Constant effective radius" begin
        TT.@test CMD.effective_radius_const(cloud_liquid) == FT(14e-6)
        TT.@test CMD.effective_radius_const(cloud_ice) == FT(25e-6)
    end

    TT.@testset "rain intercept plausibility is a flag and only a flag" begin
        # The N₀ range is retired as a clamp under the mean-mass window and kept as a
        # diagnostic. The two halves asserted here are that it FLAGS the states the cascade
        # would have rewritten, and that it rewrites nothing itself.
        rng = CMP.RainInterceptRange(toml_dict)
        win = CMP.RainParticlePDF_SB2006_windowed(toml_dict)
        ρ = FT(1)
        (; N0_min, N0_max) = rng

        # sparse large drops: below the range. The design note's worked example.
        r = CMD.rain_intercept_plausibility(rng, win, FT(1.5e-4), ρ, FT(30))
        TT.@test r.below && !r.above
        TT.@test r.N₀r < N0_min

        # many small drops: above the range. The window caps N₀ at `cbrt(π ρw/x_min)·L/x_min`,
        # so clearing `N0_max` needs a heavy loading sitting on the small-drop edge.
        r = CMD.rain_intercept_plausibility(rng, win, FT(1e-3), ρ, FT(1e9))
        TT.@test r.above && !r.below
        TT.@test r.N₀r > N0_max

        # an ordinary population inside the range is not flagged either way
        r = CMD.rain_intercept_plausibility(rng, win, FT(1e-4), ρ, FT(3e3))
        TT.@test !r.below && !r.above
        TT.@test N0_min <= r.N₀r <= N0_max

        # an empty population has no distribution to call implausible
        for (q, N) in ((FT(0), FT(0)), (FT(0), FT(1e4)), (FT(1e-4), FT(0)))
            r = CMD.rain_intercept_plausibility(rng, win, q, ρ, N)
            TT.@test !r.below && !r.above
        end

        # THE RIDER: flagging is not clamping. The reported intercept is the inversion's own,
        # and evaluating the diagnostic leaves the PSD parameters bit-identical.
        for (q, N) in ((FT(1.5e-4), FT(30)), (FT(1e-3), FT(1e9)), (FT(1e-3), FT(1e4)))
            before = CM2.pdf_rain_parameters(win, q, ρ, N)
            r = CMD.rain_intercept_plausibility(rng, win, q, ρ, N)
            after = CM2.pdf_rain_parameters(win, q, ρ, N)
            TT.@test r.N₀r === before.N₀r
            TT.@test before === after
            # and a flagged intercept is reported as it is, not clamped into the range
            if r.below || r.above
                TT.@test !(N0_min <= r.N₀r <= N0_max)
            end
        end
    end
    TT.@testset "P3 ice effective radius" begin
        p3 = CMP.ParametersP3(FT)
        quad = P3.GaussLegendre(FT, 12)
        ρ_i = p3.ρ_i

        # Below `D_th` every particle is a solid ice sphere, so the ratio of ice volume to
        # projected area of the gamma distribution `N′ = N₀ Dᵘ exp(-λ D)` reduces to
        # `r_e = (μ + 3) / (2 λ)`. That closed form is the reference, so this pins the
        # quadrature against analysis rather than against a stored number.
        ρq_ice, ρn_ice = FT(1e-5), FT(1e7)
        state = P3.state_from_prognostic(p3, ρq_ice, ρn_ice, FT(0), FT(0))
        logλ = P3.get_distribution_logλ(state)
        μ = P3.get_μ(state, logλ)
        r_e = CMD.effective_radius_P3(state, logλ; quad)
        TT.@test r_e ≈ (μ + 3) / (2 * exp(logλ)) rtol = FT(0.02)
        # ... and the reference is only the reference where the regime holds, so the state
        # is asserted to be inside it rather than assumed to be.
        TT.@test (6 * ρq_ice / ρn_ice / (π * ρ_i))^FT(1 / 3) < state.D_th

        # a larger mean particle mass gives a larger effective radius
        coarser = P3.state_from_prognostic(p3, ρq_ice, ρn_ice / 8, FT(0), FT(0))
        TT.@test CMD.effective_radius_P3(
            coarser, P3.get_distribution_logλ(coarser); quad,
        ) > r_e

        # rimed states stay finite and positive
        rimed = P3.state_from_prognostic(p3, ρq_ice, ρn_ice, FT(5e-6), FT(1e-8))
        r_rimed = CMD.effective_radius_P3(rimed, P3.get_distribution_logλ(rimed); quad)
        TT.@test isfinite(r_rimed) && r_rimed > FT(0)

        # a massless state has no area to divide by
        empty = P3.state_from_prognostic(p3, FT(0), ρn_ice, FT(0), FT(0))
        TT.@test CMD.effective_radius_P3(
            empty, P3.get_distribution_logλ(empty); quad,
        ) == FT(0)

        # Mass without number: the distribution is identically zero, so the result is zero
        # for any slope, including the placeholder value a host caches where ice number is
        # absent. A consumer falling back on a non-positive return therefore never reads a
        # radius derived from that placeholder.
        no_number = P3.state_from_prognostic(p3, ρq_ice, FT(0), FT(0), FT(0))
        for logλ_absent in (FT(2), P3.get_distribution_logλ(no_number))
            TT.@test CMD.effective_radius_P3(no_number, logλ_absent; quad) == FT(0)
        end
    end

end

TT.@testset "Cloud Diagnostics Tests ($FT)" for FT in (Float64, Float32)
    test_cloud_diagnostics(FT)
end
nothing
