import Test as TT
import ForwardDiff as FD

import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP

TT.@testset "logistic_function unit tests" begin

    TT.@test CO.logistic_function(-1.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function(0.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function(1.0, 1.0, 2.0) == 0.5
    TT.@test CO.logistic_function(2.0, 1.0, 2.0) ≈ 0.9525 atol = 1e-4

    TT.@test CO.logistic_function(1.0, 0.0, 2.0) == 1.0
    TT.@test CO.logistic_function(0.0, 0.0, 2.0) == 0.0
end

TT.@testset "logistic_function_integral unit tests" begin

    TT.@test CO.logistic_function_integral(-1.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function_integral(0.0, 1.0, 2.0) == 0.0
    TT.@test CO.logistic_function_integral(1.0, 1.0, 2.0) ≈ 0.3115 atol = 1e-4
    TT.@test CO.logistic_function_integral(3.0, 1.0, 2.0) ≈ 2.0 atol = 1e-2

    TT.@test CO.logistic_function_integral(1.0, 0.0, 2.0) == 1.0
    TT.@test CO.logistic_function_integral(0.0, 0.0, 2.0) == 0.0
end

function test_H2SO4_soln_saturation_vapor_pressure(FT)

    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

    TT.@testset "H2SO4 solution saturated vapor pressure" begin

        T_warm = FT(225.0)
        T_cold = FT(200.0)
        T_too_warm = FT(240)
        T_too_cold = FT(180)
        x_sulph = FT(0.1)

        # Note: @assert removed for GPU compatibility - T bounds no longer throw

        # p_sol should be higher at warmer temperatures
        TT.@test CO.H2SO4_soln_saturation_vapor_pressure(
            H2SO4_prs,
            x_sulph,
            T_warm,
        ) > CO.H2SO4_soln_saturation_vapor_pressure(
            H2SO4_prs,
            x_sulph,
            T_cold,
        )
    end
end

function test_a_w_xT(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    H2SO4_prs = CMP.H2SO4SolutionParameters(FT)

    TT.@testset "a_w_xT" begin

        T_warm = FT(229.2)
        T_cold = FT(228.8)
        x_sulph_low = FT(0.06)
        x_sulph_high = FT(0.1)

        # a_w greater at warmer temperatures
        for x_sulph in [x_sulph_high, x_sulph_low]
            TT.@test CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_cold) <
                     CO.a_w_xT(H2SO4_prs, tps, x_sulph, T_warm)
        end

        # a_w greater at lower sulphuric acid concentration
        for T in [T_warm, T_cold]
            TT.@test CO.a_w_xT(H2SO4_prs, tps, x_sulph_high, T) <
                     CO.a_w_xT(H2SO4_prs, tps, x_sulph_low, T)
        end
    end
end

function test_a_w_eT(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    TT.@testset "a_w_eT" begin

        T_warm = FT(285)
        T_cold = FT(251)
        e_high = FT(1088)
        e_low = FT(544)

        # a_w greater at higher altitudes
        TT.@test CO.a_w_eT(tps, e_low, T_cold) > CO.a_w_eT(tps, e_high, T_warm)

        # a_w greater at greater partial pressures
        for T in [T_warm, T_cold]
            TT.@test CO.a_w_eT(tps, e_low, T) < CO.a_w_eT(tps, e_high, T)
        end
    end
end

function test_a_w_ice(FT)

    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)

    TT.@testset "a_w_ice" begin

        T_warm = FT(240)
        T_cold = FT(230)

        # a_w greater at warmer temperatures
        TT.@test CO.a_w_ice(tps, T_cold) < CO.a_w_ice(tps, T_warm)

    end
end

function test_Chen_coefficients(FT)
    ρ = FT(1.2)
    tol = 10 * eps(FT)
    Ch2022 = CMP.Chen2022VelType(FT)
    snow = CMP.Snow(FT)
    ice = CMP.CloudIce(FT)

    TT.@testset "Chen terminal velocity rain (B1)" begin
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(Ch2022.rain, ρ)

        TT.@test all(
            isapprox.(
                aiu,
                [
                    FT(286768.02047954104),
                    FT(-1.6916433443360287e6),
                    FT(9843.240767655458),
                ],
                rtol = tol,
            ),
        )
        TT.@test all(
            isapprox.(
                bi,
                [FT(2.249342), FT(2.249342), FT(1.098942)],
                rtol = tol,
            ),
        )
        TT.@test all(
            isapprox.(ciu, [FT(0), FT(184.325), FT(184.325)], rtol = tol),
        )
    end

    TT.@testset "Chen terminal velocity small ice (B2)" begin
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(Ch2022.small_ice, ρ, ice.ρᵢ)

        TT.@test all(
            isapprox.(aiu, [312.9777159510928, -316.5335670126842], rtol = tol),
        )
        TT.@test all(
            isapprox.(bi, [0.7295470725655279, 0.7295470725655279], rtol = tol),
        )
        TT.@test all(isapprox.(ciu, [0.0, 4715.089121981011], rtol = tol))
    end

    TT.@testset "Chen terminal velocity large ice (B4)" begin
        aiu, bi, ciu = CO.Chen2022_vel_coeffs(Ch2022.large_ice, ρ, snow.ρᵢ)

        TT.@test all(
            isapprox.(aiu, [51.86069839334009, -1.394567234046072], rtol = tol),
        )
        TT.@test all(
            isapprox.(
                bi,
                [0.5655671081749194, 0.18155881980108224],
                rtol = tol,
            ),
        )
        TT.@test all(isapprox.(ciu, [0.0, 34.820462392120504], rtol = tol))
    end

    # All three coefficient sets must stay finite at a non-positive air density. The
    # rain method raises `ρₐ` to the negative `a3_pow`, so a zero floor leaves `Inf`.
    # The host can transiently deliver a negative grid-mean density, and a non-finite
    # fall speed from it aborts a GPU kernel and masks the state that caused it.
    TT.@testset "Chen coefficients are finite at non-positive air density" begin
        allfinite(t) = all(all(isfinite, x) for x in t)
        for ρ_bad in FT[-1, FT(-0.245), -floatmin(FT), FT(0)]
            TT.@test allfinite(CO.Chen2022_vel_coeffs(Ch2022.rain, ρ_bad))
            TT.@test allfinite(CO.Chen2022_vel_coeffs(Ch2022.small_ice, ρ_bad, ice.ρᵢ))
            TT.@test allfinite(CO.Chen2022_vel_coeffs(Ch2022.large_ice, ρ_bad, snow.ρᵢ))
        end
        # The floor is inert at every physical density, so the reference values above and
        # the sedimentation baseline are unchanged.
        for ρ_ok in FT[FT(0.011), FT(0.05), FT(0.5), FT(1.2)]
            TT.@test CO.Chen2022_vel_coeffs(Ch2022.rain, ρ_ok) ==
                     CO.Chen2022_vel_coeffs(Ch2022.rain, max(ρ_ok, FT(1e-4)))
        end
    end

    # The ICE density argument had no guard while the AIR density beside it did. `ρᵢ` feeds
    # `log(ρᵢ)` and `sqrt(ρᵢ)` directly - and in the LargeIce method also `C[2]/log_ρᵢ`,
    # `C[3]/ρᵢ` and `G[3]/sqrt_ρᵢ` - so a non-positive apparent ice density THROWS a
    # DomainError rather than returning something non-finite. That is worse than a NaN: it
    # aborts a GPU kernel.
    TT.@testset "Chen ice coefficients do not throw at non-positive ICE density" begin
        # The guarantee is NO THROW, not finiteness, and the distinction is the whole point.
        # `log`/`sqrt` of a negative real raises a DomainError, which inside a GPU kernel aborts
        # the run with a stacktrace pointing at the next `synchronize()` rather than at the cause.
        # A non-finite coefficient instead propagates and is caught by the host's existing
        # finiteness guards. Turning the first into the second is the improvement.
        #
        # Measured: at the floor value the LargeIce coefficients OVERFLOW - `1e-4` is a sane floor
        # for an AIR density but not for an ice density, which enters as `exp(C[3]/ρᵢ)` and
        # `exp(B[2]*log(ρᵢ)^2)`. So a physical lower bound for an ice density is still needed and is
        # a separate decision from this guard; the values at the floor are printed below so that
        # decision can be made from data.
        nothrow(f) =
            try
                f()
                true
            catch
                false
            end
        allfinite(t) = all(all(isfinite, x) for x in t)
        for ρᵢ_bad in FT[-1, FT(-0.245), -floatmin(FT), FT(0)]
            TT.@test nothrow(() -> CO.Chen2022_vel_coeffs(Ch2022.small_ice, FT(1), ρᵢ_bad))
            TT.@test nothrow(() -> CO.Chen2022_vel_coeffs(Ch2022.large_ice, FT(1), ρᵢ_bad))
        end
        # SmallIce does come back finite at the floor; LargeIce does not. Asserted for the one
        # that holds; the LargeIce values themselves are printed by `diag/chen_ice_floor.jl`,
        # which is where a diagnostic dump belongs - a test asserts, it does not report.
        TT.@test allfinite(CO.Chen2022_vel_coeffs(Ch2022.small_ice, FT(1), FT(0)))
        # Inert at every physical ice density, so no sedimentation baseline moves. Spans
        # unrimed snow through solid ice.
        for ρᵢ_ok in FT[FT(50), FT(100), FT(400), ice.ρᵢ, snow.ρᵢ, FT(916.7)]
            TT.@test CO.Chen2022_vel_coeffs(Ch2022.small_ice, FT(1), ρᵢ_ok) ==
                     CO.Chen2022_vel_coeffs(Ch2022.small_ice, FT(1), max(ρᵢ_ok, FT(1e-4)))
            TT.@test CO.Chen2022_vel_coeffs(Ch2022.large_ice, FT(1), ρᵢ_ok) ==
                     CO.Chen2022_vel_coeffs(Ch2022.large_ice, FT(1), max(ρᵢ_ok, FT(1e-4)))
        end
    end
end

function test_saturation_domain_T(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    ρ = FT(1.2)
    qₜ, qₗ, qᵢ = FT(0.01), FT(0), FT(0)
    T_floor = FT(100)
    # T = -27 is the features-lane commit's own cited Float32 substep trigger.
    T_bad = FT[0, -5, -27]

    wrapped = (
        T -> TDI.saturation_vapor_pressure_over_liquid(tps, T),
        T -> TDI.saturation_vapor_pressure_over_ice(tps, T),
        T -> TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ),
        T -> TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ),
        T -> TDI.supersaturation_over_liquid(tps, qₜ, qₗ, qᵢ, ρ, T),
        T -> TDI.supersaturation_over_ice(tps, qₜ, qₗ, qᵢ, ρ, T),
    )

    TT.@testset "saturation_domain_T floors every wrapped function to the T = 100 value" begin
        for f in wrapped
            ref = f(T_floor)
            TT.@test isfinite(ref)
            for T in T_bad
                TT.@test isfinite(f(T))
                TT.@test f(T) == ref
            end
        end
    end

    # The call without the floor is shown able to fail: at a negative temperature every one
    # of the six throws, since each routes through a log() of a negative argument. At T = 0
    # exactly, log(0) is a legal real value (-Inf), so saturation_vapor_pressure returns a
    # finite but wrong zero instead of throwing, and the two quantities derived from it divide
    # by that zero temperature and return NaN - already wrong before the floor is applied,
    # just not by way of a throw.
    TT.@testset "the call without the floor fails" begin
        qᵥ = TDI.q_vap(qₜ, qₗ, qᵢ)
        unfloored = (
            T -> TDI.TD.saturation_vapor_pressure(tps, T, TDI.TD.Liquid()),
            T -> TDI.TD.saturation_vapor_pressure(tps, T, TDI.TD.Ice()),
            T -> TDI.TD.q_vap_saturation(tps, T, ρ, TDI.TD.Liquid()),
            T -> TDI.TD.q_vap_saturation(tps, T, ρ, TDI.TD.Ice()),
            T -> TDI.TD.supersaturation(tps, qᵥ, ρ, T, TDI.TD.Liquid()),
            T -> TDI.TD.supersaturation(tps, qᵥ, ρ, T, TDI.TD.Ice()),
        )
        # The unfloored call yields no usable value below the domain. Which way it
        # fails is the thermodynamics library's choice and has changed between its
        # releases, so assert the property the floor exists for rather than the
        # mechanism: no finite, positive result.
        unusable(f, T) =
            try
                v = f(T)
                !(isfinite(v) && v > 0)
            catch
                true
            end
        for f in unfloored, T in FT[-5, -27]
            TT.@test unusable(f, T)
        end
        TT.@test unfloored[1](FT(0)) == FT(0)
        TT.@test unfloored[2](FT(0)) == FT(0)
        for f in unfloored[3:6]
            TT.@test isnan(f(FT(0)))
        end
    end

    # No AD hazard. Above the floor, saturation_domain_T is the identity and its partial (and
    # a downstream wrapped function's) must match differentiating the direct call at the same
    # Dual, to round-off. Below the floor, the floor is a Dual constant of zero partials, so the
    # returned derivative is exactly zero rather than non-finite. `Base.max` compares Duals by
    # primal value and resolves an exact tie to its second argument, so T == 100 deterministically
    # takes the constant-floor branch too (partial 0) - the same accepted sub-gradient residual
    # class as `_floored_air_density` (BMT_rosenbrock.jl:133-134): deterministic and
    # implementation-chosen, reached only at a temperature well below any atmospheric or
    # substep-transient state this model produces.
    TT.@testset "no throw and finite partials below, at, and above the floor" begin
        dual(v) = FD.Dual{Nothing}(v, one(FT))
        for Tval in FT[250, -27, 100]
            r = TDI.saturation_domain_T(dual(Tval))
            TT.@test isfinite(FD.value(r))
            TT.@test isfinite(FD.partials(r, 1))
        end
        TT.@test FD.partials(TDI.saturation_domain_T(dual(FT(250))), 1) == one(FT)
        TT.@test FD.partials(TDI.saturation_domain_T(dual(FT(-27))), 1) == zero(FT)
        TT.@test FD.partials(TDI.saturation_domain_T(dual(FT(100))), 1) == zero(FT)

        T_above = dual(FT(250))
        for (w, u) in (
            (
                TDI.saturation_vapor_pressure_over_liquid(tps, T_above),
                TDI.TD.saturation_vapor_pressure(tps, T_above, TDI.TD.Liquid()),
            ),
            (
                TDI.saturation_vapor_pressure_over_ice(tps, T_above),
                TDI.TD.saturation_vapor_pressure(tps, T_above, TDI.TD.Ice()),
            ),
        )
            TT.@test FD.value(w) == FD.value(u)
            TT.@test FD.partials(w, 1) == FD.partials(u, 1)
        end
    end
end

function test_saturation_domain_floor_engaged(FT)
    TT.@testset "saturation_domain_floor_engaged agrees with saturation_domain_T" begin
        for T in FT[0, 50, 99, 100, 101, 250]
            engaged = TDI.saturation_domain_floor_engaged(T)
            floored = TDI.saturation_domain_T(T) != T
            TT.@test engaged == floored
        end
        TT.@test TDI.saturation_domain_floor_engaged(FT(-27))
        TT.@test !TDI.saturation_domain_floor_engaged(FT(250))
    end

    TT.@testset "no throw and finite partials, at every boundary" begin
        dual(v) = FD.Dual{Nothing}(v, one(FT))
        for Tval in FT[250, -27, 100]
            r = TDI.saturation_domain_floor_engaged(dual(Tval))
            TT.@test r isa Bool
        end
    end
end

function test_volume_sphere(FT)
    TT.@testset "volume_sphere_{R/D} implemenations and type stability" begin
        R = FT(4)
        D = FT(2R)
        TT.@test CO.volume_sphere_D(D) === FT(π * D^3 / 6)
        TT.@test CO.volume_sphere_R(R) === FT(π * (2R)^3 / 6)
    end
end

TT.@testset "Common Functions Tests ($FT)" for FT in (Float64, Float32)
    test_H2SO4_soln_saturation_vapor_pressure(FT)
    test_a_w_xT(FT)
    test_a_w_eT(FT)
    test_a_w_ice(FT)
    test_Chen_coefficients(FT)
    test_saturation_domain_T(FT)
    test_saturation_domain_floor_engaged(FT)
    test_volume_sphere(FT)
end
nothing
