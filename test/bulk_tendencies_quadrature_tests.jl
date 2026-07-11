using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI

"""
Sweep test for the P3-ice quadrature order. The quadrature rule lives on
`P3IceParams`, so the sweep constructs `Microphysics2MParams(FT; with_ice = true,
quad = CM.Quadrature.GaussLegendre(FT, n))` at several orders `n` and compares the
full BMT tendency vector against a high-order reference across a set of physically
plausible column states.

Per-order relative tolerances, applied on `max(abs(a), abs(b), ε)` with
`ε = 1e-16 · mass_scale` so tendencies that are near-zero by physics are skipped:

- `n = 100`: < 2e-3
- `n = 50`:  < 6e-3
- `n = 25`:  < 5e-2
- `n = 15`:  < 2e-1
- `n = 12`:  < 5e-3
- `n = 8`:   < 5e-3
- `n = 7`:   < 1e-2
- `n = 6`:   < 2e-2

The `n ≤ 12` rows hold because the integrand breakpoints (velocity crossing,
decay scale) keep low orders accurate; the state set includes graupel/hail
cores for this reason. Non-finite (`NaN` / `Inf`) outputs fail regardless of
tolerance.

The tolerances lock in the error study in `p3_quadrature_error_study.jl`
(included below for the shared state set); rerun that study when changing the
quadrature rule, integral bounds, or P3 process integrands.
"""

include("p3_quadrature_error_study.jl")

function compare_with_reference(tendencies_ref, tendencies_n, tol;
    mass_scale = 1e-12, report = false)
    fail_count = 0
    for field in keys(tendencies_ref)
        a = getfield(tendencies_ref, field)
        b = getfield(tendencies_n, field)
        FT = typeof(a)
        if !isfinite(a) || !isfinite(b)
            fail_count += 1
            continue
        end
        # Absolute floor: tendencies below `mass_scale` (kg/kg/s or similar)
        # are physically negligible; they often include cross-process
        # cancellation that amplifies relative error without mattering for
        # the column budget. The floor keeps the test honest: a 1e-16
        # tendency that flips sign at n=100 vs n=200 is noise, not a bug.
        scale = max(abs(a), abs(b), FT(mass_scale))
        rel = abs(a - b) / scale
        if rel > tol
            fail_count += 1
            if report
                @info "field=$field a=$a b=$b rel=$rel tol=$tol"
            end
        end
    end
    return fail_count
end

function test_quadrature_order_sweep(FT)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    # The quadrature rule lives on `P3IceParams`; build one `Microphysics2MParams`
    # per order to sweep it.
    make_mp(n) = CMP.Microphysics2MParams(FT; with_ice = true, quad = CM.Quadrature.GaussLegendre(FT, n))

    orders_and_tol = [
        (100, FT(2e-3)),
        (50, FT(6e-3)),
        (25, FT(5e-2)),
        (15, FT(2e-1)),
        # With the velocity-crossing and decay-scale breakpoints, low orders
        # stay accurate down to the default order (8), including in the
        # large-mean-size states below. See #741.
        (12, FT(5e-3)),
        (8, FT(5e-3)),
        (7, FT(1e-2)),
        (6, FT(2e-2)),
    ]
    reference_n = 200

    mp_ref = make_mp(reference_n)
    states = generate_column_states(FT)
    @test length(states) >= 10
    append!(
        states,
        hail_core_states(FT, ((0.8, 600.0, 100.0), (0.8, 600.0, 20.0), (0.95, 800.0, 50.0))),
    )

    @testset "Quadrature order sweep" begin
        for (idx, s) in enumerate(states)
            # Use a realistic logλ. For ice-bearing cells solve it; for
            # others, logλ is unused by the ice branches.
            logλ = if s.q_ice > FT(0) && s.n_ice > FT(0)
                L_ice = s.q_ice * s.ρ
                N_ice = s.n_ice * s.ρ
                F_rim = iszero(s.q_ice) ? FT(0) : s.q_rim / max(s.q_ice, eps(FT))
                ρ_rim = iszero(s.b_rim) ? FT(0) : s.q_rim / max(s.b_rim, eps(FT))
                F_rim = min(F_rim, FT(0.99))
                ρ_rim = clamp(ρ_rim, FT(0), FT(0.8) * mp_ref.ice.scheme.ρ_l)
                state = CM.P3Scheme.P3State(mp_ref.ice.scheme, L_ice, N_ice, F_rim, ρ_rim)
                CM.P3Scheme.get_distribution_logλ(state)
            else
                FT(0)
            end

            ref = BMT.bulk_microphysics_tendencies(
                BMT.Microphysics2Moment(), mp_ref, tps,
                s.ρ, s.T, s.q_tot, s.q_lcl, s.n_lcl, s.q_rai, s.n_rai,
                s.q_ice, s.n_ice, s.q_rim, s.b_rim, logλ,
            )
            for (n, tol) in orders_and_tol
                mp_n = make_mp(n)
                t_n = BMT.bulk_microphysics_tendencies(
                    BMT.Microphysics2Moment(), mp_n, tps,
                    s.ρ, s.T, s.q_tot, s.q_lcl, s.n_lcl, s.q_rai, s.n_rai,
                    s.q_ice, s.n_ice, s.q_rim, s.b_rim, logλ,
                )
                fails = compare_with_reference(ref, t_n, tol; report = true)
                if fails != 0
                    @info "quadrature sweep: state=$idx n=$n fails=$fails tol=$tol"
                end
                @test fails == 0
            end
        end
    end
end

@testset "BMT quadrature-order kwarg (Float64)" begin
    test_quadrature_order_sweep(Float64)
end
nothing
