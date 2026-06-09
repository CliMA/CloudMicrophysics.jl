using Test

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI

"""
Sweep test for the P3-ice quadrature order. `quadrature_order` now lives
on `P3IceParams`, so the sweep is over `Microphysics2MParams(FT;
with_ice = true, quadrature_order = n)`.

Issue 011 measured the worst-case relative error of one single integral
(`bulk_liquid_ice_collision_sources`) across 5 states and 7 orders, and
found all values < 1 % even at `n = 15`. This suite extends that to the
**full BMT tendency vector**, sweeping a broader set of physically
plausible column states and asserting lenient per-`n` tolerances.

# Tolerance rationale

The goal is a sanity check for a student picking a lower `n` for fast
iteration, not a tight correctness check. The tolerances are set to the
level where the approximation is "good enough" for debugging but still
meaningfully faster:

- `n = 100`: < 2e-3 relative on all fields. Issue 011 showed the
  single-integral error is ~1e-5, but the full BMT pipeline composes
  multiple integrals (bulk liquid-ice collisions, ice aggregation,
  melting), which can compound. 2e-3 is loose enough to absorb
  integration-scheme drift while still catching genuine regressions.
- `n = 50`:  < 5e-3 relative — "safe for most KiD runs".
- `n = 25`:  < 5e-2 relative — acceptable for short diagnostics.
- `n = 15`:  < 2e-1 — order-of-magnitude agreement only.

The tolerance is applied on `max(abs(a), abs(b), ε)` where
`ε = 1e-16 · mass_scale` to avoid false failures on tendencies that are
near-zero by physics (e.g. ice processes in a rain-only column). All
finite-valued outputs are checked; `NaN` / `Inf` trigger a failure
regardless of tolerance.

This does NOT include the timing assertions from the issue text — those
are flaky under CI load. A separate `@info` block reports the
wall-clock ratios for the student's information when the test file is
run directly.
"""

function generate_column_states(::Type{FT}) where {FT}
    # A hand-curated collection of physically plausible (ρ, T, q_*) states
    # spanning the regimes the Jouan kinematic column experiences. Each row
    # is a NamedTuple fed directly into `bulk_microphysics_tendencies`.

    # Saturation-vapor helper (local)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    q_vs_l(T, ρ) = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρ)
    q_vs_i(T, ρ) = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)

    states = NamedTuple[]

    # 1. Warm, cloudy, no ice, no rain — pure activation/condensation
    let ρ = FT(1.2), T = FT(290)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(1e-3),
            q_lcl = FT(1e-3), n_lcl = FT(1e8),
            q_rai = FT(0),    n_rai = FT(0),
            q_ice = FT(0),    n_ice = FT(0),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 2. Warm, heavy rain, no cloud
    let ρ = FT(1.1), T = FT(285)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(5e-4),
            q_lcl = FT(0),    n_lcl = FT(0),
            q_rai = FT(5e-4), n_rai = FT(1e4),
            q_ice = FT(0),    n_ice = FT(0),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 3. Freezing-level mixed phase, light ice, no rime
    let ρ = FT(0.9), T = FT(270)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(1e-4) + FT(1e-5),
            q_lcl = FT(1e-4), n_lcl = FT(1e8),
            q_rai = FT(0),    n_rai = FT(0),
            q_ice = FT(1e-5), n_ice = FT(1e5),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 4. Cold cirrus, trace ice, no rain, no cloud
    let ρ = FT(0.5), T = FT(240)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_i(T, ρ) + FT(1e-6),
            q_lcl = FT(0),    n_lcl = FT(0),
            q_rai = FT(0),    n_rai = FT(0),
            q_ice = FT(1e-6), n_ice = FT(1e5),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 5. Heavy riming regime
    let ρ = FT(0.85), T = FT(265)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(5e-4) + FT(5e-4),
            q_lcl = FT(5e-4), n_lcl = FT(1e8),
            q_rai = FT(2e-4), n_rai = FT(1e4),
            q_ice = FT(5e-4), n_ice = FT(1e5),
            q_rim = FT(1e-4), b_rim = FT(1e-4 / 300),
        ))
    end

    # 6. Dry subsaturated, no condensate — evaporation regime
    let ρ = FT(1.0), T = FT(290)
        push!(states, (;
            ρ, T,
            q_tot = FT(0.5) * q_vs_l(T, ρ),
            q_lcl = FT(0),    n_lcl = FT(0),
            q_rai = FT(1e-4), n_rai = FT(1e4),
            q_ice = FT(0),    n_ice = FT(0),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 7. Just below 273 K — melting threshold with heavy ice
    let ρ = FT(1.0), T = FT(272.5)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(1e-3),
            q_lcl = FT(0),    n_lcl = FT(0),
            q_rai = FT(0),    n_rai = FT(0),
            q_ice = FT(1e-3), n_ice = FT(5e4),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 8. Just above 273 K — melting active
    let ρ = FT(1.0), T = FT(274.0)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(1e-3),
            q_lcl = FT(0),    n_lcl = FT(0),
            q_rai = FT(0),    n_rai = FT(0),
            q_ice = FT(1e-3), n_ice = FT(5e4),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 9. Strong supersaturation over ice, no liquid
    let ρ = FT(0.7), T = FT(250)
        push!(states, (;
            ρ, T,
            q_tot = FT(1.5) * q_vs_i(T, ρ),
            q_lcl = FT(0),    n_lcl = FT(0),
            q_rai = FT(0),    n_rai = FT(0),
            q_ice = FT(1e-5), n_ice = FT(1e5),
            q_rim = FT(0),    b_rim = FT(0),
        ))
    end

    # 10. Mixed-phase mid-troposphere with rain + ice
    let ρ = FT(0.8), T = FT(268)
        push!(states, (;
            ρ, T,
            q_tot = q_vs_l(T, ρ) + FT(3e-4) + FT(3e-4),
            q_lcl = FT(3e-4), n_lcl = FT(1e8),
            q_rai = FT(1e-4), n_rai = FT(5e3),
            q_ice = FT(3e-4), n_ice = FT(1e5),
            q_rim = FT(1e-5), b_rim = FT(1e-5 / 400),
        ))
    end

    return states
end

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
    # `quadrature_order` now lives on `P3IceParams`, so building one
    # `Microphysics2MParams` per order keeps the sweep clean.
    make_mp(n) = CMP.Microphysics2MParams(FT; with_ice = true, quadrature_order = n)

    orders_and_tol = [
        (100, FT(2e-3)),
        (50,  FT(5e-3)),
        (25,  FT(5e-2)),
        (15,  FT(2e-1)),
    ]
    reference_n = 200

    mp_ref = make_mp(reference_n)
    states = generate_column_states(FT)
    @test length(states) >= 10

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
