"""
Quadrature error study for the P3-ice integrals.

Measure the relative error of every quadrature-consuming P3 quantity against a
high-order Gauss-Legendre reference, across a set of physically plausible
column states, and attribute it to one of five levels:

- `transport`: the sedimentation velocities and melt rate.
- `collision_efficiency`: the pooled components proportional to the collision
  efficiency (`E = 1` assumed): the liquid-ice collision sources and ice
  self-collection. Their parameterization uncertainty exceeds the quadrature
  error, so they tolerate a coarser rule than the transport components.
- `collision`: the liquid-ice collision sources alone, a subset of
  `collision_efficiency` reported separately since it dominates the
  runtime cost of the two.
- `selfcol`: ice self-collection alone, the other `collision_efficiency`
  component.
- `bulk`: the full `bulk_microphysics_tendencies` vector.

The states are `generate_column_states` plus graupel/hail cores
(`hail_core_states`), whose large mean particle sizes stress the
size-distribution tail. The error metric is `|a - b| / max(|a|, |b|, floor)`,
with a per-component floor of `1e-9` times the largest reference magnitude
across states, so components that are zero by physics do not register as
relative error.

Run from the repository root:

    julia --project=test -e 'include("test/p3_quadrature_error_study.jl");
                             run_quadrature_error_study()'

Keyword arguments of `run_quadrature_error_study`:

- `FT`: float type. By default, `Float64`.
- `orders`: Gauss-Legendre orders to evaluate. By default, `(5, 6, 7, 8, 10, 12)`.
- `reference_order`: order of the reference rule. By default, `128`.
- `include_timing`: also benchmark `bulk_microphysics_tendencies` per order on
  a hail-core state. By default, `false`.

Return a vector of `(; order, level, median, p95, max)` rows.

The default-order regression test in `test/bulk_tendencies_quadrature_tests.jl`
locks in the outcome of this study; rerun the study when changing the
quadrature rule, the integral bounds or breakpoints, or the P3 process
integrands, and revisit the default order recorded in
`src/parameters/Microphysics2MParams.jl`. See #741 for the study behind the
current default.
"""

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
import CloudMicrophysics.ThermodynamicsInterface as TDI
import BenchmarkTools as BT
import Statistics: median, quantile

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
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(1e-3),
                q_lcl = FT(1e-3), n_lcl = FT(1e8),
                q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(0), n_ice = FT(0),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 2. Warm, heavy rain, no cloud
    let ρ = FT(1.1), T = FT(285)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(5e-4),
                q_lcl = FT(0), n_lcl = FT(0),
                q_rai = FT(5e-4), n_rai = FT(1e4),
                q_ice = FT(0), n_ice = FT(0),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 3. Freezing-level mixed phase, light ice, no rime
    let ρ = FT(0.9), T = FT(270)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(1e-4) + FT(1e-5),
                q_lcl = FT(1e-4), n_lcl = FT(1e8),
                q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(1e-5), n_ice = FT(1e5),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 4. Cold cirrus, trace ice, no rain, no cloud
    let ρ = FT(0.5), T = FT(240)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_i(T, ρ) + FT(1e-6),
                q_lcl = FT(0), n_lcl = FT(0),
                q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(1e-6), n_ice = FT(1e5),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 5. Heavy riming regime
    let ρ = FT(0.85), T = FT(265)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(5e-4) + FT(5e-4),
                q_lcl = FT(5e-4), n_lcl = FT(1e8),
                q_rai = FT(2e-4), n_rai = FT(1e4),
                q_ice = FT(5e-4), n_ice = FT(1e5),
                q_rim = FT(1e-4), b_rim = FT(1e-4 / 300),
            ),
        )
    end

    # 6. Dry subsaturated, no condensate — evaporation regime
    let ρ = FT(1.0), T = FT(290)
        push!(
            states,
            (;
                ρ, T,
                q_tot = FT(0.5) * q_vs_l(T, ρ),
                q_lcl = FT(0), n_lcl = FT(0),
                q_rai = FT(1e-4), n_rai = FT(1e4),
                q_ice = FT(0), n_ice = FT(0),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 7. Just below 273 K — melting threshold with heavy ice
    let ρ = FT(1.0), T = FT(272.5)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(1e-3),
                q_lcl = FT(0), n_lcl = FT(0),
                q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(1e-3), n_ice = FT(5e4),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 8. Just above 273 K — melting active
    let ρ = FT(1.0), T = FT(274.0)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(1e-3),
                q_lcl = FT(0), n_lcl = FT(0),
                q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(1e-3), n_ice = FT(5e4),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 9. Strong supersaturation over ice, no liquid
    let ρ = FT(0.7), T = FT(250)
        push!(
            states,
            (;
                ρ, T,
                q_tot = FT(1.5) * q_vs_i(T, ρ),
                q_lcl = FT(0), n_lcl = FT(0),
                q_rai = FT(0), n_rai = FT(0),
                q_ice = FT(1e-5), n_ice = FT(1e5),
                q_rim = FT(0), b_rim = FT(0),
            ),
        )
    end

    # 10. Mixed-phase mid-troposphere with rain + ice
    let ρ = FT(0.8), T = FT(268)
        push!(
            states,
            (;
                ρ, T,
                q_tot = q_vs_l(T, ρ) + FT(3e-4) + FT(3e-4),
                q_lcl = FT(3e-4), n_lcl = FT(1e8),
                q_rai = FT(1e-4), n_rai = FT(5e3),
                q_ice = FT(3e-4), n_ice = FT(1e5),
                q_rim = FT(1e-5), b_rim = FT(1e-5 / 400),
            ),
        )
    end

    return states
end

"""
    hail_core_states(FT, specs)

Mixed-phase states with large mean particle sizes, one per `(F_rim, ρ_rim,
n_ice)` row of `specs`. The size-distribution tail of these states spans
several decades of particle diameter.
"""
function hail_core_states(::Type{FT}, specs) where {FT}
    states = NamedTuple[]
    for (F_rim, ρ_rim, n_ice) in specs
        q_ice = 1e-3
        q_rim = F_rim * q_ice
        b_rim = q_rim / ρ_rim
        push!(
            states,
            (;
                ρ = FT(0.9), T = FT(262.0), q_tot = FT(4e-3),
                q_lcl = FT(2e-4), n_lcl = FT(5e7),
                q_rai = FT(2e-4), n_rai = FT(2e4),
                q_ice = FT(q_ice), n_ice = FT(n_ice),
                q_rim = FT(q_rim), b_rim = FT(b_rim),
            ),
        )
    end
    return states
end

const STUDY_HAIL_CORES = (
    (0.5, 400.0, 1e5), (0.95, 400.0, 1e4), (0.95, 800.0, 1e4),
    (0.8, 600.0, 1e3), (0.8, 600.0, 100.0), (0.8, 600.0, 20.0),
    (0.95, 800.0, 50.0),
)

"""
    evaluate_quadrature_levels(mp, tps, s)

Evaluate every quadrature-consuming P3 quantity at state `s` with the
quadrature rule stored in `mp.ice`. Return `(out, f_shd)`: `out` maps each
RATE quantity to its component vector as before (ice quantities absent for
ice-free states); `f_shd` (card #17, `nothing` where ice is absent) is
`bulk_liquid_ice_collision_sources`'s bulk wet-growth fraction, tracked
SEPARATELY from `out` rather than folded into `out[:collision]`.

`f_shd` is a DIAGNOSTIC FRACTION, not a rate, and its own numerator
(the collection rate over the wet set) is already on this study's exclusion list -
"the wet-growth-gated shed components ... `p95`/`max` inherit the onset-search
limitation ... and are not useful regression signals". Pooling it into
`out[:collision]` would contaminate the LOCKED `collision`/`collision_efficiency`
p95 with that same onset-search noise (measured: two real regression-test
failures, `p95` 0.036/0.041 against a 0.03 lock, the two levels built from
`out[:collision]`, while every one of the twelve original rate fields stayed
bit-identical - the lock itself was never wrong, only what got pooled into it).
Tracked here instead of dropped, so it is reported, not silently uncovered;
see `run_quadrature_error_study`'s own `f_shd` block for where it prints.
"""
function evaluate_quadrature_levels(mp, tps, s)
    (; quad, terminal_velocity, cloud_pdf, rain_pdf) = mp.ice
    aps = mp.warm_rain.air_properties
    (; ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = s
    FT = typeof(ρ)
    out = Dict{Symbol, Vector{FT}}()
    has_ice = q_ice > 0 && n_ice > 0
    state = has_ice ?
            P3.state_from_prognostic(mp.ice.scheme, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim) :
            nothing
    logλ = has_ice ? P3.get_distribution_logλ(state) : FT(0)
    t = BMT.bulk_microphysics_tendencies(
        BMT.Microphysics2Moment(), mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    )
    out[:bulk] = collect(FT, values(t))
    f_shd = nothing
    if has_ice
        coll = P3.bulk_liquid_ice_collision_sources(
            state, logλ, cloud_pdf, rain_pdf, ρ * q_lcl, n_lcl, ρ * q_rai, n_rai,
            aps, tps, terminal_velocity, ρ, T; B_rim = ρ * b_rim, quad,
        )
        # Explicit twelve-field list (not `values(coll)`) so the pooled rate-error vector this
        # lock depends on cannot silently change shape again if another diagnostic field is
        # ever added to `bulk_liquid_ice_collision_sources`'s return value.
        out[:collision] = FT[
            coll.∂ₜq_c, coll.∂ₜq_r, coll.∂ₜN_c, coll.∂ₜN_r, coll.∂ₜL_rim, coll.∂ₜL_ice,
            coll.∂ₜB_rim, coll.∂ₜq_c_frz, coll.∂ₜq_c_shd, coll.∂ₜq_r_frz, coll.∂ₜB_rim_c,
            coll.∂ₜB_rim_r,
        ]
        f_shd = coll.f_shd
        sc = P3.ice_self_collection(state, logλ, terminal_velocity, ρ; quad)
        out[:selfcol] = [sc.dNdt]
        out[:vN] = [P3.ice_terminal_velocity_number_weighted(terminal_velocity, ρ, state, logλ; quad)]
        out[:vM] = [P3.ice_terminal_velocity_mass_weighted(terminal_velocity, ρ, state, logλ; quad)]
        melt = P3.ice_melt(terminal_velocity, aps, tps, T, ρ, state, logλ; quad)
        out[:melt] = collect(FT, values(melt))
    end
    return out, f_shd
end

error_study_level(k) =
    k == :bulk ? :bulk :
    (k in (:collision, :selfcol) ? :collision_efficiency : :transport)

function run_quadrature_error_study(;
    FT = Float64,
    orders = (5, 6, 7, 8, 10, 12),
    reference_order = 128,
    include_timing = false,
    states = vcat(generate_column_states(FT), hail_core_states(FT, STUDY_HAIL_CORES)),
)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    make_mp(n) =
        CMP.Microphysics2MParams(FT; with_ice = true, quad = CM.Quadrature.GaussLegendre(FT, n))

    ref_pairs = [evaluate_quadrature_levels(make_mp(reference_order), tps, s) for s in states]
    refs = first.(ref_pairs)
    f_shd_refs = last.(ref_pairs)
    floors = Dict(
        k =>
            FT(1e-9) .* max.(
                [
                    maximum(abs(r[k][i]) for r in refs if haskey(r, k); init = FT(0)) for
                    i in 1:maximum(length(r[k]) for r in refs if haskey(r, k))
                ],
                eps(FT),
            ) for k in (:bulk, :vN, :vM, :melt, :selfcol, :collision)
    )
    # `f_shd` (card #17) gets its own floor, on the SAME `1e-9 x max-reference-magnitude` recipe
    # as everything else - kept separate from `floors` above because it is not one of the keys
    # `out` carries (see `evaluate_quadrature_levels`'s docstring for why it is tracked apart).
    f_shd_floor = FT(1e-9) * max(maximum(abs, filter(!isnothing, f_shd_refs); init = FT(0)), eps(FT))

    results = NamedTuple[]
    println("order | level | median | p95 | max")
    for n in orders
        mp = make_mp(n)
        errs = Dict(
            :bulk => FT[], :collision_efficiency => FT[], :transport => FT[],
            :collision => FT[], :selfcol => FT[],
        )
        worst = Tuple{FT, Int, Int}[]
        f_shd_errs = FT[]
        for (si, (s, r, f_shd_ref)) in enumerate(zip(states, refs, f_shd_refs))
            e, f_shd_e = evaluate_quadrature_levels(mp, tps, s)
            for k in keys(r)
                for (i, (a, b, f)) in enumerate(zip(e[k], r[k], floors[k]))
                    err = abs(a - b) / max(abs(a), abs(b), f)
                    push!(errs[error_study_level(k)], err)
                    k in (:collision, :selfcol) && push!(errs[k], err)
                    k == :bulk && push!(worst, (err, si, i))
                end
            end
            if f_shd_ref !== nothing && f_shd_e !== nothing
                push!(f_shd_errs, abs(f_shd_e - f_shd_ref) / max(abs(f_shd_e), abs(f_shd_ref), f_shd_floor))
            end
        end
        for level in (:transport, :collision_efficiency, :collision, :selfcol, :bulk)
            v = errs[level]
            row = (;
                order = n, level,
                median = median(v), p95 = quantile(v, 0.95), max = maximum(v),
            )
            push!(results, row)
            println(rpad("GL($n)", 7), " | ", rpad(level, 20), " | ",
                rpad(round(row.median, sigdigits = 3), 10), " | ",
                rpad(round(row.p95, sigdigits = 3), 10), " | ",
                round(row.max, sigdigits = 3))
        end
        sort!(worst; rev = true)
        println("  worst bulk components (err, state, component): ",
            join(["($(round(e, sigdigits = 3)), $si, $i)" for (e, si, i) in worst[1:5]], " "))
        # REPORTED, NOT LOCKED (card #17): f_shd inherits the wet-growth onset-search
        # sensitivity this study's own docstring already documents for the removed indicator and its
        # siblings - printed so its convergence behaviour stays visible rather than dropped, but
        # deliberately outside `results`/`errs`, so no p95/median lock applies to it here.
        if !isempty(f_shd_errs)
            println("  f_shd (reported, not locked - onset-search-limited): median=",
                round(median(f_shd_errs), sigdigits = 3), " p95=",
                round(quantile(f_shd_errs, 0.95), sigdigits = 3), " max=",
                round(maximum(f_shd_errs), sigdigits = 3), " n=", length(f_shd_errs))
        end
    end

    if include_timing
        s = last(states)
        println("\norder | bulk tendency min time (μs)")
        for n in (orders..., reference_order)
            mp = make_mp(n)
            (; ρ, T, q_tot, q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = s
            st = P3.state_from_prognostic(mp.ice.scheme, ρ * q_ice, ρ * n_ice, ρ * q_rim, ρ * b_rim)
            logλ = P3.get_distribution_logλ(st)
            b = BT.@benchmark $(BMT.bulk_microphysics_tendencies)(
                $(BMT.Microphysics2Moment()), $mp, $tps, $ρ, $T, $q_tot,
                $q_lcl, $n_lcl, $q_rai, $n_rai, $q_ice, $n_ice, $q_rim, $b_rim, $logλ,
            ) samples = 100 evals = 5
            println(rpad("GL($n)", 7), " | ", round(BT.minimum(b).time / 1000, digits = 1))
        end
    end

    return results
end
