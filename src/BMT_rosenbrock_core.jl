#####
##### Rosenbrock-Euler substepping: the scheme-agnostic core
#####

# The substep driver below is generic over the length `N` of the substep state. Every
# quantity it forms is built from the state vector, the tendency, the Jacobian and the
# species mask, so the same driver serves the eight-species 2M+P3 state, the
# temperature-coupled nine-vector, and the four-species 1M state. What is specific to a
# state - which species may be masked, how its condensate splits by phase, and how its
# components are floored - reaches the driver through the hooks declared here, whose
# methods live with their state types.

"""
    AIR_DENSITY_FLOOR

Lower bound applied to the air density before it is used in a microphysics tendency, as
`max(ρ, AIR_DENSITY_FLOOR)`.

The host can hand the microphysics a non-positive grid-mean air density transiently, and
the scheme divides by it and takes logs and negative fractional powers of it. A
non-negative clamp is *not* enough there: at exactly zero the condensation timescale
returns `NaN` and the whole manual Jacobian follows, which merely trades a `DomainError`
for a non-finite Jacobian and sends the substep down its explicit fallback. A positive
floor keeps the rates and the Jacobian finite, with the condensation timescale of order
`1e5` s, that is, no condensation, which is the right degenerate behavior.

Zero remains the correct floor for the *masses* and *numbers*, which vanish physically;
only the density needs a positive one. The same value and reasoning apply to the Chen
2022 velocity coefficients in `Common.jl`.

The value is far below any atmospheric density (a box minimum is of order `1e-2`), so the
floor is inert on physical states and the reference trajectory does not move.
"""
const AIR_DENSITY_FLOOR = 1e-4

@inline _floored_air_density(ρ) = max(ρ, oftype(ρ, AIR_DENSITY_FLOOR))

# ---- State hooks: what the driver needs from a substep state ----

"""
    _rosenbrock_species_mask(x)

Diagonal of the species projection matrix `P` used by [`_rosenbrock_update`](@ref): 1 for
an active species, 0 for a near-empty one. A masked species takes the forward-Euler
update while the active species stay implicit.

Generic declaration only. Each substep state type supplies its own method, because which
species count as near-empty, and which are never masked, is a property of the state rather
than of the solver.

Reached only by the [`ExactJacobian`](@ref) and [`ImplicitGrowth`](@ref) combination;
every other combination resolves to [`_full_species_mask`](@ref) through
[`_species_mask`](@ref).
"""
function _rosenbrock_species_mask end

"""
    _full_species_mask(x)

The all-ones species projection: every species stays in the implicit solve.
"""
@inline _full_species_mask(x::SA.StaticVector{N, FT}) where {N, FT} =
    ones(SA.SVector{N, FT})

"""
    _species_mask(jacobian, growth)

The species projection `x -> z` for a [`Jacobian`](@ref) and [`GrowthTreatment`](@ref)
pair. The donor-based matrices are bounded by their rate flooring, so every species stays
implicit ([`_full_species_mask`](@ref)). An explicit growth diagonal removes the unbounded
growth of the exact Jacobian, so its species stay implicit too; the exact Jacobian with
implicit growth uses the near-empty mask ([`_rosenbrock_species_mask`](@ref)). The manual
Jacobians drop the unbounded quadrature growth couplings, so they stay on the full mask
under either growth treatment.
"""
@inline _species_mask(::DonorJacobian, ::GrowthTreatment) = _full_species_mask
@inline _species_mask(::CoupledDonorJacobian, ::GrowthTreatment) = _full_species_mask
@inline _species_mask(::ExactJacobian, ::ExplicitGrowthDiagonal) = _full_species_mask
@inline _species_mask(::ExactJacobian, ::ImplicitGrowth) = _rosenbrock_species_mask
@inline _species_mask(::ManualJacobian, ::GrowthTreatment) = _full_species_mask
@inline _species_mask(::TemperatureCoupledJacobian, ::GrowthTreatment) = _full_species_mask

"""
    _condensate_phases(x)

The condensed water of a substep state split by phase, as `(q_liq, q_ice)` in kg/kg:
everything that freezes on cooling first, then everything already frozen.

Generic declaration only; each substep state type supplies its own method, since the split
runs over that state's own species. Rime mass is a *fraction* of the ice content rather
than an addition to it, so a state carrying it does not add it here; including it would
double-count the rimed part of the ice mass.

The split is the whole of what the water bound ([`_condensate_total`](@ref)) and the
end-state saturation limiter ([`_apply_limiter`](@ref)) need from a state, so supplying
this one method gives a state both.
"""
function _condensate_phases end

"""
    _per_process_rates(g, x)

The per-process breakdown of one substep's tendency evaluation `g(x)`, in the state's own
space, for a [`VerboseSink`](@ref) to attribute the increment against.

Generic declaration only; each substep tendency callable supplies its own method, because
the breakdown comes from that callable's own process-rate function and, for a state that
carries a temperature component, that state's own latent-heating map.

Reached only through [`_recorded_processes`](@ref), so it costs a second evaluation of the
process rates beyond whatever [`_tendency_and_jacobian`](@ref) already needed, and only on
a [`VerboseSink`](@ref) run - never on the production path, where `_recorded_processes`
dispatches the request away without calling this at all. This is what keeps the attribution
free of any knowledge of the [`Jacobian`](@ref) option: [`ExactJacobian`](@ref)'s method
never forms `pp` at all, and [`ManualJacobian`](@ref)'s and
[`TemperatureCoupledJacobian`](@ref)'s discard the `pp` their own Jacobian construction
already computed, so this is one place either mode reaches for it, keyed on the tendency
callable rather than on the Jacobian.
"""
function _per_process_rates end

"""
    _condensate_total(x)

Condensed-water content of a substep state, the sum of its two
[`_condensate_phases`](@ref).
"""
@inline function _condensate_total(x)
    q_liq, q_ice = _condensate_phases(x)
    return q_liq + q_ice
end

"""
    _apply_positivity_floor(x, Δx, ρ_min, ρ_max)

The accepted end state of a substep increment `Δx` taken from state `x`: the
component-wise positivity floor `max.(x .+ Δx, 0)`.

This is the generic method, for a state whose components are independent. A state whose
components are not supplies its own, so the treatment follows from the state type rather
than from a mode switch: the 2M+P3 state projects its rime mass/volume pair onto the
admissible density cone `[ρ_min, ρ_max]` instead of flooring the two independently, and a
state carrying a temperature component exempts that component from the floor.

`ρ_min` and `ρ_max` are the rime density bounds. Only a state that projects a rime pair
reads them, so a caller without one passes `nothing` and the bounds are inert by
construction.
"""
@inline _apply_positivity_floor(x, Δx, ρ_min, ρ_max) = max.(x .+ Δx, 0)

# ---- Jacobian assembly ----

"""
    _apply_growth(growth, J)

Apply a [`GrowthTreatment`](@ref) to the Jacobian `J`. [`ImplicitGrowth`](@ref) returns
`J` unchanged; [`ExplicitGrowthDiagonal`](@ref) removes the positive diagonal entries,
leaving the off-diagonals and the negative diagonals.
"""
@inline _apply_growth(::ImplicitGrowth, J) = J
@inline function _apply_growth(::ExplicitGrowthDiagonal, J::SA.SMatrix{N, N, FT}) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    return J - max.(Iₙ .* J, zero(FT))
end

"""
    _tendency_and_jacobian(jacobian, g, x)

The raw substep tendency `f = g(x)` and the substep Jacobian before the growth treatment,
for a [`Jacobian`](@ref) option, returned as `(f, J)`.

For [`ExactJacobian`](@ref) the primal `f` and the `N×N` Jacobian both come out of the one
`ForwardDiff` pass. A hand-built matrix produces no tendency by-product of its own, so its
method evaluates `f = g(x)` separately; those methods live with the states they are built
for.
"""
@inline function _tendency_and_jacobian(::ExactJacobian, g, x::SA.FieldVector{N, FT}) where {N, FT}
    Tag = typeof(FD.Tag(g, FT))
    dx = SA.SVector(
        ntuple(i -> FD.Dual{Tag}(x[i], ntuple(s -> ifelse(s == i, one(FT), zero(FT)), Val(N))...), Val(N)),
    )
    y = g(dx)
    f = typeof(x)(ntuple(i -> @inbounds(FD.value(y[i])), Val(N))...)
    J = SA.SMatrix{N, N, FT}(
        ntuple(k -> @inbounds(FD.partials(y[(k - 1) % N + 1], (k - 1) ÷ N + 1)), Val(N * N)),
    )
    return f, J
end

# ---- The substep increment ----

"""
    _euler_update(x, f, h)

Forward-Euler substep, floored at zero.
"""
@inline _euler_update(x, f, h) = max.(x .+ h .* f, 0)

"""
    _water_bounded_increment_diag(x, d, q_tot)

[`_water_bounded_increment`](@ref)'s full computation, returning `(σ .* d, σ)`: the
scaled increment and the scale factor `σ` itself, the water bound's contribution to a
[`VerboseSink`](@ref)'s `α_w`. `_water_bounded_increment` is `first ∘
_water_bounded_increment_diag`, so recording `σ` costs no second pass over the condensate
totals.

The increment `d` is scaled so the condensate it produces stays within the cell's total
water `q_tot`, which the vapor budget cannot exceed.

Used on the substep's two explicit fallback branches. Both span the whole substep width, so
a stiff phase-change rate can convert far more mass than the cell holds, and the resulting
latent release then carries the substep temperature outside the thermodynamic domain - the
fallback produces states worse than the ones it was reached to rescue. The condensate total
is linear in the increment, so the scale factor is exact rather than iterative, and it is
one whenever the step is already admissible.
"""
@inline function _water_bounded_increment_diag(x::SA.StaticVector{N, FT}, d, q_tot) where {N, FT}
    c₀ = _condensate_total(x)
    c₁ = _condensate_total(x .+ d)
    over = (c₁ > FT(q_tot)) & (c₁ > c₀)
    σ = ifelse(over, clamp((FT(q_tot) - c₀) / (c₁ - c₀), zero(FT), one(FT)), one(FT))
    return σ .* d, σ
end

"""
    _water_bounded_increment(x, d, q_tot)

The scaled increment of [`_water_bounded_increment_diag`](@ref), without its scale factor.
"""
@inline _water_bounded_increment(x, d, q_tot) = first(_water_bounded_increment_diag(x, d, q_tot))

"""
    _bounded_explicit_step_diag(x, d, q_tot)

[`_bounded_explicit_step`](@ref)'s full computation, returning `(d′, α_w)`. No bound where
the caller supplies no water budget - the 1M substep and the temperature-coupled entry keep
the bare explicit step, at `α_w = 1` - otherwise [`_water_bounded_increment_diag`](@ref).
"""
@inline _bounded_explicit_step_diag(x, d, ::Nothing) = (d, one(eltype(d)))
@inline _bounded_explicit_step_diag(x, d, q_tot) = _water_bounded_increment_diag(x, d, q_tot)

"""
    _bounded_explicit_step(x, d, q_tot)

The scaled increment of [`_bounded_explicit_step_diag`](@ref), without its scale factor.
"""
@inline _bounded_explicit_step(x, d, q_tot) = first(_bounded_explicit_step_diag(x, d, q_tot))

"""
    _rosenbrock_system(x, f, J, z, h)

Build the equilibrated linear system of one linearized-implicit (Rosenbrock-Euler) substep
at state `x` with raw tendency `f`, Jacobian `J`, species mask `z`, and substep `h`.
Returns `(S, S⁻¹, A)`, the equilibration matrix `S = Diagonal(|x| + h |f| + ϵ)`, its
inverse, and the equilibrated system matrix `A = I/h - S⁻¹ B S`, where `B = P J P` is the
masked Jacobian and `P = Diagonal(z)` is the species projection built from the per-scheme
species mask `z` (e.g. [`_rosenbrock_species_mask`](@ref)).

`P` and `S` are diagonal, so `P J P` and `S⁻¹ B S` are row and column scalings of `J` and
are formed as such rather than as matrix products.

The system build is separated from the solve so the full-step update and the per-process
attribution can reuse one factorization, solving against the same `S`, `S⁻¹`, `A`.
"""
@inline function _rosenbrock_system(
    x::SA.StaticVector{N, FT}, f, J, z, h,
) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    s = abs.(x) .+ h .* abs.(f) .+ eps(FT)
    s⁻¹ = inv.(s)
    S = Iₙ .* s'
    S⁻¹ = Iₙ .* s⁻¹'
    B = z .* J .* z'
    A = Iₙ / h - s⁻¹ .* B .* s'
    return S, S⁻¹, A
end

"""
    _rosenbrock_solve(S, S⁻¹, A, v)

Solve the equilibrated Rosenbrock system from [`_rosenbrock_system`](@ref) for the
unclamped increment of right-hand side `v`: `Δ = S (A \\ (S⁻¹ v))`, the equilibrated form
of `(I/h - P J P)⁻¹ v`. Linear in `v`, so per-process increments sum to the full-step
increment.
"""
@inline _rosenbrock_solve(S, S⁻¹, A, v) = S * (A \ (S⁻¹ * v))

# Acceptance bound for the linearized-implicit increment relative to the
# explicit-step scale. A well-conditioned solve of `(I/h - B) Δx = f` satisfies
# `‖Δx‖∞ = O(h ‖f‖∞)`; increments far beyond that scale indicate a
# near-singular system matrix.
const ROSENBROCK_INCREMENT_LIMIT = 10

"""
    _solve_increment_acceptable(d, f, h)

Whether the Rosenbrock increment `d` is consistent with a well-conditioned solve,
`‖d‖∞ ≤ $(ROSENBROCK_INCREMENT_LIMIT) h ‖f‖∞`, with both vectors in the equilibrated units
of [`_rosenbrock_system`](@ref) (`S⁻¹ Δx` and `S⁻¹ f`), so the bound is relative to each
component's own scale. A rejected increment is replaced with the explicit update `h f`.
Non-finite entries in `d` fail the bound.
"""
@inline function _solve_increment_acceptable(d, f, h)
    FT = eltype(d)
    return maximum(abs, Tuple(d)) <= FT(ROSENBROCK_INCREMENT_LIMIT) * h * maximum(abs, Tuple(f))
end

"""
    _rosenbrock_update_diag(x, f, J, z, h, q_tot = nothing,
        ρ_min = nothing, ρ_max = nothing)

The substep's full computation, returning `(x_new, A, W, α_w)`:

- `A`: the equilibrated system matrix [`_rosenbrock_system`](@ref) built and
  [`_rosenbrock_solve`](@ref) inverted for this update, valid whenever this function is
  reached (see [`_rosenbrock_substep_diag`](@ref) for when that is).
- `W`: the [`SubstepOperator`](@ref) whose solve actually produced the accepted increment -
  [`EquilibratedSolve`](@ref) wrapping this same system when the solved increment passed
  [`_solve_increment_acceptable`](@ref), [`ExplicitStep`](@ref) when it did not and the
  explicit `h f` replaced it. `A` is still returned in the second case (a near-singular
  system is still the system that was built), but `W` is not `EquilibratedSolve` there,
  because production did not solve against it for the increment it accepted; a sink that
  attributed against `A` regardless would be attributing an increment the substep never
  took.
- `α_w`: the water bound's rescale of whichever increment `W` names, from
  [`_bounded_explicit_step_diag`](@ref).

Production takes `first` of this, so recording `A`, `W` and `α_w` for a diagnostic caller
costs no second call to `_rosenbrock_system`: production and any diagnostic run the
identical single build-and-solve.

A caller supplies the rime density bounds explicitly. The rime-pair projection inside
[`_apply_positivity_floor`](@ref) divides by them, so a convenience wrapper that defaulted
them to `nothing` would be unusable for any state carrying rime.
"""
@inline function _rosenbrock_update_diag(
    x::SA.StaticVector{N, FT}, f, J, z, h, q_tot = nothing,
    ρ_min = nothing, ρ_max = nothing,
) where {N, FT}
    S, S⁻¹, A = _rosenbrock_system(x, f, J, z, h)
    Δx = _rosenbrock_solve(S, S⁻¹, A, f)
    solved = _solve_increment_acceptable(S⁻¹ * Δx, S⁻¹ * f, h)
    Δx_raw = solved ? Δx : h .* f
    # Bound BOTH branches, not only the rejected one. The accepted increment is the path
    # the damage takes: `_rosenbrock_system` equilibrates, so a near-singular system's
    # diagonal spike scales out of `_solve_increment_acceptable` and the huge increment is
    # ACCEPTED, minting condensate orders beyond the cell's total water. The positivity
    # floor below bounds sign and not magnitude, and the saturation limiter bisects on
    # supersaturation rather than on the water budget, so nothing else catches it.
    Δx_bounded, α_w = _bounded_explicit_step_diag(x, Δx_raw, q_tot)
    x_new = _apply_positivity_floor(x, Δx_bounded, ρ_min, ρ_max)
    W = solved ? EquilibratedSolve(S, S⁻¹, A) : ExplicitStep(h)
    return x_new, A, W, α_w
end

# ---- The increment limiter ----

"""
    _apply_limiter_diag(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)

[`_apply_limiter`](@ref)'s full computation, returning `(d′, α_l)`: the limited increment
and the scale factor `α_l` it was rescaled by, a [`VerboseSink`](@ref)'s other scalar
alongside [`_bounded_explicit_step_diag`](@ref)'s `α_w`. `_apply_limiter` is
`first ∘ _apply_limiter_diag`.

[`NoLimiter`](@ref) returns `(d, 1)`. [`EndStateSaturationAdjustment`](@ref), for a cell at
or above saturation over its more-supersaturated phase whose full-increment end state would
drop below it, scales `d` by `α_l ∈ [0, 1]` to keep that latent-heated end state at or above
saturation over that phase; otherwise `α_l = 1`.

The limiter reads the state only through [`_condensate_phases`](@ref), so it serves every
substep state that supplies that split.
"""
@inline _apply_limiter_diag(::NoLimiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps) =
    (d, one(eltype(d)))

"""
    _apply_limiter(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)

The limited increment of [`_apply_limiter_diag`](@ref), without its scale factor.
"""
@inline _apply_limiter(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps) =
    first(_apply_limiter_diag(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps))

"""
    _saturation_bisection_count(FT)

Number of bisection iterations to resolve a fraction in `[0, 1]` to the precision of `FT`.
"""
@inline _saturation_bisection_count(::Type{FT}) where {FT} = ceil(Int, -log2(eps(FT)))

"""
    _saturation_bisection_diag(Ssat, latent, x, d, Tsub)

[`_saturation_bisection`](@ref)'s full computation, returning `(d′, s)`: the scaled
increment and the scale factor `s` itself. `Ssat(x, T)` and `latent(d)` close over the
substep context.

Scales the increment `d` at state `x` so the latent-heated end state keeps `Ssat >= 0`, for
a state that begins with `Ssat >= 0`, at `s = 1` otherwise.
"""
@inline function _saturation_bisection_diag(
    Ssat::FS, latent::FL, x::SA.StaticVector{N, FT}, d, Tsub,
) where {FS, FL, N, FT}
    xf = max.(x .+ d, 0)
    if Ssat(x, Tsub) >= 0 && Ssat(xf, Tsub + latent(xf .- x)) < 0
        lo = zero(FT)
        hi = one(FT)
        for _ in 1:_saturation_bisection_count(FT)
            s = (lo + hi) / 2
            xs = max.(x .+ s .* d, 0)
            if Ssat(xs, Tsub + latent(xs .- x)) >= 0
                lo = s
            else
                hi = s
            end
        end
        return lo .* d, lo
    end
    return d, one(FT)
end

"""
    _saturation_bisection(Ssat, latent, x, d, Tsub)

The scaled increment of [`_saturation_bisection_diag`](@ref), without its scale factor.
"""
@inline _saturation_bisection(Ssat, latent, x, d, Tsub) =
    first(_saturation_bisection_diag(Ssat, latent, x, d, Tsub))

@inline function _apply_limiter_diag(::EndStateSaturationAdjustment,
    x::SA.StaticVector{N, FT}, d::SA.StaticVector{N, FT},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {N, FT}
    function Ssat(xx, TT)
        q_liq, q_ice = _condensate_phases(xx)
        return max(
            TDI.supersaturation_over_ice(tps, q_tot, q_liq, q_ice, ρ, TT),
            TDI.supersaturation_over_liquid(tps, q_tot, q_liq, q_ice, ρ, TT),
        )
    end
    function latent(dd)
        dq_liq, dq_ice = _condensate_phases(dd)
        return Lv_over_cp * dq_liq + Ls_over_cp * dq_ice
    end
    return _saturation_bisection_diag(Ssat, latent, x, d, Tsub)
end

# ---- One substep ----

"""
    _rosenbrock_substep_diag(mode, g, x, h, q_tot, ρ, Tsub, Lv_over_cp, Ls_over_cp, tps,
        ρ_min = nothing, ρ_max = nothing, sink = nothing)

[`_rosenbrock_substep`](@ref)'s full computation, returning `(x_new, diag)` where `diag` is
the substep's own internals: the tendency `f` and post-growth Jacobian `J` evaluated at `x`,
the equilibrated system matrix `A` [`_rosenbrock_update_diag`](@ref) inverted (accepted
branch only), whether the accepted-Jacobian solve was taken, and the limiter's engagement
(`limiter_engaged`, `limiter_correction`, the per-component `limiter_correction_vec`).
`_rosenbrock_substep` is `first ∘ _rosenbrock_substep_diag`, so the production path and any
diagnostic path share this one function and cannot diverge.

`J` and `A` are INVALID whenever `accepted` is `false` - the state-non-finite branch
attempts no Jacobian at all, and the Jacobian-non-finite branch attempts no accepted solve,
so no system is ever built. Both fields hold a concrete `NaN`-filled matrix in that case,
purely so `diag`'s type does not vary with control flow; the `NaN` is a marker, not the
fact. A consumer MUST gate on `accepted`, never on `isnan(J[i, j])` - `accepted` is
authoritative and `isnan` would go silently wrong if a finite linearization were ever
corrupted some other way, or if the sentinel value changed.

`sink` observes this substep through [`record!`](@ref): this function builds the context
([`_record_context`](@ref)) - the per-process rates ([`_recorded_processes`](@ref),
[`_per_process_rates`](@ref)), the [`SubstepOperator`](@ref) `W` whose solve actually
produced the accepted increment ([`_rosenbrock_update_diag`](@ref)'s or, on either fallback
branch, [`ExplicitStep`](@ref)), the water bound's and the limiter's scale factors, and the
accepted increment `x_new - x` - and hands it to `sink`. Production's default `sink =
nothing` resolves to the no-op `record!(::Nothing, ctx)`, so this costs production nothing
by construction: `_recorded_processes` never calls its `per_process` argument for a
non-[`VerboseSink`](@ref), and `record!`'s empty body on the no-op sinks leaves the whole
context unused, which the compiler removes along with everything built only to feed it.
"""
@inline function _rosenbrock_substep_diag(
    mode::RosenbrockAverage, g, x::SA.StaticVector{N, FT}, h, q_tot, ρ, Tsub,
    Lv_over_cp, Ls_over_cp, tps, ρ_min = nothing, ρ_max = nothing, sink = nothing,
) where {N, FT}
    A_absent = SA.SMatrix{N, N, FT}(ntuple(_ -> FT(NaN), Val(N * N)))
    if all(isfinite, x)
        f, J_raw = _tendency_and_jacobian(mode.jacobian, g, x)
        J = _apply_growth(mode.growth, J_raw)
        z = _species_mask(mode.jacobian, mode.growth)(x)
        accepted = all(isfinite, J)
        d_raw, A, W, α_w = if accepted
            x1, Amat, W1, α_w1 = _rosenbrock_update_diag(x, f, J, z, h, q_tot, ρ_min, ρ_max)
            (x1 - x, Amat, W1, α_w1)
        else
            # Bound this branch too: the linearization is unusable here, so the bare
            # explicit step spans the whole substep width and a stiff phase-change rate
            # can convert more mass than the cell holds.
            # `f` comes from the differentiated evaluation's value lane, so a poisoned dual
            # can leave it non-finite where the primal is fine; fall back to the primal
            # tendency rather than propagating NaN through the substep state.
            f_safe = all(isfinite, f) ? f : g(x)
            d1, α_w1 = _bounded_explicit_step_diag(x, _euler_update(x, f_safe, h) - x, q_tot)
            (d1, A_absent, ExplicitStep(h), α_w1)
        end
        d_limited, α_l = _apply_limiter_diag(
            mode.limiter, x, d_raw, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)
        limiter_correction_vec = d_limited .- d_raw
        limiter_correction = maximum(abs.(limiter_correction_vec))
        x_new = _apply_positivity_floor(x, d_limited, ρ_min, ρ_max)
        pp = _recorded_processes(sink, () -> _per_process_rates(g, x))
        record!(sink, _record_context(pp, W, h, α_w, α_l, x_new - x))
        diag = (; f, J, A, accepted, limiter_engaged = limiter_correction > zero(FT),
            limiter_correction, limiter_correction_vec)
        return x_new, diag
    else
        f = g(x)
        x_new = _euler_update(x, f, h)
        pp = _recorded_processes(sink, () -> _per_process_rates(g, x))
        record!(sink, _record_context(pp, ExplicitStep(h), h, one(FT), one(FT), x_new - x))
        J_absent = SA.SMatrix{N, N, FT}(ntuple(_ -> FT(NaN), Val(N * N)))
        diag = (; f, J = J_absent, A = A_absent, accepted = false, limiter_engaged = false,
            limiter_correction = zero(FT), limiter_correction_vec = zero(x))
        return x_new, diag
    end
end

"""
    _rosenbrock_substep(mode, g, x, h, q_tot, ρ, Tsub, Lv_over_cp, Ls_over_cp, tps,
        ρ_min = nothing, ρ_max = nothing, sink = nothing)

One substep of the Rosenbrock-Euler average at state `x`: the accepted-vs-rejected Jacobian
branch with its `f_safe` primal fallback, the water bound, the increment limiter, and the
positivity floor. Shared between the production entry and any diagnostic entry, so both run
the identical substep physics rather than separate implementations.
"""
@inline _rosenbrock_substep(args...) = first(_rosenbrock_substep_diag(args...))
