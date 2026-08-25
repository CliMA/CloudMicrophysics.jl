# Rosenbrock-average microphysics substepping

The [`RosenbrockAverage`](@ref CloudMicrophysics.BulkMicrophysicsTendencies.RosenbrockAverage) tendency mode
returns time-averaged microphysics tendencies over a time step `Δt` by taking `nsub` linearized-implicit
(Rosenbrock-Euler) substeps. Each substep solves

```math
\left(\frac{I}{h} - J\right)\, \Delta = f(x), \qquad x \leftarrow \max(x + \Delta,\, 0), \qquad h = \Delta t / n_\mathrm{sub},
```

where `f` is the raw pointwise tendency, `x` the species state, and `J` a matrix that approximates the tendency
Jacobian. The averaged tendency returned is `(x_final - x_initial) / Δt`. Temperature is advanced between
substeps from the latent heat of the realized increment.

## Options

`RosenbrockAverage` is parameterized by three independent option families:

- **`Jacobian`** - the matrix `J` used in the substep solve.
  - `DonorJacobian` - the donor-based linearization `M`: each transfer is linearized in its donor species,
    vapor sources enter as a constant, and rates are floored by `max(q_min, q_donor)`. This is the matrix the
    operational `LinearizedAverage` mode uses.
  - `CoupledDonorJacobian` - the donor-based matrix with the vapor-competition (Wegener–Bergeron–Findeisen)
    coupling added. The donor-based linearization keeps only donor-species slopes; restoring the dependence of
    each rate on the shared vapor specific content recovers the cross-species coupling and corrects the sign of
    the snow-from-cloud-liquid entry. The direct condensate dependence of the rates (rain ventilation, the
    availability terms) is not recovered; use `ExactJacobian` for the full derivative.
  - `ExactJacobian` - the exact tendency derivative, formed with `ForwardDiff`.
  - `ManualJacobian` - a hand-built approximation for the two-moment + P3 model, closed-form on the
    phase-change and number-adjustment couplings and donor-linearized on the remaining transfers. See
    [The folded and temperature-coupled Jacobians on the 2M+P3 model](@ref) below.
  - `TemperatureCoupledJacobian` - the two-moment + P3 model with the substep temperature promoted to a
    ninth prognostic variable, so the phase-change couplings are represented without folding the
    psychrometric feedback into the relaxation timescale. See the same section.

- **`GrowthTreatment`** - how the positive (growth) diagonal of `J` enters the implicit operator.
  - `ImplicitGrowth` - leave `J` unchanged.
  - `ExplicitGrowthDiagonal` - zero the positive diagonal entries of `J`, so a growth mode is taken explicitly
    and only the decay diagonal remains in the implicit operator.

- **`TendencyLimiter`** - a limiter applied to the realized substep increment.
  - `NoLimiter`.
  - `EndStateSaturationAdjustment` - scale the increment so the latent-heated end state does not cross saturation
    over its more-supersaturated phase, `max(S_ice, S_liq)` (see below).

Four preset configurations are supported:

| preset | Jacobian | growth | limiter |
|---|---|---|---|
| `rosenbrock_coupled()` | `CoupledDonorJacobian` | `ImplicitGrowth` | `NoLimiter` |
| `rosenbrock_exact()` | `ExactJacobian` | `ExplicitGrowthDiagonal` | `EndStateSaturationAdjustment` |
| `rosenbrock_manual()` | `ManualJacobian` | `ExplicitGrowthDiagonal` | `EndStateSaturationAdjustment` |
| `rosenbrock_manual_temperature()` | `TemperatureCoupledJacobian` | `ExplicitGrowthDiagonal` | `NoLimiter` |

The remaining combination, `DonorJacobian` with `ImplicitGrowth` and `NoLimiter`, is the operational
donor-based scheme itself, and it keeps its own name:
[`LinearizedAverage`](@ref CloudMicrophysics.BulkMicrophysicsTendencies.LinearizedAverage) is that
configuration expressed within the unified framework, so no preset duplicates it. The donor-based
matrices are not available on the two-moment + P3 model; use `rosenbrock_exact()`, `rosenbrock_manual()`,
or `rosenbrock_manual_temperature()` there.

The `Verbose(mode)` wrapper additionally returns the per-process tendencies realized by the implicit solve,
attributed through the same substep factorization so that they sum to the net of the unlimited solve.

### Extending the framework

To add a new Jacobian, define `struct MyJacobian <: Jacobian end` and the methods `_jacobian_provider(::MyJacobian)`
(returning a `(g, x) -> J` provider) and `_species_mask(::MyJacobian, ::GrowthTreatment)`. A new growth treatment
is a `GrowthTreatment` subtype plus an `_apply_growth(::MyGrowth, J)` method; a new limiter is a `TendencyLimiter`
subtype plus an `_apply_limiter(::MyLimiter, x, Δ, ...)` method. The substep driver dispatches on the option types
at compile time, so a configured mode resolves with no run-time branch.

## The folded and temperature-coupled Jacobians on the 2M+P3 model

The two-moment + P3 model exposes two closures of the same phase-change (condensation/evaporation,
deposition/sublimation) physics: the 8×8 (`ManualJacobian`, preset `rosenbrock_manual()`) and the
temperature-coupled 9×9 (`TemperatureCoupledJacobian`, preset `rosenbrock_manual_temperature()`). Both
linearize the same relaxation-to-saturation process and differ only in whether the substep temperature is
carried as an explicit prognostic variable or held frozen across the substep. This section states what each
closure computes, why a psychrometric correction on the relaxation timescale would double-count the
latent-heat feedback the coupled system already carries, and what the exact block-elimination identity
between the two systems does and does not say about the matrices the package assembles.

### Two closures of the same relaxation

The phase-change tendency in both closures is a relaxation toward saturation, with capacitance-integral
timescale `τ` and supersaturation `s`,

```math
s = q_v - q_{v,\mathrm{sat}}(T).
```

Both `Instantaneous2MP3Tendency` and `Temperature2MP3Tendency` evaluate the bare relaxation

```math
\partial_t q = \frac{s}{\tau}
```

directly (a limited branch, active when the relaxation would exceed the donor's own mass, replaces `s` with
the donor-capped `-min(-s, max(0, q_limit))`; both closures carry it identically).
`Instantaneous2MP3Tendency` freezes `T` for the duration of the substep and represents the latent-heat
feedback only through the outer time-stepping loop; `Temperature2MP3Tendency` instead carries `T` as a ninth
prognostic component and lets it respond within the substep.
[`_condevap_derivs`](@ref CloudMicrophysics.BulkMicrophysicsTendencies._condevap_derivs) supplies the
closed-form derivative for both, including the timescale's own dependence on the donor species.
[`_jacobian_2mp3t_manual`](@ref CloudMicrophysics.BulkMicrophysicsTendencies._jacobian_2mp3t_manual) builds
its species block by calling
[`_jacobian_2mp3_manual`](@ref CloudMicrophysics.BulkMicrophysicsTendencies._jacobian_2mp3_manual) on the
same state before appending the temperature row and column; the shared thermodynamic quantities (`τ`, `Γ`,
latent heats, the moist heat capacity) come from
[`_phase_relaxation_context`](@ref CloudMicrophysics.BulkMicrophysicsTendencies._phase_relaxation_context),
evaluated once per substep and consumed by both.

### Why the relaxation is posted bare, without a psychrometric fold

A natural-looking variant of the relaxation folds a psychrometric correction into the timescale,

```math
\partial_t q = \frac{s}{\tau\,\Gamma}, \qquad
\Gamma = 1 + \frac{L}{c_p}\frac{\partial q_{v,\mathrm{sat}}}{\partial T} \quad (\texttt{CMNonEq.gamma\_helper}),
```

intended to slow condensation in proportion to the warming its own latent heat produces. Write `s`'s
dependence on the condensed mass `q` along the path either closure generates: condensing an increment of `q`
draws `q_v` down directly and, through the latent-heat response of temperature (weight `c = L/c_p` on the
mass species, `0` on number and rime species - the same vector `Temperature2MP3Tendency`'s temperature
tendency uses), draws `q_v,sat` up in proportion. The total draw-down of `s` per unit `q` condensed is

```math
\frac{\mathrm{d}s}{\mathrm{d}q} = -\Bigl(1 + \frac{L}{c_p}\frac{\partial q_{v,\mathrm{sat}}}{\partial T}\Bigr) = -\Gamma,
```

the same `Γ` above. On the temperature-coupled system, `s` evolves continuously as `T` responds to the bare
relaxation `∂ₜq = s/τ`:

```math
\dot s = \frac{\mathrm{d}s}{\mathrm{d}q}\,\dot q = -\Gamma \cdot \frac{s}{\tau} = -\frac{\Gamma}{\tau}\,s,
```

an accelerated decay: the coupled state already carries the correction, through `T` itself. A closure that
additionally posts the rate `s/(τΓ)` applies the same correction a second time. Read through the same
`ds/dq = -Γ` slope, `s/(τΓ)` corresponds to `ṡ = -Γ·s/(τΓ) = -s/τ`, the bare, unaccelerated decay - `Γ` times
too slow relative to the coupled system's own rate. Both forms reach the identical equilibrium
`q_eq - q = s/Γ`, since the same `Γ` sets the equilibrium either way; only the rate of approach differs, by a
full factor of `Γ`, not a small residual. At a representative mixed-phase state (`ρ = 0.85` kg m⁻³,
`T = 265` K, cloud liquid, rain, and rimed ice all present and supersaturated), the folded form's ice-channel
Jacobian entry differs from the bare rate's by a factor of `2.86`, matching `Γ_ice² = 2.85` at that state's
`Γ_ice = 1.69` - a further factor of `Γ` beyond the rate itself, since the Jacobian differentiates the rate.

Two published treatments of this feedback exist, and the folded form is neither. Korolev & Mazin
[KorolevMazin2003](@cite) and Milbrandt & Yau (2005) place the correction in the growth coefficient itself -
the Mason-form coefficient `G = 1/(F_k + F_d)`, with `F_k` and `F_d` the particle-scale heat-conduction and
vapor-diffusion resistances - so a timescale built from `G` is already latent-heat corrected and needs no
further factor. Morrison & Grabowski [MorrisonGrabowski2008_supersat](@cite) and P3's own reference
implementation instead use the bare vapor-diffusivity coefficient and apply `Γ` - a bulk-heat-capacity
correction, not a particle-scale one - at the point the rate is used, paired with an analytic integration of
the coupled supersaturation relaxation over the full model time step.

`G_func_liquid` and `G_func_ice` in this package already carry the particle-scale correction (the full
`F_k + F_d` form), so applying `Γ` on top of them would combine both treatments' corrections on one rate, a
combination neither published scheme uses. Posting the relaxation bare leaves the growth coefficient's own
correction as the only one applied, matching Korolev & Mazin and Milbrandt & Yau exactly; no comparison
against the Morrison & Grabowski convention is implied, since it corrects a different (bare) coefficient at
a different point in the calculation and remains internally consistent on its own terms.

### The general block-elimination identity

Independent of which closure `f(q, T)` computes, the temperature-coupled system's structure admits an exact
reduction. Write the eight-species tendency as `f(q, T)` with Jacobian blocks `f_q = ∂f/∂q`, `f_T = ∂f/∂T`.
`MicroState2MP3T`'s temperature equation is exactly the latent-heat combination of the species tendency,
`Ṫ = cᵀf(q, T)` with the constant `c` above, so its own Jacobian blocks are `cᵀ` applied to the species
block: the temperature row is `cᵀf_q` and the temperature corner is `cᵀf_T`. Both hold to floating-point
precision by construction (verified below) - this is how `_temperature_2mp3_tendency` and
`_jacobian_2mp3t_manual` build the temperature row and column, not an independent physical assumption.

The 9×9 Rosenbrock substep solves `(I/h - J₉)Δy = f(y)` for `y = (q, T)`, with

```math
J_9 = \begin{pmatrix} f_q & f_T \\ c^\mathsf{T}f_q & c^\mathsf{T}f_T \end{pmatrix}.
```

Impose the ansatz `ΔT = cᵀΔq` - the temperature increment tracks the species increment through the same
weights as the tendency itself - and substitute into the species-block row:

```math
\Bigl(\frac{I}{h} - f_q - f_T c^\mathsf{T}\Bigr)\Delta q = f_q(y).
```

Left-multiplying this 8×8 equation by `cᵀ` reproduces the temperature-row equation under the same ansatz
exactly, so the temperature-row equation is satisfied automatically whenever the species-block equation is.
The 9×9 step therefore reduces exactly to an 8×8 step with effective Jacobian

```math
J_8^{\mathrm{elim}} = f_q + f_T\,c^\mathsf{T}, \qquad \Delta T = c^\mathsf{T}\Delta q,
```

for any `f` satisfying `Ṫ = cᵀf` with `c` constant. The identity needs two conditions: `c` constant (`L` and
`c_p` are evaluated once per substep and not differentiated), and `T` staying on the latent-heat manifold
`dT = cᵀdq` - true for the phase-change and freezing/melting sources here, since every temperature change in
this substep is the latent heat of a resolved species tendency, but not true of a temperature tendency from
any other source (radiative heating, advection, mixing), which the 8×8's analytic elimination has no state to
represent and the 9×9 does.

`ManualJacobian` and `_jacobian_2mp3t_manual`'s species block compute the identical `f_q` directly:
`_jacobian_2mp3t_manual` calls `_jacobian_2mp3_manual` on the same state and appends the temperature row and
column to its output, so the 9×9's species block and `ManualJacobian`'s output are the same matrix, `f_q`,
by construction.

That is a different statement from the elimination identity, and the two should not be run together.
`ManualJacobian` returns `f_q` alone: it is the frozen-temperature form, with no temperature column and no
`f_T` to contract. `J₈ᵉˡⁱᵐ = f_q + f_T cᵀ` is `f_q` plus the rank-one term that column supplies, so it is a
different matrix, and no code path in this package assembles it. The elimination is a design option - an 8×8
step that carries the temperature coupling without a ninth state, marching `T` by `ΔT = cᵀΔq` afterwards -
and the package implements the 9×9 instead. What `rosenbrock_manual()` solves with is `f_q`.

### Consequences for substep behavior

With both closures computing the identical species-block Jacobian, the difference between an 8×8 substep and
a 9×9 substep is no longer a difference in the phase-change rate. It is the difference between freezing `T`
for the duration of the substep, with the outer time-stepping loop supplying the only temperature feedback
(8×8), and letting `T` respond within the substep (9×9). Stability comparisons between the two substep
structures - outside this section's scope - find the temperature-coupled form tolerates substantially larger
outer steps before the substep's own linearization error becomes visible, since the coupled state removes
the mismatch between a rate evaluated at fixed `T` and a feedback applied outside the solve.

### Numerical verification

Verified directly against the code, at a representative mixed-phase state (`ρ = 0.85` kg m⁻³, `T = 265` K,
cloud liquid, rain, and rimed ice all present and supersaturated):

- `Ṫ = cᵀf`: `f₉[9]` and `c·f₉[1:8]` agree to `1e-14` relative.
- The temperature row and corner: `J₉` row 9 vs. `cᵀf_q`, and `J₉[9,9]` vs. `cᵀf_T`, agree to floating-point
  roundoff (`~1e-9`, `~1e-11`).
- `J₉`'s species block matches `ManualJacobian`'s output exactly at `Float64`, and to within `28` units in
  the last place at `Float32` (`~2×10⁻⁶` relative) - consistent with ordinary rounding-order differences
  between two independently compiled evaluations of the same formula rather than a residual mismatch. This
  measures `f_q` against `f_q`; it says nothing about `J₈ᵉˡⁱᵐ`, which is not assembled anywhere and so is not
  measured here.
- Across a corpus of twenty states spanning the vapor and evaporation-limited branches, `ManualJacobian`'s
  condensation self-derivative agrees with a full automatic-differentiation reference to a ratio of `1.0`, to
  four significant figures, at every state.

The invariant the
[`p3_2m_process_rates`](@ref CloudMicrophysics.BulkMicrophysicsTendencies.p3_2m_process_rates) breakdown
itself satisfies - that summing it reproduces the raw entry tendency bit for bit - holds by construction,
since the entry tendency is that sum, and a test asserts the slot names and their order. The
block-elimination identity is not checked, because there is nothing in the package to check it against: it
relates the 9×9 system to a hypothetical 8×8 one, and a test would have to assemble `f_q + f_T cᵀ` from the
9×9's own blocks and compare the step it produces with the 9×9's own.

## The coarse-step ice-growth instability

At a cold, ice-supersaturated state carrying supercooled cloud liquid, the ice-growth tendency has an
autocatalytic mode: rime mass grows by collecting cloud droplets, and denser rimed particles fall faster and
sweep out more liquid, so the rime-mass tendency increases with rime mass. The exact Jacobian carries this as a
positive diagonal in the rime-mass (`q_rim`) row. At a representative state - `ρ = 1.0` kg m⁻³, `T = 263` K,
`q_tot = 10⁻²`, `q_lcl = 2 × 10⁻³`, `n_lcl = 10⁸` m⁻³, `q_ice = 2 × 10⁻³`, `n_ice = 10⁴` m⁻³, unrimed, ice
supersaturation `S_ice ≈ 1.8` - the diagonal is `+5 × 10⁻² s⁻¹` (time scale ≈ 20 s, identical in `Float32` and
`Float64`), and it grows past `10⁻¹ s⁻¹` at colder, more liquid-rich states. The mode requires supercooled
liquid: with the same ice state but no cloud liquid the rime-mass diagonal falls to order `10⁻⁵ s⁻¹`, and the
pure-deposition diagonal is negative there (vapor depletion opposes further deposition). With the exact Jacobian
and `ImplicitGrowth`, the implicit operator `I/h − J` loses positive-definiteness once the growth eigenvalue
exceeds `1/h`, i.e. once the substep is coarse relative to the growth time scale. The single substep then
overshoots the nonlinear limit the linear operator does not see: ice is over-grown past the available condensate,
the latent heating drives a spurious temperature excursion, and the state goes non-physical. This crash is a
property of the single-column convective configuration, not of an isolated cell; the growth diagonal above is an
isolated-cell measurement, but the crash itself appears only in the coupled single-column run.

### What resolves it

The exact preset removes the growth mode from the implicit operator and bounds the now-explicit growth by the
physical saturation limit:

- **`ExplicitGrowthDiagonal`** zeros the positive diagonal, so the implicit operator carries only non-positive
  modes and is well-conditioned at any substep size. The exact off-diagonal couplings (which the donor-based
  matrix drops) are retained, so accuracy at cold, supersaturated cells is better than the donor scheme.
- **`EndStateSaturationAdjustment`** scales the substep increment by the largest `s ∈ [0, 1]` for which the
  latent-heated end state stays at or above saturation over its more-supersaturated phase (`max(S_ice, S_liq)`;
  see the next section). It acts only on cells that begin at or above saturation (a subsaturated, evaporating or
  sublimating cell cannot over-deposit, so its increment is returned unchanged). It is a no-op at fine substeps
  and engages only when the full step would cross saturation. The bisection count is set from the float
  precision.

Both pieces are required: zeroing the growth diagonal alone leaves the explicit growth unbounded at the coarsest
single-substep steps, and the saturation adjustment supplies the missing nonlinear bound. Together they make the
exact scheme robust across the resolved time-step envelope.

!!! note "Use two or more substeps for accurate climate"
    At a single substep the explicit growth is bounded only by the saturation adjustment, which over-produces
    precipitation at coarse time steps. Two or more substeps recover accurate precipitation; the saturation
    adjustment is then rarely active.

## The mixed-phase saturation criterion

A cell has two saturation thresholds, over liquid and over ice, and the condensation and deposition processes draw
on a single shared vapor reservoir. `EndStateSaturationAdjustment` limits the increment on the more-supersaturated
phase, `max(S_ice, S_liq)`: it keeps the latent-heated end state at or above the saturation of whichever phase
carries the larger supersaturation. Equivalently the vapor floor is the lower of the two saturation specific
humidities, since the smaller `q_sat` is the larger supersaturation.

The two saturation curves cross at the freezing point (panel (c) below). Below freezing `q_sat_ice < q_sat_liq`, so
ice carries the larger supersaturation and the criterion binds on ice - reducing exactly to the ice-saturation
limit that bounds the ice-growth instability, so the cold behaviour is unchanged. Above freezing the curves swap
and the criterion binds on liquid; an ice-only limit there would stop vapor depletion early and leave the cell
supersaturated over liquid, suppressing warm cloud, which binding on the more-supersaturated phase avoids.

A single 0-D parcel that exchanges vapor with cloud liquid and cloud ice only (no collection or precipitation)
isolates the mixed-phase physics behind the two thresholds. Cloud liquid is kinetically fast and cloud ice slow, so
while liquid is present it pins the vapor near water saturation (panel (a), `S_liq ≈ 0`); the parcel stays
supersaturated over ice (`S_ice > 0`) and ice deposits, drawing mass from the evaporating liquid (panel (b)) - the
Wegener–Bergeron–Findeisen transfer. Only once the liquid is exhausted does the vapor relax to ice saturation
(`S_ice → 0`). The drawdown from water saturation to ice saturation is therefore inherently a multi-step process at
the resolved time step.

```@example
include("plots/SatAdjustmentWBF_plots.jl")
```
![](SatAdjustmentWBF.svg)

Binding on `max(S_ice, S_liq)` is a single shared scalar on the whole increment, which is what keeps the limiter
stable: scaling processes independently breaks the coupling between paired transfers and the shared vapor draw. The
criterion binds the correct phase wherever a single phase grows from vapor - ice only (the cold deposition cells,
where it reduces to the ice limit) or liquid only (warm cells). When both phases are supersaturated below freezing,
a vigorous mixed-phase updraft core, the criterion binds on the ice floor (the more-supersaturated phase there) and
so condenses the co-present, still-growing liquid past its own saturation in a single step rather than transferring
it gradually. A fully per-phase floor - binding on the first growing phase to saturate and letting the slower phase
relax over subsequent steps - would represent the gradual transfer more faithfully and is a candidate refinement.
The distinction is inactive at fine substeps, where the limiter rarely engages, and does not affect the cold
instability resolution, which is set by the ice floor.

## Approaches that do not resolve the instability

Each of these addresses the linear operator, or a different error entirely, and none of them bounds the
nonlinear saturation overshoot that the crash consists of. The first three are not part of the supported
framework; the fourth is supported, for the separate reason given in its own section above.

- **Field-of-values growth clamp** (a uniform diagonal shift bringing the operator's rightmost eigenvalue to
  `α/h`). It stabilizes the linear operator but does not bound the single-step explicit overshoot of the
  nonlinear source as it approaches saturation: the crash is a saturation overshoot, not an operator
  amplification, so the shift only delays it. With `α` near one the near-singular resolvent it leaves amplifies
  the growth into a larger overshoot.
- **Diagonal growth clamp** (limit each positive diagonal to `α/h`). It shares the same limitation, and limiting
  a positive diagonal balanced by off-diagonal structure can itself destabilize an otherwise-stable step.
- **Smooth species mask** (a differentiable replacement for the near-empty species mask). At coarse single
  substeps it routes activating ice and liquid species to a forward-Euler step, which itself overshoots the fast
  growth.
- **Implicit temperature** (`TemperatureCoupledJacobian`, promoting `T` into the implicitly solved state). It
  removes the error of the operator-split between-substep temperature update, which is why it is supported,
  but that is not the error the crash is made of: the dominant brake on the growth - the nonlinear condensate
  depletion - is not linear, so a linear implicit temperature feedback does not bound the growth overshoot at
  fixed coarse substeps. Its preset therefore keeps `ExplicitGrowthDiagonal`.

## Internals referenced above

These are not part of the public interface and their signatures may change. They are listed
so that the cross-references on this page resolve: Documenter cannot link to a docstring that
appears on no page, and a broken cross-reference fails the build.

```@docs
CloudMicrophysics.BulkMicrophysicsTendencies.p3_2m_process_rates
CloudMicrophysics.BulkMicrophysicsTendencies._phase_relaxation_context
CloudMicrophysics.BulkMicrophysicsTendencies._condevap_derivs
CloudMicrophysics.BulkMicrophysicsTendencies._jacobian_2mp3_manual
CloudMicrophysics.BulkMicrophysicsTendencies._jacobian_2mp3t_manual
```
