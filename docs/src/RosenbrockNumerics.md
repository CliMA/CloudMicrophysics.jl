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

- **`Jacobian`** — the matrix `J` used in the substep solve.
  - `ExactJacobian` — the exact tendency derivative, formed with `ForwardDiff`.

- **`GrowthTreatment`** — how the positive (growth) diagonal of `J` enters the implicit operator.
  - `ImplicitGrowth` — leave `J` unchanged.
  - `ExplicitGrowthDiagonal` — zero the positive diagonal entries of `J`, so a growth mode is taken explicitly
    and only the decay diagonal remains in the implicit operator.

- **`TendencyLimiter`** — a limiter applied to the realized substep increment.
  - `NoLimiter`.
  - `EndStateSaturationAdjustment` — scale the increment so the latent-heated end state does not cross ice
    saturation (see below).

The supported preset configuration is:

| preset | Jacobian | growth | limiter |
|---|---|---|---|
| `rosenbrock_exact()` | `ExactJacobian` | `ExplicitGrowthDiagonal` | `EndStateSaturationAdjustment` |

On the two-moment + P3 model only `ExactJacobian` is available; use `rosenbrock_exact()`.

### Extending the framework

To add a new Jacobian, define `struct MyJacobian <: Jacobian end` and a `_species_mask(::MyJacobian, ::GrowthTreatment)`
method (returning a `x -> z` species projection). A new growth treatment is a `GrowthTreatment` subtype plus an
`_apply_growth(::MyGrowth, J)` method; a new limiter is a `TendencyLimiter` subtype plus an
`_apply_limiter(::MyLimiter, x, Δ, ...)` method. The substep driver dispatches on the option types at compile time,
so a configured mode resolves with no run-time branch.

## The coarse-step ice-growth instability

At a cold, ice-supersaturated state carrying supercooled cloud liquid, the ice-growth tendency has an
autocatalytic mode: rime mass grows by collecting cloud droplets, and denser rimed particles fall faster and
sweep out more liquid, so the rime-mass tendency increases with rime mass. The exact Jacobian carries this as a
positive diagonal in the rime-mass (`q_rim`) row. At a representative state — `ρ = 1.0` kg m⁻³, `T = 263` K,
`q_tot = 10⁻²`, `q_lcl = 2 × 10⁻³`, `n_lcl = 10⁸` m⁻³, `q_ice = 2 × 10⁻³`, `n_ice = 10⁴` m⁻³, unrimed, ice
supersaturation `S_ice ≈ 1.8` — the diagonal is `+5 × 10⁻² s⁻¹` (time scale ≈ 20 s, identical in `Float32` and
`Float64`), and it grows past `10⁻¹ s⁻¹` at colder, more liquid-rich states. The mode requires supercooled
liquid: with the same ice state but no cloud liquid the rime-mass diagonal falls to order `10⁻⁵ s⁻¹`, and the
pure-deposition diagonal is negative there (vapor depletion opposes further deposition). With the exact Jacobian
and `ImplicitGrowth`, the implicit operator `I/h − J` loses positive-definiteness once the growth eigenvalue
exceeds `1/h`, i.e. once the substep is coarse relative to the growth time scale. The single substep then
overshoots the nonlinear limit the linear operator does not see: ice is over-grown past the available condensate,
the latent heating drives a spurious temperature excursion, and the state goes non-physical. This crash is a
property of the single-column convective configuration, not of an isolated cell; the growth diagonal above is an
isolated-cell measurement, but the crash itself appears only in the coupled single-column run.

### What cures it

The exact preset removes the growth mode from the implicit operator and bounds the now-explicit growth by the
physical saturation limit:

- **`ExplicitGrowthDiagonal`** zeros the positive diagonal, so the implicit operator carries only non-positive
  modes and is well-conditioned at any substep size. The exact off-diagonal couplings are retained, so accuracy
  at cold, supersaturated cells is good.
- **`EndStateSaturationAdjustment`** scales the substep increment by the largest `s ∈ [0, 1]` for which the
  latent-heated end state stays at or above ice saturation. It acts only on cells that begin at or above
  saturation (a subsaturated, evaporating or sublimating cell cannot over-deposit, so its increment is returned
  unchanged). It is a no-op at fine substeps and engages only when the full step would cross saturation. The
  bisection count is set from the float precision.

Both pieces are required: zeroing the growth diagonal alone leaves the explicit growth unbounded at the coarsest
single-substep steps, and the saturation adjustment supplies the missing nonlinear bound. Together they make the
exact scheme robust across the resolved time-step envelope.

!!! note "Use two or more substeps for accurate climate"
    At a single substep the explicit growth is bounded only by the saturation adjustment, which over-produces
    precipitation at coarse time steps. Two or more substeps recover accurate precipitation; the saturation
    adjustment is then rarely active.

## Approaches that did not cure the instability

These were tried and are not part of the supported framework.

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
- **Implicit temperature** (promoting `T` into the implicitly solved state). It removes the error of the
  operator-split between-substep temperature update, but the dominant brake on the growth — the nonlinear
  condensate depletion — is not linear, so a linear implicit temperature feedback does not bound the growth
  overshoot at fixed coarse substeps.
