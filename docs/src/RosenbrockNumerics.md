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
  - `DonorJacobian` — the donor-based linearization `M`: each transfer is linearized in its donor species,
    vapor sources enter as a constant, and rates are floored by `max(q_min, q_donor)`. This is the matrix the
    operational `LinearizedAverage` mode uses.
  - `CoupledDonorJacobian` — the donor-based matrix with the vapor-competition (Wegener–Bergeron–Findeisen)
    coupling added. The donor-based linearization keeps only donor-species slopes; restoring the dependence of
    each rate on the shared vapor specific content recovers the cross-species coupling and corrects the sign of
    the snow-from-cloud-liquid entry. The direct condensate dependence of the rates (rain ventilation, the
    availability terms) is not recovered; use `ExactJacobian` for the full derivative.
  - `ExactJacobian` — the exact tendency derivative, formed with `ForwardDiff`.

- **`GrowthTreatment`** — how the positive (growth) diagonal of `J` enters the implicit operator.
  - `ImplicitGrowth` — leave `J` unchanged.
  - `ExplicitGrowthDiagonal` — zero the positive diagonal entries of `J`, so a growth mode is taken explicitly
    and only the decay diagonal remains in the implicit operator.

- **`TendencyLimiter`** — a limiter applied to the realized substep increment.
  - `NoLimiter`.
  - `EndStateSaturationAdjustment` — scale the increment so the latent-heated end state does not cross ice
    saturation (see below).

Three preset configurations are supported:

| preset | Jacobian | growth | limiter |
|---|---|---|---|
| `rosenbrock_donor()` | `DonorJacobian` | `ImplicitGrowth` | `NoLimiter` |
| `rosenbrock_coupled()` | `CoupledDonorJacobian` | `ImplicitGrowth` | `NoLimiter` |
| `rosenbrock_exact()` | `ExactJacobian` | `ExplicitGrowthDiagonal` | `EndStateSaturationAdjustment` |

`rosenbrock_donor()` reproduces `LinearizedAverage` (the operational donor-based scheme), now expressed within
the unified framework; in `Float64` the two agree to round-off. On the two-moment + P3 model only
`ExactJacobian` is available (there is no donor-based matrix there); use `rosenbrock_exact()`.

### Extending the framework

To add a new Jacobian, define `struct MyJacobian <: Jacobian end` and the methods `_jacobian_provider(::MyJacobian)`
(returning a `(g, x) -> J` provider) and `_species_mask(::MyJacobian, ::GrowthTreatment)`. A new growth treatment
is a `GrowthTreatment` subtype plus an `_apply_growth(::MyGrowth, J)` method; a new limiter is a `TendencyLimiter`
subtype plus an `_apply_limiter(::MyLimiter, x, Δ, ...)` method. The substep driver dispatches on the option types
at compile time, so a configured mode resolves with no run-time branch.

## The coarse-step deposition instability

At a cold, ice-supersaturated state the ice-deposition tendency has an autocatalytic growth mode: the snow
deposition rate increases with snow content (collector self-gain), a positive Jacobian diagonal of order
`+5 × 10⁻² s⁻¹` (time scale ≈ 20 s). With the exact Jacobian and `ImplicitGrowth`, the implicit operator
`I/h − J` loses positive-definiteness once the growth eigenvalue exceeds `1/h`, i.e. once the substep is coarse
relative to the growth time scale. The single substep then overshoots the nonlinear saturation limit (which the
linear operator does not see): snow is over-deposited far past the available vapor, the latent heating drives a
spurious temperature excursion, and the state goes non-physical. With the exact Jacobian this crash occurs at
every coarse time step in the single-column convective test.

### What resolves it

The exact preset removes the growth mode from the implicit operator and bounds the now-explicit growth by the
physical saturation limit:

- **`ExplicitGrowthDiagonal`** zeros the positive diagonal, so the implicit operator carries only non-positive
  modes and is well-conditioned at any substep size. The exact off-diagonal couplings (which the donor-based
  matrix drops) are retained, so accuracy at cold, supersaturated cells is better than the donor scheme.
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

## Approaches that did not resolve the instability

These were tried and are not part of the supported framework; they are recorded here because the failure modes
are instructive.

- **Field-of-values growth clamp** (a uniform diagonal shift bringing the operator's rightmost eigenvalue to
  `α/h`). It stabilizes the linear operator but does not bound the single-step explicit overshoot of the
  nonlinear source as it approaches saturation: the crash is a saturation overshoot, not an operator
  amplification, so the clamp only delays it. With `α` near one the near-singular resolvent it leaves actually
  amplifies the deposition into a larger overshoot.
- **Diagonal growth clamp** (cap each positive diagonal at `α/h`). A cheaper variant of the above with the same
  limitation, and capping a positive diagonal balanced by off-diagonal structure can itself destabilize an
  otherwise-stable step.
- **Smooth species mask** (a differentiable replacement for the near-empty species mask). At coarse single
  substeps it routes activating ice and liquid species to a forward-Euler step, which itself overshoots
  the fast deposition.
- **Implicit temperature** (promoting `T` into the implicitly solved state). It removes the operator-split
  ringing of the between-substep temperature update, but the dominant brake on the growth — the nonlinear vapor
  depletion — is not linear, so a linear implicit temperature feedback does not bound the deposition overshoot
  at fixed coarse substeps.
