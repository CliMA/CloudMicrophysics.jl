# Bulk Tendencies

## Linearized average tendencies

Microphysical source terms can be stiff, especially for depletion processes such as evaporation, sublimation, and melting. To improve stability and allow larger timesteps, we introduce a **linearized implicit formulation** for computing *time-averaged bulk tendencies*.

The idea is to approximate the nonlinear microphysics tendencies locally as a linear system:

```math
\frac{dq}{dt} \approx M q + e
```

where $q = (q_{\mathrm{lcl}}, q_{\mathrm{icl}}, q_{\mathrm{rai}}, q_{\mathrm{sno}})$, and the matrix $M$ and vector $e$ are constructed from the instantaneous tendencies.

### Donor-based linearization

Each microphysical process is linearized with respect to its **donor species**:

- Transfer processes (e.g. accretion, conversion):
  ```math
  S \;\rightarrow\; D \, q_{\text{donor}}, \quad D = \frac{S}{\max(\epsilon, q_{\text{donor}})}
  ```

- Vapor → condensate sources are treated as **constant sources** (added to $e$)

- Condensate → vapor sinks are treated as **linear sinks**:
  ```math
  S \;\rightarrow\; -D q
  ```

With this formulation, sink terms take the form:

```math
\frac{dq}{dt} = -D q
```

which corresponds to **exponential decay over the timestep**, providing strong numerical stability.

---

## Linearized implicit solve

For a timestep $\Delta t$, we solve the linearized system implicitly:

```math
\frac{q^\star - q^0}{\Delta t} = M q^\star + e
```

which gives:

```math
\left(I/\Delta t - M\right) q^\star = e + q^0/\Delta t
```

The average tendency is then:

```math
\overline{T} = \frac{q^\star - q^0}{\Delta t}
```

---

## Sparse 4×4 structure

The system has a fixed sparse structure:

```math
\begin{bmatrix}
a_{11} & a_{12} & 0      & 0 \\
0      & a_{22} & 0      & 0 \\
a_{31} & 0      & a_{33} & a_{34} \\
a_{41} & a_{42} & a_{43} & a_{44}
\end{bmatrix}
```

The entry $a_{12}$ is the melting of cloud ice into cloud liquid, which couples the $q_{\mathrm{lcl}}$ row to $q_{\mathrm{icl}}$; below freezing it is zero. This allows an efficient solve:

-  $q_{\mathrm{icl}}$ is a scalar solve, and $q_{\mathrm{lcl}}$ follows from it by back-substitution
-  $q_{\mathrm{rai}}$ and $q_{\mathrm{sno}}$ are solved as a $2\times2$ system

This avoids forming or inverting a full dense matrix.

---

## Substepping

A single linearization assumes the operator $M$ is constant over the timestep. To better capture nonlinear effects and regime changes (e.g. near freezing), we apply **substepping**:

- Split the timestep into `nsub` substeps
- At each substep:
  - rebuild $M$ and $e$ from the updated state
  - solve the linearized system
  - update $q$ and temperature

As `nsub` increases, the solution approaches the nonlinear evolution of the system.

The same substep structure also admits a purely explicit variant, in which each substep takes a forward-Euler step of the raw tendency rather than a linearized-implicit solve. That variant is `SubsteppedAverage`; the two are contrasted in [Implicit versus explicit averaging](@ref) below.

---

## Thermodynamic assumption

Within each timestep, we assume that **thermodynamic variables such as density and energy remain approximately constant**. As a result, temperature changes are modeled solely through latent heating:

```math
\frac{dT}{dt}
=
\frac{L_v}{c_p} \left(\dot{q}_{\mathrm{lcl}} + \dot{q}_{\mathrm{rai}}\right)
+
\frac{L_s}{c_p} \left(\dot{q}_{\mathrm{icl}} + \dot{q}_{\mathrm{sno}}\right)
```

This is consistent with the microphysics-only update and avoids coupling to a full thermodynamic solve.

---

## Tendency modes

`bulk_microphysics_tendencies` selects its output through a `TendencyMode` argument. The modes fall into three groups:

- **Instantaneous evaluation.** `Instantaneous` returns the raw pointwise tendency from a single evaluation of all processes, with no linearization and no time-averaging. `InstantaneousVerbose` additionally returns every individual source term for diagnostics. These take no ``\Delta t``.

- **Linearized-implicit averaging.** `LinearizedAverage` and `RosenbrockAverage` return tendencies averaged over ``\Delta t`` from `nsub` linearized-implicit (Rosenbrock-Euler) substeps, each solving the system in [Linearized implicit solve](@ref). `LinearizedAverage` is the donor-based scheme described above and is the configuration used operationally by ClimaAtmos; it forwards to the donor preset of `RosenbrockAverage`. `RosenbrockAverage` exposes the matrix, growth-diagonal, and limiter choices documented in [Linearized-implicit options](@ref).

- **Explicit averaging.** `SubsteppedAverage` returns tendencies averaged over ``\Delta t`` from `nsub` forward-Euler substeps of the raw tendency — the explicit counterpart of the linearized-implicit modes.

The `Verbose(mode)` wrapper returns, alongside the net tendencies, the per-process contributions realized by the implicit solve of `mode`, attributed through the same substep factorization so that they sum to the net of the unlimited solve.

## Implicit versus explicit averaging

`SubsteppedAverage` shares the substep loop, the fixed-``\log\lambda`` and conserved-``q_{\mathrm{tot}}`` assumptions, the latent-heating temperature update, and the positivity clamp of the linearized-implicit modes. The two differ only in how each substep advances the state:

- `LinearizedAverage` / `RosenbrockAverage` solve the linearized-implicit increment ``\left(I/h - J\right)\Delta = f(x)``.
- `SubsteppedAverage` takes the forward-Euler increment ``\Delta = h\, f(x)``.

The implicit solve is unconditionally stable for the linear decay modes it captures: a depletion process linearized as ``\dot q = -D q`` decays monotonically over the substep for any ``h`` (cf. [Donor-based linearization](@ref)). The forward-Euler increment of the same process, ``q \leftarrow q (1 - h D)``, oscillates and diverges once ``h D > 2``, so the explicit average requires ``h`` small relative to the fastest process time scale ``1/D``. For the stiff depletion processes that motivate averaging — evaporation, sublimation, melting — this can demand many substeps, whereas the implicit average stays stable at one. The example below shows that explicit substepping with the raw tendency does not converge even at ten substeps for the near-freezing case, while the linearized-implicit average is already accurate at `nsub = 2`.

Because the forward-Euler step has no implicit damping of its own, `SubsteppedAverage` carries a per-substep limiter to bound the otherwise-unstable increments:

- A coupled-sink limiter scales paired mass and number sinks together (per the warm-rain pairs ``(q_{\mathrm{lcl}}, n_{\mathrm{lcl}})`` and ``(q_{\mathrm{rai}}, n_{\mathrm{rai}})``) so neither member depletes its species below zero within the substep.
- When the limiter is `EndStateSaturationAdjustment`, a saturation-adjustment limiter on the net condensation/deposition prevents a single explicit step from overshooting saturation. With `NoLimiter` this limiter is omitted, which is appropriate when the increment is taken inside an outer implicit timestepper that supplies the stability.

`SubsteppedAverage` reduces to the single-shot limited tendency when ``n_{\mathrm{sub}} \le 1``.

### Consistency of the increment limiter

For the substepped average to converge to the pointwise trajectory as ``h \to 0``, its limiter must be *consistent*: the correction it applies must vanish as the substep shrinks, so that it changes only steps that would otherwise cross saturation and leaves the converged solution unbiased. The one-moment `SubsteppedAverage` and both linearized-implicit averages use the same consistent limiter: a bisection that scales the whole increment only when the full latent-heated end state would cross saturation over its more-supersaturated phase, and is inactive otherwise (see [Increment scaling (`RosenbrockAverage`, `SubsteppedAverage` one-moment)](@ref)). It vanishes as ``h \to 0``, so these modes converge in the warm-rain, cold/ice, and mixed-phase regimes.

The two-moment + P3 `SubsteppedAverage` uses a different saturation-adjustment limiter: it scales each phase's tendencies by an analytic saturation-mass ratio rather than by the bisected end-state factor (see [Analytic condensation/deposition bound (`SubsteppedAverage` two-moment + P3)](@ref)). That ratio is not proportional to ``h``, so the limiter does not vanish under refinement and biases the converged solution in mixed-phase cells: it zeros the ice-melting sink while passing the rain-number source from that melting, over-producing rain number in mixed-cold states, and it carries no per-substep positivity floor, so ice can be driven negative in mixed-warm states. Warm-rain and cold/ice two-moment + P3 states are unaffected. Making this limiter ``h``-consistent and adding a per-substep positivity floor — so two-moment + P3 explicit averaging is both inexpensive and convergent — is being addressed.

The two families are appropriate in different settings. The linearized-implicit average is preferred when stability at coarse time steps is the priority and the cost of forming and solving the substep system is acceptable; it is the operational one-moment choice. The explicit average is cheaper per substep — no Jacobian, no linear solve — but for the two-moment + P3 scheme the linearized-implicit `rosenbrock_exact` (`ExactJacobian` with `ExplicitGrowthDiagonal` and `EndStateSaturationAdjustment`) is roughly an order of magnitude more expensive per microphysics evaluation than the explicit `SubsteppedAverage`; amortized against the rest of a host-model step this difference is typically much smaller, so at the full-step level the cost difference is often modest. The numerics that motivate `rosenbrock_exact` — the coarse-step deposition instability and the mixed-phase saturation criterion — are detailed in [Bulk-tendency averaging: modes and numerics](@ref).

---

## Example figures

```@example
include("plots/BulkTendencies_plots.jl")
```

![](bulk_microphysics_linearized_convergence.svg)

The figure compares:

- a **nonlinear reference solution**, obtained using a finely substepped explicit integration
- the **linearized implicit method** with different numbers of substeps (`nsub`)
- a **single explicit update** using the instantaneous tendency at $t=0$
- **explicit updates**, using the instantaneous tendency with $10$ substeps 

### Initial conditions

-  $\rho = 1\,\mathrm{kg/m^3}$
-  $q_{\mathrm{tot}} = 13\,\mathrm{g/kg}$
-  $q_{\mathrm{lcl}} = q_{\mathrm{rai}} = 1\,\mathrm{g/kg}$
-  $q_{\mathrm{icl}} = q_{\mathrm{sno}} = 0.5\,\mathrm{g/kg}$
-  $T = 278\,\mathrm{K}$

These conditions activate multiple processes simultaneously (liquid, ice, rain, and snow interactions) and are close to freezing, making the case strongly nonlinear.

### Interpretation

- `nsub = 1` corresponds to a **single linearization over the full step**, which is the least accurate but cheapest approximation.
- Increasing `nsub` improves the solution by updating the linearization more frequently.
- For sufficiently large `nsub`, the solution approaches the nonlinear reference trajectory. Even `nsub = 2` agrees well with the nonlinear solution.
- The dashed line (instantaneous tendency) shows a simple explicit Euler step, which can significantly deviate from the true evolution.
- The yellow dash-dotted line shows an integration using instantaneous tendencies with 10 substeps and exhibits significant instabilities. Thus, without the linearized model, even 10 substeps do not converge.

This demonstrates that the linearized implicit substepping method provides a controllable trade-off between **cost and accuracy**, while maintaining stability.

---

## Current limitations

- The donor-based linearization, and hence `LinearizedAverage`, is implemented only for the one-moment scheme. The one-moment scheme also supports the donor, coupled-donor, and exact configurations of `RosenbrockAverage`, and the explicit `SubsteppedAverage`.
- On the two-moment + P3 model there is no donor-based matrix, so averaged tendencies use either the exact-Jacobian `RosenbrockAverage` (`rosenbrock_exact`) or the explicit `SubsteppedAverage`. The zero-moment scheme provides instantaneous tendencies only.

## Rosenbrock-averaged tendencies (2M+P3)

For the two-moment + P3 configuration the linearized-implicit average uses the exact Jacobian of the fused tendency, obtained by forward-mode automatic differentiation, because there is no donor-based matrix there. This is the `rosenbrock_exact` preset, `RosenbrockAverage` with `ExactJacobian`, `ExplicitGrowthDiagonal`, and `EndStateSaturationAdjustment`. The interval ``\Delta t`` is divided into `nsub` substeps of length ``h``, and each substep performs one linearized-implicit (Rosenbrock-Euler) update of the eight prognostic species ``x = (q_{\mathrm{lcl}}, n_{\mathrm{lcl}}, q_{\mathrm{rai}}, n_{\mathrm{rai}}, q_{\mathrm{ice}}, n_{\mathrm{ice}}, q_{\mathrm{rim}}, b_{\mathrm{rim}})`` of the form ``\left(I/h - P J P\right) \Delta x = f(x)`` followed by a positivity clamp, where ``f`` is the raw instantaneous tendency and ``J`` is its exact ``8\times8`` `ForwardDiff` Jacobian.

The end-state saturation-adjustment limiter is part of this preset and is required, not optional. The two-moment + P3 ice-deposition tendency carries an autocatalytic growth mode whose positive Jacobian diagonal makes the implicit operator lose positive-definiteness once the substep is coarse relative to the growth time scale; the single substep then overshoots the nonlinear saturation limit, which the linear operator does not see, and the state goes non-physical. `ExplicitGrowthDiagonal` removes that growth mode from the implicit operator and `EndStateSaturationAdjustment` supplies the missing nonlinear bound. Both are needed: zeroing the growth diagonal alone leaves the now-explicit growth unbounded at the coarsest single substep, so the saturation adjustment is what keeps it physical. This is why the L-stability of the one-stage linearized-implicit step does not by itself remove the limiter: that stability covers only the linear modes the operator captures, not the nonlinear vapor-depletion shutoff.

The full treatment — the substep semantics, the equilibrated linear solve, the species projection ``P``, the growth-diagonal and limiter options, the coarse-step deposition instability, the mixed-phase saturation criterion (``\max(S_{\mathrm{ice}}, S_{\mathrm{liq}})``), and the convergence behaviour under substep refinement — is in [Bulk-tendency averaging: modes and numerics](@ref). At a single substep the explicit growth is bounded only by the saturation adjustment, which over-produces precipitation at coarse time steps; two or more substeps recover accurate precipitation and the saturation adjustment is then rarely active.
