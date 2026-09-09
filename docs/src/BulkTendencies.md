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

- Vapor ↔ cloud condensate phase changes (condensation/evaporation of cloud
  liquid, deposition/sublimation of cloud ice) are treated as **implicit
  relaxations** toward their equilibrium. The non-equilibrium schemes compute
  $S = (q^\star - q)/\tau$ for a relaxation timescale $\tau$ (returned by
  `τ_vap_to_q_lcl` / `τ_vap_to_q_icl`; $q^\star = q + S\tau$ already includes
  the latent-heat factor $\Gamma$ and the available-condensate bound), so we write
  ```math
  \frac{dq}{dt} = S - \frac{q^{new} - q}{\tau}
  ```
  i.e. $-1/\tau$ on the diagonal of $M$, the drive $S$ (either sign) in $e$, and
  a hold term $h = q/\tau$ on the right-hand side. The substep result
  $\Delta q = S\,\Delta t/(1 + \Delta t/\tau)$ can never overshoot $q^\star$,
  for any $\Delta t/\tau$, and reduces to a plain explicit source when
  $\Delta t \ll \tau$. This matters when $\tau$ is a few seconds (e.g.
  `PrescribedIceNumber` with a large prescribed ice number concentration), where
  treating the source as a constant produced a deposition/sublimation flip-flop.

- Vapor → snow deposition is treated as a **constant source** (added to $e$)

- Other condensate sinks (snow sublimation, rain evaporation) are treated as
  **linear sinks**:
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
\frac{q^\star - q^0}{\Delta t} = M q^\star + e + h
```

which gives:

```math
\left(I/\Delta t - M\right) q^\star = e + h + q^0/\Delta t
```

The average tendency is then:

```math
\overline{T} = \frac{q^\star - q^0}{\Delta t}
```

---

## Vapor-budget cap on vapor → condensate sources

Vapor → condensate processes (condensation on cloud liquid, deposition on
cloud ice, deposition on snow — the positive contributions to `e_1`, `e_2`,
`e_4`) together consume vapor. If their combined rate is fast enough, an
unlimited substep can drive `q_v` below saturation or even negative. To
prevent this, the positive drives are uniformly scaled by

```math
\alpha = \min\!\left(1,\; \frac{\max(0,\, q_v - q^\star_{\min})}
                                 {d_1 + d_2 + d_4}\right),
\qquad
q^\star_{\min} = \min\!\bigl(q^\star_{\mathrm{liq}}, q^\star_{\mathrm{ice}}\bigr),
```

where $d_i$ is the vapor each source would consume over the substep as
realized by the solver: $d_i = e_i^+/(1/\Delta t + 1/\tau_i)$ for the implicitly
relaxed cloud terms (which equals $e_i^+ \Delta t$ when $\tau_i \to \infty$) and
$d_4 = e_4 \Delta t$ for snow. This keeps `q_v` from being driven below
`q^\star_{\min}` over one substep while preserving the relative rates of the
three processes. Negative drives (evaporation/sublimation), the hold terms
`h`, and the sinks (`M` blocks) are unaffected.

- Below freezing: `q^\star_{\min} = q^\star_{\mathrm{ice}}`, the natural
  ice-saturation floor (permits the Bergeron process to drive `q_v` below
  liquid saturation).
- Above freezing: `q^\star_{\min} = q^\star_{\mathrm{liq}}`, so the liquid
  floor is enforced (mathematical `q^\star_{\mathrm{ice}}` is unphysical
  there).

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

This allows an efficient solve:

-  $q_{\mathrm{lcl}}$ and $q_{\mathrm{icl}}$ are solved as an upper-triangular
   **2×2 system** (coupled through cloud ice melt via $a_{12}$)
-  $q_{\mathrm{rai}}$ and $q_{\mathrm{sno}}$ are solved as a **2×2 system**

This avoids forming or inverting a full dense matrix and is efficient on both CPU and GPU.

---

## Substepping

A single linearization assumes the operator $M$ is constant over the timestep. To better capture nonlinear effects and regime changes (e.g. near freezing), we apply **substepping**:

- Split the timestep into `nsub` substeps
- At each substep:
  - rebuild $M$ and $e$ from the updated state
  - solve the linearized system
  - update $q$ and temperature

As `nsub` increases, the solution approaches the nonlinear evolution of the system.

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

Here ``L_v`` and ``L_s`` are the constant reference latent heats
  (at the thermodynamic reference temperature) and ``c_p`` is the
  dry-air specific heat capacity.
This is consistent with the microphysics-only update and avoids coupling to a full thermodynamic solve.

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
-  $T = 278.15\,\mathrm{K}$

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

- Average (implicit) bulk tendencies are currently implemented **only for the one-moment microphysics scheme**.
- For other microphysics schemes, only **instantaneous bulk tendencies** are available at present.
