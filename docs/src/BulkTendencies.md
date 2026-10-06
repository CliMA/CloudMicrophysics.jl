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

- Vapor ↔ cloud condensate phase changes (condensation/evaporation of cloud liquid, deposition/sublimation of cloud ice)
  are relaxations toward equilibrium. The scheme computes $S = (q^\star - q)/\tau$ for a relaxation timescale $\tau$.
  $q^\star = q + S\tau$ already includes the latent-heat factor $\Gamma$ and the available-condensate bound.
  Their transfer over the substep is the time average of the relaxation,
    see [MorrisonMilbrandt2015](@cite) Appendix C.
  ```math
  \Delta q = S\,\tau\,\bigl(1 - e^{-\Delta t/\tau}\bigr)
           = S\,\Delta t\,\varphi(\Delta t/\tau), \qquad
  \varphi(x) = \frac{1 - e^{-x}}{x},
  ```
  The substep never crosses $q^\star$ for any $\Delta t/\tau$ and the
  instantaneous rate is recovered for $\Delta t \ll \tau$. A source
  ($\Delta q > 0$) is added to $e$ as a non-negative constant. A sink
  ($\Delta q < 0$) is added to $M$ as an implicit decay $-D\,q$ with
  $D = |\Delta q| / \bigl(\max(q + \Delta q, q_{\min})\,\Delta t\bigr)$, which
  removes exactly $|\Delta q|$ when acting alone (the $q_{\min}$ floor keeps $D$
  finite when the whole pool sublimates), keeps $q \ge 0$ when combined with
  the other sinks, and, unlike a plain $S/q$ decay, does not remove all the
  condensate when $\tau \ll \Delta t$. This matters when $\tau$ is a few
  seconds (e.g. `PrescribedIceNumber` with a large prescribed ice number
  concentration), where treating the instantaneous rate as a constant over the
  substep produced a deposition/sublimation flip-flop.

- Vapor → snow deposition is treated as a constant source (added to $e$)

- Other condensate sinks (snow sublimation, rain evaporation) are treated as
  linear sinks:
  ```math
  S \;\rightarrow\; -D q
  ```

With this formulation, sink terms take the form:

```math
\frac{dq}{dt} = -D q
```

which corresponds to exponential decay over the timestep, providing strong numerical stability.

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

## Vapor and latent heating limiters

Two quantities are not part of the linear solve. Vapor is not one of its four
unknowns, so the vapor budget of a substep is only known afterwards from
$q_v = q_{tot} - \sum_k q_k$, and so is the latent heating of the substep. Each
substep therefore does two linear solves. The first, with the transfers above,
gives the realized transfers $\Delta q^1_k$ of the four species and the heating

```math
\Delta T_1 = \frac{L_v \left(\Delta q^1_{lcl} + \Delta q^1_{rai}\right) + L_s \left(\Delta q^1_{icl} + \Delta q^1_{sno}\right)}{c_p},
```

with the reference latent heats and the dry-air specific heat of the substep
temperature update, so that the fusion heat of freezing and melting is included.
Two factors are derived from it, $\alpha$ for the vapor and $f$ for the heating,
and the system is solved a second time with them.

### Vapor limiter

A vapor source (condensation, deposition) relaxes the vapor toward the
saturation of its phase and stops there. With consistent transfers the sources
can therefore at most bring the vapor down to the lower of the two saturations,
evaluated at the end-of-substep temperature. The solved vapor can end below it
when a source was sized for vapor that is not there: the transfers are computed
from the start-of-substep state, and the evaporation or sublimation that would
have supplied vapor during the substep is clamped to its pool, or shared with
other sinks in the implicit solve, while the deposition computed from the same
state is not reduced.

With $q^\star_{\min}$ the lower of the two saturations, $\lambda_{\min} = dq^\star_{\min}/dT$
its slope, $\Gamma_{\min} = 1 + (L/c_p)\,\lambda_{\min}$ with the latent heat of that phase,
$e_k$ the vapor sources of the linear system and $q^1_v$ the vapor after the first solve,
the vapor shortfall and the vapor taken by the sources are

```math
\Delta q_{gap} = q^\star_{\min} + \lambda_{\min}\,\Delta T_1 - q^1_v, \qquad
\Delta q_{src} = \sum_k e_k\,\Delta t,
```

and all vapor sources are scaled by

```math
\alpha = \max\!\left(0,\; 1 - \frac{\Delta q_{gap}}{\Gamma_{\min}\,\Delta q_{src}}\right)
\quad \text{if } \Delta q_{gap} > 0 \text{ and } \Delta q_{src} > 0, \text{ else } \alpha = 1 .
```

Removing the fraction $1 - \alpha$ of the sources returns $(1 - \alpha)\,\Delta q_{src}$ of
vapor and lowers the saturation by $(\Gamma_{\min} - 1)(1 - \alpha)\,\Delta q_{src}$ through
the heating it removes, so with this $\alpha$ the vapor ends exactly on the lower
saturation at the corrected temperature. This replaces the cap of earlier versions,
which compared the sources with the excess at the start-of-substep temperature and
therefore lacked $\Gamma$. The limiter never adds vapor: with no sources there is
nothing to scale and subsaturated air is left to the evaporation and sublimation
sinks. Vapor above the upper saturation is not limited either; it is removed by the
sources of the next substep and poses no positivity risk.

### Latent heating limiter

With $\Delta T_\alpha$ the heating of the first solve after the vapor limiter,

```math
\Delta T_\alpha = \Delta T_1 - (1 - \alpha)\,\frac{L_v \left(e_{lcl} + e_{rai}\right) + L_s \left(e_{icl} + e_{sno}\right)}{c_p}\,\Delta t,
```

and $\Delta T_{\max}$ = `max_latent_heating_rate` × $\Delta t$ (ClimaParams
`microphysics_max_latent_heating_rate`, K/s, default 2 K per minute, `inf`
disables, must be positive),

```math
f = \min\!\left(1, \frac{\Delta T_{\max}}{|\Delta T_\alpha|}\right), \qquad
s_k = \frac{f\,(1 + D^{c}_k \Delta t)}{1 + D^{c}_k \Delta t + (1 - f)\,D^{p}_k \Delta t} ,
```

where $D^{p}_k$ is the sum of the decays of donor $k$ that change phase
(evaporation or sublimation, freezing, melting, riming) and $D^{c}_k$ the sum of
its collision and conversion decays. In the second solve the vapor sources are
scaled by $\alpha f$ and the phase-change decays of each donor by $s_k$, which
scales the realized transfer of a donor that is not refilled by another phase
change by exactly $f$; collision and conversion transfers are not scaled.

After the second solve, with $\Delta q^2_k$ its transfers, $\Delta T_2$ the
corresponding heating, $\Delta q_{cond} = \sum_k \Delta q^2_k$ the condensate growth
and $q^0_v = \max(0, q_{tot} - \sum_k q^0_k)$ the initial vapor,

```math
g_T = \min\!\left(1, \frac{\Delta T_{\max}}{|\Delta T_2|}\right), \qquad
g_v = \frac{q^0_v}{\Delta q_{cond}} \text{ if } \Delta q_{cond} > q^0_v, \text{ else } 1, \qquad
g = \min(g_T, g_v), \qquad
\frac{dq_k}{dt} = g\,\frac{\Delta q^2_k}{\Delta t} .
```

The uniform factor $g$ covers the cases the per-donor factors miss (pools that
refill each other, such as rain freezing on cloud ice while snow melts) and
bounds the condensate growth by the vapor available; heating and the vapor
budget are linear in the tendencies, so both bounds are met exactly. With the
implicit decays keeping every pool non-negative, no tracer can become negative
in a substep however stiff the processes are. Because the limiter bounds a
rate, it acts per substep: with fine substeps it can engage during the fast
initial part of a relaxation that is below the bound when averaged over a
longer substep. Condensate inputs are clamped to non-negative values before
the solve.

---

## Sparse 4×4 structure

The system has a fixed sparse structure:

```math
\begin{bmatrix}
a_{11} & a_{12} & 0      & 0 \\
a_{21} & a_{22} & 0      & 0 \\
a_{31} & 0      & a_{33} & a_{34} \\
a_{41} & a_{42} & a_{43} & a_{44}
\end{bmatrix}
```

This allows an efficient solve:

-  $q_{\mathrm{lcl}}$ and $q_{\mathrm{icl}}$ are solved as a coupled **2×2 system**
   (cloud ice melt via $a_{12}$, cloud liquid freezing via $a_{21}$)
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
