# Bulk Tendencies

## Linearized average tendencies

Microphysical source terms can be stiff, especially for depletion processes
such as evaporation, sublimation, and melting. To improve stability and
allow larger timesteps, we introduce a linearized implicit formulation
for computing time-averaged bulk tendencies.

The idea is to approximate the nonlinear microphysics tendencies as a linear system:

```math
\frac{dq}{dt} \approx M q + e
```

where $q = (q_{\mathrm{lcl}}, q_{\mathrm{icl}}, q_{\mathrm{rai}}, q_{\mathrm{sno}})$,
and the matrix $M$ and vector $e$ are constructed from the instantaneous tendencies.

### Donor-based linearization

Each microphysical process is linearized with respect to its donor species:

- Vapor-driven phase changes (condensation/evaporation of cloud liquid,
  evaporation of rain, deposition/sublimation of cloud ice and of snow) are
  relaxations of the vapor excess:
  ```math
  \delta = q_v - q^\star
  ```
  where $q^\star$ is the saturation specific humidity over liquid or ice.
  Instantaneous rate of vapor transfer can be written as
  ```math
  S = c\,\delta
  ```
  For the cloud formation processes
  ``math
  c = \frac{1}{\tau \, \Gamma}
  ```
  where
  $\tau$ is the process timescale, and
  $\Gamma = 1 + (L/c_p)\,dq^\star/dT$ accounts for the latent heat moving the saturation.
  For rain and snow processes
  ```math
  c = \frac{S}{\delta}
  ```
  whose 1-moment rates are linear in the vapor excess.

  For one process acting alone the vapor excess decays as
  ```math
  d\frac{\delta}{dt} = -\Gamma c\,\delta = -\frac{\delta}{\tau}.
  ```
  The transfer over the substep is the time average of that relaxation
  ([MorrisonMilbrandt2015](@cite), Appendix C):
  ```math
  \Delta q = S\,\tau\,\bigl(1 - e^{-\Delta t/\tau}\bigr)
           = S\,\Delta t\,\varphi(\Delta t/\tau), \qquad
  \varphi(x) = \frac{1 - e^{-x}}{x}.
  ```
  This formulation ensures that the solution does not cross the equilibrium for any $\Delta t/\tau$
  and equals $S\,\Delta t$ for $\Delta t \ll \tau$.

  With several processes acting on the same vapor reservoir,
  each process changes the excess the others see, through the vapor it
  takes and through its latent heat.
  The excess of the primary phase $p$ then
  follows a relaxation with a combined rate and a constant drive,
  ```math
  \frac{d\delta_p}{dt} = A - \frac{\delta_p}{\tau}, \qquad
  \bar\delta_p = A\tau + (\delta_{p,0} - A\tau)\,\varphi(\Delta t/\tau), \qquad
  \Delta q_k = c_k\,\bar\delta_k\,\Delta t,
  ```
  with $\tau$, $A$ and the excess of the other phase given in the section on
  the joint relaxation below; for a single process $A = 0$ and
  $1/\tau = \Gamma c$, which is the first formula. Sinks are clamped to their
  pools. A source ($\Delta q > 0$) is added to $e$ as a constant. A sink
  ($\Delta q < 0$) is added to $M$ as an implicit decay $-D\,q$ with
  $D = |\Delta q| / \bigl(\max(q + \Delta q, q_{\min})\,\Delta t\bigr)$, which
  removes exactly $|\Delta q|$ when acting alone, keeps $q \ge 0$ together
  with the other sinks of the pool and shares the pool among them, and, unlike
  a plain $S/q$ decay, does not empty the pool when $\tau \ll \Delta t$.

- Transfer processes (e.g. accretion, conversion):
  ```math
  S \;\rightarrow\; D \, q_{\text{donor}}, \quad D = \frac{S}{\max(\epsilon, q_{\text{donor}})}
  ```

- Fusion transfers (freezing and melting of cloud condensate, riming, the
  freeze/melt arms of the accretion processes, melting of snow) are donor decays
  like the collision/conversion transfers above; their latent heat enters through
  the substep temperature update and the latent-heating limiter.

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

## Joint Γ-consistent relaxation

Condensation on cloud liquid, evaporation of rain and deposition on cloud ice
and snow draw on the same vapor, and the latent heat of each moves the
saturation the others relax toward. Relaxing each process from the initial
excess on its own over-deposits: in the glaciating updraft that motivated this
change, three fast processes deposited about twice the consistent amount and
sublimated it back the next substep.

Let $\delta_l = q_v - q^\star_l$, $\delta_i = q_v - q^\star_i$,
$\Delta s = q^\star_l - q^\star_i$, $C_l = c_{lcl} + c_{rai}$,
$C_i = c_{icl} + c_{sno}$, $\Gamma_p = 1 + (L_p/c_p)\,dq^\star_p/dT$ and
$\kappa_{sp} = 1 + (L_s/c_p)\,dq^\star_p/dT$, the latent heat of phase $s$
acting on the saturation of phase $p$. The vapor uptake of all four processes
and their heating give

```math
\frac{d\delta_l}{dt} = -\Gamma_l C_l\,\delta_l - \kappa_{il} C_i\,\delta_i, \qquad
\frac{d\delta_i}{dt} = -\Gamma_i C_i\,\delta_i - \kappa_{li} C_l\,\delta_l ,
```

with $\Delta s$, $\Gamma$, $\kappa$ and the coefficients held at their
start-of-substep values, so that $\delta_i = \delta_l + \Delta s$ throughout
and one equation suffices. The primary phase $p$ is the one with the larger
own decay rate $\Gamma_p C_p$ (a single active process is then integrated
exactly), $s$ is the other, $\sigma = +1$ if liquid is primary and $-1$ if
ice is. Substituting $\delta_s = \delta_p + \sigma\Delta s$,

```math
\frac{d\delta_p}{dt} = A - \frac{\delta_p}{\tau}, \qquad
\frac{1}{\tau} = \Gamma_p C_p + \kappa_{sp} C_s, \qquad
A = -\sigma\,\Delta s\,\kappa_{sp} C_s ,
```

where $A$ is the Wegener-Bergeron-Findeisen drive (liquid evaporating while
ice deposits). The exact average over the substep,
$\bar\delta_p = A\tau + (\delta_{p,0} - A\tau)\,\varphi(\Delta t/\tau)$, gives
$\bar\delta_s = \bar\delta_p + \sigma\Delta s$ and the transfers
$\Delta q_k = c_k \bar\delta_k \Delta t$, sinks clamped to their pools. The
average never crosses the equilibrium $A\tau$, so the transfers cannot
overshoot the joint saturation for any $\Delta t/\tau$.

Pools are not tracked within the substep (as in P3): the implicit solve shares
each pool among its sinks and the vapor check below removes deposition that an
exhausted pool could not feed. Consequences: the residue left after a complete
evaporation still counts as a small vapor supply, so the air can stay a few
per cent subsaturated over ice while the residue persists; a small cloud in dry
air evaporates on the excess relaxation time rather than at the pool-bounded
rate of the instantaneous scheme; and the split of an uptake between cloud ice
and snow within one substep can differ from the finely resolved one while
vapor and heating agree.

## Vapor check and latent heating limiter

Each substep does two linear solves. The first, with the transfers above,
gives the realized transfers $\Delta q^1$ and the heating
$\Delta T_1 = (L_v\,\Delta q^1_{\mathrm{liq}} + L_s\,\Delta q^1_{\mathrm{ice}})/c_p$,
with the same latent-heat factors as the substep temperature update, so the
fusion heat of freezing and melting is included. Two factors are derived from
it and applied in the second solve.

Vapor check. With $q^\star_{\min}$ the lower of the two saturations,
$\lambda_{\min} = dq^\star_{\min}/dT$ and $\Sigma e = \sum_k e_k \Delta t$ the
vapor taken by the sources,

```math
\mathrm{gap} = q^\star_{\min} + \lambda_{\min}\Delta T_1 - q^1_v, \qquad
\alpha = \max\!\left(0,\; 1 - \frac{\mathrm{gap}}{\Gamma_{\min}\,\Sigma e}\right)
\quad \text{if } \mathrm{gap} > 0, \text{ else } \alpha = 1 .
```

Scaling the vapor sources by $\alpha$ raises the vapor by $(1-\alpha)\Sigma e$
and lowers the floor by $(\Gamma_{\min} - 1)(1-\alpha)\Sigma e$, so this
$\alpha$ is the closed-form solution of $q_v(\alpha) = \text{floor}(\alpha)$.
It removes the deposition that a pool exhausted within the substep (for
example liquid taken by freezing) could not feed, and replaces the former
raw-excess cap, which lacked $\Gamma$.

Latent heating limiter. With $\Delta T_\alpha$ the heating after the vapor
check and $\Delta T_{\max}$ = `max_latent_heating_rate` × $\Delta t$
(ClimaParams `microphysics_max_latent_heating_rate`, K/s, default 2 K per
minute, `inf` disables, must be positive),

```math
f = \min\!\left(1, \frac{\Delta T_{\max}}{|\Delta T_\alpha|}\right), \qquad
s_k = \frac{f\,(1 + D^{c}_k \Delta t)}{1 + D^{c}_k \Delta t + (1 - f)\,D^{p}_k \Delta t} ,
```

where $D^{p}_k$ is the sum of the phase-change decays of donor $k$ and
$D^{c}_k$ the sum of its collision/conversion decays. The vapor sources are
scaled by $f$ and each donor's phase-change decays by $s_k$, which scales the
realized transfer of a donor that is not refilled by another phase change by
exactly $f$; collision and conversion transfers are not scaled.

After the second solve, with $\Delta q^2$ its transfers, $\Delta T_2$ the
corresponding heating, $\Delta q_{\mathrm{cond}} = \sum_k \Delta q^2_k$ and
$q^0_v$ the initial vapor,

```math
g = \min\!\left(1, \frac{\Delta T_{\max}}{|\Delta T_2|}, \frac{q^0_v}{\Delta q_{\mathrm{cond}}}\right), \qquad
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
longer substep.

For comparison runs, `Microphysics1MParams.joint_vapor_relaxation = false`
(optional ClimaParams entry `microphysics_joint_vapor_relaxation`, type
`bool`, `true` when absent) replaces the joint transfers by the single-process
formula for each process ($S\,\Delta t\,\varphi(\Delta t/\tau)$ for the cloud
condensates, $S\,\Delta t$ for rain and snow) under the same vapor check,
limiter and positivity guard; with several fast processes on the same excess
that branch shows the substep alternation described above.

Accuracy with `PrescribedIceNumber` and $N_0 = 5\times10^8$
($\tau_{\mathrm{ice}} \ll \Delta t$): with three substeps a glaciating
updraft, a Wegener-Bergeron-Findeisen cloud and a liquid pool that evaporates,
freezes and feeds deposition end on ice saturation within 1 % with the step
heating within 1 % of a finely substepped reference; a single 40-180 s
substep stays within 3 % of saturation, and when a pool is exhausted within
it part of the heating is deferred to the next substep. When a coefficient
changes during the substep (cloud ice growing from a very small pool while the
liquid holds the vapor at liquid saturation) the error is that of holding it
at its initial value, about 10 % of the step's deposition per 60 s substep in
the test case. The default limiter clips the first 40 s substep of the
glaciating updraft slightly (factor 0.97); the deferred heat is realized in
the second substep and the step total is unchanged to 1e-4 K.

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
