# Bulk Tendencies

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

## Linearization

Each microphysical process is linearized with respect to its donor species:

- Vapor-driven phase changes (condensation/evaporation of cloud liquid,
  evaporation of rain, deposition/sublimation of cloud ice and of snow) are
  based on ([MorrisonMilbrandt2015](@cite), Appendix C).
  They are written as relaxations of the vapor excess:
  ```math
  \delta = q_v - q^\star
  ```
  where $q^\star$ is the saturation specific humidity over liquid or ice.

  For a single process acting alone, the vapor excess decays as
  ```math
  \frac{d \delta}{dt} = \frac{dq_v}{dt} − \left(\frac{dq^\star}{dT} \right) \frac{dT}{dt}
  ```

  Taking $dq_v/dt = −S$ and $dT/dt = (L/c_p) S$ results in
  ```math
  \frac{d \delta}{dt} = −\left(1 + \frac{L}{c_p} \frac{d q^\star}{dT} \right) S = − \Gamma S = - \frac{\delta}{\tau}.
  ```
  where $L$ is the latent heat and $c_{p}$ is the specific heat.
  For cloud liquid and ice to vapor transfers
  the functional form of the rate is $S = \frac{\delta}{\tau \Gamma}$
  We can define $c$ as the process coefficient:
  ```math
  c_{lcl} = \frac{1}{\tau_{lcl} \Gamma_l}, \qquad
  c_{icl} = \frac{1}{\tau_{icl} \Gamma_i}.
  ```
  For rain and snow $\tau = \frac{1}{\Gamma c}$ is the timescale implied by their rate,
  ```math
  c_{rai} = \frac{S_{rai}}{\delta_l}, \qquad
  c_{sno} = \frac{S_{sno}}{\delta_i},
  ```
  and is evaluated once per substep.

  The solution of the differential equation for vapor excess is
  ```math
  \delta(t) = \delta_0 e^{−t/\tau}.
  ```
  The average saturation excess over the substep $\bar\delta$ and, for example,
  the mass transfer of cloud liquid water $\Delta q_{lcl}$ are
  ```math
  \bar\delta = \frac{1}{\Delta t}\int_0^{\Delta t} \delta(t)\,dt = \delta_0\,\varphi(\Delta t/\tau), \qquad
  \Delta q_{lcl} = c_{lcl}\,\bar\delta\,\Delta t, \qquad
  \varphi(x) = \frac{1 - e^{-x}}{x}.
  ```
  This formulation ensures that the solution does not cross the equilibrium for any $\Delta t/\tau$
    and equals $S\,\Delta t$ for $\Delta t \ll \tau$.

  When several processes are acting on the same vapor reservoir,
  each process changes the excess the others see, through both vapor it takes and latent heating.
  Because saturation vapor pressures over ice and liquid are different,
  there are two vapor excesses and two states the model is relaxing towards.
  ```math
  \frac{d \delta_l}{dt} = −\Gamma_l c_l \delta_l − \Gamma_{il} c_i \delta_i, \qquad
  \frac{d \delta_i}{dt} = −\Gamma_i c_i \delta_i − \Gamma_{li} c_l \delta_l
  ```
  where $\Gamma_{il} = 1 + (L_s/c_p) dq^\star_l / dT$ represents the
  latent heating of the deposition/sublimation acting on the saturation over liquid,
  $\Gamma_{li} = 1 + (L_v/c_p) dq^\star_i / dT$ represents the
  latent heating of vaporization acting on the saturation over ice, and
  $c_l = c_{lcl} + c_{rai}$ and $c_i = c_{icl} + c_{sno}$ represent
  the sum of all liquid and ice process coefficients.

  The two saturation excesses are separated by
  ```math
  \Delta s = q^\star_l − q^\star_i
  ```
  The $\Delta s$ is positive below the triple point and negative above.
  ```math
  \delta_i = \delta_l + \Delta s
  ```

  We can reorder this as
  ```math
  \frac{d \delta_l}{dt} = − \Gamma_{il} c_i \Delta s −\left(\Gamma_l c_l + \Gamma_{il} c_i \right) \delta_l = A - \frac{\delta_l}{\tau}
  ```
  where $A$ represents the constant forcing that the ice phase processes exert on the vapor excess over liquid,
  and $\tau$ now represents the multi-process timescale.
  The average saturation excess over liquid during the substep $\bar\delta_{l}$ and,
  for example, the mass transfer $\Delta q_{lcl}$ are
  ```math
  \bar\delta_l = A\tau + (\delta_{l,0} - A\tau)\,\varphi(\Delta t/\tau), \qquad
  \Delta q_{lcl} = c_{lcl}\,\bar\delta_l\,\Delta t,
  ```
  A similar equation can be written for the evolution of vapor excess over ice,
  and other microphysics tracers.
  ```math
  \bar\delta_i = \bar\delta_l + \Delta s, \qquad
  \Delta q_{rai} = c_{rai} \bar\delta_l \Delta t, \qquad
  \Delta q_{icl} = c_{icl} \bar\delta_i \Delta t, \qquad
  \Delta q_{sno} = c_{sno} \bar\delta_i \Delta t.
  ```
  With one process $A = 0$ and $1/\tau = \Gamma c$, and the $\bar\delta_l$
  equation reduces to pure liquid process described above.
  The CloudMicrophysics solver integrates the equation for the phase with the larger $\Gamma c$,
  and computes the other one based on $\Delta s$.

  The four vapor transfers ($\Delta q_{lcl}$, $\Delta q_{icl}$, $\Delta q_{rai}$, $\Delta q_{sno}$)
  enter the linear system as follows:
  A source ($\Delta q > 0$) is added to $e$ as a constant.
  A sink ($\Delta q < 0$) is limited by the available tracer amount, and
  is added to $M$ as an implicit decay $-D\,q$ with
  $D = |\Delta q| / \bigl(\max(q + \Delta q, q_{\min})\,\Delta t\bigr)$.
  This removes exactly $|\Delta q|$ when acting alone, keeps $q \ge 0$ together
  with the other sinks, and, unlike
  a plain $S/q$ decay, does not empty the available tracer pool when $\tau \ll \Delta t$.

- Transfer processes (accretion, conversion) and fusion transfers
  (freezing, melting, riming, the freeze/melt parts of accretion, snow melt)
  are linearized in their donor species
  ```math
  S \;\rightarrow\; D\, q_{\text{donor}}, \qquad D = \frac{S}{\max(\epsilon, q_{\text{donor}})}
  ```
  The linearization supplies the coefficient to the implicit solve.
  The implicit solve then gives the transfer as $\Delta q = −D q^{n+1} \Delta t$
  with $q^{n+1}$ the end-of-substep donor content.
  For a process acting alone on its pool,
  ```math
  \Delta q = -\frac{S\,\Delta t}{1 + S\,\Delta t / q^n}
  ```
  which never removes more than the pool $q^n}$ and tends to $−S \Delta t$ when $S \Delta t \ll q^n$.
  This is a backward-Euler decay $1 / (1 + D \Delta t)$,
  not the exponential decay $e^{−D\Delta t}$ of the vapor relaxation.
  Therefre the two $\Delta q$ formulas differ in form, but are both bounded and monotone.
  Several sinks on one pool share it in proportion to their D. The receiving species gains $\Delta q$.
  Fusion latent heat enters through the substep temperature update.

In short, the CloudMicrophysics library provides $S$ rates of individual processes.
The actual realized transfers are not those rates times $\Delta t$.
The $\varphi$-averaging of the vapor processes, and the implicit solve
for the donor decays replace $S \Delta t$ by bounded quantities.

---

## Linearized implicit solve

For a timestep $\Delta t$, we solve the linearized system implicitly:

```math
\frac{q^{n+1} - q^n}{\Delta t} = M q^{n+1} + e
```
where $q^n$ and $q^{n+1}$ are tracer values at the beginning and end of the substep.
This results in

```math
\left(I/\Delta t - M\right) q^{n+1} = e + q^n/\Delta t
```

The average tendency is then:

```math
\overline{T} = \frac{q^{n+1} - q^n}{\Delta t}
```

---

## Vapor and latent heating limiters

Additional considerations:

- Vapor is not one of the prognostic variables.
  When a vapor source (for example from evaporating cloud) exceeds
  the available donor pool (i.e. the available cloud water),
  the sink of cloud water is clamped, but the implied water vapor source is not.
  Within a substep the joint relaxation knows the coefficients $c$
  and the initial excesses, but not how much water each donor pool holds while the substep runs.
  Cloud liquid evaporating is treated as able to supply vapor at the rate $c_{lcl} \delta_l$
  for the whole substep, whether or not there is enough liquid.
  The pool size enters only at the end, as the clamp of the transfer to minus the pool for the donor,
  and in the implicit solve, which shares each pool among its sinks.
  This leads to an imbalance between how much vapor the solver though was available and provided
  to the other phase changes, and how much was actually depleted from the donor.
  This imbalance residue is then split between the remaining phase changes.

- Additionally, for the host model stability, one may want to limit the total
  amount of heating the microphysics can provide.

To address those issues each substep does two linear solves.
The first solve gives the transfers of mass $\Delta q^1$ and the heating
$\Delta T_1 = (L_v\,\Delta q^1_{\mathrm{liq}} + L_s\,\Delta q^1_{\mathrm{ice}})/c_p$
Two factors are derived from it and applied in the second solve
$\alpha$ - to address the vapor inconsistency and `f` to allow for heating limiters.

### Vapor limiter

With $q^\star_{\min}$ the lower of the two saturations,
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

### Latent heating limiter

With $\Delta T_\alpha$ the heating after the vapor
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
$q^n_v$ the initial vapor,

```math
g = \min\!\left(1, \frac{\Delta T_{\max}}{|\Delta T_2|}, \frac{q^0_v}{\Delta q_{\mathrm{cond}}}\right), \qquad
\frac{dq_k}{dt} = g\,\frac{\Delta q^2_k}{\Delta t} .
```

- In the limiter section Δq^1 and Δq^2 index the two solves. One sentence there, "q^{n+1} is the state after the second solve scaled by g", ties the two notations together.

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
