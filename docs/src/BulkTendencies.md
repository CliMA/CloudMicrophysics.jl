# Bulk Tendencies

## Linearized average tendencies

The microphysics tendency of the condensate species $q = (q_{\mathrm{lcl}}, q_{\mathrm{icl}}, q_{\mathrm{rai}}, q_{\mathrm{sno}})$ is

```math
\frac{dq}{dt} = f(q),
```

where $f$ is the sum of the instantaneous rates $S_p$ of the individual processes $p$.

Some of these rates are stiff, especially those of depletion processes such as evaporation, sublimation and melting.
To improve stability and allow larger time steps, the tendency is averaged over the model time step with a linearized implicit formulation, which approximates $f(q)$ locally by a linearized tendency:

```math
\frac{dq}{dt} \approx M q + e,
```

with the matrix $M$ and the vector $e$ constructed from the rates $S_p$.

The model time step is divided into `nsub` equal substeps.
Each substep:

1. computes the rates $S_p$ at the current state and builds $M$ and $e$ from them by the donor-based linearization;
2. scales the vapor sources in $e$ so that the substep does not reduce the vapor below saturation;
3. takes a backward Euler step of the linearized tendency for $q$;
4. updates the temperature from the latent heat of the change in $q$.

The average tendency is the net change of $q$ over the model time step divided by its length.
The following sections describe each part.

### Donor-based linearization

Each process moves mass from a donor species, the species that loses mass, to a receiving species at the rate $S_p$.
The linearization writes each rate according to its donor and receiver:

| Donor      | Receiver   | Linearized rate       | Entries                                                                            |
| ---------- | ---------- | --------------------- | ---------------------------------------------------------------------------------- |
| condensate | condensate | $D\,q_{\text{donor}}$ | $-D$ in $M_{\text{donor},\text{donor}}$, $D$ in $M_{\text{receiver},\text{donor}}$ |
| condensate | vapor      | $D\,q_{\text{donor}}$ | $-D$ in $M_{\text{donor},\text{donor}}$                                             |
| vapor      | condensate | $S_p$, constant       | $S_p$ in $e_{\text{receiver}}$                                                      |

Here $D \ge 0$ is the decay coefficient of the donor.
The vapor is not part of $q$, so a loss to the vapor has no entry for its receiver, and a gain from the vapor enters $e$.

A sink of a species therefore takes the form

```math
\frac{dq}{dt} = -D q,
```

which corresponds to an exponential decay over the substep and keeps the species non-negative in the implicit step.

The decay coefficients $D$ and the sources in $e$ are computed from the state at the start of each substep and are held constant over it.

### Kinds of processes

The 1-moment processes are of three kinds.

- Transfers between two condensate species, such as accretion, autoconversion, melting and freezing.
  A transfer from $q_{\text{donor}}$ to $q_{\text{receiver}}$ at the rate $S \ge 0$ contributes
  ```math
  \frac{dq_{\text{donor}}}{dt} = -S, \qquad \frac{dq_{\text{receiver}}}{dt} = S,
  ```
  and its donor is $q_{\text{donor}}$, with
  ```math
  D = \frac{S}{\max(q_{\min}, q_{\text{donor}})}.
  ```

- Exchanges between the vapor and a condensate species $q_{\text{cond}}$ in either direction, such as rain evaporation and snow deposition and sublimation.
  An exchange at the rate $S$, positive from the vapor to $q_{\text{cond}}$, contributes $dq_{\text{cond}}/dt = S$, and the vapor changes by $-S$.
  A gain ($S \ge 0$) has the vapor as its donor and adds $S$ to $e_{\text{cond}}$.
  A loss has $q_{\text{cond}}$ as its donor, with
  ```math
  D = \frac{-S}{\max(q_{\min}, q_{\text{cond}})}.
  ```

- Relaxations of the cloud condensate toward equilibrium: the condensation and evaporation of cloud liquid, and the deposition and sublimation of cloud ice.
  A relaxation contributes the same tendency as an exchange, but enters the linearized tendency through its transfer over the substep, as described in the next section.

Here $q_{\min}$ is a small positive floor that keeps $D$ finite for a vanishing donor.

### Relaxation of the cloud condensate

The scheme computes the rate of a relaxation as $S = (q^\star - q_{\text{cond}})/\tau$, with the relaxation timescale $\tau$, where $q^\star = q_{\text{cond}} + S\tau$ is the condensate that the relaxation approaches.
The rate $S$ includes two corrections, so $q^\star$ includes them as well: the latent heat of the phase change reduces the supersaturation that relaxes, through the factor $\Gamma = 1 + (L/c_p)\,dq_{\mathrm{sat}}/dT$, and evaporation or sublimation cannot remove more condensate than exists, so $q^\star \ge 0$.

The transfer over the substep is the time average of the relaxation, see [MorrisonMilbrandt2015](@cite) Appendix C:

```math
\Delta q = S\,\tau\,\bigl(1 - e^{-\Delta t/\tau}\bigr)
         = S\,\Delta t\,\varphi(\Delta t/\tau), \qquad
\varphi(x) = \frac{1 - e^{-x}}{x}.
```

The substep never crosses $q^\star$ for any $\Delta t/\tau$, and the instantaneous rate is recovered for $\Delta t \ll \tau$.
A gain ($\Delta q \ge 0$) adds $\Delta q / \Delta t$ to $e_{\text{cond}}$.
A loss ($\Delta q < 0$) adds $-D$ to $M_{\text{cond},\text{cond}}$, with

```math
D = \frac{|\Delta q|}{\max(q_{\text{cond}} + \Delta q, q_{\min})\,\Delta t}.
```

This decay removes exactly $|\Delta q|$ when acting alone (the $q_{\min}$ floor keeps $D$ finite when the whole pool sublimates), keeps $q \ge 0$ when combined with the other sinks, and, unlike a plain $S/q$ decay, does not remove all the condensate when $\tau \ll \Delta t$.
This matters when $\tau$ is a few seconds (e.g. `PrescribedIceNumber` with a large prescribed ice number concentration), where treating the instantaneous rate as a constant over the substep produced a deposition/sublimation flip-flop.

---

## Linearized implicit solve

A substep of width $\Delta t$ starts from the species $q^n$ and solves the linearized tendency implicitly for the species $q^{n+1}$ at its end:

```math
\frac{q^{n+1} - q^n}{\Delta t} = M q^{n+1} + \alpha e,
```

where $\alpha \le 1$ scales the vapor sources, as described in the next section.
This gives

```math
\left(I/\Delta t - M\right) q^{n+1} = \alpha e + q^n/\Delta t,
```

and the average tendency over the substep is

```math
\overline{f} = \frac{q^{n+1} - q^n}{\Delta t}.
```

Over the substep, a transfer moves mass at the rate $D\,q^{n+1}_{\text{donor}}$, and an exchange or a relaxation at the rate $\alpha s - D\,q^{n+1}_{\text{cond}}$, where $s$ is its entry in $e$.
These rates sum to the average tendency, and `LinearizedAverageVerbose` returns their average over the substeps.

---

## Vapor-budget limit on the vapor sources

The vapor sources in $e$, from condensation on cloud liquid and deposition on cloud ice and on snow, together consume vapor over the substep.
If their combined rate is large enough, an unlimited substep can reduce the vapor $q_v$ below saturation, or even below zero.
To prevent this, all three sources are scaled by the same factor

```math
\alpha = \min\!\left(1,\; \frac{\max(0,\, q_v - q_{\mathrm{sat},\min})}
                                 {\Delta t\;(e_{\mathrm{lcl}} + e_{\mathrm{icl}} + e_{\mathrm{sno}})}\right),
\qquad
q_{\mathrm{sat},\min} = \min\!\bigl(q_{\mathrm{sat,liq}}, q_{\mathrm{sat,ice}}\bigr),
```

where $q_{\mathrm{sat,liq}}$ and $q_{\mathrm{sat,ice}}$ are the saturation specific contents over liquid and over ice.
Then $q_v$ does not fall below $q_{\mathrm{sat},\min}$ over one substep.
The common factor keeps the relative rates of the three processes unchanged, and the sinks in $M$ are unaffected.

- Below freezing, $q_{\mathrm{sat},\min} = q_{\mathrm{sat,ice}}$, which allows the Bergeron process to reduce $q_v$ below liquid saturation.
- Above freezing, $q_{\mathrm{sat},\min} = q_{\mathrm{sat,liq}}$, so the limit is liquid saturation.


---

## Sparse 4×4 structure

The matrix $A = I/\Delta t - M$ of the implicit solve has a fixed sparse structure,

```math
A =
\left[
\begin{array}{cc|cc}
A_{\mathrm{lcl},\mathrm{lcl}} & A_{\mathrm{lcl},\mathrm{icl}} & 0 & 0 \\
A_{\mathrm{icl},\mathrm{lcl}} & A_{\mathrm{icl},\mathrm{icl}} & 0 & 0 \\
\hline
A_{\mathrm{rai},\mathrm{lcl}} & 0 & A_{\mathrm{rai},\mathrm{rai}} & A_{\mathrm{rai},\mathrm{sno}} \\
A_{\mathrm{sno},\mathrm{lcl}} & A_{\mathrm{sno},\mathrm{icl}} & A_{\mathrm{sno},\mathrm{rai}} & A_{\mathrm{sno},\mathrm{sno}}
\end{array}
\right].
```

Melting and freezing couple the two cloud species, and the two precipitation species.
Autoconversion and accretion move mass from the cloud to the precipitation species.
No process transfers mass from a precipitation species to a cloud species, so the upper right block is zero.
This allows an efficient solve:

-  $q_{\mathrm{lcl}}$ and $q_{\mathrm{icl}}$ are solved as a coupled **2×2 system**
-  $q_{\mathrm{rai}}$ and $q_{\mathrm{sno}}$ are then solved as a **2×2 system**, with the new cloud species as sources

This avoids forming or inverting a full dense matrix and is efficient on both CPU and GPU.

---

## Substepping

A single linearization assumes the operator $M$ is constant over the model time step.
Substeps rebuild $M$ and $e$ from the updated state, which captures nonlinear effects and regime changes, for example near freezing.
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

## Processes in the code

`_microphysics_source_terms` returns each 1-moment process as a term of its kind:

| Kind                    | Constructor                           | Example |
| ----------------------- | ------------------------------------- | ------- |
| Transfer                | `Transfer(:Donor => :Receiver, S)`    | `Transfer(:q_lcl => :q_rai, S)`, autoconversion |
| Exchange with the vapor | `VaporExchange(:Condensate, S)`       | `VaporExchange(:q_sno, S)`, deposition on snow  |
| Relaxation              | `VaporRelaxation(:Condensate, S, τ)`  | `VaporRelaxation(:q_lcl, S, τ)`, condensation   |

To add a process to the 1-moment scheme, compute its rate in `Microphysics1M` and add its term to `_microphysics_source_terms`.
The instantaneous tendencies, the entries of $M$ and $e$, and the rate of the process in the verbose modes follow from the term.
A transfer from a precipitation species to a cloud species is not supported, because the block solve assumes that the upper right block of $M$ is zero.

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
