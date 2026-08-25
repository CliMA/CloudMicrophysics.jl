# Rain number as a relaxation

[SeifertBeheng2006](@cite) writes raindrop self-collection and breakup as a sink and a
correction to it. This page shows that the pair is *already* a relaxation of the raindrop
number toward a collisional equilibrium, derives the timescale that makes that explicit, and
shows that rewriting it changes no value.

The claim is an identity, not a model. Everything below is algebra applied to the two rates the
scheme already computes; no parameter is introduced, no closure is approximated, and no default
moves. What the rewriting buys is the *timescale*, which is the object the implicit solver's
Jacobian needs and which the sink form does not expose.

Throughout, ``L ≡ ρ q_{rai}`` [kg m``^{-3}``] is the rain mass density and ``N ≡ N_{rai}``
[m``^{-3}``] the number density. The mean drop mass is ``\overline{x}_r = L/N`` and the
mean-volume diameter is

```math
\overline{D}_r = \left( \frac{6 \overline{x}_r}{π ρ_w} \right)^{1/3},
```

a strictly decreasing function of ``N`` at frozen ``L``. That monotonicity is used repeatedly
below and is the reason a statement about drop size can be turned into a statement about number.

The derivation is stated for the windowed size distribution of
[the raindrop mean-mass window](@ref "The raindrop mean-mass window"). That is not a
presentational choice: it is the variant the result requires. Where the alternative clamp cascade
binds, ``\overline{D}_r`` is no longer the honest mean-volume diameter of the pair ``(L, N)``, the
sign argument in step 1 loses its hypothesis, and ``τ_{eff}`` loses its guarantee.

## 1. The pair is one rate with one zero

Self-collection ([SeifertBeheng2006](@cite) Eq. 11) is

```math
\left. \frac{∂N}{∂t} \right|_{sc}
  = -k_{rr} \, N \, L \left(1 + \frac{κ_{rr}}{B_r}\right)^{d} \left(\frac{ρ_0}{ρ}\right)^{1/2}
  \; ≡ \; -|sc| \; < \; 0 ,
```

strictly negative on any populated state. Breakup (Eq. 13) is coupled to it as
``\partial_t N|_{br} = -[Φ_{br} + 1] \, \partial_t N|_{sc}``. The pair's net number tendency is
therefore

```math
f(N; L) \; ≡ \; \left.\frac{∂N}{∂t}\right|_{sc} + \left.\frac{∂N}{∂t}\right|_{br}
  \; = \; -Φ_{br} \left.\frac{∂N}{∂t}\right|_{sc}
  \; = \; Φ_{br}(\overline{D}_r) \, |sc| .
\tag{1}
```

The pair is the self-collection *magnitude* modulated by ``Φ_{br}``. Every statement about the
pair's sign is therefore a statement about ``Φ_{br}`` alone, since ``|sc| > 0``.

**The zero.** ``f = 0`` exactly where ``Φ_{br} = 0``. On both branches that admit a zero,

```math
k_{br} Δ\overline{D}_r = 0
\quad\text{and}\quad
\exp(κ_{br} Δ\overline{D}_r) - 1 = 0
\qquad ⟺ \qquad \overline{D}_r = \overline{D}_{eq},
```

and the two branches meet continuously there. Writing ``x_{eq} = \tfrac{π}{6} ρ_w
\overline{D}_{eq}^3`` for the drop mass at the equilibrium diameter, the zero in ``N`` at frozen
``L`` is

```math
N_{eq}(L) \; = \; \frac{L}{x_{eq}} .
\tag{2}
```

``N_{eq}`` is linear in ``L``, defined for every ``L ≥ 0``, and exactly zero on an empty state, so
a relaxation toward it cannot manufacture number without mass and needs no special empty-state
arm.

**The sign, branch by branch.** Because ``\overline{D}_r`` decreases with ``N``,

```math
\overline{D}_r > \overline{D}_{eq}
\;⟺\; \overline{x}_r > x_{eq}
\;⟺\; N < N_{eq}
\;⟺\; N_{eq} - N > 0 ,
```

so it is enough to show ``\operatorname{sign} Φ_{br} = \operatorname{sign}(\overline{D}_r -
\overline{D}_{eq})`` on each of the three branches:

- ``\overline{D}_r > \overline{D}_{eq}`` (exponential): ``κ_{br} Δ\overline{D}_r > 0``, so
  ``\exp(κ_{br} Δ\overline{D}_r) > 1`` and ``Φ_{br} > 0``. Breakup dominates, number is *added*,
  and the mean size shrinks toward equilibrium.
- ``\overline{D}_{th} ≤ \overline{D}_r ≤ \overline{D}_{eq}`` (linear): ``k_{br} > 0`` and
  ``Δ\overline{D}_r ≤ 0``, so ``Φ_{br} ≤ 0``, vanishing only at ``\overline{D}_{eq}``. Collection
  dominates, number is removed, and the mean size grows toward equilibrium.
- ``\overline{D}_r < \overline{D}_{th}`` (pure collection): ``Φ_{br} ≡ -1``. This branch is worth
  stating explicitly, because breakup is switched **off** entirely and the sign nevertheless still
  agrees. Since ``\overline{D}_{th} < \overline{D}_{eq}`` (``0.35`` mm against ``0.9`` mm),
  ``\overline{D}_r < \overline{D}_{th}`` implies ``\overline{D}_r < \overline{D}_{eq}``, hence
  ``N > N_{eq}``, and ``f = -|sc| < 0`` removes number, moving ``N`` *down* toward ``N_{eq}``. The
  pure-collection branch is not an exception to the relaxation reading; it is the relaxation
  running with its restoring factor saturated at ``-1``.

Hence

```math
\operatorname{sign} f(N; L) \; = \; \operatorname{sign}\left(N_{eq}(L) - N\right)
\qquad\text{on every branch,}
\tag{3}
```

with both sides vanishing together at ``N = N_{eq}``. Equation (3) is the entire content of the
recast; everything after it is bookkeeping.

!!! note "``Φ_{br}`` jumps at ``\overline{D}_{th}``, and the sign argument survives it"
    ``Φ_{br}`` is *discontinuous* at the threshold diameter: approaching from above it tends to
    ``k_{br}(\overline{D}_{th} - \overline{D}_{eq}) = -0.55``, while below it is ``-1``. The jump
    is real and inherited from the fit. It does not disturb (3), because both one-sided values are
    negative and the branch lies entirely below ``\overline{D}_{eq}``. What jumps with it is the
    *magnitude* of ``τ_{eff}``, not its sign.

## 2. The timescale, and the removable singularity at equilibrium

Define

```math
τ_{eff}(N; L) \; ≡ \; \frac{N_{eq} - N}{f(N; L)} .
\tag{4}
```

By (3) the numerator and denominator share sign, so ``τ_{eff} > 0`` wherever the population is
present. It is nonzero and finite away from equilibrium: ``f`` vanishes only at ``N = N_{eq}``,
where the numerator vanishes too, and ``|Φ_{br}|`` is bounded because the mean-mass window bounds
``\overline{x}_r`` and hence ``\overline{D}_r``, so ``f`` cannot run away and ``τ_{eff}`` cannot
collapse.

At ``N = N_{eq}`` equation (4) is ``0/0``. The singularity is removable, and the limit is the
linearization of the pair about its own fixed point. By L'Hôpital in ``N`` at frozen ``L``, the
numerator's derivative being ``-1``,

```math
\frac{1}{τ_{relax}} \; = \; \lim_{N → N_{eq}} \frac{f}{N_{eq} - N}
  \; = \; -\left.\frac{∂f}{∂N}\right|_{N_{eq}} .
```

Differentiating (1) gives ``∂_N f = (∂_N Φ_{br})|sc| + Φ_{br} \, ∂_N |sc|``, and ``Φ_{br} = 0`` at
equilibrium, so **only ``Φ_{br}`` carries a derivative there** — the ``|sc|`` factor's own
``N`` dependence drops out. With ``\overline{D}_r ∝ N^{-1/3}`` at frozen ``L``, so that
``∂\overline{D}_r/∂N = -\overline{D}_r/(3N)``, and writing ``κ`` for the slope
``∂Φ_{br}/∂\overline{D}_r`` at ``Δ\overline{D}_r = 0``,

```math
\frac{1}{τ_{relax}(L)}
  \; = \; \frac{κ \overline{D}_{eq}}{3} \, \frac{|sc|}{N_{eq}}
  \; = \; \frac{κ \overline{D}_{eq}}{3} \, k_{rr} \, L
          \left(1 + \frac{κ_{rr}}{B_{r,eq}}\right)^{d}
          \left(\frac{ρ_0}{ρ}\right)^{1/2} ,
\qquad B_{r,eq} = \left(\frac{6}{x_{eq}}\right)^{1/3} .
\tag{5}
```

``N_{eq}`` cancels, because ``|sc|/N`` carries no ``N``: **the limit is a function of ``L``
alone.** That is what makes it usable as a relaxation rate — the damping at equilibrium depends
on how much rain there is, not on how the number happens to be distributed.

**The limit is two-sided.** [SeifertBeheng2006](@cite) fits different slopes on the two sides of
``\overline{D}_{eq}``, so ``κ = κ_{br}`` approaching from ``N < N_{eq}`` (that is, from
``\overline{D}_r > \overline{D}_{eq}``) and ``κ = k_{br}`` from ``N > N_{eq}``. ``f`` is
continuous at ``N_{eq}`` but not differentiable there, and ``1/τ`` has a finite jump of ratio
``κ_{br}/k_{br} = 2.3``. This is a property of the fit, not of the recast.

With the shipped constants, and ``L`` in kg m``^{-3}``:

```math
\frac{1}{τ_{relax}} \; = \; C \, L \left(\frac{ρ_0}{ρ}\right)^{1/2} \ \mathrm{s^{-1}},
\qquad
C = \frac{κ \overline{D}_{eq} k_{rr}}{3}\left(1 + \frac{κ_{rr}}{B_{r,eq}}\right)^{d} .
```

| quantity | value |
|---|---|
| ``x_{eq} = \tfrac{π}{6} ρ_w \overline{D}_{eq}^3`` | ``3.8170 \times 10^{-7}`` kg |
| ``B_{r,eq} = (6/x_{eq})^{1/3}`` | ``250.50`` kg``^{-1/3}`` |
| ``(1 + κ_{rr}/B_{r,eq})^{d}`` | ``0.33794`` |
| ``C`` approached from ``\overline{D}_r > \overline{D}_{eq}`` (``κ = κ_{br}``) | ``1.6602`` m``^3`` kg``^{-1}`` s``^{-1}`` |
| ``C`` approached from ``\overline{D}_r < \overline{D}_{eq}`` (``κ = k_{br}``) | ``0.72184`` m``^3`` kg``^{-1}`` s``^{-1}`` |

At ``L = 1`` g m``^{-3}`` and ``ρ = ρ_0`` these are ``τ_{relax} = 602`` s and ``1385`` s. The
units check as ``[k_{rr}][L] = \mathrm{m^3\,kg^{-1}\,s^{-1}} \cdot \mathrm{kg\,m^{-3}} =
\mathrm{s^{-1}}``, with ``κ \overline{D}_{eq}`` dimensionless.

## 3. The identity

```math
-\frac{N - N_{eq}}{τ_{eff}}
  \; = \; -(N - N_{eq}) \, \frac{f}{N_{eq} - N}
  \; = \; f .
\tag{6}
```

This is algebra. ``τ_{eff}`` is *defined* by (4) as the quotient that makes (6) hold, so whatever
``sc + br`` evaluates to, the relaxation form evaluates to the same thing — to the bit, not to a
tolerance. The recast introduces no approximation, adds no parameter, and changes no default. It
supplies one quantity the sink form does not: ``τ_{eff}`` itself.

!!! note "At the fixed point the identity is exact in ``ℝ`` and roundoff-level in floating point"
    At ``N = N_{eq}`` both sides of (6) are zero. In floating point they are zero to *roundoff*
    rather than exactly, because ``\overline{D}_r(N_{eq})`` round-trips through a cube root and
    lands ``Φ_{br}`` at ``\mathcal{O}(10^{-16})`` instead of ``0``. This is why the implementation
    selects the linearized branch on a *relative* deviation ``|N_{eq} - N| ≤ \sqrt{ε}\,N_{eq}``
    rather than testing ``N = N_{eq}``, and why tests at the fixed point must compare against a
    measured scale — the rate a factor two away — instead of a relative tolerance on a zero.

## 4. The Jacobian entry

The implicit substep linearizes at frozen distribution shape, holding ``τ_{eff}`` and ``N_{eq}``
fixed. Differentiating (6) under that convention,

```math
\left.\frac{∂f}{∂N}\right|_{\text{frozen shape}} \; = \; -\frac{1}{τ_{eff}} \; ≤ \; 0 ,
\tag{7}
```

non-positive wherever the population is present, by step 2. It is moreover **exact at the fixed
point**: the term dropped by freezing the shape is ``Φ_{br}\,∂_N|sc|``, which carries the factor
``Φ_{br} = 0`` there. So (7) is not merely sign-correct at equilibrium — it is the correct
linearization.

This is what the recast is for. Applying the scheme's donor-linearization recipe to the *sink*
form instead gives a diagonal entry ``(sc + br)/N``, which by (3) is **positive** wherever
breakup dominates, that is wherever ``N < N_{eq}``. That anti-damping entry reflects no physical
instability; it is the artifact of applying a donor recipe to a rate whose donor is not the
quantity it relaxes. Writing the rate in the variable it actually relaxes removes it.

!!! note "Where the mean-mass window binds, the entry is exactly zero"
    On states where the window clamp is active, the rate reaches ``N`` only through the clamp, so
    the chain rule gives exactly zero rather than something small. The row's damping there belongs
    to the number adjustment, which is the designated restorer for that regime.

## Provenance of ``\overline{D}_{eq}``

The equilibrium mean-volume diameter is not fitted here. It is the coalescence/breakup equilibrium
of the Low–List collision–breakup kernel ([LowList1982a](@cite), [LowList1982b](@cite)), whose
self-similar equilibrium distribution [SeifertBeheng2006](@cite) parameterizes through
``Φ_{br}``; ``\overline{D}_{eq} = 0.9`` mm is the value that equilibrium selects. The relaxation
target ``N_{eq} = L/x_{eq}`` therefore inherits its meaning from that lineage rather than from the
recast, which only rewrites where the scheme was already heading.

!!! warning "A documentation/implementation discrepancy in the exponential branch"
    The ``Φ_{br}`` definition rendered under "Raindrops breakup" in
    [the 2-moment page](@ref "Microphysics 2M") carries a factor two above
    ``\overline{D}_{eq}``, ``2[\exp(κ_{br} Δ\overline{D}_r) - 1]``, while the implementation uses
    ``\exp(κ_{br} Δ\overline{D}_r) - 1``. Both are continuous at ``Δ\overline{D}_r = 0`` and both
    are defensible readings of an equation the source paper prints ambiguously (see the note on
    that page). They are **not** interchangeable here: the factor propagates directly into (5), so
    the doubled form would give ``C = 3.3205`` rather than ``1.6602`` above ``\overline{D}_{eq}``
    and a two-sided ratio of ``4.6`` rather than ``2.3``. Everything on this page describes the
    implementation. Reconciling the two is a physics decision and is deliberately left open here.
