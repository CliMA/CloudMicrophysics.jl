# P3 Scheme

The `P3Scheme.jl` module implements the predicted particle properties (P3) scheme for ice-phase microphysics developed by 
  [MorrisonMilbrandt2015](@cite). The P3 scheme is a 2-moment, bulk scheme involving a single ice-phase category with 4 
  degrees of freedom: total ice content, rime content, rime volume, and number concentration. Traditionally, cloud ice 
  microphysics schemes use various predefined categories (such as ice, graupel, or hail) to represent ice modes, but the 
  P3 scheme sidesteps the problem of prescribing transitions between ice categories by adopting a single ice category and 
  evolving its properties. This simplification aids in attempts to constrain the scheme's free parameters.

The prognostic variables are:

| Variable            | Units         | Description                                        |
|:--------------------|:--------------|:---------------------------------------------------|
| $N_\mathrm{ice}$    | m$^{-3}$      | number concentration                               |
| $L_\mathrm{ice}$    | kg m$^{-3}$   | total ice particle mass content                    |
| $L_\mathrm{rim}$    | kg m$^{-3}$   | rime mass content                                  |
| $B_\mathrm{rim}$    | m$^3$ m$^{-3}$ | rime volume (volume of rime per total air volume) |

From these, we derive the following rime quantities:
 -  $F_{rim} = \frac{L_{rim}}{L_{ice}}$, the rime mass fraction,
 -  $ρ_{rim} = \frac{L_{rim}}{B_{rim}}$, the rime density.

!!! note "Change of symbol convention"
    The original paper [MorrisonMilbrandt2015](@cite) uses the symbol ``q`` to denote the mass of a tracer per volume of air 
    (named mass mixing ratio). In our documentation of the 1-moment and 2-moment schemes we used ``q`` to denote the mass 
    of a tracer per mass of air (specific content). To keep the notation consistent between the 1,2-moment schemes 
    and P3, and to highlight the difference between normalizing by air mass or volume, we denote the mass of a tracer per 
    volume of air as ``L`` (named content).

## [Assumed particle size relationships](@id P3-assumed-particle-size-relationships)

The mass ``m`` and projected area ``A`` of particles as a function of maximum particle dimension ``D`` are piecewise 
  functions with variable thresholds described by the following table,

| Particle properties                  | Size condition        | Rime condition | m(D) relation                         |    A(D) relation                                 |
|:-------------------------------------|:----------------------|:------------------------|:-----------------------------|:-------------------------------------------------|
|small, spherical ice                  | $D < D_{th}$          |                | $\frac{π}{6} ρ_i D^3$                 | $\frac{π}{4} D^2$                                |
|large, unrimed ice                    | $D_{th} < D$          | $L_{rim} = 0$  | $α_{va} D^{β_{va}}$                   | $γ D^{σ}$                                        |
|dense nonspherical ice                | $D_{th} < D < D_{gr}$ | $L_{rim} > 0$  | $α_{va} D^{β_{va}}$                   | $γ D^{σ}$                                        |
|graupel (completely rimed, spherical) | $D_{gr} < D < D_{cr}$ | $L_{rim} > 0$  | $\frac{π}{6} ρ_g D^3$                 | $\frac{π}{4} D^2$                                |
|partially rimed ice                   | $D_{cr} < D$          | $L_{rim} > 0$  | $\frac{α_{va}}{1-F_{rim}} D^{β_{va}}$ | $F_{rim} \frac{π}{4} D^2 + (1-F_{rim})γ \ D^{σ}$ |

where $D_{th}$, $D_{gr}$, $D_{cr}$, and $ρ_g$ are determined by relations below; $ρ_i$, $β_{va}$, $α_{va}$, $γ$, and $σ$ are fixed.

!!! details "Symbol definitions"
    | Symbol    | Value                                   | Units                    | Description                                                                                                                                                        |
    |:----------|:----------------------------------------|:-------------------------|:-------------------------------------------------------------------------------------------------------------------------------------------------------------------|
    | $D_{th}$  |                                         | $\text{m}$               | small, spherical ice threshold                                                                                                                                     |
    | $D_{gr}$  |                                         | $\text{m}$               | graupel threshold                                                                                                                                                  |
    | $D_{cr}$  |                                         | $\text{m}$               | partially rimed ice threshold                                                                                                                                      |
    | $ρ_i$     | 916.7                                   | $\text{kg m}^{-3}$       | cloud ice density                                                                                                                                                  |
    | $ρ_g$     |                                         | $\text{kg m}^{-3}$       | bulk density of graupel                                                                                                                                            |
    | $ρ_{rim}$ |                                         | $\text{kg m}^{-3}$       | rime density                                                                                                                                                       |
    | $ρ_d$     |                                         | $\text{kg m}^{-3}$       | bulk density of unrimed ice                                                                                                                                        |
    | $F_{rim}$ |                                         | --                       | rime mass fraction                                                                                                                                                 |
    | $β_{va}$  | 1.9                                     | --                       | dimensionless parameter from [BrownFrancis1995](@cite) (based on measurements of vapor diffusion and aggregation in midlatitude cirrus)                            |
    | $α_{va}$  | $7.38 \; 10^{-11} \; 10^{6 β_{va} - 3}$ | $\text{kg m}^{-β_{va}}$  | parameter modified for units from [BrownFrancis1995](@cite) in base SI units (also based on measurements of vapor diffusion and aggregation in midlatitude cirrus) |
    | $γ$       | 0.2285                                  | $\text{m}^{2-σ}$         | parameter fit from aggregates of side planes, columns, bullets, and planar polycrystals; see Table 1 in [Mitchell1996](@cite)                                      |
    | $σ$       | 1.88                                    | --                       | dimensionless parameter fit from aggregates of side planes, columns, bullets, and planar polycrystals; see Table 1 in [Mitchell1996](@cite)                        |


The first threshold, $D_{th}$, is solely determined by the free parameters:
```math
D_{th} = \left( \frac{π ρ_i}{6α_{va}} \right)^{\frac{1}{β_{va} - 3}}
```
The remaining thresholds: $D_{gr}$, $D_{cr}$, as well as the bulk density of graupel $ρ_g$, and the bulk density of the unrimed part $ρ_d$ form a system of equations:

```math
\begin{align*}
  D_{gr}  &= \left( \frac{6α_{va}}{π ρ_g} \right)^{\frac{1}{3 - β_{va}}} \\
  D_{cr}  &= \left[ \left(\frac{1}{1-F_{rim}} \right) \frac{6 α_{va}}{π ρ_g} \right]^{\frac{1}{3 - β_{va}}} \\
  ρ_g     &= ρ_{rim} F_{rim} + (1 - F_{rim}) ρ_d \\
  ρ_d     &= \frac{6α_{va}(D_{cr}^{β_{va} - 2} - D_{gr}^{β_{va} - 2})}{π \ (β_{va} - 2)(D_{cr}-D_{gr})}
\end{align*}
```

Unlike (REFS), which solve the system of equations with iterative methods, we here derive an analytical solution for $ρ_d$ in terms of $ρ_{rim}$ and $F_{rim}$.

!!! details "Click here to see the derivation"
    To derive the expression for $ρ_d$ we first observe that $D_{cr}$ is proportional to $D_{gr}$,

    ```math
    D_{cr} 
      = (1 - F_{rim})^{-1/(3 - β_{va})} \left( \frac{6 α_{va}}{π ρ_g} \right)^{1/(3 - β_{va})}
      = (1 - F_{rim})^{-1/(3 - β_{va})} D_{gr} 
      ≡ k D_{gr},
    ```
    where we define $k ≡ (1 - F_{rim})^{-1/(3 - β_{va})}$. This implies $k^{β_{va} - 2} = (1 - F_{rim})k$. 
    It can also be shown that
    ```math
    \frac{6α_{va}}{π}D_{gr}^{β_{va} - 2} = ρ_g D_{gr}.
    ```
    We now rewrite terms in the numerator and denominator of the expression for $ρ_d$,

    ```math
    \begin{align*}
    D_{cr} - D_{gr} &= (k - 1) D_{gr} \\
    % 
    \frac{6α_{va}}{π} \left[ D_{cr}^{β_{va} - 2} - D_{gr}^{β_{va} - 2} \right]
      &= \frac{6α_{va}}{π} \left[ (k^{β_{va} - 2} - 1) D_{gr}^{β_{va} - 2} \right]
       = ((1 - F_{rim})k - 1) ρ_g D_{gr}.
    \end{align*}
    ```
    Then, we substitute these expressions into the equation for $ρ_d$ (above), to obtain,
    ```math
    \begin{align*}
    ρ_d &= \frac{((1 - F_{rim}) k - 1) ρ_g D_{gr}}{(β_{va} - 2)(k - 1)D_{gr}} \\
      &= \frac{((1 - F_{rim}) k - 1)}{(β_{va} - 2)(k - 1)} ( ρ_{rim} F_{rim} + (1 - F_{rim}) ρ_d )
    \end{align*}
    ```
    and rearrange to obtain the expression for $ρ_d$ below.

We obtain the following expression for $ρ_d$
```math
ρ_d = \frac{
  ρ_{rim} F_{rim}
}{
  \frac{(β_{va} - 2)(k - 1)}{(1 - F_{rim}) k - 1} - (1 - F_{rim})
}
```

!!! details "Click here to see a numerically stable form"
    Evaluated directly, this expression loses accuracy as ``F_{rim} → 0``. In that limit ``k → 1`` and
    ``(1 - F_{rim}) k → 1``, so the ratio ``(k - 1) / ((1 - F_{rim}) k - 1)`` is ``0/0``. In `Float32` the
    resulting cancellation drives the computed ``ρ_d``, and therefore ``ρ_g``, negative. We rewrite the
    expression so that the cancellation is removed.

    Let ``F_u = 1 - F_{rim}`` be the unrimed fraction, let ``L = \log F_u``, and let ``p = 1 / (3 - β_{va})``,
    so that ``k = F_u^{-p}`` and ``(1 - F_{rim}) k = F_u^{1 - p}``. Every power of the form ``F_u^a - 1`` is
    equal to ``e^{aL} - 1``. We write these through the relative exponential functions

    ```math
    φ(a) = \frac{e^{aL} - 1}{aL} = \mathrm{exprel}_1(aL), \qquad
    ψ(a) = \frac{e^{aL} - 1 - aL}{(aL)^2} = \mathrm{exprel}_2(aL),
    ```

    both of which stay finite and accurate as ``L → 0``. Using ``β_{va} - 2 = -(1 - p) / p``, the denominator
    of ``ρ_d`` becomes ``φ(-p) / φ(1 - p) - F_u``. Both the numerator ``ρ_{rim} F_{rim}`` and the denominator
    are proportional to ``L`` as ``F_{rim} → 0``, so we factor ``L`` out of each. With ``F_{rim} = -L\,φ(1)``
    and ``φ(-p) - φ(1 - p) = L\,[-p\,ψ(-p) - (1 - p)\,ψ(1 - p)]``, the denominator divided by ``L`` is

    ```math
    G = \big[-p\,ψ(-p) - (1 - p)\,ψ(1 - p)\big] - φ(1 - p)\,φ(1),
    ```

    which tends to ``-3/2`` as ``F_{rim} → 0``. The factor ``L`` cancels, and we are left with

    ```math
    ρ_d = -\frac{ρ_{rim}\,φ(1)\,φ(1 - p)}{G}.
    ```

    This form has no subtraction of nearly equal numbers, so it stays accurate in `Float32`, and ``ρ_g``
    stays positive for every physical input without any clipping. The functions ``\mathrm{exprel}_1`` and
    ``\mathrm{exprel}_2`` (implemented by the internal `exprel` function) use a Taylor series at small
    arguments, where the closed forms would themselves cancel.

    The figure below shows the relative error of ``ρ_g`` in `Float32`, against a reference computed in
    `BigFloat`, for both the direct evaluation and this rewrite. The direct form is order-one wrong at small
    rime fractions and only reaches `Float32` precision near ``F_{rim} = 0.35``, while the rewrite stays at
    the precision floor across the whole range.

    ```@example
    include("plots/P3RhoDStability.jl")

    nothing # hide
    ```
    ![](P3RhoDStability.svg)

Given $ρ_d$, we can obtain $ρ_g$, $D_{gr}$, and $D_{cr}$ using the expressions above. Depending on the value of $ρ_{rim}$ and $F_{rim}$, 
  these thresholds and densities obtain a range of values, as shown in the plot below.

```@example
include("plots/P3Thresholds.jl")

nothing # hide
```
![](P3Thresholds.svg)

Below we show the ``m(D)`` and ``A(D)`` regimes replicating Figures 1 (a) and (b) from [MorrisonMilbrandt2015](@cite).
We also show the density as a function of ``D``.
Note that because graupel is completely filled with rime, the density (``ρ_g``) is independent of ``D`` between ``D_{gr}`` and ``D_{cr}``.
Following [MorrisonMilbrandt2015](@cite), for nonspherical particles ``ρ_{ice}`` is assumed to be equal to the mass of the particle
  divided by the volume of a sphere with the same particle size.

```@example
include("plots/P3SchemePlots.jl")

nothing # hide
```
![](P3Scheme_relations.svg)

## Assumed particle size distribution

Following [MorrisonMilbrandt2015](@cite), the scheme assumes a gamma distribution for the concentration of ice particles 
 per unit volume based on particle size measurements obtained by [Heymsfield2003](@cite) in tropical and midlatitude 
 ice clouds and implemented by [MorrisonGrabowski2008](@cite):

```math
N'(D) = N_{0} D^μ e^{-λ D}
```
where:
 -  $N'$   [m$^{-4}$]   $~~~~~$ is the number concentration,
 -  $D~~$  [m]      $~~~~~~~~~$ is the maximum particle dimension,
 -  $N_0$  [m$^{-4 - μ}$]   $~$ is the intercept parameter,
 -  $μ~~~$ [--]     $~~~~~~~~~$ is the shape parameter,
 -  $λ~~~$ [m$^{-1}$]   $~~~~~$ is the slope parameter.

The model predicted ice number concentration, $N_\mathrm{ice}$, and ice content, $L_\mathrm{ice}$, are defined as
```math
\begin{align*}
N_\mathrm{ice} &= \int_{0}^{∞}      N'(D) \mathrm{d}D \\
L_\mathrm{ice} &= \int_{0}^{∞} m(D) N'(D) \mathrm{d}D
\end{align*}
```

To close the system of equations, it has been customary (e.g. [MorrisonMilbrandt2015](@cite)) to assume $μ$ to be a function of $λ$, $μ = μ(λ)$.

## Calculating shape parameters

Given the system

```math
\begin{cases}
    N_\mathrm{ice} &= N_0 ∫_{0}^{∞}      D^μ e^{-λD}\ \mathrm{d}D \\
    L_\mathrm{ice} &= N_0 ∫_{0}^{∞} m(D) D^μ e^{-λD}\ \mathrm{d}D \\
    % μ &= 0.00191 λ^{0.8} - 2 \\    
\end{cases}
```

we seek to solve for $N_0$ and $λ$. This is most conveniently achieved by first solving for $λ$, then obtain $N_0$ by 
  rewriting the first equation. For numerical stability, we compute $\log \left( L_\mathrm{ice} \big/ N_\mathrm{ice} \right)$ 
  and solve for $λ$ using a root-finding algorithm for the expression:

```math
0
% = \log\left(\frac
%     {∫_{0}^{∞} m(D) D^μ e^{-λD}\ \mathrm{d}D}
%     {∫_{0}^{∞} D^μ e^{-λD}\ \mathrm{d}D}
% \right)
% - \log\left(\frac{L_\mathrm{ice}}{N_\mathrm{ice}}\right) 
= \log\left( ∫_{0}^{∞} m(D) D^μ e^{-λD}\ \mathrm{d}D \right)
- \log\left( ∫_{0}^{∞} D^μ e^{-λD}\ \mathrm{d}D \right)
- \log\left(\frac{L_\mathrm{ice}}{N_\mathrm{ice}}\right) 
```

The slope is a derived quantity of the ice state.
Before the solve, the distribution construction bounds the ice number so that the mean particle mass $L_\mathrm{ice} / N_\mathrm{ice}$ lies between the mass of a newly nucleated crystal (a solid-ice sphere of diameter $D_\mathrm{nuc}$) and a largest mean mass of $10^{-4}$ kg.
The root then lies in the search interval of the solve, from $\log λ = 2$ to the slope at which the mean particle mass equals the nucleation mass.
Some of the gamma functions used to evaluate $\log λ$ do not have analytical derivatives, so $\log λ$ is treated as a constant with respect to the state when using automatic differentiation.
$\log λ$ and the other distribution parameters are stored in the [`CloudMicrophysics.P3Scheme.P3Distribution`](@ref) type.
The [`CloudMicrophysics.P3Scheme.P3State`](@ref) holds the prognostic moments, and constructs the distribution from them and from the [`CloudMicrophysics.P3Scheme.P3Particles`](@ref), the mass, area, and size thresholds of a particle.
The process rates read the distribution from the state to compute the integrals over the distribution.

!!! details "Computing Distribution Integrals Using Gamma Functions"
    The integrals presented above can be evaluated analytically using the Gamma function family. Here we present the key 
      formulas and their numerical implementation.

    ##### Core Definitions
    The Gamma function is defined as
    ```math
    Γ(a) = ∫_0^∞ t^{a-1} e^{-t}\ \mathrm{d}t.
    ```
    For partial intervals, we use the upper incomplete Gamma function,
    ```math
    Γ(a,x) = ∫_x^∞ t^{a-1} e^{-t}\ \mathrm{d}t
    ```

    ##### Analytical Solutions
    Given a power law mass-diameter relationship $m(D)$, the key integrals resolve to
    ```math
    \begin{align*}
    ∫_0^∞ D^μ e^{-λD}\ \mathrm{d}D &= \frac{Γ(μ+1)}{λ^{μ+1}} \\
    ∫_0^∞ (aD^b) D^μ e^{-λD}\ \mathrm{d}D &= a \frac{Γ(b+μ+1)}{λ^{b+μ+1}} \\
    ∫_{D_1}^{D_2} (aD^b) D^μ e^{-λD}\ \mathrm{d}D &= \frac{a}{λ^{b+μ+1}} \Big( Γ(b+μ+1,λD_1) - Γ(b+μ+1,λD_2) \Big)
    \end{align*}
    ```

    ##### Numerical Implementation
    For numerical stability, we compute these integrals in log space. The regularized Gamma functions 
      $p(a,x)=\frac{1}{Γ(a)}\int_0^x t^{a-1} e^{-t}\ \mathrm{d}t$ and $q(a,x)=1-p(a,x)$ help us evaluate the incomplete Gamma functions:
    ```math
    \begin{align*}
    \log\left( ∫_{0}^{∞} D^μ e^{-λD}\ \mathrm{d}D \right) &= \log\Big( Γ(μ+1,λx) \Big) - (μ+1)\log(λ) \\
    \log\left( ∫_{D_1}^{D_2} (aD^b) D^μ e^{-λD}\ \mathrm{d}D \right) &= \log(a) - (b+μ+1)\log(λ) + \log\Big(Γ(b+μ+1)\Big) \\
    &\quad + \log\Big( q(b+μ+1,λD_1) - q(b+μ+1,λD_2) \Big)
    \end{align*}
    ```
    We define this latter integral as
    ```math
    G(D_1,D_2,a,b,μ,λ) ≡ \log\left( ∫_{D_1}^{D_2} (aD^b) D^μ e^{-λD}\ \mathrm{d}D \right),
    ```
    
    The integral of $m(D) N'(D)$ can be computed using the following function definitions, $G_i$:

    | Particle properties                  | Size condition        | Rime condition | m(D) relation                         | $G_i$ definition                                                    |
    |:-------------------------------------|:----------------------|:---------------|:--------------------------------------|:--------------------------------------------------------------------|
    |small, spherical ice                  | $D < D_{th}$          |                | $\frac{π}{6} ρ_i D^3$                 | $G_1 = G\left(0,D_{th},\frac{π}{6} ρ_i,3,μ,λ\right)$                |
    |large, unrimed ice                    | $D_{th} < D$          | $L_{rim} = 0$  | $α_{va} D^{β_{va}}$                   | $G_2 = G\left(D_{th},∞,α_{va},β_{va},μ,λ\right)$                    |
    |dense nonspherical ice                | $D_{th} < D < D_{gr}$ | $L_{rim} > 0$  | $α_{va} D^{β_{va}}$                   | $G_3 = G\left(D_{th},D_{gr},α_{va},β_{va},μ,λ\right)$               |
    |graupel (completely rimed, spherical) | $D_{gr} < D < D_{cr}$ | $L_{rim} > 0$  | $\frac{π}{6} ρ_g D^3$                 | $G_4 = G\left(D_{gr},D_{cr},\frac{π}{6} ρ_g,3,μ,λ\right)$           |
    |partially rimed ice                   | $D_{cr} < D$          | $L_{rim} > 0$  | $\frac{α_{va}}{1-F_{rim}} D^{β_{va}}$ | $G_5 = G\left(D_{cr},∞,\frac{α_{va}}{1-F_{rim}},β_{va},μ,λ\right)$  |

    Let LSE denote the log-sum-exp operation: $\mathrm{LSE}(x_1,\ldots,x_n) = \log\sum_{i=1}^n \exp(x_i)$. Then,

    - If $L_\mathrm{rim} = 0$:

    ```math
    \log ∫_{0}^{∞} m(D) D^μ e^{-λD}\ \mathrm{d}D = \mathrm{LSE}(G_1, G_2)
    ```

    -  If $L_\mathrm{rim} > 0$:

    ```math
    \log ∫_{0}^{∞} m(D) D^μ e^{-λD}\ \mathrm{d}D = \mathrm{LSE}(G_1, G_3, G_4, G_5)
    ```

    The LSE operation is computed in a numerically stable manner using the `logsumexp` function from the `LogExpFunctions.jl` package.


### Parameterizations for the slope parameter $μ$

#### $μ$ as a power law in $λ$

The de-facto standard parameterization for $μ$ as a function of $λ$ is a power law on the form

```math
μ = a λ^b - c
```

typically with lower and upper limiters. In [MorrisonMilbrandt2015](@cite), the following coefficients and limiting behavior are used:

```math
μ = \begin{cases}
    0 & \text{if } μ < 0 \\
    0.00191 λ^{0.8} - 2 & \text{if } 0 ≤ μ ≤ 6 \\
    6 & \text{if } μ > 6
\end{cases}
```
This hard-limited form is compared here with the smooth form, `SmoothSlopePowerLaw`, which `ParametersP3` uses by default:

```@example
include("plots/P3SlopeParameterizations.jl")

nothing # hide
```
![](P3SlopeParameterizations_slope_laws.svg)

The vertical lines mark the transitions in the hard-limited relation.

For the ice state shown here, the hard limits give three solutions for $\lambda$ at the same $\log(L/N)$.
The smooth form gives one solution.
The left panel marks the roots of each shape residual, and the right panel shows the corresponding number distributions $N'(D)$.

![](P3SlopeParameterizations_multiple_solutions.svg)

The smooth form replaces the hard limits with smooth transitions, because the lower limit can make $\log(L/N)$ nonmonotonic in $\lambda$.

#### $μ$ as a constant

An alternative parameterization for $μ$ is a constant value:

```math
μ = μ_{const}
```

## Terminal Velocity

We use the Chen et al. [Chen2022](@cite) velocity parameterization,
  see [here](https://clima.github.io/CloudMicrophysics.jl/dev/TerminalVelocity.html#Chen-et.-al.-2022) for details.

### Aspect ratio

The ice particle is modeled as an oblate spheroid of solid material with maximum
dimension $D$, projected area $a_i = $ `ice_area(state, D)`, and mass
$m_i = $ `ice_mass(state, D)`. With equatorial diameter $a = 2\sqrt{a_i/π}$ and polar
diameter $c = ϕ\,a$, the spheroid encloses the volume $(π/6)\,a^2 c$ and holds the mass
$ρ\,(π/6)\,a^2 c$, where $ρ$ is the material density `ϕ_material_density`: the density
of the solid the particle consists of in its mass regime ($ρ_i$ for ice, $ρ_g$ for
graupel). Equating this mass to $m_i$ gives the oblate aspect ratio

```math
ϕ = \frac{3 \sqrt{π}\, m_i}{4\, ρ\, a_i^{3/2}}.
```

The aspect ratio is a statement about particle shape: it measures how far the solid
material filling the spheroid falls short of a sphere of the same projected area. The
density in this relation is therefore the density of the solid, which is constant
within each mass regime, and $ϕ$ tracks $m_i / a_i^{3/2}$. In the spherical regimes
($D < D_{th}$, and graupel $D_{gr} ≤ D < D_{cr}$) the mass and area laws are exactly
spherical and $ϕ = 1$. In the dense-nonspherical regime $ϕ < 1$ (oblate) and decreases
with $D$, slowing large unrimed ice through the $\sqrt[3]{ϕ}$ factor in
`ice_particle_terminal_velocity`.

!!! note "Residual $ϕ > 1$ band just above $D_{th}$"
    $ϕ$ can slightly exceed 1 (peak $≈ 1.2$) in a narrow size band immediately above
    $D_{th}$. This is not a sign error: it is caused by the area discontinuity at
    $D_{th}$, where the projected area drops from the spherical law $(π/4) D^2$ to the
    nonspherical law $γ D^σ$ while the mass stays continuous. The smaller area inflates
    $m_i / a_i^{3/2}$ just past the threshold. $\sqrt[3]{ϕ}$ remains finite and bounded
    there ($\sqrt[3]{1.2} ≈ 1.06$), so the value is used as-is. Making the area law
    continuous at $D_{th}$, which would remove the band, is a property of the area
    power law itself.

The figure below shows the implied aspect ratio for different particle size regimes.
Note that $ϕ = 1$ corresponds to spherical particles (small spherical ice ($D < D_{th}$) and graupel ($D_{gr} < D < D_{cr}$)),
and that $ϕ$ exceeds 1 in a narrow band just above $D_{th}$ (the area discontinuity noted above; a follow-on item).
```@example
include("plots/P3AspectRatioPlot.jl")

nothing # hide
```
![](P3Scheme_aspect_ratio.svg)

The mass-weighted fall speed (``V_m``) and the number-weighted fall speed (``V_n``) are calculated as
```math
\begin{align*}
V_m &= \frac{\int_0^∞ V(D) m(D) N'(D) \mathrm{d}D}{\int_0^∞ m(D) N'(D) \mathrm{d}D} \\
V_n &= \frac{\int_0^∞ V(D)      N'(D) \mathrm{d}D}{\int_0^∞      N'(D) \mathrm{d}D}
\end{align*}
```

We also plot the mass-weighted mean particle size ``D_m`` which is given by:
```math
D_m = \frac{\int_0^∞ D m(D) N'(D) \mathrm{d}D}{\int_0^∞ m(D) N'(D) \mathrm{d}D}
```

Below we provide plots of these relationships for small, medium, and large ``D_m``:
  the first row highlights the particle size regime,
  the second displays ``D_m`` of the particles,
  the third shows the aspect ratio ``ϕ (D_m)``,
  the fourth shows ``V_m`` without using aspect ratio in the computation (i.e. ``ϕ = 1``),
  and the final row exhibits ``V_m``, computed using ``ϕ(D)``.
They can be compared with Figure 2 from [MorrisonMilbrandt2015](@cite).
```@example
include("plots/P3TerminalVelocityPlots.jl")

nothing # hide
```
![](MorrisonandMilbrandtFig2.svg)

## Liquid Fraction

!!! todo "TODO: Update code to include liquid fraction"
  Currently, the code doesn't work with liquid fraction. So ignore the following section.

To allow for the modeling of mixed-phase particles with P3, a new prognostic variable can be introduced:
  ``L_{liq}``, the content of liquid on mixed-phase particles. As described in [Choletteetal2019](@cite),
  this addition to the framework of P3 opens the door to tracking the gradual melting (and refreezing) of a
  particle population. Here, we describe the characteristics of the scheme with the addition of ``L_{liq}``.

Liquid fraction, analogous to ``F_{rim}`` for rime, is defined ``F_{liq} = \frac{L_{liq}}{L_{p3, tot}}``
  where ``L_{p3, tot} = L_{ice} + L_{liq}`` and where ``L_{ice}`` is the solid ice mass content which includes
  rime mass content and mass grown by vapor deposition. This is important notably for
  ``F_{rim} = \frac{L_{rim}}{L_{ice}}``, which differs from ``F_{rim} = \frac{L_{rim}}{L_{p3, tot}}``
  because we want to normalize only by solid ice mass content for ``F_{rim}``.

Based on Fig. 1 from [Choletteetal2019](@cite), we can expect the accumulation of liquid on an ice core to increase
  velocities of small particles for all ``F_{rim}`` values. Below, we reproduce this figure with ``\rho_{r} = 900 kg m^{-3}``,
  and, notably, we continue to use the terminal velocity parameterizations from [Chen2022](@cite)
  (described [here](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics2M/#Terminal-Velocity) for rain), 
  whereas [Choletteetal2019](@cite) uses [MitchellHeymsfield2005](@cite) for snow and ice and [Simmeletal2002](@cite) 
  for rain. Despite these different choices, we reproduce similar behavior with, and we include a dashed line for the 
  velocity computed without aspect ratio. This dashed line mimics velocity with ``F_{rim} = 1``, since both ``\phi = 1`` 
  and ``F_{rim} = 1`` shift us into a spherical particle regime.

```@example
nothing
# include("plots/Cholette2019_fig1.jl")
# ![](Choletteetal2019_fig1.svg)
```

The addition of the liquid fraction does not change the thresholds ``D_{th}``, ``D_{gr}``, ``D_{cr}``,
  since the threshold regime depends only on ice core properties.

However, the assumed particle properties become ``F_{liq}``-weighted averages of particles' solid and liquid components:

```math
\begin{align*}
m(D, F_{liq}) &= (1 - F_{liq}) m(D, F_{liq} = 0) + F_{liq} m_{liq}(D) \\
A(D, F_{liq}) &= (1 - F_liq) A(D, F_{liq} = 0) + F_{liq} A_{liq}(D)
\end{align*}
```

where ``m_{liq}(D) = \frac{π}{6} ρ_{liq} D^3`` and ``A_{liq}(D) = \frac{π}{4} D^2``.

When calculating shape parameters and integrating over the particle size distribution (PSD), it is important to
  keep in mind whether the desired moment of the PSD is tied only to ice (in which case we concern ourselves with the
  ice core diameter ``D_{core}``) or to the whole mixed-phase particle (in which case we need ``D_p``). If calculating
  the ice core parameters, ``F_{liq} = 0`` is passed into the solver framework as indicated in the code.

For the above particle properties and for terminal velocity, we use the PSD corresponding to the whole mixed-phase particle,
  so our terminal velocity of a mixed-phase particle is:

```math
V(D_{p}, F_{liq}) = (1 - F_{liq}) V_{ice}(D_{p}) + F_{liq} V_{rain}(D_{p})
```

To continue with the same plotting format as we see above for terminal velocity, below, we show terminal velocity with
  ``F_{liq} = 0.0, 0.5, 0.9`` using the mass-weighted terminal velocity with aspect ratio. Clearly, the velocity increases
  with ``F_{liq}`` and becomes decreasingly dependent on ``F_{rim}`` and ``\rho_{r}`` as we shrink the size of the ice core.
  Of note relating to our calculation of terminal velocity with liquid fraction, we let ``F_{liq} = 0`` in our calculation of
  ``\phi``. This is because we want ``V_{ice}(D_{p})`` (velocity of a mixed-phase particle treating it as an ice particle
  with the same D).

Visualizing mass-weighted terminal velocity as a function of ``F_{liq}``, ``F_{rim}`` with ``ρ_{rim} = 900 kg m^{-3}`` for small,
  medium, and large particles, we have mostly graupel (``D_{gr} < D < D_{cr}``) for small ``D_m`` and mostly partially rimed ice
  (``D > D_{cr}``) for medium and large ``D_m``. Thus, we can attribute the non-monotonic behavior of velocity with ``F_liq`` in the
  medium and large ``D_m`` plots to the variations in ``ϕ`` caused by nonspherical particle shape, whereas the small ``D_m`` plot
  confirms a more monotonic change in ``V_m`` for spherical ice. The ``L`` and ``N`` values used to generate small, medium, and large ``D_{m}``
  are the same as in the plot above.

When modifying process rates, we now need to consider whether they are concerned with the ice core or the whole particle, 
  in addition to whether they become sources and sinks of different prognostic variables with the inclusion of ``F_{liq}``.
  With the addition of liquid fraction, too, come new process rates.

## Microphysical Process Rates

At a high level, we can categorize the microphysical process rates that affect the P3 prognostic variables into the following categories:

| Symbol | Description |
|:-------|:------------|
| $\textcolor{magenta}{\textsf{NUC}}$    | Ice nucleation, including homogeneous and heterogeneous freezing |
| $\textcolor{brown}{\textsf{COL}}$    | Collision/collection between liquid and ice, including wet growth and shedding |
| $\textcolor{lime}{\textsf{DEP/SUB}}$    | Ice deposition/sublimation |
| $\textcolor{pink}{\textsf{MLT}}$    | Melting of ice |
| $\textcolor{orange}{\textsf{SLF}}$    | Self-collection of ice |

The scheme is coupled to the 2-moment Seifert & Beheng (2006) warm rain microphysics scheme with non-equilibrium moisture.
Some of the "P3 processes" naturally affect the 2-moment prognostic variables, and will be described in the following sections.
At a high level, the processes enter as sources and sinks of the prognostic variables, as follows:

Liquid phase:
```math
S_x = \textcolor{pink}{    \text{MLT}} 
      \textcolor{brown}{ - \text{COL}},
      \qquad x ∈ \{q_c, q_r, N_c, N_r\}
```

Ice phase:
```math
\begin{alignat*}{6}
  S_{L_{rim}} &= % NUC + COL + F_rim (- SUB - MLT)
                              \textcolor{magenta}{\text{NUC}}
    &&                        \textcolor{brown}{ + \text{COL}}
    &&         + F_{rim} (    
                              \textcolor{lime}{  - \text{SUB}}
    &&                        \textcolor{pink}{  - \text{MLT}}
                         )
    && 
    \\
  S_{L_{ice}} &= % NUC + COL - SUB - MLT + DEP
                              \textcolor{magenta}{\text{NUC}}
    &&                        \textcolor{brown}{+ \text{COL}}
    &&\phantom{ + F_{rim} ( }
                              \textcolor{lime}{  - \text{SUB}}
    &&                        \textcolor{pink}{ - \text{MLT}}
                          %)
    &&                        \textcolor{lime}{+ \text{DEP}}
    \\ 
  S_{N_{ice}} &= % NUC - SUB - MLT - SLF
                              \textcolor{magenta}{\text{NUC}} 
    &&\phantom{               \textcolor{brown}{+ \text{COL}} }
    &&\phantom{ + F_{rim} ( } 
                              \textcolor{lime}{ - \text{SUB}}
    &&                        \textcolor{pink}{ - \text{MLT}}
                              \textcolor{orange}{ - \text{SLF}}
    \\ 
  S_{B_{rim}} &=
                              \textcolor{magenta}{\text{NUC}}
    &&                        \textcolor{brown}{ + \text{COL}}
    &&         + F_{rim} (    
                              \textcolor{lime}{  - \text{SUB}}
    &&                        \textcolor{pink}{  - \text{MLT}}
                         )
    && 
    \\
              % NUC + COL + F_rim (- SUB - MLT)
\end{alignat*}
```

!!! note "Differences with Morrison & Milbrandt (2015)"
    There are numerous differences between the P3 scheme in Morrison & Milbrandt (2015) and the P3 scheme in this package:
    - the maximum freezing rate (wet growth limit) is computed at the particle level, instead of as a bulk process
    - changes in rime volume due to (dry) collisions with cloud and rain are computed at the particle level, with a local rime density evaluated as a function of both the liquid and ice particle diameters
    - in the wet growth regime, the densification process is computed as a rapidly adjusting bulk process, instead of instantenous particle-level densification. We also scale by the fraction of the collected liquid mass that is shed.


Below, we describe the different processes in more detail.

### Collisions with liquid droplets (``\textcolor{brown}{\textsf{COL}}``)

At a high level, collisions between liquid droplets (cloud or rain) and ice particles lead to sources and sinks of the form:

Liquid phase:
```math
\begin{align*}
  S_{q_c} &= \frac{\textcolor{brown}{ - \text{QCFRZ} - \text{QCSHD}}}{ρ_a}
          \\
  S_{q_r} &= \frac{\textcolor{brown}{ - \text{QRFRZ} + \text{QCSHD}}}{ρ_a}
          \\
  S_{N_c} &= \textcolor{brown}{ - \text{NCCOL}} 
          \\
  S_{N_r} &= \textcolor{brown}{ - \text{NRCOL} + \text{NRSHD}} 
          \\
\end{align*}
```
Ice phase:
```math
\begin{align*}
  S_{L_{rim}} &= \textcolor{brown}{\text{QCFRZ} + \text{QRFRZ} + \text{QIWET}}
              \\
  S_{L_{ice}} &= \textcolor{brown}{\text{QCFRZ} + \text{QRFRZ}}
              \\
  S_{N_{ice}} &= 0
              \\
  S_{B_{rim}} &= \textcolor{brown}{\text{BCCOL} + \text{BRCOL} + \text{BIWET}} 
              \\
\end{align*}
```

where $\textcolor{brown}{\text{COL}}$ is the collision rate between liquid droplets and ice particles, $\textcolor{brown}{\text{FRZ}}$ is the component of the collisions that freezes, and $\textcolor{brown}{\text{SHD}}$ is the component of the collisions that lead to droplet shedding (shedding + freezing = collisions). The first letter of the prefix ``Q``, ``N``, and ``B`` specifies that the process rates affect the mass, number, and bulk rime volume of the ice particle, respectively. The second letter, ``C``, ``R``, and ``I``, specify that the process rates affect the cloud, rain, and ice particles, respectively. For multispecies interactions, the reduced species determine the second letter.

Collisions between liquid droplets (cloud or rain) and ice particles are parameterized 
 by integrating the volumetric collision rate $∂_t\mathcal{V}$ [$\text{m}^3/\text{s}$] 
 over the particle size distributions. The volumetric collision rate is on the form
```math
∂_t\mathcal{V}_l(D_i, D_l) = E(D_i, D_l) K(D_i, D_l) |v_i(D_i) - v_l(D_l)|
```
where $l ∈ \{c, r\}$ denotes cloud or rain, 
- ``E(D_i, D_l)`` is the collision efficiency (which we assume to be 1), 
- ``K(D_i, D_l) = π (r_i + D_l/2)^2`` is the collision cross section, 
  - ``r_i`` is the effective radius of the ice particle, which is the radius of a circle with the same cross-sectional area as the (in general, non-spherical) ice particle,
- ``v_i(D_i)`` and ``v_l(D_l)`` are the terminal velocities of the ice particle and the liquid droplet, respectively.

!!! note "Number of collisions per unit time"
    If we multiply $∂_t\mathcal{V}_l(D_i, D_l)$ by the liquid droplet size distribution $N_l(D_l)$, we obtain the kernel
    ```math
    ∂_t\mathcal{V}_l(D_i, D_l) N_l(D_l),
    ```
    which represents the number of liquid droplets of size ``D_l`` colliding with an ice particle of size ``D_i`` per unit time. If we know, say, the mass of each liquid droplet, and integrate the mass-kernel over all liquid drop sizes, we get the total mass of the collected liquid droplets by an ice particle of size ``D_i``.

The quantity $∂_t\mathcal{N}_{\text{col},l}(D_i)$ [$\text{s}^{-1}$],
```math
∂_t\mathcal{N}_{\text{col},l}(D_i) = ∫_0^∞ ∂_t\mathcal{V}_l(D_i, D_l) N_l(D_l) \mathrm{d}D_l
```
quantifies how many liquid droplets collide with an ice particle of size $D_i$ per unit time.

!!! note "Bulk loss of liquid droplets"
    The bulk loss of liquid droplets is then given by
    ```math
    \textcolor{brown}{\text{NCCOL}} = ∫_0^∞ ∂_t\mathcal{N}_{\text{col},l}(D_i) N'_i(D_i) \mathrm{d}D_i,
    \quad
    \textcolor{brown}{\text{NRCOL}} = ∫_0^∞ ∂_t\mathcal{N}_{\text{col},r}(D_i) N'_i(D_i) \mathrm{d}D_i,
    ```

The total mass of the collected liquid droplets $∂_t\mathcal{M}_{col,l}(D_i)$ [$\text{kg/s}$] can be expressed as
```math
∂_t\mathcal{M}_{\text{col},l}(D_i) 
  = ∫_0^∞ ∂_t\mathcal{V}_l(D_i, D_l) N_l(D_l) m_l(D_l) \mathrm{d}D_l
```

When liquid droplets collide with ice particles, we assume that they can either freeze or be shed as
 rain particles. Above freezing, all particles are shed. In subfreezing temperatures, normally all particles are frozen. However,
 close to the freezing temperature, there may not be enough time for all particles to freeze.
 We refer to this as the wet growth regime. In this case, the particles that are not frozen are shed as rain.
 The wet growth regime is determined by comparing the freezing rate to the collection rate.

The maximum freezing rate is based on Musil (1970) [Musil1970](@cite),
 which considers the heat transfer rate for a spherical wet hailstone.
 The maximum freezing rate [kg/s] is computed as
```math
∂_t\mathcal{M}_\text{max}(D_i) 
  = 2π D_i F_v(D_i) \frac{K_t ΔT + L_v D_v Δρ_{v,\text{sat}}}{L_f - C_p ΔT}
```
where the first term in the numerator is the heat transfer rate to the air surrounding the hailstone,
 and the second term is evaporative cooling. The denominator is the heat that must be dissipated by the hailstone.
The latent heats ``L_v`` and ``L_f`` are evaluated at the freezing temperature, the temperature of the particle surface.

!!! details "Symbol definitions"
    | Symbol     | Units | Description |
    |:-----------|:------|:------------|
    | $D_i$      | $\text{m}$       | diameter of the ice particle |
    | $F_v(D_i)$ | $-$              | ventilation factor |
    | $K_t$      | $\text{W/(m K)}$ | thermal conductivity |
    | $ΔT$       | $\text{K}$       | temperature difference between the surface of the ice particle and the surrounding air |
    | $L_v$      | $\text{J/kg}$    | latent heat of vaporization |
    | $D_v$      | $\text{m}^2/\text{s}$    | vapor diffusivity |
    | $Δρ_{v,\text{sat}}$  | $\text{kg}/\text{m}^3$ | difference in saturation vapor density over ice, between the surface of the ice particle and the surrounding air |
    | $L_f$      | $\text{J/kg}$     | latent heat of fusion |
    | $C_p$      | $\text{J/(kg K)}$ | specific heat capacity of water |

For collisions with ice particles of size $D_i$, the freezing (riming) rate $∂_t\mathcal{M}_\text{frz}$ [$\text{kg/s}$] is thus
```math
∂_t\mathcal{M}_\text{frz} = \min \left( ∂_t\mathcal{M}_\text{col}, ∂_t\mathcal{M}_\text{max} \right), 
\quad \textsf{where} \quad
∂_t\mathcal{M}_\text{col} = ∂_t\mathcal{M}_\text{col,c} + ∂_t\mathcal{M}_\text{col,r},
```
where $∂_t\mathcal{M}_{col}$ is the total mass of the collected liquid droplets. The fraction of the collected liquid that freezes is given by
```math
f_\text{frz}(D_i) 
  = \frac{∂_t \mathcal{M}_\text{frz}(D_i)}{∂_t \mathcal{M}_\text{col}(D_i)}
```

!!! note "Bulk freezing rates"
    The bulk freezing rates are obtained by integrating the per-particle freezing rates over the ice particle size distribution,
    ```math
    \textcolor{brown}{\text{QCFRZ}} 
      = ∫_0^∞ ∂_t\mathcal{M}_{\text{col},c}(D_i) f_\text{frz}(D_i) N'_i(D_i) \mathrm{d}D_i,
    \quad
    \textcolor{brown}{\text{QRFRZ}} 
      = ∫_0^∞ ∂_t\mathcal{M}_{\text{col},r}(D_i) f_\text{frz}(D_i) N'_i(D_i) \mathrm{d}D_i, 
    ```
    which can be combined to give the bulk liquid freezing rate
    ```math
    \textcolor{brown}{\text{QCFRZ} + \text{QRFRZ}} 
      = ∫_0^∞ ∂_t \mathcal{M}_\text{frz}(D_i)                   N'_i(D_i) \mathrm{d}D_i 
      = ∫_0^∞ ∂_t \mathcal{M}_\text{col}(D_i) f_\text{frz}(D_i) N'_i(D_i) \mathrm{d}D_i
    ```

Any excess mass $∂_t\mathcal{M}_{shd}$ [$\text{kg/s}$] is shed as rain, where
```math
∂_t\mathcal{M}_\text{shd} = ∂_t\mathcal{M}_{col} - ∂_t\mathcal{M}_\text{frz}
```
We assume that all shed droplets are shed at some fixed diameter $D_\text{shd}$. 
  Following [MorrisonMilbrandt2015](@cite), we set $D_\text{shd} = 1 \text{ mm}$. 
  This implies that the number of shed droplets is
```math
∂_t\mathcal{N}_\text{shd} = \frac{∂_t\mathcal{M}_\text{shd}}{m(D_\text{shd})}
```

!!! note "Bulk source to rain mass and number"
    The corresponding bulk source to rain mass and number is
    ```math
    \begin{align*}
    \textcolor{brown}{\text{QCSHD}}
    &= ∫_0^∞ (1 - f_{\text{frz}}(D_i)) ∂_t\mathcal{M}_\text{c,col}(D_i) N'_i(D_i) \mathrm{d}D_i
    \\
    \textcolor{brown}{\text{QRSHD}}
    &= ∫_0^∞ (1 - f_{\text{frz}}(D_i)) ∂_t\mathcal{M}_\text{r,col}(D_i) N'_i(D_i) \mathrm{d}D_i
    \\
    \textcolor{brown}{\text{NRSHD}}
    &= \frac{\textcolor{brown}{\text{QRSHD}}}{m(D_\text{shd})},
    \end{align*}
    ```
    Because we assume that any shed droplets are rain, the shed cloud water
    ``\text{QCSHD}`` is a source of rain mass, while the shed rain water
    ``\text{QRSHD}`` returns to the rain category and does not change the
    rain mass.
    The shed drops re-enter rain at the shedding diameter ``D_\text{shd}``,
    which gives the rain number source ``\text{NRSHD}``
    (the corresponding number source from shed cloud water is currently
    neglected).

Finally, how does this affect the rime volume? 

Collisions with rain, cloud, and shedded wet growth all contribute to the rime volume.

The change in rime volume is given by the change in mass divided by some rime density.

!!! note "Local rime density parameterization"
    Following [MorrisonMilbrandt2015](@cite), we apply the Cober & List (1993) [CoberList1993](@cite) parameterization to calculate the local rime density for each collision pair. 
    For an ice particle of size $D_i$ colliding with a liquid droplet of size $D_l$ (where $l ∈ \{c, r\}$ for cloud or rain), the rime density is given by
    ```math
    \begin{align}
    V_\text{term}(D_i, D_l) &= |v_i(D_i) - v_l(D_l)|, 
    \\
    R_i(D_i, D_l) &= \frac{D_l \, V_\text{term}(D_i, D_l)}{2 |T - T_\text{freeze}|},
    \\
    ρ'_{rim}(R_i) &= 
    \begin{cases}
      a_\text{CL} + b_\text{CL} R_i + c_\text{CL} R_i^2         & \text{if }\, 1 ≤ R_i ≤ 8, 
      \\
      (1 - \frac{R_i-8}{12-8}) ρ'_{rim}(8) + \frac{R_i-8}{12-8} ρ^* & \text{if }\, 8 < R_i ≤ 12,
    \end{cases}
    \end{align}
    ```
    where $ρ^* = 916.7$ kg/m³ is the density of solid bulk ice, and $a_\text{CL}, b_\text{CL}, c_\text{CL}$ are the coefficients for the Cober & List (1993) parameterization [CoberList1993](@cite).
    The $R_i$ quantity is limited to the range $1 ≤ R_i ≤ 12$.
    Their values are $a_\text{CL} = 51$, $b_\text{CL} = 114$, and $c_\text{CL} = -5.5$.

The change in rime volume at some diameter $D_i$ is then given by
```math
∂_t\mathcal{B}_{\text{rim},l}(D_i) 
  = ∫_0^∞ \frac{∂_t\mathcal{V}_l(D_i, D_l) m_l(D_l)}{ρ'_{rim}(R_i(D_i, D_l))} N_l(D_l) \mathrm{d}D_l
```

!!! note "Bulk rime volume change from collisions"
    The bulk rate of rime volume change is limited by the fraction of mass that freezes, and is then given by
    ```math
    \textcolor{brown}{\text{BCCOL} + \text{BRCOL}} 
      = \sum_{l ∈ \{c, r\}} ∫_0^∞ ∂_t\mathcal{B}_{\text{rim},l}(D_i) f_\text{frz}(D_i) N'_i(D_i) \mathrm{d}D_i
    ```

In the wet growth regime, the rime compacts, rapidly relaxing towards the solid bulk ice density of $916.7$ kg/m³.
Because the P3 scheme only tracks the bulk rime volume $B_\text{rim}$, the densification rate is applied to the bulk rime volume.
The intensity of wet growth is measured by the fraction of the collected liquid mass that is shed,
```math
f_\text{shd}
= \frac{\textcolor{brown}{\text{QCSHD} + \text{QRSHD}}}
{\textcolor{brown}{\text{QCFRZ} + \text{QCSHD} + \text{QRFRZ} + \text{QRSHD}}},
```
which is zero when all collected liquid freezes, and approaches one when little of it freezes.
Densification requires freezing, so $f_\text{shd}$ is zero at and above the freezing temperature.

!!! note "Bulk rime volume change from wet growth"
    The bulk rate of rime volume change is then given by
    ```math
    \begin{align*}
    \textcolor{brown}{\text{QIWET}}
      &= f_\text{shd} ⋅ \frac{1}{τ_\text{wet}} (L_\text{ice} - L_\text{rim})
      = f_\text{shd} ⋅ \frac{1}{τ_\text{wet}}  L_\text{ice} (1 - F_{rim})
    \\
    \textcolor{brown}{\text{BIWET}}
      &= f_\text{shd} ⋅ \frac{1}{τ_\text{wet}} \left(\frac{L_\text{ice}}{ρ^*} - B_\text{rim}\right) \\
    \end{align*}
    ```

### Heterogeneous Freezing

Immersion freezing is parameterized based on water activity and follows the ABIFM
  parameterization from [KnopfAlpert2013](@cite).
See also the derivation notes about different
  [ice nucleation parameterizations](https://clima.github.io/CloudMicrophysics.jl/dev/IceNucleation/).
The immersion freezing nucleation rates for ice number and ice content are
  computed from the available cloud droplet number concentration and
  cloud water content:
```math
\frac{dN}{dt} = J_\textrm{ABIFM} \, A_a \, N_{lcl}
```
```math
\frac{dL}{dt} = J_\textrm{ABIFM} \, A_a \, q_{lcl} \, \rho_a
```
where
- ``J_\textrm{ABIFM}`` - is the immersion freezing nucleation rate coefficient,
- ``A_a`` - is the assumed surface area of insoluble ice nucleating particles
  per droplet (currently a constant ``10^{-10}`` m``^2``),
- ``N_{lcl}`` and ``q_{lcl}`` - cloud droplet number concentration and
  cloud water specific content,
- ``\rho_a`` - air density.

```@example
include("plots/P3ImmersionFreezing.jl")
```
![](P3_het_ice_nucleation.svg)

### Melting

Melting follows from heat conduction to the particles, as in the
  [1-moment scheme](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/#Snow-melt).
A particle of diameter ``D`` and capacitance ``C = D/2`` melts at
```math
\frac{dm}{dt} = \frac{2 \pi D \, K_\mathrm{thermo}}{L_f} \left(T - T_\mathrm{freeze}\right) F_v(D),
```
with ``L_f`` evaluated at ``T_\mathrm{freeze}``, the temperature of the particle surface,
  and the ventilation factor of [SeifertBeheng2006](@cite) and the terminal velocity of [Chen2022](@cite).
Mass and number melt at the same fractional rate, that of the particle size distribution:
```math
f = \frac{\int_{0}^{\infty} \frac{dm}{dt} N(D) \, \mathrm{d}D}{\int_{0}^{\infty} m(D) N(D) \, \mathrm{d}D}, \qquad
\left. \frac{dL}{dt} \right|_\mathrm{melt} = L f, \qquad
\left. \frac{dN}{dt} \right|_\mathrm{melt} = N f.
```
The size distribution solves the shape problem, so its mass equals ``L``, and ``L f`` is the integral of ``dm/dt`` over the distribution.
``N`` is the prognostic ice number, so that a step removes at most the number that exists, while the size distribution has the number bounded by the mean particle mass limits.
The number tendency of ice self-collection scales with the prognostic number in the same way.

```@example
include("plots/P3Melting.jl")
```
![](P3_ice_melt.svg)

## [Lookup tables](@id P3-lookup-tables)

The size-distribution integrals and the slope solve can be read from lookup tables instead of evaluated by quadrature.
[`CloudMicrophysics.P3Scheme.with_tables`](@ref) builds the tables for a complete set of parameters and returns a copy of the parameters whose `integration` is a [`CloudMicrophysics.P3Scheme.P3Tables`](@ref):

```julia
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

mp = CMP.Microphysics2MParams(Float32; with_ice = true)
mp = P3.with_tables(mp)
```

A [`CloudMicrophysics.P3Scheme.P3State`](@ref) built from these parameters reads every tabulated integral from the tables, so the process rates need no other change.
Quadrature stays the default, and it is the generator of the tables.
The tables depend on the mass, area and slope laws, the terminal velocities and the cloud and rain size distributions, so a change to any of them needs new tables.
The tables store a fingerprint of these parameters, and [`CloudMicrophysics.P3Scheme.tables_are_current`](@ref) returns `false` for parameters that differ from the ones the tables were built from.
The constructors of the parameters from a parameter dictionary throw an `ArgumentError` for tables built from other parameters.
The tables support the cloud size distribution `CloudParticlePDF_SB2006` and the limited rain size distribution `RainParticlePDF_SB2006_limited`, and the construction throws an `ArgumentError` for any other.

### What is tabulated

| Tables | Values | Axes |
|:-------|:-------|:-----|
| slope | ``\log(L/N)`` of the size distribution | ``\log λ``, ``w``, ``ρ_{rim}`` |
| ice | ventilation moment ``M_v``, self-collection rate, ``v_n``, ``v_m``, per unit number | ``\log λ``, ``w``, ``ρ_{rim}``, ``\log ρ_a`` |
| liquid-ice | collected number and mass of cloud or rain per unit ice number and per unit liquid number or mass, impact scale ``\log \bar{κ}`` | the ice axes and ``\log \bar{x}`` of the liquid |
| rime density | mean local rime density of the collected mass | ``\log R``, the ice axes, ``\log \bar{x}`` |
| wet set | mass one particle collects per unit liquid mass and its impact scale, ``D^{3/2} v^{1/2}``, mean local rime density of the mass one particle collects, fractions of the integrals below a diameter | ``\log D`` or the ice axes, with ``\log \bar{x}``, ``\log R``, the number logit and ``μ`` |

The rime coordinate ``w = 1 - \sqrt{1 - F_{rim}}`` spreads the rime mass fractions near one, where the particle properties change fastest.
The slope ``\log λ`` is the key of the integral tables, because the shape parameter ``μ`` changes steeply with the mean mass between ``\log \bar{x} = -20.2`` and ``-19.8``.
The slope comes from the inverse of the slope table, by a bisection with a fixed number of steps and one Newton step on the closed form.
The first moment ``M_1`` of the ventilation has a closed form, and so does the integral of ``D`` over the wet set.
The local rime density depends on the temperature only through ``R = κ / |T|``, so the rime density tables cover every temperature below freezing on ``\log R``.
No table holds an air property other than the air density: the ventilation coefficients, the viscosity, the diffusivity and the latent heats enter at run time.
The freezing of liquid by Gauss-Laguerre quadrature and the effective radius are not tabulated.

As in the quadrature, each ice particle freezes the fraction ``\min(m / (u + v), 1)`` of the liquid mass ``u + v`` it collects from cloud and rain, with ``m`` its maximum freezing rate, and sheds the rest.
The tables locate the wet set, where a particle sheds, from the mass one particle collects and its maximum freezing rate.
The fraction tables give the total frozen mass, the mass collected outside the wet set plus the maximum freezing rate on it, and the total shed mass.
The division of the frozen and the shed mass between cloud and rain, and the mean rime density of the frozen mass, follow from a Gauss-Legendre rule of 8 points in ``\log D`` on each interval between the ends of the wet set.

The slope ``\log λ`` is held under differentiation, as for the quadrature.
The other coordinates of the tables, the rime coordinate, the rime density, the air density, the liquid mean masses and ``\log R``, keep their derivatives through the multilinear interpolation, so the derivatives with respect to the state follow those of the quadrature to the accuracy of the tables.
On the default grid, at a partly frozen state, the elasticities of the collected masses and numbers with respect to the state differ from those of the quadrature by less than 0.03, and those of the frozen mass and the rime volume by less than 3%.
A coordinate outside a table reads its nearest edge, in particular an air density outside ``0.005`` to ``1.7`` kg/m³, and the derivative along that axis is zero there.

### Accuracy

The default grid, [`CloudMicrophysics.P3Scheme.P3TableGrid`](@ref), was chosen from accuracy measurements against Gauss-Legendre quadrature of order 16 on 2000 seeded states between the grid points.
The states have:

- the mean ice mass log-uniform between the limits of the ice number, and the rime density uniform between its bounds,
- the rime mass fraction zero for a fifth of the states, log-uniform from ``10^{-8}`` to ``10^{-2}`` for a fifth, uniform from 0.01 to 0.9 for two fifths, and with ``1 - F_{rim}`` log-uniform from ``10^{-3}`` to ``10^{-1}`` for a fifth,
- the air density log-uniform from 0.005 to 1.7 kg/m³,
- for the collisions, the ice number log-uniform from ``10^2`` to ``10^6`` m⁻³, the cloud mean mass log-uniform from ``5 \cdot 10^{-15}`` to ``2.5 \cdot 10^{-10}`` kg, the cloud number log-uniform from ``10^7`` to ``10^9`` m⁻³, the rain content log-uniform from ``10^{-6}`` to ``10^{-3}`` kg/m³, the rain number log-uniform from ``10^2`` to ``10^6`` m⁻³, and the temperature uniform from 235 to 278 K.

The relative errors in `Float64` are

| Quantity | Median | 99th percentile | Maximum |
|:---------|:-------|:----------------|:--------|
| ``\log λ`` (absolute) | ``8 \cdot 10^{-10}`` | ``4 \cdot 10^{-4}`` | ``1.5 \cdot 10^{-3}`` |
| Ventilation integral | 0.003% | 0.24% | 0.57% |
| Self-collection | 0.026% | 1.0% | 3.2% |
| ``v_n`` | 0.030% | 0.79% | 1.3% |
| ``v_m`` | 0.027% | 0.88% | 1.8% |
| Melting | 0.003% | 0.23% | 0.57% |
| Collected cloud mass | 0.12% | 3.0% | 4.7% |
| Collected rain mass | 0.11% | 2.5% | 7.4% |
| Collected cloud number | 0.13% | 2.8% | 5.1% |
| Collected rain number | 0.086% | 1.8% | 4.5% |
| Frozen mass, ice mass tendency | 0.042% | 2.0% | 6.5% |
| Rime volume | 0.14% | 2.6% | 8.2% |
| Rime volume of cloud | 0.69% | 7.1% | 22% |
| Rime volume of rain | 0.084% | 3.6% | 17% |
| Rime volume tendency | 0.14% | 2.6% | 16% |
| Shed fraction (absolute) | ``5 \cdot 10^{-6}`` | ``5 \cdot 10^{-3}`` | 0.058 |

On the 990 states that freeze and shed more than 0.1% of the collected mass each, the errors of the frozen and the shed mass of each species, relative to the mass the species gives to the ice, are

| Quantity | Median | 99th percentile | Maximum |
|:---------|:-------|:----------------|:--------|
| Frozen cloud | 0.068% | 2.9% | 7.9% |
| Shed cloud | 0.061% | 2.8% | 7.3% |
| Frozen rain | 0.020% | 2.1% | 4.3% |
| Shed rain | 0.082% | 1.6% | 5.0% |

The errors in `Float32` agree with these to two digits.
The rime volume of cloud has the largest errors.
The mass one particle collects from cloud has a sharp minimum at the diameter where the particle falls as fast as most of the droplets, and the per-particle tables interpolate across this minimum.
The tendencies sum the rime volumes of cloud and rain, so their errors are those of the total rime volume.

### Building the tables

The tables are built at initialization, from the parameters of the run.
[`CloudMicrophysics.P3Scheme.with_tables`](@ref) fills them in two layers, and both are public.
The building blocks are:

1. [`CloudMicrophysics.P3Scheme.P3TableGrid`](@ref), the axes of every table, with bounds derived from the parameters,
2. [`CloudMicrophysics.P3Scheme.table_point_functions`](@ref), one point function per table group, which maps the coordinates of one grid point to the values of its tables by quadrature in `Float64`,
3. [`CloudMicrophysics.P3Scheme.table_axes`](@ref), the axes of the point function of each group,
4. [`CloudMicrophysics.P3Scheme.tabulate`](@ref), which evaluates a point function at the points of its axes in an array of a given type,
5. [`CloudMicrophysics.P3Scheme.P3Tables`](@ref) and the method `with_tables(mp, grid, values)`, which assemble the arrays of values into the tables in the float type of the parameters, with the fingerprint of the parameters.

The point functions are `isbits` callable objects, and `tabulate` computes the coordinates of a point from its index, so no array of coordinates is stored.
`with_tables(mp, ArrayType)` tabulates each group in an array of type `ArrayType` and stores it as a table before it fills the next group, so the tables are built where `ArrayType` places its data and they are stored there.
For an `Array`, the default, the points are evaluated on all threads.
A host with GPU arrays, such as ClimaAtmos with `ClimaComms.array_type(device)`, builds the tables on the device without a dependency of CloudMicrophysics on a GPU package:

```julia
mp = P3.with_tables(mp, ClimaComms.array_type(device))
```

A host that maps the point functions itself uses the building blocks:

```julia
grid = P3.P3TableGrid(mp)
points = P3.table_point_functions(mp.ice, grid)
values = map((f, axes) -> P3.tabulate(f, axes, ArrayType), points, P3.table_axes(grid))
mp = P3.with_tables(mp, grid, values)
@assert P3.tables_are_current(mp)
```

On the default grid, the tables take 176 MB in `Float64` and 88 MB in `Float32`.
The build is dominated by the liquid-ice collision tables, a double integral at each of 4.7 million points, and by the ice tables, whose self-collection integral is a double integral at each of 0.9 million points:

| Group | Points | Build, 16 threads | Build, 1 thread |
|:------|-------:|------------------:|----------------:|
| slope | 32768 | 0.01 s | 0.03 s |
| ice | 912384 | 9.7 s | 109 s |
| cloud collisions | 3133440 | 13.4 s | 154 s |
| rain collisions | 1566720 | 12.2 s | 143 s |
| rime density | 22400 | 0.8 s | 6.2 s |
| wet set | 422464 | 4.1 s | 41 s |
| all, with the assembly | | 40 s | 455 s |

These are times on an Apple M4 Max with 12 performance and 4 efficiency cores, the same in `Float32` and `Float64`.
The tables evaluate the rates 5 to 12000 times faster than the quadrature on one CPU thread:

| Function | Quadrature | Tables |
|:---------|-----------:|-------:|
| `P3State`, with the slope | 2.9 μs | 0.53 μs |
| `ice_terminal_velocities` | 6.5 μs | 0.03 μs |
| `ice_self_collection` | 343 μs | 0.03 μs |
| `ice_melt` | 6.8 μs | 0.31 μs |
| `∫liquid_ice_collisions`, below freezing | 210 to 217 μs | 9.3 μs |
| `∫liquid_ice_collisions`, above freezing | 172 μs | 0.30 μs |

These are `Float32` times with the state built, except for `P3State`.

## Acknowledgments

The reference Fortran implementation of the P3 scheme, developed by the authors
  of [MorrisonMilbrandt2015](@cite) and others, is available in the
  [P3-microphysics repository](https://github.com/P3-microphysics/P3-microphysics).
