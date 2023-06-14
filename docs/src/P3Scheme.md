# P3 Scheme

The `P3Scheme.jl` module implements the predicted particle properties
 (P3) scheme for ice-phase microphysics developed by [cite articles].

The P3 scheme is a 2-moment, bulk scheme involving a
 single ice-phase category with 4 degrees of freedom: total mass,
 rime mass, rime volume, and number mixing ratios.

Traditionally, cloud ice microphysics models use various predefined
 categories to represent ice modes, but the P3 scheme sidesteps the
 problem of prescribing transitions between ice categories by adopting
 a single ice category. This simplification might
 aid in attempts to constrain scheme parameters.

## Assumed particle size distribution (PSD)

Following [Morrison and Milbrandt 2015], the scheme assumes a
 gamma distribution for the concentration of particles per unit volume:

```math
N'(D) = N_{0} D^\mu \, e^{-\lambda \, D}
```

where:
 - ``N'`` is the number concentration,
 - ``D`` is the maximum particle dimension (m),
 - ``N_0`` is the intercept parameter,
 - ``\mu `` is the shape parameter,
 - ``\lambda `` is the slope parameter.

and ``\mu \ = 0.00191 \lambda \ ^{0.8} - 2`` for ``\mu \ \in (0,6)``.

## Assumed particle mass relationships

The mass (m) of particles as a function of maximum particle dimension (D)
 is a piecewise function with variable breakpoints described
 by the following table.

| particle properties |      condition(s)     |    m(D) relation      |
|:--------------------|:----------------------|:----------------------|
|small, spherical ice | ``D < D_{th}`` | ``\frac{\pi}{6} \rho_i \ D^3`` |
|large, unrimed ice   | ``q_{rim} = 0`` and ``D > D_{th}`` | ``\alpha_{va} \ D^{\beta_{va}}`` |
|dense nonspherical ice | ``q_{rim} > 0`` and ``D_{th} < D < D_{gr}`` | ``\alpha_{va} \ D^{\beta_{va}}`` |
|partially rimed ice | ``q_{rim} > 0`` and ``D < D_{cr}`` | ``\frac{\alpha_{va}}{1-F_r} D^{\beta_{va}}`` |
|graupel (completely rimed, spherical)| ``q_{rim} > 0``, ``D < D_{cr}``, and ``D > D_{gr}`` | ``\frac{\pi}{6} \rho_g \ D^3`` |

> **_NOTE:_**  unsure about D_cr: should graupel be defined by D > D_cr? and partially rimed ice?

where:
 - ``\rho_i \ = 917 kgm^{-3}`` (density of bulk ice)
 - ``D_{th} = \frac{\pi \ \rho_i\}{6\alpha_{va}}^{\frac{1}{\beta_{va}\ - 3}}``,
  the threshold between spherical and nonspherical ice;
 - ``q_{rim}`` is the rime mass mixing ratio for ice;
  (should this be specific humidity for us?)
 - ``\alpha_{va} = 7.38*10^{-11}``, a parameter from [cite:Brown and Francis 1995],
  derived from measurements of mass grown
  by vapor diffusion and aggregation in midlatitude cirrus;
 - ``\beta_{va} = 1.9``, another parameter from [cite:Brown and Francis 1995];
 - ``D_{gr} = \frac{6\alpha_{va}}{\pi \\rho_g \}^{\frac{1}{3 - \beta_{va}\}}``,
  a threshold defined to ensure continuity
  and reasonable masses for smaller D;
 - ``D_{cr} = [(\frac{1}{1-F_r}\frac{6\alpha_{va}}{\pi \\rho_g \}]^{\frac{1}{3 - \beta_{va}\}}``,
  the threshold separating partially rimed ice from graupel;
 - ``\rho_g \ = \rho_r \ F_r + (1 - F_r)\rho_d \``,
  the rime mass fraction (``F_r = \frac{q_{rim}}{q_i,tot}``)
  weighted average of the predicted rime density
  ``\rho_r \ = \frac{q_{rim}}{B_{rim}}``
  and the density of the unrimed part
  ``\rho_d \ = \frac{6\alpha_{va}(D_{cr}^{\beta{va} \ - 2} - D_{gr}^{\beta{va} \ - 2})}{\pi \ (\beta_{va} \ - 2)(D_{cr}-D_{gr})}``

> **_NOTE:_**  unsure about what B_rim is
## Assumed projected area relationships

The projected area (A) of particles as a function of maximum particle dimension (D)
 is another piecewise function with variable breakpoints described
 by the following table.

| particle properties |      condition(s)     |    A(D) relation      |
|:--------------------|:----------------------|:----------------------|
|small, spherical ice | ``D < D_{th}``        | ``\frac{\pi}{4} D^2`` |
|graupel (completely rimed, spherical)| ``q_{rim} > 0`` and ``D_{gr} < D < D_{cr}`` | ``\frac{\pi}{4} D^2`` |
|large, unrimed ice   | ``q_{rim} = 0`` and ``D > D_{th}`` | see below |
|dense nonspherical ice | ``q_{rim} > 0`` and ``D_{th} < D < D_{gr}`` | see p7 of mitchell 1996 and do aggregates |
|partially rimed ice | ``q_{rim} > 0`` and ``D < D_{cr}`` | ...todo |
## Velocity

> **_TODO:_**  finish documentation and add more sections