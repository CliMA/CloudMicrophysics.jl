# P3 Scheme
TO DO: This is where P3 documentation will be.

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
 - ``\lamda `` is the slope parameter.

and ``\mu \ = 0.00191 \lambda \ ^{0.8} - 2`` for ``\mu \ \in (0,6)``.

## Assumed particle mass relationships

The mass (m) of particles as a function of maximum particle dimension (D)
 is a piecewise function with variable breakpoints described
 by the following table.

| particle properties |      conditions       |    m(D) relation      |
|---------------------|-----------------------|-----------------------|
|small, spherical ice | ``D < D_th`` | ``\frac{\pi}{6} \rho_i \ D^3`` |
|large, unrimed ice   | ``q_i > 0``, ``q_rim = 0``, and ``D > D_th`` | ``\alpha_va \ D^{\beta_va}`` |
|dense nonspherical ice | ``q_rim > 0`` and ``D_th < D < D_gr`` | ^ |
|partially rimed ice | ``q_rim > 0`` and ``D < D_cr`` | ``\frac{\alpha_va}{1-F_r} D^{\beta_va}`` |
|graupel | ``q_rim > 0`` and ``D < D_cr`` | ``\frac{\pi}{6} \rho_g \ D^3`` |

## A(D)

## Velocity
