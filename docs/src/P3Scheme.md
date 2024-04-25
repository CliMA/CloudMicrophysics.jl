# P3 Scheme

The `P3Scheme.jl` module implements the predicted particle properties
 (P3) scheme for ice-phase microphysics developed by [MorrisonMilbrandt2015](@cite).
The P3 scheme is a 2-moment, bulk scheme involving a
 single ice-phase category with 4 degrees of freedom: total mass,
 rime mass, rime volume, and number mixing ratios.
Traditionally, cloud ice microphysics schemes use various predefined
 categories (such as ice, graupel, or hail) to represent ice modes, but the P3 scheme sidesteps the
 problem of prescribing transitions between ice categories by adopting
 a single ice category and evolving its properties.
This simplification aids in attempts to constrain the scheme's free parameters.

The prognostic variables are:
 - ``N_{ice}`` - number concentration 1/m3
 - ``q_{ice}`` - ice mass density kg/m3
 - ``q_{rim}`` - rime mass density kg/m3
 - ``B_{rim}`` - rime volume - (volume of rime per total air volume: dimensionless)

!!! note
    TODO - At some point we should switch to specific humidities...

## Assumed particle size relationships

The mass ``m`` and projected area ``A`` of particles
  as a function of maximum particle dimension ``D``
  are piecewise functions with variable thresholds described
  by the following table.

| particle properties                  |      condition(s)                            |    m(D) relation                             |    A(D) relation                                           |
|:-------------------------------------|:---------------------------------------------|:---------------------------------------------|:-----------------------------------------------------------|
|small, spherical ice                  | ``D < D_{th}``                               | ``\frac{\pi}{6} \rho_i \ D^3``               | ``\frac{\pi}{4} D^2``                                      |
|large, unrimed ice                    | ``q_{rim} = 0`` and ``D > D_{th}``           | ``\alpha_{va} \ D^{\beta_{va}}``             | ``\gamma \ D^{\sigma}``                                    |
|dense nonspherical ice                | ``q_{rim} > 0`` and ``D_{gr} > D > D_{th}``  | ``\alpha_{va} \ D^{\beta_{va}}``             | ``\gamma \ D^{\sigma}``                                    |
|graupel (completely rimed, spherical) | ``q_{rim} > 0`` and ``D_{cr} > D > D_{gr}``  | ``\frac{\pi}{6} \rho_g \ D^3``               | ``\frac{\pi}{4} D^2``                                      |
|partially rimed ice                   | ``q_{rim} > 0`` and ``D > D_{cr}``           | ``\frac{\alpha_{va}}{1-F_r} D^{\beta_{va}}`` | ``F_{r} \frac{\pi}{4} D^2 + (1-F_{r})\gamma \ D^{\sigma}`` |

where:
 - ``D_{th}``, ``D_{gr}``, ``D_{cr}`` are particle size thresholds in ``m``,
 - ``\rho_i`` is cloud ice density in ``kg m^{-3}``,
 - ``\beta_{va} = 1.9`` is a dimensionless parameter from [BrownFrancis1995](@cite) (based on measurements of vapor diffusion and aggregation in midlatitude cirrus),
 - ``\alpha_{va} = 7.38 \; 10^{-11} \; 10^{6 \beta_{va} - 3}`` in ``kg \; m^{-β_{va}}`` is a parameter modified for units from [BrownFrancis1995](@cite) in base SI units (also based on measurements of vapor diffusion and aggregation in midlatitude cirrus),
 - ``\rho_g`` is the bulk density of graupel in ``kg \; m^{-3}``
 - ``\gamma = 0.2285`` (``m^{2 - \sigma}``) where
 - ``\sigma = 1.88`` (dimensionless), both from the aggregates of side planes, columns, bullets, and planar polycrystals in [Mitchell1996](@cite).

The first threshold is solely determined by the free parameters:
  ``D_{th} = (\frac{\pi \rho_i}{6\alpha_{va}})^{\frac{1}{\beta_{va} - 3}}``.
The remaining thresholds: ``D_{gr}``, ``D_{cr}``, as well as the
  bulk density of graupel ``\rho_{g}``,
  and the bulk density of the unrimed part ``\rho_d``
  form a nonlinear system:
 - ``D_{gr} = (\frac{6\alpha_{va}}{\pi \rho_g})^{\frac{1}{3 - \beta_{va}}}``
 - ``D_{cr} = [ (\frac{1}{1-F_r}) \frac{6 \alpha_{va}}{\pi \rho_g} ]^{\frac{1}{3 - \beta_{va}}}``
 - ``\rho_g = \rho_r F_r + (1 - F_r) \rho_d``
 - ``\rho_d = \frac{6\alpha_{va}(D_{cr}^{\beta{va} \ - 2} - D_{gr}^{\beta{va} \ - 2})}{\pi \ (\beta_{va} \ - 2)(D_{cr}-D_{gr})}``
where
 - ``F_r = \frac{q_{rim}}{q_{ice}}`` is the rime mass fraction,
 - ``\rho_{r} = \frac{q_{rim}}{B_{rim}}`` is the predicted rime density.

!!! note
    TODO - Check units, see in [issue #151](https://github.com/CliMA/CloudMicrophysics.jl/issues/151)

Below we show the m(D) and a(D) regimes replicating Figures 1 (a) and (b) from [MorrisonMilbrandt2015](@cite).
```@example
include("plots/P3SchemePlots.jl")
```
![](P3Scheme_relations.svg)

## Assumed particle size distribution

Following [MorrisonMilbrandt2015](@cite), the scheme assumes a
 gamma distribution for the concentration of ice particles per unit volume
 based on particle size measurements obtained by [Heymsfield2003](@cite)
 in tropical and midlatitude ice clouds and implemented by
 [MorrisonGrabowski2008](@cite):

```math
N'(D) = N_{0} D^\mu \, e^{-\lambda \, D}
```
where:
 - ``N'`` is the number concentration in ``m^{-4}``
 - ``D`` is the maximum particle dimension in ``m``,
 - ``N_0`` is the intercept parameter in ``m^{-5 - \mu }``,
 - ``\mu`` is the shape parameter (dimensionless),
 - ``\lambda`` is the slope parameter in ``m^{-1}``.

We assume ``\mu \ = 0.00191 \; \lambda \ ^{0.8} - 2``.
Following [MorrisonGrabowski2008](@cite) we limit ``\mu \ \in (0,6)``.
A negative ``\mu`` can occur only for very small mean particle sizes``\frac{1}{\lambda} < ~0.17 mm``.

The model predicted ice number concentration and ice mass density are defined as
```math
N_{ice} = \int_{0}^{\infty} \! N'(D) \mathrm{d}D
```
```math
q_{ice} = \int_{0}^{\infty} \! m(D) N'(D) \mathrm{d}D
```

## Calculating shape parameters

As a next step we need to find the mapping between
  predicted moments of the size distribution ``N_{ice}`` and ``q_{ice}``
  and the shape parameters ``\lambda`` and ``N_0``.
Solving for ``N_{ice}`` is relatively straightforward:
```math
N_{ice} = \int_{0}^{\infty} \! N'(D) \mathrm{d}D = \int_{0}^{\infty} \! N_{0} D^\mu \, e^{-\lambda \, D} \mathrm{d}D = N_{0} (\lambda \,^{-(\mu \, + 1)} \Gamma \,(1 + \mu \,))
```

``q_{ice}`` depends on the variable mass-size relation ``m(D)`` defined above.
We solve for ``q_{ice}`` in a piece-wise fashion defined by the same thresholds as ``m(D)``.
As a result ``q_{ice}`` can be expressed as a sum of inclomplete gamma functions,
  and the shape parameters are found using iterative solver.

|      condition(s)                            |    ``q_{ice} = \int \! m(D) N'(D) \mathrm{d}D``                                          |         gamma representation          |
|:---------------------------------------------|:-----------------------------------------------------------------------------------------|:---------------------------------------------|
| ``D < D_{th}``                               | ``\int_{0}^{D_{th}} \! \frac{\pi}{6} \rho_i \ D^3 N'(D) \mathrm{d}D``                    | ``\frac{\pi}{6} \rho_i N_0 \lambda \,^{-(\mu \, + 4)} (\Gamma \,(\mu \, + 4) - \Gamma \,(\mu \, + 4, \lambda \,D_{th}))``|
| ``q_{rim} = 0`` and ``D > D_{th}``           | ``\int_{D_{th}}^{\infty} \! \alpha_{va} \ D^{\beta_{va}} N'(D) \mathrm{d}D``             | ``\alpha_{va} \ N_0 \lambda \,^{-(\mu \, + \beta_{va} \, + 1)} (\Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{th}))`` |
| ``q_{rim} > 0`` and ``D_{gr} > D > D_{th}``  | ``\int_{D_{th}}^{D_{gr}} \! \alpha_{va} \ D^{\beta_{va}} N'(D) \mathrm{d}D``             | ``\alpha_{va} \ N_0 \lambda \,^{-(\mu \, + \beta_{va} \, + 1)} (\Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{th}) - \Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{gr}))`` |
| ``q_{rim} > 0`` and ``D_{cr} > D > D_{gr}``  | ``\int_{D_{gr}}^{D_{cr}} \! \frac{\pi}{6} \rho_g \ D^3 N'(D) \mathrm{d}D``               | ``\frac{\pi}{6} \rho_g N_0 \lambda \,^{-(\mu \, + 4)} (\Gamma \,(\mu \, + 4, \lambda \,D_{gr}) - \Gamma \,(\mu \, + 4, \lambda \,D_{cr}))`` |
| ``q_{rim} > 0`` and ``D > D_{cr}``           | ``\int_{D_{cr}}^{\infty} \! \frac{\alpha_{va}}{1-F_r} D^{\beta_{va}} N'(D) \mathrm{d}D`` | ``\frac{\alpha_{va}}{1-F_r} N_0 \lambda \,^{-(\mu \, + \beta_{va} \, + 1)} (\Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{cr}))``  |

where ``\Gamma \,(a, z) = \int_{z}^{\infty} \! t^{a - 1} e^{-t} \mathrm{d}D``
  and ``\Gamma \,(a) = \Gamma \,(a, 0)`` for simplicity.

In our solver, we approximate ``\mu`` from ``q/N`` and keep it constant throughout the solving step.
We approximate ``\mu`` by an exponential function given by the ``q/N`` points
  corresponding to ``\mu = 6`` and ``\mu = 0``.
This is shown below as well as how this affects the solvers ``\lambda`` solutions.

```@example
include("plots/P3LambdaErrorPlots.jl")
```
![](MuApprox.svg)

An initial guess for the non-linear solver is found by approximating the gamma functions
  as a simple linear function from ``x = \log{(q/N)}`` to ``y = \log{(\lambda)}``.
The equation is given by ``(x - x_1) = A \; (y - y_1)`` where
  ``A = \frac{x_1 - x_2}{y_1 - y_2}``.
``y_1`` and ``y_2`` define ``log(\lambda)`` values for three ``q/N`` ranges

|             q/N                |         ``y_1``        |      ``y_2``         |
|:-------------------------------|:-----------------------|:---------------------|
|        ``q/N >= 10^-8``        |         ``1``          |     ``6 * 10^3``     |
|   ``2 * 10^9 <= q/N < 10^-8``  |     ``6 * 10^3``       |     ``3 * 10^4``     |
|        ``q/N < 2 * 10^9``      |     ``4 * 10^4``       |       ``10^6``       |

We use this approximation to calculate initial guess for the shape parameter
  ``\lambda_g = \lambda_1 (\frac{q}{q_1})^{(\frac{y_1 - y_2}{x_1 - x_2})}``.

```@example
include("plots/P3ShapeSolverPlots.jl")
```
![](SolverInitialGuess.svg)

Using this approach we get the following relative errors for ``\lambda``

```@example
include("plots/P3LambdaErrorPlots.jl")
```
![](P3LambdaHeatmap.png)

## Terminal Velocity

We use the [Chen2022](@cite) velocity parametrization:
```math
V(D) = \phi^{\kappa} \sum_{i=1}^{j} \; a_i D^{b_i} e^{-c_i \; D}
```
where ``\phi = (16 \rho_{ice}^2 A(D)^3) / (9 \pi m(D)^2)`` is the aspect ratio,
and ``\kappa``, ``a_i``, ``b_i`` and ``c_i`` are the free parameters.

The mass-weighted fall speed (``V_m``) and the number-weighted fall speed (``V_n``) are calculated as
```math
V_m = \frac{\int_{0}^{\infty} \! V(D) m(D) N'(D) \mathrm{d}D}{\int_{0}^{\infty} \! m(D) N'(D) \mathrm{d}D}
```
```math
V_n = \frac{\int_{0}^{\infty} \! V(D) N'(D) \mathrm{d}D}{\int_{0}^{\infty} \! N'(D) \mathrm{d}D}
```

We also plot the mass-weighted mean particle size ``D_m`` which is given by:
```math
D_m = \frac{\int_{0}^{\infty} \! D m(D) N'(D) \mathrm{d}D}{\int_{0}^{\infty} \! m(D) N'(D) \mathrm{d}D}
```

Below w show these relationships for small, medium, and large ``D_m``
They can be compared with Figure 2 from [MorrisonMilbrandt2015](@cite).
```@example
include("plots/P3TerminalVelocityPlots.jl")
```
![](MorrisonandMilbrandtFig2.svg)
