# P3 Scheme

The `P3Scheme.jl` module implements the predicted particle properties
 (P3) scheme for ice-phase microphysics developed by [MorrisonMilbrandt2015](@cite).
The P3 scheme is a 2-moment, bulk scheme involving a
 single ice-phase category with 4 degrees of freedom: total ice content,
 rime content, rime volume, and number mixing ratios.
Traditionally, cloud ice microphysics schemes use various predefined
 categories (such as ice, graupel, or hail) to represent ice modes, but the P3 scheme sidesteps the
 problem of prescribing transitions between ice categories by adopting
 a single ice category and evolving its properties.
This simplification aids in attempts to constrain the scheme's free parameters.

The prognostic variables are:
 - ``N_{ice}`` - number concentration 1/m3
 - ``L_{ice}`` - ice mass content kg/m3
 - ``L_{rim}`` - rime mass content kg/m3
 - ``B_{rim}`` - rime volume - (volume of rime per total air volume: dimensionless)

!!! note
    The original paper [MorrisonMilbrandt2015](@cite) uses symbol ``q``
    to denote the mass of a tracer per volume of air (named mass mixing ratio).
    In our documentation of the 1-moment and 2-moment schemes we used ``q``
    to denote the mass of a tracer per mass of air (named specific humidity).
    To keep the notation consistent between the 1,2-moment schemes and P3,
    and to highlight the difference between normalizing by air mass or volume,
    we denote the mass of a tracer per volume of air as ``L`` (named content).

## Assumed particle size relationships

The mass ``m`` and projected area ``A`` of particles
  as a function of maximum particle dimension ``D``
  are piecewise functions with variable thresholds described
  by the following table.

| particle properties                  |      condition(s)                            |    m(D) relation                             |    A(D) relation                                           |
|:-------------------------------------|:---------------------------------------------|:---------------------------------------------|:-----------------------------------------------------------|
|small, spherical ice                  | ``D < D_{th}``                               | ``\frac{\pi}{6} \rho_i \ D^3``               | ``\frac{\pi}{4} D^2``                                      |
|large, unrimed ice                    | ``L_{rim} = 0`` and ``D > D_{th}``           | ``\alpha_{va} \ D^{\beta_{va}}``             | ``\gamma \ D^{\sigma}``                                    |
|dense nonspherical ice                | ``L_{rim} > 0`` and ``D_{gr} > D > D_{th}``  | ``\alpha_{va} \ D^{\beta_{va}}``             | ``\gamma \ D^{\sigma}``                                    |
|graupel (completely rimed, spherical) | ``L_{rim} > 0`` and ``D_{cr} > D > D_{gr}``  | ``\frac{\pi}{6} \rho_g \ D^3``               | ``\frac{\pi}{4} D^2``                                      |
|partially rimed ice                   | ``L_{rim} > 0`` and ``D > D_{cr}``           | ``\frac{\alpha_{va}}{1-F_rim} D^{\beta_{va}}`` | ``F_{rim} \frac{\pi}{4} D^2 + (1-F_{rim})\gamma \ D^{\sigma}`` |

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
 - ``D_{cr} = [ (\frac{1}{1-F_rim}) \frac{6 \alpha_{va}}{\pi \rho_g} ]^{\frac{1}{3 - \beta_{va}}}``
 - ``\rho_g = \rho_r F_rim + (1 - F_rim) \rho_d``
 - ``\rho_d = \frac{6\alpha_{va}(D_{cr}^{\beta{va} \ - 2} - D_{gr}^{\beta{va} \ - 2})}{\pi \ (\beta_{va} \ - 2)(D_{cr}-D_{gr})}``
where
 - ``F_rim = \frac{L_{rim}}{L_{ice}}`` is the rime mass fraction,
 - ``\rho_{r} = \frac{L_{rim}}{B_{rim}}`` is the predicted rime density.

!!! note
    TODO - Check units, see in [issue #151](https://github.com/CliMA/CloudMicrophysics.jl/issues/151)

Below we show the m(D) and a(D) regimes replicating Figures 1 (a) and (b) from [MorrisonMilbrandt2015](@cite).
We also show the density as a function of D.
Note that because graupel is completely filled with rime,
  the density (``\rho_{g}``) is independent of ``D`` between ``D_{gr}`` and ``D_{cr}``.
Following [MorrisonMilbrandt2015](@cite), for nonspherical particles
  ``\rho_{ice}``is assumed to be equal to the mass of the particle
  divided by the volume of a sphere with the same particle size
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

The model predicted ice number concentration and ice content are defined as
```math
N_{ice} = \int_{0}^{\infty} \! N'(D) \mathrm{d}D
```
```math
L_{ice} = \int_{0}^{\infty} \! m(D) N'(D) \mathrm{d}D
```

## Calculating shape parameters

As a next step we need to find the mapping between
  predicted moments of the size distribution ``N_{ice}`` and ``L_{ice}``
  and the shape parameters ``\lambda`` and ``N_0``.
Solving for ``N_{ice}`` is relatively straightforward:
```math
N_{ice} = \int_{0}^{\infty} \! N'(D) \mathrm{d}D = \int_{0}^{\infty} \! N_{0} D^\mu \, e^{-\lambda \, D} \mathrm{d}D = N_{0} (\lambda \,^{-(\mu \, + 1)} \Gamma \,(1 + \mu \,))
```

``L_{ice}`` depends on the variable mass-size relation ``m(D)`` defined above.
We solve for ``L_{ice}`` in a piece-wise fashion defined by the same thresholds as ``m(D)``.
As a result ``L_{ice}`` can be expressed as a sum of incomplete gamma functions,
  and the shape parameters are found using iterative solver.

|      condition(s)                            |    ``L_{ice} = \int \! m(D) N'(D) \mathrm{d}D``                                          |         gamma representation          |
|:---------------------------------------------|:-----------------------------------------------------------------------------------------|:---------------------------------------------|
| ``D < D_{th}``                               | ``\int_{0}^{D_{th}} \! \frac{\pi}{6} \rho_i \ D^3 N'(D) \mathrm{d}D``                    | ``\frac{\pi}{6} \rho_i N_0 \lambda \,^{-(\mu \, + 4)} (\Gamma \,(\mu \, + 4) - \Gamma \,(\mu \, + 4, \lambda \,D_{th}))``|
| ``L_{rim} = 0`` and ``D > D_{th}``           | ``\int_{D_{th}}^{\infty} \! \alpha_{va} \ D^{\beta_{va}} N'(D) \mathrm{d}D``             | ``\alpha_{va} \ N_0 \lambda \,^{-(\mu \, + \beta_{va} \, + 1)} (\Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{th}))`` |
| ``L_{rim} > 0`` and ``D_{gr} > D > D_{th}``  | ``\int_{D_{th}}^{D_{gr}} \! \alpha_{va} \ D^{\beta_{va}} N'(D) \mathrm{d}D``             | ``\alpha_{va} \ N_0 \lambda \,^{-(\mu \, + \beta_{va} \, + 1)} (\Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{th}) - \Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{gr}))`` |
| ``L_{rim} > 0`` and ``D_{cr} > D > D_{gr}``  | ``\int_{D_{gr}}^{D_{cr}} \! \frac{\pi}{6} \rho_g \ D^3 N'(D) \mathrm{d}D``               | ``\frac{\pi}{6} \rho_g N_0 \lambda \,^{-(\mu \, + 4)} (\Gamma \,(\mu \, + 4, \lambda \,D_{gr}) - \Gamma \,(\mu \, + 4, \lambda \,D_{cr}))`` |
| ``L_{rim} > 0`` and ``D > D_{cr}``           | ``\int_{D_{cr}}^{\infty} \! \frac{\alpha_{va}}{1-F_rim} D^{\beta_{va}} N'(D) \mathrm{d}D`` | ``\frac{\alpha_{va}}{1-F_rim} N_0 \lambda \,^{-(\mu \, + \beta_{va} \, + 1)} (\Gamma \,(\mu \, + \beta_{va} \, + 1, \lambda \,D_{cr}))``  |

where ``\Gamma \,(a, z) = \int_{z}^{\infty} \! t^{a - 1} e^{-t} \mathrm{d}D``
  and ``\Gamma \,(a) = \Gamma \,(a, 0)`` for simplicity.

In our solver, we approximate ``\mu`` from ``L/N`` and keep it constant throughout the solving step.
We approximate ``\mu`` by an exponential function given by the ``L/N`` points
  corresponding to ``\mu = 6`` and ``\mu = 0``.
This is shown below as well as how this affects the solvers ``\lambda`` solutions.

```@example
include("plots/P3LambdaErrorPlots.jl")
```
![](MuApprox.svg)

An initial guess for the non-linear solver is found by approximating the gamma functions
  as a simple linear function from ``x = \log{(L/N)}`` to ``y = \log{(\lambda)}``.
The equation is given by ``(x - x_1) = A \; (y - y_1)`` where
  ``A = \frac{x_1 - x_2}{y_1 - y_2}``.
``y_1`` and ``y_2`` define ``log(\lambda)`` values for three ``L/N`` ranges

|             L/N                |         ``y_1``        |      ``y_2``         |
|:-------------------------------|:-----------------------|:---------------------|
|        ``L/N >= 10^-8``        |         ``1``          |     ``6 * 10^3``     |
|   ``2 * 10^9 <= L/N < 10^-8``  |     ``6 * 10^3``       |     ``3 * 10^4``     |
|        ``L/N < 2 * 10^9``      |     ``4 * 10^4``       |       ``10^6``       |

We use this approximation to calculate initial guess for the shape parameter
  ``\lambda_g = \lambda_1 (\frac{L}{L_1})^{(\frac{y_1 - y_2}{x_1 - x_2})}``.

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

Note that ``\phi = 1`` corresponds to spherical particles
  (small spherical ice (``D < D_{th}``) and graupel (``D_{gr} < D < D_{cr}``)).
The plot provided below helps to visualize the transitions between spherical and nonspherical regimes.
```@example
include("plots/P3AspectRatioPlot.jl")
```
![](P3Scheme_aspect_ratio.svg)

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

Below we provide plots of these relationships for small, medium, and large ``D_m``:
  the first row highlights the particle size regime,
  the second displays ``D_m`` of the particles,
  the third shows the aspect ratio ``\phi (D_m)``,
  and the final row exhibits ``V_m``.
They can be compared with Figure 2 from [MorrisonMilbrandt2015](@cite).
```@example
include("plots/P3TerminalVelocityPlots.jl")
```
![](MorrisonandMilbrandtFig2.svg)

## Liquid Fraction

To allow for the modeling of mixed-phase particles with P3, a new prognostic variable can be introduced:
  ``L_{liq}``, the content of liquid on mixed-phase particles. As described in [Choletteetal2019](@cite),
  this addition to the framework of P3 opens the door to tracking the gradual melting (and refreezing) of a
  particle population with hydrometeor categories such as wet snow. Here, we describe the characteristics
  of the scheme with the addition of ``L_{liq}``.

Replicating Fig. 1 from [Choletteetal2019](@cite), we can expect more liquid to increase velocities of small particles for all ``F_rim`` values, while exhibiting more complex behavior for larger particles:

```@example
include("plots/Cholette2019_fig1.jl")
```
![](Choletteetal2019_fig1.svg)

Liquid fraction, analogous to ``F_rim`` for rime, is defined ``F_{liq} = \frac{L_{liq}}{L_{ice}}``.
  Importantly, the introduction of ``L_{liq}`` changes the formulation of ``F_rim``: whereas above,
  we have ``F_rim = \frac{L_{rim}}{L_{ice}}``, we now need ``F_rim = \frac{L_{rim}}{L_{ice} - L_{liq}}``
  since the total mass of ice particles ``L_{ice} = L_{rim} + L_{dep} + L_{liq}`` now includes liquid mass.

The addition of the liquid fraction does not change the thresholds ``D_{th}``, ``D_{gr}``, ``D_{cr}``,
  since the threshold regime depends only on ice core properties.

However, the assumed particle properties become ``F_{liq}``-weighted averages of particles' solid and liquid
  components:

```math
m(D, F_{liq}) = (1 - F_liq) * m(D, F_{liq} = 0) + F_liq * m_{liq}(D)
```
```math
A(D, F_{liq}) = (1 - F_liq) * A(D, F_{liq} = 0) + F_liq * A_{liq}(D)
```

where ``m_{liq}(D) = \frac{\pi}{6} \rho_liq * D^3`` and ``A_{liq}(D) = \frac{\pi}{4} D^2``.

When calculating shape parameters and integrating over the particle size distribution (PSD), it is important to
  keep in mind whether the desired moment of the PSD is tied only to ice (in which case we concern ourselves with the
  ice core diameter ``D_{core}``) or to the whole mixed-phase particle (in which case we need ``D_p``). If calculating
  the ice core parameters, ``F_{liq} = 0`` is passed into the solver framework as indicated.

For the above particle properties and for terminal velocity, we use the PSD corresponding to the whole mixed-phase particle,
  so our terminal velocity of a mixed-phase particle is:

```math
V(D_{p}, F_{liq}) = (1 - F_{liq}) * V_{ice}(D_{p}) + F_{liq} * V_{rain}(D_{p})
```

We continue to use the terminal velocity parameterizations from [Chen2022](@cite) for rain, which
  is described [here](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics2M/#Terminal-Velocity).

Visualizing mass-weighted terminal velocity as a function of ``F_{liq}``, ``F_{rim}`` with ``\rho_{r} = 900`` for small, medium, and large particles, we observe mostly higher fall speeds for all ``F_{rim}``  with increasing ``F_{liq}``. In all the plots we observe partially rimed ice (``D > D_{cr}``) except for a small band of graupel for very high ``F_{rim}`` in medium and large particles.

```@example
include("plots/P3TerminalVelocity_F_liq_rim.jl")
```
![](P3TerminalVelocity_F_liq_rim.svg)

When modifying process rates, we now need to consider whether they are concerned with
  the ice core or the whole particle, in addition to whether they become sources
  and sinks of different prognostic variables. With the addition of liquid fraction, too,
  come new process rates.

!!! note
    TODO - Implement process rates, complete docs.

## Acknowledgments

Click on the P3 mascot duck to be taken to the repository
  in which the authors of [MorrisonMilbrandt2015](@cite) and others
  have implemented the P3 scheme in Fortran!

[![P3 mascot](assets/p3_mascot.png)](https://github.com/P3-microphysics/P3-microphysics)
