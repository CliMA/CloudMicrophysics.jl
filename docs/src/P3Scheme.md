# P3 Scheme

The `P3Scheme.jl` module implements the predicted particle properties
 (P3) scheme for ice-phase microphysics developed by [MorrisonMilbrandt2015](@cite)
The P3 scheme is a 2-moment, bulk scheme involving a
 single ice-phase category with 4 degrees of freedom: total mass,
 rime mass, rime volume, and number mixing ratios.
Traditionally, cloud ice microphysics schemes use various predefined
 categories (such as ice, graupel, or hail) to represent ice modes, but the P3 scheme sidesteps the
 problem of prescribing transitions between ice categories by adopting
 a single ice category and evolving its properties. This simplification
 aids in attempts to constrain the scheme's free parameters.

The prognostic variables are:
 - ``N_{ice}`` - number concentration 1/m3
 - ``q_{ice}`` - ice mass density kg/m3
 - ``q_{rim}`` - rime mass density kg/m3
 - ``B_{rim}`` - rime volume - (volume of rime per total air volume: dimensionless)

!!! note
    TODO - At some point we should switch to specific humidities...

## Assumed particle size distribution (PSD)

Following [MorrisonMilbrandt2015](@cite), the scheme assumes a
 gamma distribution for the concentration of particles per unit volume
 based on particle size measurements obtained by [Heymsfield2003](@cite)
 in tropical and midlatitude ice clouds and implemented by
 [MorrisonGrabowski2008](@cite):

```math
N'(D) = N_{0} D^\mu \, e^{-\lambda \, D}
```
where:
 - ``N'`` is the number concentration in ``m^{-4}``
 - ``D`` is the maximum particle dimension in ``m``,
 - ``N_0`` is the intercept parameter in ``m^{-4}``,
 - ``\mu`` is the shape parameter (dimensionless),
 - ``\lambda`` is the slope parameter in ``m^{-1}``.

We assume ``\mu \ = 0.00191 \; \lambda \ ^{0.8} - 2``.
Following [MorrisonGrabowski2008](@cite) we limit ``\mu \ \in (0,6)``.
A negative ``\mu`` can occur only for very small mean particle sizes``\frac{1}{\lambda} < ~0.17 mm``.
``N_0`` and ``\lambda`` can be found using different moments of the PSD,
namely the total number concentration ``N`` and mass mixing ratio ``q``, where

```math
N = \int_{0}^{\infty} \! N'(D) \mathrm{d}D
```

```math
q = \int_{0}^{\infty} \! m(D) N'(D) \mathrm{d}D
```

For liquid droplets, these equations are solved without issue, but for ice, the third moment of the size distribution listed above (i.e. ``\int_{0}^{\infty} \! m(D) N'(D) \mathrm{d}D``) varies as the mass relation varies across the PSD (see below for the mass regime documentation). On the other hand, the first moment of the PSD, the number concentration, does not vary across the PSD and yields ``N = \frac{N_0}{\lambda}``.

!!! note
    TODO - The scheme uses a mean particle size value ``D_m`` for each time step to determine which mass relation to employ. In other words, ``N_0`` and ``\lambda`` must be calculated for the five different mass relations below to accommodate ranges of ``D_m`` corresponding to each mass relation. For mean particle sizes that employ the mass relations characterized by graupel and by partially rimed ice, the mass relations are time-dependent due to the presence of ``\rho_g`` and ``F_r``. This complicates the scheme's use of the PSD, and as a result, writing analytical formulas for the PSD parameters is challenging.

## Assumed particle mass relationships

The mass ``m`` of particles as a function of maximum particle dimension ``D``
 is a piecewise function with variable thresholds described
 by the following table.

| particle properties                  |      condition(s)                            |    m(D) relation                             |
|:-------------------------------------|:---------------------------------------------|:---------------------------------------------|
|small, spherical ice                  | ``D < D_{th}``                               | ``\frac{\pi}{6} \rho_i \ D^3``               |
|large, unrimed ice                    | ``q_{rim} = 0`` and ``D > D_{th}``           | ``\alpha_{va} \ D^{\beta_{va}}``             |
|dense nonspherical ice                | ``q_{rim} > 0`` and ``D_{th} < D < D_{gr}``  | ``\alpha_{va} \ D^{\beta_{va}}``             |
|graupel (completely rimed, spherical) | ``q_{rim} > 0`` and ``D_{gr} < D < D_{cr}``  | ``\frac{\pi}{6} \rho_g \ D^3``               |
|partially rimed ice                   | ``q_{rim} > 0`` and ``D > D_{cr}``           | ``\frac{\alpha_{va}}{1-F_r} D^{\beta_{va}}`` |

where:
 - ``D_{th}``, ``D_{gr}``, ``D_{cr}`` are particle size thresholds in ``m``,
 - ``\rho_i`` is cloud ice density in ``kg m^{-3}``,
 - ``\beta_{va} = 1.9`` is a dimensionless parameter from [BrownFrancis1995](@cite) (based on measurements of vapor diffusion and aggregation in midlatitude cirrus),
 - ``\alpha_{va} = 7.38 \; 10^{-11} \; 10^{6 \beta_{va} - 3}`` in ``kg m^{-β_{va}}`` is a parameter from [BrownFrancis1995](@cite) in base SI units,
 - ``\rho_g`` is the bulk density of graupel in ``kg m^{-3}``.

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
The system is solved using [`NonlinearSolve.jl`](https://docs.sciml.ai/NonlinearSolve/stable/).

!!! note
    TODO - The use of NonlinearSolve.jl is not ideal because of its runtime and memory allocation requirements. Currently, there is also a look-up table NetCDF file which could be used to look up values of the quantities which form the nonlinear system. However, the look-up table is not GPU-compatible and would require too much memory in an Earth System Model. The current approach may be of use for testing and for visualization of the system, but other options, such as using RootSolvers.jl or using a simpler fit that approximates the solver output, are more suitable long term solutions which do not require outside packages which employ auto-differentiation or use memory, both of which do not suit the needs of CliMA.

Below we show the m(D) regime, replicating Figures 1 (a) and (b) from [MorrisonMilbrandt2015](@cite).

```@example
include("P3SchemePlots.jl")
p3_m_plot1(["cyan2", "cyan4", "midnightblue"], ["hotpink", "magenta3", "purple4"])
p3_m_plot2(["cyan2", "cyan4", "midnightblue"], ["hotpink", "magenta3", "purple4"])
```

![](MorrisonandMilbrandtFig1a.svg)

![](MorrisonandMilbrandtFig1b.svg)

### Threshold solver and look-up table

The look-up table environment involves three functions: one to generate the table using [`NonlinearSolve.jl`](https://docs.sciml.ai/NonlinearSolve/stable/), another to load the table as Julia Matrix objects, and a last function that can be called at each time step to generate look-up quantities. The table itself is a NetCDF file with dimensions of predicted rime density and rime mass fraction, which form domains over which look-up values for four variables (``D_{cr}``, ``D_{gr}``, ``\rho_{g}``, ``\rho_{d}``) are stored. Rather than call upon [`NCDatasets.jl`](https://github.com/Alexander-Barth/NCDatasets.jl) to open and close a NetCDF file for each time step, passing a 4-vector containing matrices corresponding to each quantity optimizes performance. The look-up function indexes into the aforementioned matrices to generate values corresponding to specific predicted rime density and rime mass fraction inputs, and in the event that the inputs do not align perfectly with axis points, a linear interpolation is used between the values closest to the inputs to approximate the result given by the solver.

However, the look-up table is not GPU-compatible and would require too much memory in an Earth System Model. The current approach may be of use for testing and for visualization of the system, but other options such as using [`RootSolvers.jl`](https://github.com/CliMA/RootSolvers.jl) or using a simpler fit that approximates the solver output are more suitable long-term solutions which do not require outside packages which employ auto-differentiation or use memory, both of which do not suit the needs of CliMA.

Although the table will ultimately be futile, we have provided the following heatmaps for testing and visualization of the system.

```@example
using CloudMicrophysics
CloudMicrophysics.P3Scheme.generate_threshold_table()
include("P3SchemePlots.jl")
p3_heatmap("D_cr", "../build/p3_lookup.nc")
p3_heatmap("D_gr", "../build/p3_lookup.nc")
p3_heatmap("ρ_g", "../build/p3_lookup.nc")
p3_heatmap("ρ_d", "../build/p3_lookup.nc")
```

![](p3_D_cr_lookup.svg)

![](p3_D_gr_lookup.svg)

![](p3_ρ_g_lookup.svg)

![](p3_ρ_d_lookup.svg)
