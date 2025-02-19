# Diagnostics

`CloudMicrophysics.jl` offers a couple of options to compute cloud and precipitation
radiative properties based on different available parameterizations and their
underlying assumptions about the size distribution and properties of particles.

Available diagnostics are:
 - Radar reflectivity
 - Effective radius

## Radar reflectivity

The radar reflectivity factor ``Z`` is used to measure the power returned
  by a radar signal when it encounters particles (cloud, rain droplets, etc),
  and is defined as the sixth moment of the particles distributions ``n(r)``:

```math
\begin{equation}
Z = {\int_0^\infty r^{6} \, n(r) \, dr}.
\label{eq:Z}
\end{equation}
```

``Z`` is typically normalized by radar reflectivity factor ``Z_0``
of a drop of radius ``1 mm`` in a volume of ``1 m^3``, and is reported
as a decimal logarithm to obtain the normalized
logarithmic rain radar reflectivity ``L_Z``.
```math
\begin{equation}
L_Z = {10 \, \log_{10} \left( \frac{Z}{Z_0} \right)}.
\end{equation}
```
The resulting logarithmic dimensionless unit is decibel relative to ``Z_0``, or ``dBZ``.

### 1-moment microphysics

For the [1-moment scheme](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/)
  we only consider the rain drop size distribution.
Integrating over the assumed Marshall-Palmer distribution leads to
```math
\begin{equation}
Z = {\frac{6! \, n_{0}^{rai}}{\lambda^7}},
\end{equation}
```
where:
 - ``n_{0}^{rai}`` and ``\lambda`` - rain drop size distribution parameters.

### 2-moment microphysics

For the [2-moment scheme](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics2M/)
  we take into consideration the effect of both cloud and rain droplets.
Integrating over the assumed cloud droplets Gamma distribution leads to
```math
\begin{equation}
Z_c = A_c C^{\nu_c+1} \frac{ (B_c C^{\mu_c})^{-\frac{3+\nu_c}{\mu_c}} \, \Gamma \left(\frac{3+\nu_c}{\mu_c}\right)}{\mu_c},
\end{equation}
```
where:
 - ``\Gamma \,(x) = \int_{0}^{\infty} \! t^{x - 1} e^{-t} \mathrm{d}t`` is the gamma function,
 - ``C = \frac{4}{3} \pi \rho_w``.

Similar for rain drop exponential distribution
```math
\begin{equation}
Z_r = A_r C^{\nu_r+1} \frac{ (B_r C^{\mu_r})^{-\frac{3+\nu_r}{\mu_r}} \, \Gamma \left(\frac{3+\nu_r}{\mu_r}\right)}{\mu_r},
\end{equation}
```
The final radar reflectivity factor is a sum of ``Z_c`` and ``Z_r``.

## Effective radius

The effective radius of hydrometeors (``r_{eff}``) is defined as
the area weighted radius of the population of particles.
It can be computed as the ratio of the third to the second moment
of the size distribution.

### 2-moment microphysics

We compute the total third and second moment as a sum of cloud condensate and
precipitation moments:
```math
\begin{equation}
r_{eff} = \frac{M_{3}^c + M_{3}^r}{M_{2}^c + M_{2}^r} = \frac{{\int_0^\infty r^{3} \, (n_c(r) + n_r(r)) \, dr}}{{\int_0^\infty r^{2} \, (n_c(r) + n_r(r)) \, dr}}.
\label{eq:reff}
\end{equation}
```
After integrating we obtain
```math
\begin{equation}
M_{3}^c + M_{3}^r = A_c C^{\nu_c+1} \frac{ (B_c C^{\mu_c})^{-\frac{2+\nu_c}{\mu_c}} \, \Gamma \left(\frac{2+\nu_c}{\mu_c}\right)}{\mu_c} + A_r C^{\nu_r+1} \frac{ (B_r C^{\mu_r})^{-\frac{2+\nu_r}{\mu_r}} \, \Gamma \left(\frac{2+\nu_r}{\mu_r}\right)}{\mu_r}.
\end{equation}
```
```math
\begin{equation}
M_{2}^c + M_{2}^r = A_c C^{\nu_c+1} \frac{ (B_c C^{\mu_c})^{-\frac{5+3\nu_c}{3\mu_c}} \, \Gamma \left(\frac{5+3\nu_c}{3\mu_c}\right)}{\mu_c} + A_r C^{\nu_r+1} \frac{ (B_r C^{\mu_r})^{-\frac{5+3\nu_r}{3\mu_r}} \, \Gamma \left(\frac{5+3\nu_r}{3\mu_r}\right)}{\mu_r}.
\end{equation}
```

### Liu and Hallett 1997

For 1-moment microphysics scheme the effective radius
is parameterized following [Liu1997](@cite):

```math
\begin{equation}
  r_{eff} = \frac{r_{vol}}{k^{\frac{1}{3}}},
\end{equation}
```
where:
  - ``r_{vol}`` represents the volume-averaged radius,
  - ``k = 0.8``.

Where the volume-averaged radius is computed using
```math
\begin{equation}
  r_{vol} = \left(\frac{3}{4 \pi \rho_w}\right)^{\frac{1}{3}} \, \left(\frac{L}{N}\right)^{\frac{1}{3}} = \left(\frac{3  \rho (q_{liq} + q_{rai})}{4 \pi \rho_w (N_{liq} + N_{rai})}\right)^{\frac{1}{3}},
\end{equation}
```
where:
  - ``L = \rho q``, is the liquid water content,
  - ``N = N_{liq} + N_{rai}``.

By default for the 1-moment scheme we don't consider precipitation and assume
a constant cloud droplet number concentration of 100 $cm^{-3}$.

### Constant

For testing purposes we also provide a constant effective radius option.
The default values are 14 $\mu m$ for liquid clouds and 25 $\mu m$ for ice clouds,
and can be easily overwritten via `ClimaParams.jl`.

## Example figures

Below we show effective radius and radar reflecivity plots as a function of cloud water and rain water.
The effective radius is computed assuming a constant cloud droplet number concentration
  of 100 or 1000 per cubic centimeter.
The radar reflectivity is computed assuming a constant rain drop number concentration
  of 10 or 100 per cubic centimeter.
Note the effect of using the limiters in SB2006 scheme on the radar reflectivity.

```@example
include("plots/CloudDiagnostics.jl")
```
![](CloudDiagnostics.svg)
