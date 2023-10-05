# Precipitation Susceptibility

The `PrecipitationSusceptibility.jl` module contains functions for determining
the precipitation susceptibility rates due to various processes for 2-moment
microphysics schemes, as described in [Glassmeier2016](@cite).

The precipitation susceptibility describes the rate of change of precipitation
production due to some moment of the particle size distribution:

For some partial moment $X$, such as $q_{liq}$, $q_{rai}$, $N_{liq}$, or $N_{rai}$
and some process rate of precipitation production $r$, such as autoconversion,
accretion, etc, this is defined as:
```math
\frac{\partial \ln r}{\partial \ln X}
```

The total precipitation susceptibility, then, is
```math
\frac{\partial \ln PP}{\partial \ln X} = \sum_{r \in R}
    \frac{r}{PP} \frac{\partial \ln r}{\partial \ln X}
```

Where $R$ is the set of all process rates that contribute to precipitation production.


Example reproducing Fig. 2 From [Glassmeier2016](@cite):

```@example
include("plots/PrecipitationSusceptibilityPlots.jl")
```
![](Glassmeier-Lohmann_Fig2.svg)
