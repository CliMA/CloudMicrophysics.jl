# Ice Nucleation

The `IceNucleation.jl` module currently contains parameterization
  of activation of aerosol particles into ice crystals via heterogeneous freezing.
Future work includes addition of parameterization for differential activated fraction, homogeneous freezing, and activated fraction for competing aerosol size distributions.
The parameterization computes the activated fraction and ice saturation ratio.

The modules are an implementation of the parameterization
  from [Murray et al. (2012)](@cite) and
  [Vali (1971)](@cite).


## Assumed ice nucleating aerosol composition and dependencies

There are two major descriptions for ice nucleation data: *stochastic* and *singular*. The stochastic model accounts for the time dependency of the probability of nucleation. The singular model assumes that the dependence on the distribution of ice nuclei types is of much greater relative importance than that of time in nucleation of ice. The modified singular model by Vali in which time is a contributing factor through cooling rate is also included.

Within stochastic modeling, there is the *single component stochastic (SCS)* and the *multiple component stochastic (MSC)*. SCS assumes each droplet contains the same ice nucleating materials and thus same probability of freezing. In MSC, droplets may contain varying amounts and types of particles and, therefore, have differing freezing probabilities. Similarly, the singular models hold for droplets with varying composition of aerosol types.

Particles are assumed to be internally mixed.
The following table lists the parameters defining the droplet
  physical and chemical properties.

|    variable name     |         definition                                | units                |
|----------------------|---------------------------------------------------|----------------------|
|``J``                 | nucleation rate coefficient                       | ``cm^{-2} \, s^{-1}``|
|``s``                 | nucleant surface area                             | ``cm^{2}``           |
|``N      ``           | number of liquid droplets                         | ``-``                |
|``n_{ice}``           | number of frozen droplets                         | ``-``                |
|``\alpha``            | temperature offset                                | ``°C``               |
|``r``                 | cooling rate                               | ``-``               |
|``\beta``             | beta parameter                                    | ``-``                |
|``T_c``               | characteristic temperature                        | ``°C``               |
|``n_s``               | active site density                        | ``kg \, m^{-3}``    |
|``k``                 | density of nucleation sites                | ``kg \, m^{-3}``    |


## Single Component Stochastic (SCS)
The SCS model depends on the heterogeneous nucleation rate coefficient, ``J_{het}``. Depending on the mode of nucleation, ``J_{het}`` will be dependent on different parameters. For deposition mode, ``J_{het}`` is a funciton of relative humidity and temperature. For immersion mode, ``J_{het}`` is a function of temperature only. The fraction of droplets that freeze within a time duration is defined in eq. (16) of Murray (2012):

```math
\frac{{\Delta}n_{ice}}{N_1} = 1 - exp(-J_{het}s{\Delta}t)
```
$$\frac{{\Delta}n_{ice}}{N_1} = 1 - exp(-J_{het}s{\Delta}t)$$

where:
  - ``{\Delta}t = t_2 - t_1`` is the time duration of interest,
  - ``{\Delta}n_{ice}`` is the newly frozen droplets within ``{\Delta}t``,
  - ``N_1`` is the number of liquid droplets at ``time t_1``,
  - ``J_{het}`` is the heterogeneous nucleation rate coefficient,
  - ``s`` is the particle's surface area.

## Multiple Component Stochastic (MSC)
The MSC model accounts for a mixture of particle types, each with their own ice nucleating potential. This is accounted for in the temperature dependent nucelation rate coefficient ``J_i`` for nuclei of type ``i``.

```math
\frac{{\Delta}n_{ice}}{N} = 1 - exp(-\sum_{i} J_{i}s_{i}{\Delta}t)
```
$$\frac{{\Delta}n_{ice}}{N} = 1 - exp(-\sum_{i} J_{i}s_{i}{\Delta}t)$$

where:
  - ``N`` is number of liquid droplets,
  - ``{\Delta}n_{ice}`` is the newly frozen droplets within ``{\Delta}t``,
  - ``J_i`` is the nucleation rate coefficient for nucleus of type ``i``,
  - ``s_i`` is the total surface area of nucleus type ``i``.

## Singular Model
Singular models revolve around the distribution of particle types and their differing freezing abilities rather than time dependency.
### Surface Area Based
```math
f_{ice}(T) = \frac{{\Delta}n_{ice}(T)}{N_{tot}} = 1 - exp(-n_{s}(T)s)
```
$$f_{ice}(T) = \frac{{\Delta}n_{ice}(T)}{N_{tot}} = 1 - exp(-n_{s}(T)s)$$

where:
  - ``f_{ice}(T)`` is the temperature-dependent ice activated fraction,
  - ``n_{ice}(T)`` is the cumulative number of frozen droplets at temperature ``T``,
  - ``n_{s}(T)`` is the active site density (cumulative number of nucleation sites per surface area) at temperature ``T``.
 
```math
n_s(T) = - \int_{T_0}^T k(T){d}T
```
$$n_s(T) = - \int_{T_0}^T k(T){d}T$$

### Volume Based
```math
f_{ice}(T) = \frac{{\Delta}n_{ice}(T)}{N_{tot}} = 1 - exp(-K(T)V)
```
$$f_{ice}(T) = \frac{{\Delta}n_{ice}(T)}{N_{tot}} = 1 - exp(-K(T)V)$$

where:
  - ``K(T)`` is the density of nucleation sites per unit volume,
  - ``V`` is nucleus volume.

## Modified Singular Model
The modified singular model is similar to the singular model except the active site density, ``n_{s}``, now accounts for the offset in freezing temperature due to the stochastic characteristic of freezing.

```math
f_{ice}(T) = \frac{{\Delta}n_{ice}(T)}{N} = 1 - exp(-n_{s}(T-\alpha)s)
```
$$f_{ice}(T) = \frac{{\Delta}n_{ice}(T)}{N} = 1 - exp(-n_{s}(T-\alpha)s)$$

where:
  - ``n_{s}(T)`` is the active site density (cumulative number of nucleation sites per surface area) at temperature ``T``,
  - ``\alpha`` is the temperature offset.

```math
\alpha = \beta \, log(|r|)
``` 
$$\alpha = \beta \, log(|r|)$$

where:
  - ``r`` is cooling rate,
  - ``\beta`` is an empirical parameter.

