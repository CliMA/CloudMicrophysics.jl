# Ice Nucleation

The `IceNucleation.jl` module includes
  the parameterization of activation of dust aerosol particles into ice crystals
  via deposition of water vapor.
This is a heterogeneous ice nucleation process.
The parameterization is an implementation of
  the empirical formulae from [Mohler2006](@cite)
  and is valid for two types of dust particles:
  Arizona Test Dust and desert dust from Sahara.

!!! note

    Future work includes adding parameterizations
    for other nucleation paths such as
    heterogeneous immersion freezing or homogeneous freezing
    and modeling the competition between them.


## Activated fraction for deposition freezing on dust
The parameterization models the activated fraction
  as an empirical function of ice saturation ratio,
  see eq. (3) in [Mohler2006](@cite).
```math
\begin{equation}
f_i(S_i) = exp[a(S_i - S_0)] - 1
\end{equation}
```
where:
  - ``f_i`` is the activated fraction
      (the ratio of aerosol particles acting as ice nuclei to the total number of aerosol particles),
  - ``S_i`` is the ice saturation ratio
      (the ratio of water vapor partial pressure and the water vapor partial pressure at saturation over ice),
  - ``S_0`` is the threshold ice saturation ratio,
  - ``a`` is a scaling parameter dependent on aerosol properties and temperature.

Limited experimental values for the free parameters ``S_0`` and ``a`` can be found in [Mohler2006](@cite).
Both parameters are dependent on aerosol properties and temperature.

!!! note

    For a ``f_i`` values above 0.08 or ``S_i`` between 1.35 and 1.5,
    freezing occurs in a different ice nucleation mode
    (either a second deposition or other condensation type mode).
