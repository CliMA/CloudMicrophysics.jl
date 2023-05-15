# CloudMicrophysics.jl

The CloudMicrophysics.jl is a library of bulk cloud microphysics and aerosol schemes.

The goal of a cloud microphysics scheme is to represent the micro-scale processes
  leading to the formation of clouds and precipitation.
Bulk microphysics schemes represent overall properties of cloud and precipitation particles,
  instead of solving for the evolution of individual cloud droplets or ice crystals.
Bulk schemes typically consider different somewhat arbitrary water categories
  such as cloud water, cloud ice, rain and snow.
They then predict total mass and (optionally) number concentration of particles in each category.
A scheme that predicts the total mass of particles per category is called 1-moment,
  because the total mass of particles is proportional to the 3rd moment of the particle size distribution.
A scheme that predicts the total mass and number concentration of particles per category is called 2-moment,
  because the predicted quantities are proportional to the 3rd and 0th moment of the particle size distribution.
Aerosol particles serve as nuclei for forming cloud droplets and ice crystals.
Additional schemes are needed to predict how many cloud droplets or ice crystals
  form for a given population of aerosol particles, when using a 2-moment microphysics scheme.

So far CloudMicrophysics.jl includes:
  - 0-moment scheme that instantly removes the precipitable cloud condensate,
  - 1-moment scheme for warm rain and mixed-phase clouds (cloud water and ice, ran and snow (aggregate)),
  - 2-moment scheme for warm rain clouds (cloud water and rain),
  - collection of different 2-moment autoconversion and accretion functions,
  - experimental non-equilibrium cloud formation scheme,
  - a collection of logistic functions for smooth transitions at thresholds,
  - aerosol activation scheme,
  - ice nucleation scheme through water vapor deposition on dust aerosol.

This documentation provides some use examples for different available schemes,
  along with derivation notes and links to the literature.
The CI tests include unit tests and some very simple performance benchmarks and GPU tests.

## Authors

`CloudMicrophysics.jl` is being developed by the
  [Climate Modeling Alliance](https://clima.caltech.edu/).
![Clima logo](assets/Clima_logo.png)
