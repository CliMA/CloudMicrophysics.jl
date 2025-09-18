# Thermodynamics Interface

The `CloudMicrophysics.jl` library is designed to be flexible and agnostic
  to the specific thermodynamics of the model it is coupled with.
It requires certain functions and free parameters to be provided by the host model,
  such as, for example, the formula to compute the saturation vapor pressure.

## Assumptions about water categories

The different parameterizations available within the `CloudMicrophysics.jl` library
  do make different assumptions about the partitioning of cloud condensate.
It is crucial to understand those assumptions when choosing a parameterization
  best suited for your model, and coupling it in a thermodynamics-consistent way.
  - 1-moment scheme splits condensed water into
    cloud liquid `q_lcl`, cloud ice `q_icl`, rain `q_rai` and snow `q_sno`.
  - 2-moment scheme only considers liquid phase and divides the condensed water into
    cloud liquid `q_lcl` and rain `q_rai`.
  - P3 scheme only considers ice phase and has one category for all condensed species `q_ice`

## Additional functions and parameters

Needed by the `CloudMIcrophysics.jl` library:
  - Gas constants for dry air and water vapor
  - Gas constant for moist air
  - Latent heats of vaporization, fusion and sublimation
  - Specific heat capacity under constant pressure for moist air
  - Saturation vapor pressure and supersaturation over liquid and ice
  - Water vapor specific content
  - Specific content from partial pressure and partial pressure from specific content

Needed in tests and examples shown in the documentation
  - Gravitational constant (only used in the adiabatic parcel example)
  - Air density (only used in tests)
