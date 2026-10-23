# Microphysics 0M

The `Microphysics0M.jl` module defines a 0-moment bulk parameterization of
  the moisture sink due to precipitation.
It offers a simplified way of removing the excess water
  without assuming anything about the size distributions of cloud
  or precipitation particles.

The ``q_{tot}`` (total water specific content) sink due to precipitation
  is obtained by relaxation with a constant timescale
  to a state with condensate exceeding a threshold value removed.
The threshold for removing excess ``q_{tot}`` is defined either by the
  condensate specific content or supersaturation.
The thresholds and the relaxation timescale are defined in
  `ClimaParams.jl`.

!!! note

    To remove precipitation instantly, the relaxation timescale should be
    equal to the timestep length.

## Moisture sink due to precipitation

If based on maximum condensate specific content, the sink is defined as:
``` math
\begin{equation}
  \left. \frac{d \, q_{tot}}{dt} \right|_{precip} =-
    \frac{\max(0, q_{lcl} + q_{icl} - q_{c0})}{\tau_{precip}}
\end{equation}
```
where:
  - ``q_{lcl}``, ``q_{icl}`` are cloud liquid water and cloud ice specific contents,
  - ``q_{c0}`` is the condensate specific content threshold above which water is removed,
  - ``\tau_{precip}`` is the relaxation timescale.

If based on saturation excess, the sink is defined as:
```math
\begin{equation}
  \left. \frac{d \, q_{tot}}{dt} \right|_{precip} =-
    \frac{\max(0, q_{lcl} + q_{icl} - S_{0} \, q_{vap}^{sat})}{\tau_{precip}}
\end{equation}
```
where:
  - ``q_{lcl}``, ``q_{icl}`` are cloud liquid water and cloud ice specific contents,
  - ``S_{0}`` is the supersaturation threshold above which water is removed,
  - ``q_{vap}^{sat}`` is the saturation specific humidity,
  - ``\tau_{precip}`` is the relaxation timescale.

## Moisture sink derivative

We provide the derivative of the precipitation removal rate with respect to the total water
  specific humidity ``q_{tot}``, to allow implicit time integration.
The derivative is equal to ``-1/\tau_{precip}`` if the cloud condensate
  (liquid water plus ice) exceeds the threshold, and 0 otherwise.


## Updraft precipitation parameters

`Microphysics0MParams` also stores a second set of `Parameters0M` in
  `precip_updraft`, for condensate in convective updrafts
  (e.g., EDMF updrafts in ClimaAtmos).
It uses the same thresholds as `precip` and the relaxation timescale
  `precipitation_timescale_updraft_ratio` ``\times`` `precipitation_timescale`
  from `ClimaParams.jl` (the default ratio of 1 gives the updrafts the same
  timescale as `precip`).

## Evaporation of precipitation

The precipitation removed by the 0-moment scheme can evaporate
  (or sublimate) as it falls through sub-saturated air.
Following [Kessler1969](@cite), as used in the IFS
  ([Tiedtke1993](@cite), [IFS2023](@cite)), the evaporation rate is
```math
\begin{equation}
  E = a_p \, k_E \, \max(0, RH_c \, q_{vap}^{sat} - q_{vap})
    \left[ \frac{\sqrt{p/p_s} \, F / a_p}{\alpha_2} \right]^{\alpha_3}
\end{equation}
```
where:
  - ``F`` is the grid-mean precipitation flux (positive downward),
  - ``a_p`` is the fraction of the grid box covered by precipitation,
  - ``k_E`` is the evaporation coefficient (``k_E = 0`` disables evaporation),
  - ``RH_c`` is the critical relative humidity above which precipitation does not evaporate,
  - ``q_{vap}^{sat}`` is the saturation specific humidity
    (over liquid for rain and over ice for snow),
  - ``p/p_s`` is the ratio of the air pressure to the surface pressure,
  - ``\alpha_2``, ``\alpha_3`` are the flux scale and exponent.
The parameters are stored in `Evaporation0M` and defined in `ClimaParams.jl`.
The default ``k_E = 0`` turns evaporation off; the IFS uses
  ``k_E = 5.44 \times 10^{-4} \, s^{-1}``.

The rate is not limited.
Time-stepping schemes should limit it so that no more precipitation
  evaporates than falls into a grid cell, and so that the air is not
  moistened beyond ``RH_c``.
