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
