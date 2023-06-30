# 0-dimensional ice nucleation model

## Equations

This is a 0-dimensional adiabatic parcel model for testing nucleation schemes.
It is based on [KorolevMazin2003](@cite), as well as the cirrus box model
  [Karcher2006](@cite), [Tully2022](@cite).
The model solves for saturation in a 0-dimensional
  adiabatic parcel raising with constant velocity.

!!! note
    For now, the equations are written assuming
    that only the water vapor and cloud ice
    are present in the system.

We define ice saturation ratio ``S_i``
```math
\begin{equation}
S_i = \frac{e}{e_{si}}
\end{equation}
```
where:
- ``e`` - is the partial pressure of water vapor,
- ``e_{si}`` - is the partial pressure of water vapor at saturation over ice.

The change in saturation can be described as
```math
\begin{equation}
\frac{dS_i}{dt} = \frac{1}{e_{si}} \frac{de}{dt} - \frac{e}{e_{si}^2} \frac{de_{si}}{dt}
\end{equation}
```
From ideal gas law the partial pressure of water vapor can be written as
```math
\begin{equation}
e = q_v p \frac{R_v}{R_a}
\end{equation}
```
where:
- ``q_v`` is the water vapor specific humidity
- ``p`` is the air pressure
- ``R_v``, ``R_a`` are the gas constants for water vapor and air.

The change in partial pressure can be written as
```math
\begin{equation}
\frac{de}{dt} = \frac{dq_v}{dt} p \frac{R_v}{R_a} + q_v \frac{dp}{dt}\frac{R_v}{R_a}
\end{equation}
```

From the Clausius–Clapeyron relation
```math
\begin{equation}
\frac{de_{si}}{dt} = \frac{L_s e_{si}}{R_v T^2} \frac{dT}{dt}
\end{equation}
```
where:
- ``L_s`` is the latent heat of sublimation
- ``T`` is the temperature.

From the moist adiabatic assumption
```math
\begin{equation}
\frac{dT}{dt} = \frac{R_a T}{c_p p} \frac{dp}{dt} + \frac{L_w}{c_p} \frac{dq_l}{dt} + \frac{L_s}{c_p} \frac{d q_i}{dt}
\end{equation}
```
where:
- ``L_w`` is the latent heat of evaporation.

From hydrostatic balance and assuming constant vertical velocity:
```math
\begin{equation}
\frac{dp}{dt} = -\frac{p g}{R_a T} w
\end{equation}
```
where:
- ``g`` is the gravitational acceleration
- ``w`` is the constant vertical velocity.

Accounting for conservation of water, i.e. ``\frac{dq_v}{dt} + \frac{dq_w}{dt} + \frac{dq_i}{dt} = 0``,
and rearranging the terms
```math
\begin{equation}
\frac{dS_i}{dt} = a_1 w S_i - \left(a_2 + a_4 S_i\right) \frac{dq_w}{dt} - \left(a_2 + a_3 S_i\right) \frac{dq_i}{dt}
\end{equation}
```
where:
```math
\begin{equation}
a_1 = \frac{L_s g}{c_p T^2 R_v} - \frac{g}{R_a T}
\end{equation}
```
```math
\begin{equation}
a_2 = \frac{p}{e_{si}} \frac{R_v}{R_a}
\end{equation}
```
```math
\begin{equation}
a_3 = \frac{L_s^2}{R_v T^2 c_p}
\end{equation}
```
```math
\begin{equation}
a_4 = \frac{L_s L_w}{R_v T^2 c_p}
\end{equation}
```

The crux of the problem is modeling the ``\frac{dq_i}{dt}`` for different homogeneous and
heterogeneous ice nucleation paths. For now we only consider the water vapor deposition
on dust particles.

## Water vapour deposition on dust

For the simplest case of a spherical ice particle growing through water vapor deposition
([see discussion](https://clima.github.io/CloudMicrophysics.jl/previews/PR103/Microphysics1M/#Snow-autoconversion)),
the change of it's mass is defined as
```math
\begin{equation}
  \frac{dm}{dt} = 4 \pi \, r \, \alpha_m \, (S_i - 1) \, G(T)
  \label{eq:massratesphere}
\end{equation}
```
where:
 - ``r`` is the particle radius,
 - ``\alpha_m`` is the accommodation coefficient that takes into account the fact that not all water vapor molecules
     that arrive at the particle surface will join the growing crystal,
 - ``G(T) = \left(\frac{L_s}{KT} \left(\frac{L_s}{R_v T} - 1 \right) + \frac{R_v T}{e_{si} D} \right)^{-1}``
     combines the effects of thermal conductivity and water diffusivity,
 - ``K`` is the thermal conductivity of air,
 - ``D`` is the diffusivity of water vapor.

!!! note
    The ``r`` in eq. (\ref{eq:massratesphere})
    should be replaced by the capacitance ``C``,
    when considering non-spherical particles.
    For a sphere ``C=r``, for a circular disc ``C=2r/\pi``.

The rate of change of ice specific humidity ``q_i`` can be computed as
```math
\begin{equation}
  \frac{d \, q_i}{dt} =
  \frac{1}{\rho} \int_{0}^{\infty} \frac{dm}{dt} n(r) dr
  \label{eq:qi_rat}
\end{equation}
```

Assuming that the size distribution of nucleated and growing ice crystals is a delta function
```math
\begin{equation}
n(r) = N_{act} \delta(r-r_0)
\end{equation}
```
eq. (\ref{eq:qi_rat}) can be written as
```math
\begin{equation}
  \frac{d \, q_i}{dt} =
  \frac{1}{\rho} N_{act} \alpha_m 4 \pi \, r_0 \, (S_i - 1) \, G(T)
\end{equation}
```
where ``N_{act}`` is the number of activated ice particles.

``N_{act}`` can be computed from [activated fraction](https://clima.github.io/CloudMicrophysics.jl/previews/PR103/IceNucleation/#Activated-fraction-for-deposition-freezing-on-dust) ``f_i``
```math
\begin{equation}
  N_{act} = N_{tot} * f_i
\end{equation}
```

## Example figures

Here we show example simulation results from the adiabatic parcel model.
The model is run three times for 30 minutes simulation time,
  (shown by three different colors on the plot).
Between each run the water vapor specific humidity is changed,
  while keeping all other state variables the same as at the last time step
  of the previous run.
The prescribed vertical velocity is equal to 3.5 cm/s.

```@example
include("../../parcel/parcel.jl")
```
![](cirrus_box.svg)