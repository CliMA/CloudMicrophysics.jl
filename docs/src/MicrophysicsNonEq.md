# Microphysics NonEquilibrium

The `MicrophysicsNonEq.jl` module describes a bulk parameterization of
  diffusion of water vapor on cloud droplets and cloud ice crystals
  modeled as a relaxation to equilibrium.

The cloud microphysics variables are expressed as specific humidities:
  - `q_tot` - total water specific humidity,
  - `q_vap` - water vapor specific humidity,
  - `q_liq` - cloud water specific humidity,
  - `q_ice` - cloud ice specific humidity,

Parameters used in the parameterization are defined in
  [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl) package.
They consist of:

|    symbol  |         definition                             | units | default value |
|------------|------------------------------------------------|-------|---------------|
|``\tau_{l}``| cloud water condensation/evaporation timescale | ``s`` | ``10``        |
|``\tau_{i}``| cloud ice deposition/sublimation timescale     | ``s`` | ``10``        |

## Simple condensation/evaporation and deposition/sublimation

Condensation/evaporation of cloud liquid water and
deposition/sublimation of cloud ice are parameterized
  as a relaxation to equilibrium value at the current time step.
The equilibrium value is obtained based on a prescribed phase partition function
  that divides the available excess water vapor between liquid and ice
  (based on temperature).
```math
\begin{equation}
  \left. \frac{d \, q_{liq}}{dt} \right|_{cond, evap} = \frac{q^{eq}_{liq} - q_{liq}}{\tau_{l}}; \;\;\;\;\;\;\;
  \left. \frac{d \, q_{ice}}{dt} \right|_{dep, sub}   = \frac{q^{eq}_{ice} - q_{ice}}{\tau_{i}}
\end{equation}
```
where:
 - ``q^{eq}_{liq}, q^{eq}_{ice}`` - liquid and water specific humidity in equilibrium at current temperature and
   assuming some phase partition function based on temperature
 - ``q_{liq}, q_{ice}`` - current liquid water and ice specific humidity,
 - ``\tau_{l}, \tau_{i}`` - relaxation timescales.

!!! note
    Both ``\tau_{l}`` and ``\tau_{i}`` are assumed to be constant.
    It would be great to make the relaxation time a function of
    available condensation nuclei, turbulence intensity, etc.
    See works by [prof Raymond Shaw](https://www.mtu.edu/physics/department/faculty/shaw/)
    for hints.
    In particular, [Desai2019](@cite).

## Condensation/evaporation and deposition/sublimation from Morrison and Milbrandt 2015

Condensation/evaporation and deposition/sublimation rates are based on
  the difference between the water vapor specific humidity and saturation
  vapor specific humidity over liquid and ice at the current temperature.
The process is modeled as a relaxation with a constant timescale.
This formulation is derived from [MorrisonGrabowski2008_supersat](@cite)
  and [MorrisonMilbrandt2015](@cite), but without imposing exponential time integrators.

!!! note
    The [MorrisonGrabowski2008_supersat](@cite) and [MorrisonMilbrandt2015](@cite)
    papers use mass mixing ratios, not specific humidities.
    Additionally, in their formulations they consider two different categories for liquid:
    cloud water and rain. For now we only consider cloud water and use a single relaxation timescale
    ``\tau_l`` (liquid) rather than separate ``\tau_c`` (cloud) and ``\tau_r`` (rain) values.

```math
\begin{equation}
   \left. \frac{d \, q_{liq}}{dt} \right|_{cond, evap} = \frac{q_{vap} - q_{sl}}{\tau_l \Gamma_l}; \;\;\;\;\;\;\;
   \left. \frac{d \, q_{ice}}{dt} \right|_{dep, sub}   = \frac{q_{vap} - q_{si}}{\tau_i \Gamma_i}
\end{equation}
```
where:
- ``q_{vap}`` is the water vapor specific humidity
- ``q_{sl}``, ``q_{si}`` is the saturation specific humidity over liquid and ice
- ``\tau_l``, ``\tau_i`` is the liquid and ice relaxation timescale
- ``\Gamma_l``, ``\Gamma_i`` is a psychometric correction due to latent heating/cooling:

```math
\begin{equation}
    \Gamma_l = 1 + \frac{L_{v}}{c_p} \frac{dq_{sl}}{dT}; \;\;\;\;\;\;\;\;
    \Gamma_i = 1 + \frac{L_{s}}{c_p} \frac{dq_{si}}{dT}
\end{equation}
```
```math
\begin{equation}
    \frac{dq_{sl}}{dT} = q_{sl} \left(\frac{L_v}{R_v  T^2} - \frac{1}{T} \right); \;\;\;\;\;\;\;\;\;\;
    \frac{dq_{si}}{dT} = q_{si} \left(\frac{L_s}{R_v  T^2} - \frac{1}{T} \right)
\end{equation}
```
where:
- ``T`` is the temperature,
- ``c_p`` is the specific heat of air at constant pressure,
- ``R_v`` is the gas constant of water vapor,
- ``L_v`` and ``L_s`` is the latent heat of vaporization and sublimation.
 
    Note that these forms of condensation/sublimation and deposition/sublimation are equivalent to those described in the adiabatic parcel model with some rearrangements and assumptions. It is just necessary to use the definitions of ```\tau```, ```q_{sat}```, and the thermal diffusivity ```D_v```:

    ```math
    \begin{equation}
      \tau = 4 \pi N_{tot} \bar{r}, \;\;\;\;\;\;\;\;
      q_{sat} = \frac{e_{sat}}{\rho R_v T}, \;\;\;\;\;\;\;\;
      D_v = \frac{K}{\rho c_p}
    \end{equation}
    ```
    and if we assume that the supersaturation S can be approximated by specific humidities (this is only exactly true for mass mixing ratios):
    ```math
    \begin{equation}
        S = \frac{q_{vap}}{q_{sat}}
    \end{equation}
    ```
    then we can write:
    ```math
    \begin{equation}
      q_{vap} - q_{sat} = q_{sat}(S - 1)
    \end{equation}
    ```
    and ```Gamma``` is equivalent to the ```G``` function used in our parcel and parameterizations.


### Time integrated version

Finally, we can also use the time integrated formulation of condensation/evaporation and deposition/sublimation described by [MorrisonMilbrandt2015](@cite). We assume that the condensation/evaporation and deposition/sublimation over the course of the timestep can be calculated using the average difference between the $q_v$ value and those of saturation. Ie, they write:

```math
\begin{equation}
   \left( \frac{d q_l}{dt} \right) _{cond} = \frac{\bar{\delta_l}}{\tau_l \Gamma_l}
\end{equation}
```
where $\delta$ is defined as
```math
\begin{equation}
    \delta = q_v - q_{sl}
\end{equation}
```

The same can be done for deposition/sublimation, with an addition of the correction for the difference between the specific humidities of supersaturation for liquid and for ice:

```math
\begin{equation}
   \left( \frac{d q_i}{dt} \right) _{cond} = \frac{\bar{\delta_l}}{\tau_i \Gamma_i} + \frac{q_{sl} - q_{si}}{\tau_i \Gamma_i}
\end{equation}
```

In order to calculate $\bar{\delta}$, we consider the derivative of delta:

```math
\begin{equation}
    \frac{d \delta}{dt} = \frac{dq_v}{dt} - \frac{dq_{sl}}{dt}
\end{equation}
```

We can write the derivation for the saturated specific humidity:

```math
\begin{equation}
    \frac{d q_{s}}{dt} = \frac{dq_{s}}{dT} \frac{dT}{dt} + \frac{dq_{s}}{dp} \frac{dp}{dt} = \frac{dq_{s}}{dT} \frac{dT}{dt} + \frac{q_{s} \rho_a g w}{p - e}
\end{equation}
```
where
- $\rho_a$ is the air density
- $p$ is the air pressure
- $e$ is the saturation vapor pressure
- This makes sense given the assumptions that $\frac{dp}{dt} = -\rho_a gw$ from hydrostatic balance, and $q_s = \frac{e_s}{p-e_s}$  for water vapor, where $e_s$ is the saturation water vapor pressure. Then $\frac{d q_s}{dp} = - \frac{e_s}{(p-e_s)^2} = - \frac{q_s}{p-e_s}$ .

and they write the change in temperature with time as:

```math
\begin{equation}
    \frac{dT}{dt} = - \frac{g w}{c_p} + \frac{L_v}{c_p} \frac{\delta}{\tau_l \Gamma_l} + \frac{L_s}{c_p} \left(\frac{\delta}{\tau_i \Gamma_i} + \frac{q_{sl} - q_{si}}{\tau_i \Gamma_i} \right)
\end{equation}
```

Then

```math
\begin{equation}
       \frac{d q_{sl}}{dt} = \frac{dq_{sl}}{dT} \left[- \frac{g w}{c_p} + \frac{L_v}{c_p} \frac{\delta}{\tau_l \Gamma_l} + \frac{L_s}{c_p} \left(\frac{\delta}{\tau_i \Gamma_i} + \frac{q_{sl} - q_{si}}{\tau_i \Gamma_i} \right) \right] + \frac{q_{sl} \rho_a g w}{p - e}
\end{equation}
```

!!! note
    we neglect terms due to radiation and mixing (included in the [morrison_modeling_2008](@cite) and [MorrisonMilbrandt2015](@cite) versions) in the parcel case.

Putting these back into the $\delta$ equations and rearranging based on the definitions of $\Gamma_l$ and $\Gamma_i$, we get:

```math
\begin{equation}
    \frac{d \delta}{dt} = \frac{- q_{sl} \rho_a g w}{p - e_{sl}} + \frac{gw}{c_p} \frac{dq_{sl}}{dT} - \frac{\delta}{\tau}
\end{equation}
```

where $\tau$ is

```math
\begin{equation}
    \tau = \left( \tau_{l}^{-1} + \left( 1 + \frac{L_{s}}{c_p} \frac{dq_{sl}}{dT} \right) \frac{\tau_i^{-1}}{\Gamma_i} \right)^{-1}
\end{equation}
```

where

- $\tau_l$ is the relaxation timescale for liquid droplets
- $\tau_i$ is the relaxation timescale for ice droplets
- $L_s$ is the latent heat of sublimation
- $c_p$ is the specific heat of air at constant pressure
- $T$ is temperature

If we assume that changes in $\frac{dq_s}{dT}$ and $\tau$ are small in comparison to their magnitudes, then we can approximate them as constant and both $\delta$ equations are linear differential equations with a solution of the form:

```math
\begin{equation}
    \delta (t) = A_c \tau + (\delta_{t=0} - A_c \tau) e^{-t/\tau}
\end{equation}
```

Then the $A_c$ values can be calculated based on the previous derivatives.

Condensation Ac:
```math
\begin{equation}
        A_c = - \frac{q_{sl} \rho g w}{p - e_s} + \frac{dq_{sl}}{dT} \frac{wg}{c_p} - \frac{(q_{sl} - q_{si})}{\tau_i \Gamma_i} \left( 1 + \frac{L_s}{c_p} \frac{dq_{sl}}{dT} \right)
\end{equation}
```

Finally, to get the values for condensation/evaporation and deposition/sublimation over the timestep, we calculate the time average of $\delta$ by dividing by $\Delta t$ (the timestep) and then by the relaxation time as described in the above equations.

```math
\begin{equation}
    \left( \frac{d q_l}{dt} \right)_{cond} = \frac{A_c \tau}{\tau_i \Gamma_l} + (\delta_{t=0} - A_c \tau) \frac{\tau}{\Delta t \tau_i \Gamma_l} (1-e^{-\Delta t / \tau} )
\end{equation}
```

```math
\begin{equation}
    \left( \frac{d q_{si}}{dt} \right)_{dep} = A_c \frac{\tau}{\tau_i \Gamma_i} + (\delta_{t=0} - A_c \tau) \frac{\tau}{\Delta t \tau_i \Gamma_i} (1-e^{-\Delta t / \tau} ) + \frac{(q_{sl} - q_{si})}{\tau_i \Gamma_i}
\end{equation}
```

!!! note
    Need to add limiters