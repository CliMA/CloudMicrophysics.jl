# Microphysics NonEquilibrium

The `MicrophysicsNonEq.jl` module describes a bulk parameterization of
  diffusion of water vapour on cloud droplets and cloud ice crystals
  modeled as a relaxation to equilibrium.

The cloud microphysics variables are expressed as specific humidities:
  - `q_tot` - total water specific humidity,
  - `q_vap` - water vapor specific humidity,
  - `q_liq` - cloud water specific humidity,
  - `q_ice` - cloud ice specific humidity,

Parameters used in the parameterization are defined in
  [CLIMAParameters.jl](https://github.com/CliMA/CLIMAParameters.jl) package.
They consist of:

|    symbol                  |         definition                                        | units                    | default value          | reference |
|----------------------------|-----------------------------------------------------------|--------------------------|------------------------|-----------|
|``\tau_{cond\_evap}``       | cloud water condensation/evaporation timescale            | ``s``                    | ``10``                 |           |
|``\tau_{dep\_sub}``         | cloud ice deposition/sublimation timescale                | ``s``                    | ``10``                 |           |


## Cloud water condensation/evaporation

Condensation and evaporation of cloud liquid water is parameterized
  as a relaxation to equilibrium value at the current time step.
```math
\begin{equation}
  \left. \frac{d \, q_{liq}}{dt} \right|_{cond, evap} =
    \frac{q^{eq}_{liq} - q_{liq}}{\tau_{cond\_evap}}
\end{equation}
```
where:
 - ``q^{eq}_{liq}`` - liquid water specific humidity in equilibrium,
 - ``q_{liq}`` - liquid water specific humidity,
 - ``\tau_{cond\_evap}`` - relaxation timescale.

## Cloud ice deposition/sublimation

Deposition and sublimation of cloud ice is parameterized as
  a relaxation to equilibrium value at the current time step.
```math
\begin{equation}
  \left. \frac{d \, q_{ice}}{dt} \right|_{dep, sub} =
    \frac{q^{eq}_{ice} - q_{ice}}{\tau_{dep\_sub}}
\end{equation}
```
where:
 - ``q^{eq}_{ice}`` - ice specific humidity in equilibrium,
 - ``q_{ice}`` - ice specific humidity,
 - ``\tau_{dep\_sub}`` - relaxation timescale.

!!! note
    Both ``\tau_{cond\_evap}`` and ``\tau_{dep\_sub}`` are
    assumed constant here.
    It would be great to make the relaxation time a function of
    available condensation nuclei, turbulence intensity, etc.
    See works by [prof Raymond Shaw](https://www.mtu.edu/physics/department/faculty/shaw/)
    for hints.
    In particular, [Desai2019](@cite).

