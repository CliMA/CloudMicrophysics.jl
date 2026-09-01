# Microphysics NonEquilibrium

The `MicrophysicsNonEq.jl` module describes a bulk parameterization of
  diffusion of water vapor on cloud droplets and cloud ice crystals
  modeled as a relaxation to equilibrium.

The cloud microphysics variables are expressed as specific contents:
  - `q_tot` - total water specific content,
  - `q_vap` - water vapor specific content (i.e., specific humidity),
  - `q_lcl` - cloud water specific content,
  - `q_icl` - cloud ice specific content,

Parameters used in the parameterization are defined in
  [ClimaParams.jl](https://github.com/CliMA/ClimaParams.jl) package.
They consist of:

|    symbol  |         definition                             | units | default value |
|------------|------------------------------------------------|-------|---------------|
|``\tau_{l}``| cloud water condensation/evaporation timescale | ``s`` | ``10``        |
|``\tau_{i}``| cloud ice deposition/sublimation timescale     | ``s`` | ``10``        |

!!! note
    For ice, the deposition timescale ``\tau_{dep}`` can optionally be computed
    from the prescribed cloud-ice number concentration ``N_0`` (see
    [below](@ref ice-relaxation-timescale-prescribed-ice-number)) or from the
    [Frostenberg2023](@cite) INP parameterization (see
    [below](@ref ice-relaxation-timescale-frostenberg-et-al-2023)).
    With the Frostenberg option, the sublimation timescale ``\tau_{sub}`` remains
    at the constant ``\tau_i``.

## Condensation/evaporation and deposition/sublimation from Morrison and Milbrandt 2015

Condensation/evaporation and deposition/sublimation rates are based on
  the difference between the specific humidity and the
  specific humidity at saturation over liquid and ice at the current temperature.
The process is modeled as a relaxation with a constant timescale.
This formulation is derived from [MorrisonGrabowski2008_supersat](@cite)
  and [MorrisonMilbrandt2015](@cite), but without imposing exponential time integrators.

!!! note
    The [MorrisonGrabowski2008_supersat](@cite) and [MorrisonMilbrandt2015](@cite)
    papers use mass mixing ratios, not specific contents.
    Additionally, in their formulations they consider two different categories for liquid:
    cloud water and rain. For now we only consider cloud water and use a single relaxation timescale
    ``\tau_l`` (liquid) rather than separate ``\tau_c`` (cloud) and ``\tau_r`` (rain) values.

```math
\begin{equation}
   \left. \frac{d \, q_{lcl}}{dt} \right|_{cond, evap} = \frac{q_{vap} - q_{sl}}{\tau_l \Gamma_l}; \;\;\;\;\;\;\;
   \left. \frac{d \, q_{icl}}{dt} \right|_{dep, sub}   = \frac{q_{vap} - q_{si}}{\tau_i \Gamma_i}
\end{equation}
```
where:
- ``q_{vap}`` is the specific humidity
- ``q_{sl}``, ``q_{si}`` is the specific humidity at saturation over liquid and ice
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

When the air is subsaturated, the evaporation and sublimation rates are
  additionally limited by the available condensate:
```math
\begin{equation}
   \left. \frac{d \, q_{k}}{dt} \right|_{evap, sub} =
   - \frac{\min(q_{sk} - q_{vap}, \; \max(0, q_{k}))}{\tau_k \Gamma_k},
   \;\;\;\;\;\;\; k \in \{lcl, icl\}.
\end{equation}
```
Two temperature-based limiters are applied to the tendencies:

1. **INP limiter (ice)**: Ice deposition (positive tendency) is suppressed
   at temperatures above freezing (``T > T_{freeze}``), where no ice
   nucleating particles (INPs) are available.
   Ice sublimation (negative tendency) is unaffected.

2. **Homogeneous limiter (liquid)**: Liquid condensation (positive tendency)
   is suppressed at temperatures below the homogeneous nucleation
   threshold (``T < T_{hom} \approx 233\,\text{K}``), where all liquid
   water would freeze instantaneously.
   Liquid evaporation (negative tendency) is unaffected, allowing
   any remaining cloud water to dissipate.

Both limiters apply to all relaxation timescale options.

!!! note
    Both ``\tau_{l}`` and ``\tau_{i}`` are assumed to be constant by default.
    Making the relaxation timescales functions of available condensation
    nuclei, turbulence intensity, and other environmental conditions is
    left for future work; see for example [Desai2019](@cite).

Note that these forms of condensation/sublimation and deposition/sublimation
  are equivalent to those described in the adiabatic parcel model with some rearrangements and assumptions.
To see this, it is necessary to use the definitions of ``\tau``, ``q_{sl}``, and the thermal diffusivity ``D_v``:

```math
\begin{equation}
  \tau = 4 \pi N_{tot} \bar{r}, \;\;\;\;\;\;\;\;
  q_{sl} = \frac{e_{sl}}{\rho R_v T}, \;\;\;\;\;\;\;\;
  D_v = \frac{K}{\rho c_p}.
\end{equation}
```
If we then assume that the supersaturation ``S`` can be approximated by the specific contents (this is only exactly true for mass mixing ratios):
```math
\begin{equation}
    S_l = \frac{q_{vap}}{q_{sl}},
\end{equation}
```
we can write
```math
\begin{equation}
  q_{vap} - q_{sl} = q_{sl}(S_l - 1).
\end{equation}
```
``\Gamma_l`` and ``\Gamma_i`` then are equivalent to the ``G`` function used in our parcel model and parameterizations.

```@example
include("plots/NonEqCondEvapRate.jl")
```
![](condensation_evaporation_ql_z.svg)
![](condensation_evaporation_ql_T.svg)

## [Ice relaxation timescale — prescribed ice number](@id ice-relaxation-timescale-prescribed-ice-number)

`PrescribedIceNumber` option derives the relaxation timescale ``\tau_i``
  based on the prescribed cloud-ice number concentration ``N_0``.

```math
\begin{equation}
  \tau_{i} = \frac{1}{4 \pi \, D_v \, N_0 \, r}
\end{equation}
```
where ``D_v`` is the water vapor diffusivity (``\text{m}^2/\text{s}``)
  and ``r`` is the ice crystal radius.

Assuming a mono-disperse distribution and spherical ice crystals we estimate
  the radius as
```math
\begin{equation}
  r = max \left(\left(\frac{3 \, q_{icl} \, \rho }{4 \pi \, N_0 \, \rho_i}\right)^{1/3} \, , \, r_0 \right)
\end{equation}
```
where ``\rho`` is the air density,
      ``q_{icl} is cloud ice specific humidity,
      and ``\rho_i`` is the ice density.
A minimum radius ``r_0 = 1\,\mu m`` is enforced.

Figure below show the relaxation timescale for different assumed number concentrations and
  cloud ice specific humidities. The constant timescale at low cloud ice values is the
  result of the minimum assumed particle radius.

```@example
include("plots/plotting_tau_relax_prescribed_N0.jl")
```
![](tau_relax_prescribed_N0.svg)

## [Ice relaxation timescale — Frostenberg et al. (2023)](@id ice-relaxation-timescale-frostenberg-et-al-2023)

The constant ice relaxation timescale ``\tau_i`` can be replaced with
  a temperature-dependent timescale derived from the
  Frostenberg et al. [Frostenberg2023](@cite)
  parameterization of ice nucleating particle (INP) concentrations.
This makes the deposition timescale physically dependent on the
  number of available INPs and the ice crystal size.
This could be a good approximation at cold temperatures, but should not be
  used as a sole timescale due to the importance of ice-multiplication processes
  in warmer temperatures.

The INP number concentration is estimated as a function of temperature:
```math
\begin{equation}
  N_{icl} = \exp\!\left(\overline{\ln(\mathrm{INPC})}(T)\right)
\end{equation}
```
where ``\overline{\ln(\mathrm{INPC})}(T)`` is the mean log-INP concentration
from the Frostenberg et al. (2023) parameterization.

Given ``N_{icl}`` and the cloud ice specific content ``q_{icl}``,
  the mean crystal radius is computed assuming spherical ice particles
  with a monodisperse size distribution:
```math
\begin{equation}
  r = \left(\frac{3 \, q_{icl}}{4 \pi \, N_{icl} \, \rho_i}\right)^{1/3}
\end{equation}
```
A minimum radius ``r_0 = 1\,\mu m`` is enforced to avoid singularities
  when ``q_{icl} = 0`` or ``N_{icl} = 0``.

The Frostenberg deposition timescale is then:
```math
\begin{equation}
  \tau_{dep} = \frac{1}{4 \pi \, D_v \, N_{icl} \, r_{safe}}
\end{equation}
```
where ``D_v`` is the water vapor diffusivity and ``r_{safe} = \max(r, r_0)``.

### Asymmetric deposition and sublimation timescales

The ice tendency uses **different timescales** for deposition and sublimation:

```math
\begin{equation}
  \left. \frac{d \, q_{icl}}{dt} \right|_{dep, sub} =
  \begin{cases}
    \dfrac{q_{vap} - q_{si}}{\tau_{dep} \, \Gamma_i} & \text{if } q_{vap} > q_{si} \text{ (deposition)} \\[10pt]
    \dfrac{-\min(-\Delta q,\; q_{icl})}{\tau_{sub} \, \Gamma_i} & \text{if } q_{vap} \le q_{si} \text{ (sublimation)}
  \end{cases}
\end{equation}
```
where:
- ``\tau_{dep}`` is the Frostenberg INP-based timescale (temperature-dependent),
- ``\Gamma_i`` is the psychometric correction factor (applied to both deposition and sublimation),
- ``\tau_{sub}`` is a constant sublimation timescale (default: ``\tau_i = 10\,s``),
- ``\Delta q = q_{vap} - q_{si}`` is the saturation excess.

```@example
include("plots/plotting_tau_relax_frostenberg.jl")
```
![](tau_relax_frostenberg.svg)

## Cloud condensate sedimentation

We use the analytical Stokes-regime terminal velocity for cloud liquid droplets
  and the Chen et al. [Chen2022](@cite) parameterization
  (with the coefficients for small ice particles)
  for cloud ice sedimentation velocities.
In the 1-moment precipitation scheme, we assume that cloud condensate is a continuous field
  and doesn't introduce an explicit particle size distribution.
For simplicity, we assume a monodisperse size distribution
  and compute the group terminal velocity based on the volume diameter
  and prescribed number concentration:

```math
\begin{equation}
  D_{vol} = \left( \frac{6 \, \rho_{air} \, q}{\pi \, N \, \rho} \right)^{1/3}
\end{equation}
```
where:
 - ``\rho_{air}`` is the air density,
 - ``q`` is the cloud liquid or cloud ice specific content,
 - ``N`` is the prescribed number concentration (``500/cm^3`` by default),
 - ``\rho`` is the cloud water or cloud ice density.

The sedimentation velocity then is
```math
\begin{equation}
  v_t = v_{term}(D_{vol}).
\end{equation}
```

!!! note
    The Chen et al. [Chen2022](@cite) coefficients used for cloud ice
    were fitted for small ice particles (below 625 μm).
    The Stokes-regime formula used for cloud liquid is valid for small
    droplets, where the drag is dominated by viscous forces.
