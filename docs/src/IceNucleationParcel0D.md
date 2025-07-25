# 0-dimensional adiabatic parcel model

The model solves for saturation ratio in a 0-dimensional
  adiabatic parcel raising with constant velocity.
It is based on [KorolevMazin2003](@cite), as well as the cirrus box model
  [Karcher2006](@cite), [Tully2023](@cite) and
  is meant as a sandbox for testing different nucleation schemes.

We define liquid saturation ratio ``S_l``
```math
\begin{equation}
S_l = \frac{e}{e_{sl}}
\end{equation}
```
where:
- ``e`` - is the partial pressure of water vapor,
- ``e_{sl}`` - is the partial pressure of water vapor at saturation over liquid water.
Note that some models, for example [KorolevMazin2003](@cite), denote ``S_l``
  to be supersaturation ``\frac{e}{e_{sl}} - 1``.

The change in saturation ratio can be described as
```math
\begin{equation}
\frac{dS_l}{dt} = \frac{1}{e_{sl}} \frac{de}{dt} - \frac{e}{e_{sl}^2} \frac{de_{sl}}{dt}
\end{equation}
```
From ideal gas law the partial pressure of water vapor can be written as
```math
\begin{equation}
e = q_v p \frac{R_v}{R_a}
\end{equation}
```
where:
- ``q_v`` is the specific humidity
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
\frac{de_{sl}}{dt} = \frac{L_v e_{sl}}{R_v T^2} \frac{dT}{dt}
\end{equation}
```
where:
- ``L_v`` is the latent heat of vaporization,
- ``T`` is the temperature.

From the moist adiabatic assumption
```math
\begin{equation}
\frac{dT}{dt} = \frac{R_a T}{c_p p} \frac{dp}{dt} + \frac{L_v}{c_p} \frac{dq_{l,vap}}{dt} + \frac{L_s}{c_p} \frac{d q_{i,subl}}{dt} + \frac{L_f}{c_p} \frac{d q_{i,fus}}{dt}
\end{equation}
```
where:
- ``q_{l,vap}`` is the cloud liquid water specific content from vaporization/condensation,
- ``q_{i,subl}`` is the cloud ice specific content from sublimation/deposition,
- ``q_{i,fus}`` is the cloud ice specific content from melting/freezing,
- ``L_s`` is the latent heat of sublimation,
- ``L_f`` is the latent heat of fusion.

From hydrostatic balance and assuming constant vertical velocity:
```math
\begin{equation}
\frac{dp}{dt} = -\frac{p g}{R_a T} w
\end{equation}
```
where:
- ``g`` is the gravitational acceleration
- ``w`` is the constant vertical velocity.

Accounting for conservation of water, i.e. ``\frac{dq_v}{dt} = - \frac{dq_{l,vap}}{dt} - \frac{dq_{i,subl}}{dt}``,
and rearranging the terms
```math
\begin{equation}
\frac{dS_l}{dt} = a_1 w S_l - \left(a_2 + a_3 \right)S_l  \frac{dq_{l,vap}}{dt} - \left(a_2 + a_4 \right) S_l \frac{dq_{i,subl}}{dt} - a_5 S_l \frac{dq_{i,fus}}{dt}
\end{equation}
```
where:
```math
\begin{equation}
a_1 = \frac{L_v g}{c_p T^2 R_v} - \frac{g}{R_a T}
\end{equation}
```
```math
\begin{equation}
a_2 = \frac{1}{q_v}
\end{equation}
```
```math
\begin{equation}
a_3 = \frac{L_v^2}{R_v T^2 c_p}
\end{equation}
```
```math
\begin{equation}
a_4 = \frac{L_v L_s}{R_v T^2 c_p}
\end{equation}
```
```math
\begin{equation}
a_5 = \frac{L_v L_f}{R_v T^2 c_p}
\end{equation}
```

Saturation ratio over ice can then be related to ``S_l`` by the relation
```math
\begin{equation}
  S_i = \xi S_l
  \label{eq:supersat_relation}
\end{equation}
```
where ``\xi = \frac{e_{sl}}{e_{si}}`` is the ratio of saturation vapor pressure over liquid and saturation vapor pressure over ice.

The crux of the problem is modeling the ``\frac{dq_l}{dt}`` and ``\frac{dq_i}{dt}``
  for different homogeneous and heterogeneous ice nucleation paths.

## Supported size distributions
Currently, the parcel model supports monodisperse and gamma size distributions of cloud droplets and ice crystals, and solves prognostic equations for the cloud water and cloud ice specific contents (`q_l`, `q_i`) and number concentrations (`N_l`, `N_i`). Additionally, a `monodisperseMix` option is now supported for cloud liquid droplets. This allows representing the droplet population as a mixture of two monodisperse modes: one corresponding to an initial set of preexisting droplets, and the other representing droplets formed through activation during the simulation. This feature enables more realistic treatment of scenarios involving preexisting hydrometeors.

For a monodisperse size distribution of cloud droplets or ice crystals
```math
\begin{equation}
  n(r) = N \delta(r-\bar{r})
\end{equation}
```
and assuming spherical shape, the mean radius is defined as
```math
\begin{equation}
  \bar{r} = \left( \frac{3 \rho_a q}{4 \pi \rho N} \right)^{1/3}
\end{equation}
```
where:
 - ``\rho_a`` is the air density,
 - ``\rho``  is the cloud water or cloud ice density.

For a gamma distribution of cloud droplets
```math
\begin{equation}
  n(r) = N_0 \; r \; \exp(-\lambda r)
\end{equation}
```
where:
 - ``\lambda`` and ``N_0`` are the free parameters.
```math
\begin{equation}
  N = {\int_0^\infty \, n(r) \, dr}, \;\;\;\;\;
  q = \frac{1}{\rho_a}{\int_0^\infty \, \frac{4}{3} \pi \rho r^3 n(r) \, dr}
\end{equation}
```
As a result
```math
\begin{equation}
  \lambda = \left(\frac{32 \pi N \rho}{q \rho_a}\right)^{1/3}, \;\;\;\;\;
  N_0 = N \lambda^2
\end{equation}
```

## Supported source terms
### Aerosol Activation
Aerosol activation is described by ([see discussion](https://clima.github.io/CloudMicrophysics.jl/dev/AerosolActivation/#Number-and-mass-of-activated-particles)).
It is inherently assumed that the aerosols have a lognormal size distribution. For simplicity, the parcel accepts one mode and one aerosol type at a time, therefore, internal mixing is not needed. The maxiumum supersaturation as described in the above mentioned documentation is replaced by the liquid supersaturation in the parcel as it evolves over time.

!!! note
    The standard deviation and $r_{\text{mean}}$ of the aerosol size distribution represent the initial (pre-activation) aerosol population and are assumed to remain fixed throughout the simulation. They do not evolve as aerosols activate.

### Condensation growth

There are both equilibrium and non-equilibrium options for condensation. The non-equilibrium equations are described in the non-equilibrium documentation ([here](https://clima.github.io/CloudMicrophysics.jl/dev/MicrophysicsNonEq)).

For the equilibrium formulation, the diffusional growth of individual cloud droplet is described by
([see discussion](https://clima.github.io/CloudMicrophysics.jl/dev/Microphysics1M/#Snow-autoconversion)),
```math
\begin{equation}
  r_l \frac{dr_l}{dt} = \frac{1}{\rho_l} \, (S_l - 1) \, G_l(T)
\end{equation}
```
where:
- ``r_l`` is the droplet radius
- ``\rho_l`` is the water density
- ``G_l(T) = \left(\frac{L_v}{KT} \left(\frac{L_v}{R_v T} - 1 \right) + \frac{R_v T}{e_{sl} D} \right)^{-1}``
    combines the effects of thermal conductivity and water diffusivity,
- ``K`` is the thermal conductivity of air,
- ``D`` is the diffusivity of water vapor.

!!! note
    Some formulations don't consider the ``-1`` term in the denominator of G(T)

Assuming spherical water droplets, the change in droplet mass ``m_l`` can be described as
```math
\begin{equation}
  \frac{dm_l}{dt} = 4 \pi \rho_l r^2 \frac{dr_l}{dt} = 4 \pi r_l (S_l - 1) G_l(T)
\end{equation}
```
Integrating over the assumed size distribution of droplets
```math
\begin{equation}
   \frac{dq_l}{dt} = \frac{1}{V} \frac{1}{\rho_a} 4 \pi (S_l - 1) G_l(T) \int_{0}^{\infty} n(r) r dr =
                                 \frac{1}{\rho_a} 4 \pi (S_l - 1) G_l(T) N_{tot} \bar{r}
\end{equation}
```
where:
- ``N_{tot}`` is the total number concentration of droplets per volume of air ``V``
- ``\bar{r}`` is their mean radius
- ``\rho_a`` is the density of air.

### Deposition growth

Again, there are both equilibrium and non-equilibrium options for deposition. The non-equilibrium equations are described in the non-equilibrium documentation ([here](https://clima.github.io/CloudMicrophysics.jl/dev/MicrophysicsNonEq)).

For the equilibrium formulation, we use a similar process as with condensation.
For a case of a spherical ice particle growing through water vapor deposition
```math
\begin{equation}
  \frac{dm_i}{dt} = 4 \pi \, r_i \, \alpha_m \, (S_i - 1) \, G_i(T)
  \label{eq:massratesphere}
\end{equation}
```
where:
 - ``r_i`` is the ice particle radius,
 - ``\alpha_m`` is the accommodation coefficient that takes into account the fact that not all water vapor molecules
     that arrive at the particle surface will join the growing crystal,
 - ``G_i(T) = \left(\frac{L_s}{KT} \left(\frac{L_s}{R_v T} - 1 \right) + \frac{R_v T}{e_{si} D} \right)^{-1}``
     combines the effects of thermal conductivity and water diffusivity,

!!! note
    The ``r_i`` in eq. (\ref{eq:massratesphere})
    should be replaced by the capacitance ``C``,
    when considering non-spherical particles.
    For a sphere ``C=r``, for a circular disc ``C=2r/\pi``.

It follows that
```math
\begin{equation}
  \frac{dq_i}{dt} = \frac{1}{\rho_a} \alpha_m 4 \pi (S_i - 1) G_i(T) N_{act} \bar{r}
\end{equation}
```
where:
- ``N_{act}`` is the number of activated ice particles.

### Deposition Nucleation on dust particles
There are multiple ways of running deposition nucleation in the parcel.
  `"MohlerAF_Deposition"` will trigger an activated fraction approach
  from [Mohler2006](@cite). `"MohlerRate_Deposition"` will trigger a
  nucleation rate approach from [Mohler2006](@cite). For both approaches,
  there is no nucleation if saturation over ice exceeds 1.35 as conditions
  above this value will result in nucleation in a different mode.
  `"ActivityBasedDeposition"` will trigger a water activity based approach
  from [Alpert2022](@cite). In this approach, ice production rate ``P_{ice, depo}``
  is calculated from
```math
\begin{equation}
  P_{ice, depo} = \left[ \frac{dN_i}{dt} \right]_{depo} = J_{depo}\;A_{aero}\;N_{aero}
  \label{eq:ActivityBasedDeposition_P_ice}
\end{equation}
```
where ``N_{areo}`` is total number of unactiviated ice nucleating particles and
  ``A_{aero}`` is surface area of those INP.
The deposition nucleation methods are parameterized as described in
  [Ice Nucleation](https://clima.github.io/CloudMicrophysics.jl/dev/IceNucleation/).

### Immersion Freezing
Following the water activity based immersion freezing model (ABIFM), the ABIFM derived
  nucleation rate coefficient, ``J_{immer}``, can be determined. The ice production rate,``P_{ice, immer}``,
  per second via immersion freezing can then be calculating using
```math
\begin{equation}
  P_{ice, immer} = \left[ \frac{dN_i}{dt} \right]_{immer} = J_{immer}\;A_{aero}(N_{liq})
  \label{eq:ABIFM_P_ice}
\end{equation}
```
where ``N_{liq}`` is total number of ice nuclei containing droplets and
  ``A_{aero}`` is surface area of the ice nucleating particle.

### Homogeneous Freezing
Homogeneous freezing follows the water-activity based model described in the
  [Ice Nucleation](https://clima.github.io/CloudMicrophysics.jl/dev/IceNucleation/) section which gives a nucleation rate coefficient of
  units ``cm^{-3}s^{-1}``.
The ice production rate from homogeneous freezing can then be determined:
```math
\begin{equation}
  P_{ice, hom} = \left[ \frac{dN_i}{dt} \right]_{hom} = J_{hom}V(N_{liq})
  \label{eq:hom_P_ice}
\end{equation}
```
where ``N_{liq}`` is total number of ice nuclei containing droplets and
  ``V`` is the volume of those droplets.

### P3 Ice Nucleation Parameterizations
The parcel also includes ice nucleation parameterizations used in
  the P3 scheme as described in [MorrisonMilbrandt2015](@cite).
  Deposition nucleation is based on the ice crystal number parameterization
  from Cooper (1986). The heterogeneous freezing parameterization, which
  follows Bigg(1953) with parameters from Barklie aand Gokhale (1959), is
  treated as immersion freezing in the parcel. Homogeneous freezing happens
  instantaneously at 233.15K.

### Frostenberg et al. 2023 Immersion Freezing
An immersion freezing parameterization based on [Frostenberg2023](@cite) is also available.
The concentration of ice nucleating particles (INPC) depends only on air temperature,
  and is based on a lognormal relative frequency distribution.
New ice crystals are created if the INPC exceeds the existing concentration of ice crystals,
  provided there are sufficient numbers of cloud liquid droplets to freeze.

Three different implementations of this parametrization are used in the parcel model:
- `mean` - in which INPC is equal to its mean value defined in [Frostenberg2023](@cite).
- `random` - in which INPC is sampled randomly from the distribution defined in [Frostenberg2023](@cite).
  The number of model time steps between sampling is set by `sampling_interval`.
- `stochastic` - in which INPC is solved for as a stochastic process,
  with the mean and standard deviation defined in [Frostenberg2023](@cite).
  The inverse timescale of the process is set by ``\gamma``.

The stochastic implementation is based on the [Ornstein-Uhlenbeck process](https://en.wikipedia.org/wiki/Ornstein–Uhlenbeck_process),
in which the variable ``x`` is a mean-reverting process perturbed by Gaussian random noise (i.e. increments of the Wiener process ``W``):
```math
\begin{equation}
  dx = - \gamma(x - \mu)dt + \sqrt{2\gamma} \sigma dW; \quad\quad   dW \sim N(0, dt),
\end{equation}
```
where ``N`` is a zero-mean normal distribution with variance ``dt``.
For constant ``\gamma`` and ``\sigma``, and given some initial condition ``x(0)=x_0``, ``x`` has the analytical solution:
```math
\begin{equation}
  x(t) =
      x_0 e^{-\gamma t} + \mu (1 - e^{-\gamma t})
    + \sqrt{2\gamma} \sigma \int_0^t e^{-\gamma(t-s)} dW,
\end{equation}
```
where ``\tau \equiv 1 / \gamma`` is the assumed timescale of the process. The process mean is ``x_0 e^{-\gamma t} + \mu (1 - e^{-\gamma t})``. We can calculate the variance ``\mathbb{V}(t)`` as,
```math
\begin{equation}
  \mathbb{V}(t)
  = 2\gamma \sigma^2 \int_0^t e^{-2\gamma(t-s)} ds
  = \frac{g^2}{2\gamma} \left( 1 - e^{-2\gamma t} \right).
\end{equation}
```

We use this process to model ``x=\log(\text{INPC})``, which tends toward a temperature-dependent mean value ``\mu(T)``.
The equation for ``\log(\text{INPC})`` is then,
```math
\begin{equation}
  d\log(\text{INPC}) =
    - \frac{\log(\text{INPC}) - μ}{\tau} dt
    + \sigma \sqrt{\frac{2}{\tau}} dW
\end{equation}
```

This equation is currently implemented with the simple [Euler-Maruyama method](https://en.wikipedia.org/wiki/Euler–Maruyama_method), which is the stochastic analog of the forward Euler method for (deterministic) ordinary differential equations, so that
```math
\begin{equation}
  \log(\text{INPC})_{t+dt} = \log(\text{INPC})_{t}
    - \gamma\left(\log(\text{INPC})_t - μ(T_t)\right) dt
    + \sigma \sqrt{2\gamma dt} z_t
\end{equation}
```
where ``z_t \sim N(0,1)`` is a standard normal random variable.

## Example figures

Here we show various example simulation results from the adiabatic parcel
  model. This includes examples with deposition nucleation on dust,
  liquid processes only, immersion freezing with condensation and deposition growth,
  and homogeneous freezing with deposition growth.

First, we check that aerosol activation works reasonably within the parcel.

```@example
    include("../../parcel/Example_AerosolActivation.jl")
```
![](Parcel_Aerosol_Activation.svg)

The following examples show ice nucleation, starting with deposition
  freezing on dust.
The model is run three times using the `"MohlerAF_Deposition"` approach
  for 30 minutes simulation time, (shown by three different colors on the plot).
Between each run the specific humidity is changed,
  while keeping all other state variables the same as at the last time step
  of the previous run.
The prescribed vertical velocity is equal to 3.5 cm/s.

Supersaturation is plotted for both liquid (solid lines) and ice (dashed lines).
The pale blue line uses the `"MohlerRate_Deposition"` approach.
  We only run it for the first GCM timestep because the rate approach requires
  the change in ice saturation over time. With the discontinuous jump in saturation,
  the parameterization is unable to determine a proper nucleation rate. When we force
  the initial ice crystal number concentration for this simulation to match
  that in the `"MohlerAF_Deposition"` approach, we obtain the same results as
  in the `"MohlerAF_Deposition"` approach for the first GCM timestep.

```@example
using Suppressor: @suppress #hide
  @suppress begin #hide
     include("../../parcel/Example_Tully_et_al_2023.jl")
end # hide
```
![](cirrus_box.svg)

The water activity based parameterization for deposition nucleation shows
  similar outcomes when compared to the `"MohlerRate_Deposition"` approach.
  Here, we run the parcel for 100 secs for all available aerosol types. The
  solid lines correspond to the `"MohlerRate_Deposition"` approach while the
  dashed lines correspond to `"ActivityBasedDeposition"`. Note that there
  is no common aerosol type between the two parameterizations.

```@example
  include("../../parcel/Example_Deposition_Nucleation.jl")
```
![](deposition_nucleation.svg)

In the plots below, the parcel model is ran with only condensation (no ice or freezing)
  assuming either a monodisperse or a gamma distribution of droplets.
It is compared to [Rogers1975](@cite). This is to ensure the parcel can run
with liquid phase and different size distributions without issues.

```@example
  include("../../parcel/Example_Liquid_only.jl")
```
![](liquid_only_parcel.svg)

The following plots show the parcel model running with homogeneous freezing and
  depositional growth assuming a lognormal distribution of aerosols.
  It is compared against [Jensen2022](@cite). Note that running with the initial
  conditions described in [Jensen2022](@cite) results in a ``\Delta a_w`` smaller
  than the minimum valid value for the ``J_{hom}`` parameterization. We have forced
  the ``\Delta a_w`` to be equal to the minimum valid value in this example only
  for demonstrative purposes.

```@example
  include("../../parcel/Example_Jensen_et_al_2022.jl")
```
![](Jensen_et_al_2022.svg)

Shown below are three separate parcel simulations for deposition nucleation,
  immersion freezing, and homogeneous freezing as parameterized in the P3 scheme.
  Note that initial temperature varies for each run. The deposition nucleation
  run does not conserve INP number, while the other two freezing modes do.
  Updraft velocity is set to 0.5 m/s.
```@example
  include("../../parcel/Example_P3_ice_nuc.jl")
```
![](P3_ice_nuc.svg)

The results of the immersion freezing parameterization from Frostenberg et al. 2023 are shown below.
The following plot shows resuls of the parcel model with the `mean` (black line), `random` (dotted lines)
and `stochastic` (solid lines) parameterization options. We show results for two sampling intervals
``\Delta t`` (random), two process time scales ``\tau`` (stochastic), and two model time steps `dt`.

```@example
  include("../../parcel/Example_Frostenberg_Immersion_Freezing.jl")
```
![](frostenberg_immersion_freezing.svg)

Below, we show how the non-equilibrium formulation is able to represent the
Wegener-Bergeron-Findeisen (WBF) regime in the parcel model.
These plots show the liquid supersaturation, ice supersaturation,
temperature, specific humidity `q_{vap}`, specific content of liquid `q_{liq}`,
and specific content of ice `q_{ice}` over time.
When the supersaturation is negative evaporation/sublimation occurs,
and when it is positive condensation/deposition occurs.
In the example below, the water vapor pressure is smaller than the saturation water vapor pressure of liquid
but larger than that of ice, and as a result liquid water evaporates while ice grows.

```@example
  include("../../parcel/Example_NonEq.jl")
```
![](noneq_parcel.svg)
