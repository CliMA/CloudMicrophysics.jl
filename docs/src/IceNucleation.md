# Ice Nucleation

Ice nucleation describes the process of forming ice crystals
  from aerosol particles and/or liquid droplets.
The `IceNucleation.jl` module includes:
  - the parameterization of activation of dust aerosol particles into ice crystals
    via deposition of water vapor,
  - water activity based parameterization of immersion freezing,
  - water activity based parameterization of homogeneous freezing,
  - parameterization of temperature-dependent aerosol-independent ice nucleating particles concentration for immersion freezing.

The parameterization for deposition on dust particles is an implementation of
  the empirical formulae from [Mohler2006](@cite)
  and is valid for two types of dust particles:
  Arizona Test Dust and desert dust from Sahara.
  The parameterization for immersion freezing is an implementation of [KnopfAlpert2013](@cite)
  and is valid for droplets containing sulphuric acid.
  The parameterization for homogeneous freezing is an implementation of [Koop2000](@cite).
  The parameterization for ice nucleating particles concentration for immersion freezing is an implementation of [Frostenberg2023](@cite).

!!! note

    Future work includes refining the homogeneous freezing
    parameterization and modeling the competition between
    freezing modes.

## Activated fraction for deposition freezing on dust
There are 2 parameterizations from [Mohler2006](@cite) available: one
  which calculates the activated fraction and one which calculates nucleation
  rate.
The activated fraction parameterization follows eq. (3) in the paper.
```math
\begin{equation}
f_i(S_i) = exp[a(S_i - S_0)] - 1
\end{equation}
```
where
  - ``f_i`` is the activated fraction
     (the ratio of aerosol particles acting as ice nuclei to the total number of aerosol particles),
  - ``a`` is a scaling parameter dependent on aerosol properties and temperature,
  - ``S_i`` is the ice saturation ratio,
  - ``S_0`` is an empirically derived threshold ice saturation ratio.
The other parameterization models the nucleation rate of ice
  as an empirical function of ice saturation ratio,
  see eq. (5) in [Mohler2006](@cite).
```math
\begin{equation}
\frac{dn_{icl}}{dt} = N_{aer} a \frac{dS_i}{dt}
\end{equation}
```
where:
  - ``N_{aer}`` is the number of available aerosol/ice nuclei,
  - ``a`` is a scaling parameter dependent on aerosol properties and temperature,
  - ``S_i`` is the ice saturation ratio
    (the ratio of water vapor partial pressure and the water vapor partial pressure at saturation over ice).

Limited experimental values for the free parameters ``a`` and ``S_0`` can be found in [Mohler2006](@cite). These
  free parameters are strongly dependent on aerosol properties and temperature.

!!! note

    For ``f_i`` values above 0.08 or ``S_i`` between 1.35 and 1.5,
    freezing occurs in a different ice nucleation mode
    (either a second deposition or other condensation type mode).

### Water Activity Based Deposition Nucleation
The water activity based deposition nucleation model is analagous to ABIFM
  for immersion freezing (see [ABIFM for Sulphuric Acid Containing Droplets](https://clima.github.io/CloudMicrophysics.jl/dev/IceNucleation/#ABIFM-for-Sulphuric-Acid-Containing-Droplets)
  section below). It calculates a nucleation rate coefficient, ``J``, which
  describes the number of ice nuclei formed per unit area of INP per unit time
  dependent on the water activity criterion, ``\Delta  a_w``, and aerosol type.
  The form of this empirical parameterization is taken from [KnopfAlpert2013](@cite).
  Currently, we have parameters for kaolinite, feldspar, and ferrihydrite derived
  from [China2017](@cite) and [Alpert2022](@cite), as well as for illite,
  Arizona Test Dust, Saharan dust, Asian dust, and a generic dust type.

```math
\begin{equation}
  log_{10}J_{deposition} = m \Delta a_w + c
\end{equation}
```
where ``J`` is in units of ``cm^{-2}s^{-1}``. Note that our source code returns
  ``J`` in SI units. ``m`` and ``c`` are aerosol dependent coefficients. They will
  have different values than those for ABIFM.

### Water activity based deposition nucleation plot
The following plot shows ``J`` as a function of ``\Delta a_w`` as compared to
  figure 6 in Alpert et al 2013 and figure S5 in supplementary material of China et al 2017. Intent of this
  plot is to prove that ``J`` is correctly parameterized as a function
  of ``\Delta a_w``, with no implications of whether ``\Delta a_w`` is properly
  parameterized. For more on water activity, please see above section.
```@example
include("plots/activity_based_deposition.jl")
```
![](water_activity_depo_nuc.svg)


### Cooper (1986) Ice Crystal Number
The P3 scheme as described in [MorrisonMilbrandt2015](@cite) follows
  the implementation of [Thompson2004](@cite) where a parameterization from
  Cooper 1986 was used. Number of ice crystals formed is derived from
  ambient temperature.
```math
\begin{equation}
  N_i = 0.005 \, exp[0.304(T_0 - T)]
\end{equation}
```
Where ``T_0 = 273.15 \, K``, ``T`` is the ambient temperature, and ``N_i``
  is the number of ice crystals per liter of air
  (the implementation returns the value converted to ``m^{-3}``
  and returns zero for ``T \geq T_0``).
Because the parameterization is prone to overpredict at low temperatures,
  ``N_i = N_i(T = 233 \, K)`` for all values of ``N_i`` at which ``T < 233 \, K``.

## Immersion Freezing
### ABIFM for Sulphuric Acid Containing Droplets
Water Activity-Based Immersion Freezing Model (ABFIM)
  is a method of parameterizing immersion freezing inspired by the time-dependent
  classical nucleation theory (CNT). More on CNT can be found in [Karthika2016](@cite).
  The nucleation rate coefficient, ``J``, describes the number of ice nuclei formed per unit area
  per unit time and can be determined by the water activity, ``a_w``. This parameterization follows
  [KnopfAlpert2013](@cite). In this model, aerosols are assumed to contain an insoluble and
  soluble material. When immersed in water, the soluble material diffuses into the liquid water
  to create a sulphuric acid solution.

Using empirical coefficients, ``m`` and ``c``, from [KnopfAlpert2013](@cite)
  and related studies (depending on the aerosol type),
  the heterogeneous nucleation rate coefficient in units of ``cm^{-2}s^{-1}``
  can be determined by the linear equation
```math
\begin{equation}
  log_{10}J_{ABIFM} = m \Delta a_w + c
\end{equation}
```
The water activity functions needed to compute ``\Delta a_w`` can be found
  in `Common.jl`. More information on
  them can be found in the `Water Activity` section. ``m`` and ``c`` here are different
  from the ``m`` and ``c`` parameters for deposition nucleation.

!!! note

    Our source code for the nucleation rate coefficient returns
    ``J`` in base SI units.

Once ``J_{ABIFM}`` is calculated, it can be used to determine the ice production rate, ``P_{ice}``,
per second via immersion freezing.
```math
\begin{equation}
  P_{ice} = J_{ABIFM}A(N_{aer} - N_{icl})
\end{equation}
```
where ``A`` is surface area of an individual ice nuclei, ``N_{aer}`` is total number
  of ice nuclei, and ``N_{icl}`` is number of ice crystals already in the system.

### ABIFM Example Figures
The following plot shows ``J`` as a function of ``\Delta a_w`` as compared to
  figure 1 in Knopf & Alpert 2013. Solution droplets were assumed to contain
  a constant 10% wt. sulphuric acid. Changing the concentration will simply
  shift the line, following Knopf & Alpert's parameterization. As such, this
  plot is just to prove that ``J`` is correctly parameterized as a function
  of ``\Delta a_w``, with no implications of whether ``\Delta a_w`` is properly
  parameterized. For more on water activity, please see above section.
```@example
include("plots/KnopfAlpert2013_fig1.jl")
```
![](Knopf_Alpert_fig_1.svg)

The following plot shows J as a function of temperature as compared to figure 5a in Knopf & Alpert 2013.

```@example
include("plots/KnopfAlpert2013_fig5.jl")
```
![](KnopfAlpert2013_fig5.svg)
Note that water activity of the droplet was assumed equal to relative humidity so that:
```math
\begin{equation}
  a_{w} = \frac{p_{sat}(T = T_{dew})}{p_{sat}(T)}
\end{equation}
```
where `T_dew` is the dewpoint (in this example, it is constant at -45C).

!!! note

    For the same figure in Knopf & Alpert,
    the mixed-phase cloud uses T_{dew} = T.
    We are unsure when to use constant T_{dew}
    or set it equal to T.

It is also important to note that this plot is reflective of cirrus clouds
  and shows only a very small temperature range. The two curves are slightly
  off because of small differences in parameterizations for vapor pressures.

### Bigg (1953) Volume and Time Dependent Heterogeneous Freezing
Heterogeneous freezing in the P3 scheme as described in [MorrisonMilbrandt2015](@cite)
  follows the parameterization from [Bigg1953](@cite) with parameters from
  [BarklieGokhale1959](@cite). The number of ice nucleated in a timestep via
  heterogeneous freezing is determined by
```math
\begin{equation}
  N_{icl} = N_{lcl} \left[ 1 - exp(-B V_{l} \Delta t exp(aT_s)) \right]
\end{equation}
```
where `a` and `B` are parameters taken from [BarklieGokhale1959](@cite) for
  rainwater as ``0.65`` and ``2 \times 10^{-4} \, cm^{-3}``
  (``200 \, m^{-3}`` in the implementation), respectively.
  `T_s` is the difference between
  273.15 K and ambient temperature. `V_l` is the volume of droplets to be
  frozen and `\Delta t` is timestep in which freezing occurs.

## Homogeneous Freezing
### Homogeneous Freezing for Sulphuric Acid Containing Droplets
Homogeneous freezing occurs when supercooled liquid droplets freeze on their own.
  Closely based on [Koop2000](@cite), this parameterization determines a homogeneous nucleation
  rate coefficient, ``J_{hom}``, using water activity. The water activity functions
  needed to compute ``\Delta a_w(c,T,P)`` can be found in `Common.jl` and are described in the
  `Water Activity` section. It is then used to empirically calculate ``J_{hom}(\Delta a_w)``.
  The fit below is in units of ``cm^{-3}s^{-1}``;
  the implementation converts the result to ``m^{-3}s^{-1}``.

The nucleation rate coefficient is determined with the cubic function from [Koop2000](@cite)
```math
\begin{equation}
  logJ_{hom} = -906.7 + 8502 \Delta a_w - 26924(\Delta a_w)^2 + 29180(\Delta a_w)^3
\end{equation}
```
This parameterization is valid only when ``0.26 < \Delta a_w < 0.34`` and ``185 \, K < T < 235 \, K``.
The implementation throws a domain error when ``\Delta a_w`` is outside this range.

### Homogeneous Freezing Example Figures
Here is a comparison of our parameterization of ``J_{hom}`` compared to Koop 2000 as
  plotted in figure 1 of [Spichtinger2023](@cite). Our parameterization differs in the calculation
  of ``\Delta a_w``. We define water activity to be a ratio of saturated vapor pressures whereas
  Koop 2000 uses the difference in chemical potential.

```@example
include("plots/HomFreezingPlots.jl")
```
![](HomFreezingPlots.svg)

It should be noted that the Koop 2000
  parameterization is only valid for temperatures up to 240 K and a temperature-dependent max
  pressure. The max valid pressure becomes negative around 237 K, so the Koop 2000 parameterization
  should not be valid beyond 237 K. For this reason, we limit the curve from [Spichtinger2023](@cite)
  to 237 K.

Multiple sulphuric acid concentrations, ``x``,
  are plotted since the actual concentration used in literature values is unspecified.

!!! note

    Spichtinger plot may be under the condition that x = 0 (pure liquid droplets).
    The current parameterization in CloudMicrophysics.jl is not valid for \Delta a_w
    values that are obtained from pure water droplets. Though CliMA lines look far
    from the Spichtinger 2023 line, the lines seem to move closer as x approaches 0.

### Linear fit homogeneous freezing
To avoid complications in extrapolating the Koop 2000 parameterization outside of the valide temperature range, a linear fit was created.

```@example
include("plots/linear_HOM_J.jl")
```
![](linear_HOM_J.svg)

## Deposition Nucleation in the Two-Moment P3 Tendencies
The two-moment P3 tendencies take deposition nucleation from two parameter choices of [`CloudMicrophysics.Parameters.P3IceParams`](@ref).
The target spectrum, a subtype of [`CloudMicrophysics.Parameters.INPTargetSpectrum`](@ref) selected by the `ice_nucleation` keyword, sets the target concentration of ice nucleating particles (INPs), ``N_t(T)``, the conditions under which nucleation is active, and the rate ``k`` at which nucleation relaxes toward the target.
The depletion model, selected by the `inp_depletion_model` keyword, sets the specific concentration of activated INPs, ``n_{act}``.
The default depletion model, [`CloudMicrophysics.Parameters.NIceProxyDepletion`](@ref), takes ``n_{act} = n_{ice}``, the specific ice number of the cell.
[`CloudMicrophysics.HetIceNucleation.deposition_rate`](@ref) creates crystals at the rates
```math
\begin{equation}
  \partial_t n_{frz} = k \max\left(\frac{N_t(T)}{\rho} - n_{act}, 0\right), \qquad
  \partial_t q_{frz} = m_{nuc} \, \partial_t n_{frz},
\end{equation}
```
where ``m_{nuc}`` is the mass of a newly nucleated crystal, a solid ice sphere of the P3 diameter ``D_{nuc}``.

### Frostenberg 2023 spectrum
The default target spectrum, [`CloudMicrophysics.Parameters.Frostenberg2023`](@ref), takes the median of the climatology of [Frostenberg2023](@cite) as the target, ``N_t(T) = \exp(\mu(T))``, with ``\mu`` the mean log concentration of the [INP Concentration Frequency](@ref) section.
Nucleation relaxes toward it on the timescale ``\tau`` of 300 s, ``k = 1 / \tau``, below ``T_{thr} = 258.15`` K at an ice supersaturation of at least ``S_{thr} = 0.05``.
The mass rate is at most half the vapor excess over ice saturation per relaxation time, ``\partial_t q_{frz} \leq \max(q_v - q_{sat,ice}, 0) / (2 \tau)``.
The same relaxation bounds the number rate of cloud droplet immersion freezing, through [`CloudMicrophysics.HetIceNucleation.immersion_limit_rate`](@ref).

### Cooper spectrum
[`CloudMicrophysics.Parameters.ExponentialSupercoolingINP`](@ref) is the Cooper (1986) spectrum in the form of [Thompson2004](@cite), as in Appendix C of [MorrisonMilbrandt2015](@cite),
```math
\begin{equation}
  N_t(T) = \min\left(a \exp(b (T_0 - T)), N_{max}\right),
\end{equation}
```
active below ``T_{thr} = 258.15`` K at ``S_i \geq 0.05``, with the ceiling ``N_{max} = 10^5`` m⁻³ of the P3 reference implementation.
Nucleation relaxes toward the target at the rate at which vapor deposition grows a newly nucleated crystal to its radius ``r_{nuc} = D_{nuc} / 2``.
A crystal of radius ``r`` grows as ``r \, \partial_t r = G_{ice}(T) S_i / \rho_i``, with ``G_{ice}`` from [`CloudMicrophysics.Common.G_func_ice`](@ref), so
```math
\begin{equation}
  k = \frac{2 G_{ice}(T) \max(S_i, 0)}{\rho_i r_{nuc}^2}.
\end{equation}
```
This spectrum has no INP budget for immersion freezing, as in the P3 reference implementation.

### Stochastic spectrum
[`CloudMicrophysics.Parameters.StochasticINPSpectrum`](@ref) shifts the log concentration of a base spectrum by ``s``, with a ceiling ``N_{max}`` and a lowest temperature ``T_{min}``,
```math
\begin{equation}
  N_t(T) = \min\left(\exp\left(\log N_{t,base}(T) + s\right), N_{max}\right).
\end{equation}
```
Nucleation is active where the base spectrum is active and ``T \geq T_{min}``, at the relaxation rate of the base spectrum.
On a Frostenberg 2023 base, the INP budget of immersion freezing uses the shifted target and is zero below ``T_{min}``.

### New spectra and depletion models
A new target spectrum is a callable that returns ``N_t(T)``, with methods of the two functions below.
It may define an INP budget for immersion freezing, a method of [`CloudMicrophysics.HetIceNucleation.immersion_limit_rate`](@ref).
A new depletion model defines a method of [`CloudMicrophysics.HetIceNucleation.n_active`](@ref).
```@docs
CloudMicrophysics.HetIceNucleation.is_active
CloudMicrophysics.HetIceNucleation.delivery_rate
```

## Water Activity Based vs P3 Ice Nucleation Parameterizations
The water activity based models are compared with P3 ice nucleation parameterizations as described
  in [MorrisonMilbrandt2015](@cite) using an adiabatic parcel model with depositional growth. The
  first row of plots shows deposition nucleation; the second row shows heterogeneous freezing; the
  last row shows homogeneous freezing.

```@example
include("../../parcel/Example_P3_vs_activitybased.jl")
```
![](P3_vs_activitybased.svg)

For the activity based deposition ice nucleation model (ABDINM):
 - Green corresponds to feldspar,
 - Orange corresponds to ferrihydrite,
 - Magenta corresponds to kaolinite.
We note that the P3 deposition parameterization shows a significantly greater ICNC at all times in
  the simulation. This is expected because the parameterization does not account for the number of
  available INP.

For the activity based immersion freezing model (ABIFM):
 - Green corresponds to desert dust,
 - Orange corresponds to illite,
 - Magenta corresponds to kaolinite.
The P3 heterogeneous parameterization also shows much more ICNC than any of the ABIFM runs. It
  should be noted that the P3 parameterization does not distinguish between immersion freezing,
  contact freezing, etc.

The P3 scheme allows homogeneous freezing to freeze all droplets at temperatures equal to or less than 233.15 K. No homogeneous freezing occurs at warmer temperatures.

## Freezing of Supercooled Drops in the Two-Moment P3 Tendencies
The two-moment P3 tendencies freeze cloud droplets and raindrops at any supercooling, through [`CloudMicrophysics.HetIceNucleation.cloud_freezing_rate`](@ref) and [`CloudMicrophysics.HetIceNucleation.rain_freezing_rate`](@ref).
A drop of diameter ``D`` freezes in three sequential stages: an ice embryo nucleates, the drop sheds the latent heat of fusion, and a dendritic freezing front crosses it.
The freezing rate of the drop is the reciprocal of the sum of the three timescales,
```math
\begin{equation}
  r(D) = \frac{1}{\tau_{nuc}(D) + \tau_{heat}(D) + \tau_{dend}(D)},
\end{equation}
```
and the number and mass freezing rates of a population with size distribution ``n(D)`` are
```math
\begin{equation}
  \partial_t n_{frz} = \int_0^\infty r(D) \, n(D) \, dD, \qquad
  \partial_t q_{frz} = \int_0^\infty x(D) \, r(D) \, n(D) \, dD,
\end{equation}
```
with ``x(D) = \pi \rho_w D^3 / 6`` the drop mass.
The two-moment rain distribution has a fixed shape and cannot lose its small drops alone, for example where heat loss limits freezing and the small drops freeze first.
The rain number therefore freezes at most at the fraction of the mass that freezes, ``\partial_t n_{frz} \leq (n / q) \, \partial_t q_{frz}``, so that freezing does not raise the mean drop mass.
Where nucleation limits freezing, the large drops freeze first, and the number rate is the integral above.

### Nucleation
Immersion freezing ([Bigg1953](@cite)) and homogeneous freezing ([Koop2000](@cite)) act in parallel on the same drop, so their volumetric rate coefficients add,
```math
\begin{equation}
  \frac{1}{\tau_{nuc}(D)} = \left( J_{Bigg}(T) + J_{Koop}(T) \right) \frac{\pi}{6} D^3.
\end{equation}
```
``J_{Koop}`` is evaluated for pure water, with ``\Delta a_w = 1 - e_{sat,ice}(T) / e_{sat,liq}(T)``.
It is zero where ``\Delta a_w`` is below the fitted range, and it is held at the value at the cold edge of the range where ``\Delta a_w`` is above it.
For cloud droplets, the INP budget of a target spectrum, [`CloudMicrophysics.HetIceNucleation.immersion_limit_rate`](@ref), bounds ``J_{Bigg}``.
The immersion number rate at the coefficient ``J`` is ``J q_{lcl} / \rho_w``, so the bound is ``J_{Bigg} \leq \rho_w \partial_t n_{lim} / q_{lcl}``.
Homogeneous freezing has no INP budget.

### Heat loss
The drop first warms to ``T_{freeze}`` by freezing the fraction ``c_w \Delta T / L_f`` of its mass, with ``\Delta T = T_{freeze} - T``.
The remaining latent heat leaves through the ventilated surface by conduction and evaporation, following Eq. (A7) of [Musil1970](@cite):
```math
\begin{equation}
  \tau_{heat}(D) = \frac{x(D) \max(L_f - c_w \Delta T, 0)}
    {2 \pi D F_v(D) \left( K_{therm} \Delta T + L_v D_{vapor} \Delta\rho_v \right)},
  \qquad
  \Delta\rho_v = \max(\rho_{v,sat,liq}(T_{freeze}) - \rho_a q_v, 0).
\end{equation}
```
The surface of the drop stays at ``T_{freeze}`` while it freezes, so the latent heats ``L_f`` and ``L_v`` are evaluated at ``T_{freeze}``.
``L_f - c_w \Delta T`` is then the enthalpy difference between liquid at ``T`` and ice at ``T_{freeze}``.
As in [Musil1970](@cite), one ventilation factor ``F_v``, that of vapor diffusion with the Schmidt number ``\nu_{air} / D_{vapor}``, applies to heat conduction as well.
The fall speed in the ventilation factor is that of the Stokes regime for cloud droplets and the power law of [SeifertBeheng2006](@cite) for raindrops.
``\tau_{heat}`` is zero at and above freezing, and it is zero at and below ``\Delta T^* = L_f / c_w \approx 80`` K, where the drop absorbs all of its latent heat as it warms.

### Dendrite growth
The freezing front crosses the drop at the speed ``v_{dend}`` of [`CloudMicrophysics.Parameters.DendriticFreezing`](@ref), ``\tau_{dend}(D) = D / v_{dend}``.
The default ``v_{dend} = 0.5`` m/s lies within the range of 0.1 to 1 m/s of measured dendrite tip velocities at large supercooling.
``\tau_{dend}`` bounds the freezing rate where ``\tau_{heat}`` is zero.

### Quadrature
The integrals over the size distribution use generalized Gauss-Laguerre rules, [`CloudMicrophysics.Quadrature.GaussLaguerre`](@ref), whose weight function is the normalized size distribution.
For rain, the rule integrates in ``u = D / \bar{D}_r`` over the exponential distribution of [SeifertBeheng2006](@cite).
For cloud droplets, the rule integrates in ``t = \lambda_c D^{\mu_{cD}}``, with the exponent ``\alpha = (\nu_c + 1)/\mu_c - 1`` of the generalized gamma distribution.
The rules have eight nodes by default, set by the `freezing_quadrature_order` keyword of [`CloudMicrophysics.Parameters.P3IceParams`](@ref).

## INP Concentration Frequency

With this parameterization, the concentration of ice nucleating particles (INPs) is found based on the relative frequency distribution, which depends on the temperature, but does not depend on aerosol types. It is based on [Frostenberg2023](@cite) and is derived from measurements, for marine
data sets. It is a lognormal distribution, described by

```math
\begin{equation}
  D(\mu, \sigma^2) = \frac{1}{\sqrt{2 \pi} \sigma} exp \left(- \frac{[ln(a \cdot INPC) - \mu(T)]^2}{2 \sigma^2} \right),
\end{equation}
```
where ``T`` is the temperature in degrees Celsius, ``INPC`` is the INP concentration in m``^{-3}`` and the temperature-dependent mean ``\mu(T)`` is given by
```math
\begin{equation}
  \mu(T) = ln(-(b \cdot T)^9 \times 10^{-9})
\end{equation}
```
 ``\sigma^2`` is the variance, ``a`` and ``b`` are coefficients. The parameters defined in [Frostenberg2023](@cite) for marine data sets are ``\sigma=1.37``, ``a=1`` m``^3``, ``b=1``/C.

!!! note

    Our implementation uses base SI units and takes ``T`` in Kelvin.

The following plot shows the relative frequency distribution for INPCs, as a function of temperature (the same as figure 1 in [Frostenberg2023](@cite)).

```@example
include("plots/Frostenberg_fig1.jl")
```
![](Frostenberg_fig1.svg)

The following plot shows the relative frequency distribution for INPCs at the temperature ``T=-16 C``.

![](Frostenberg_fig1_T16.svg)
