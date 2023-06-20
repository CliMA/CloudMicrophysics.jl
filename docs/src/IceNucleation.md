# Ice Nucleation

The `IceNucleation.jl` module includes
  the parameterization of activation of dust aerosol particles into ice crystals
  via deposition of water vapor and water activity based immersion freezing model for sulphuric acid containing droplets.
These are heterogeneous ice nucleation processes.
The parameterization for deposition on dust particles is an implementation of
  the empirical formulae from [Mohler2006](@cite)
  and is valid for two types of dust particles:
  Arizona Test Dust and desert dust from Sahara.

!!! note

    Future work includes adding parameterizations
    for other nucleation paths such as
    heterogeneous immersion freezing or homogeneous freezing
    and modeling the competition between them.


## Activated fraction for deposition freezing on dust
The parameterization models the activated fraction
  as an empirical function of ice saturation ratio,
  see eq. (3) in [Mohler2006](@cite).
```math
\begin{equation}
f_i(S_i) = exp[a(S_i - S_0)] - 1
\end{equation}
```
where:
  - ``f_i`` is the activated fraction
      (the ratio of aerosol particles acting as ice nuclei to the total number of aerosol particles),
  - ``S_i`` is the ice saturation ratio
      (the ratio of water vapor partial pressure and the water vapor partial pressure at saturation over ice),
  - ``S_0`` is the threshold ice saturation ratio,
  - ``a`` is a scaling parameter dependent on aerosol properties and temperature.

Limited experimental values for the free parameters ``S_0`` and ``a`` can be found in [Mohler2006](@cite).
Both parameters are dependent on aerosol properties and temperature.

!!! note

    For a ``f_i`` values above 0.08 or ``S_i`` between 1.35 and 1.5,
    freezing occurs in a different ice nucleation mode
    (either a second deposition or other condensation type mode).

## ABIFM for Sulphuric Acid Containing Droplets
Water Activity-Based Immersion Freezing Model (ABFIM) is a method of parameterizing immersion freezing inspired by the time-dependent classical nucleation theory. The nucleation rate coefficient, ``J``, can be determined by the water activity, ``a_w``. This parameterization follows [KnopfAlpert2013](@cite), [Koop2002](@cite), [MurphyKoop2005](@cite), and [Luo1995](@cite). In this model, aerosols are assumed to contain an insoluble and soluble material. When immersed in water, the soluble material diffuses into the liquid water to create a sulphuric acid solution.

Using empirical coefficients, ``m`` and ``c``, from [KnopfAlpert2013](@cite), the heterogeneous nucleation rate coefficient in units of ``cm^{-2}s^{-1}`` can be determined by the linear equation
```math
\begin{equation}
  log_{10}J_{het} = m \Delta a_w + c
\end{equation}
```

``\Delta a_w``is the difference between the water activity of the droplet, ``a_w``, and the water activity of ice at the same temperature, ``a_{w,ice}(T)``. From [Koop2002](@cite), 
```math
\begin{equation}
  a_w = \frac{p_{sol}}{p_{sat}}
\end{equation}
```
```math
\begin{equation}
  a_{w,ice} = \frac{p_{i,sat}}{p_{sat}}
\end{equation}
```

where ``p`` is saturated vapor pressure of water above solution, pure liquid water, and ice for subscripts ``sol``, ``sat``, and ``i,sat`` respectively. ``p_{sol}`` is determined in mbar using a parameterization for supercooled, binary ``H_2SO_4/H_2O`` solution from [Luo1995](@cite) which is valid for 185K < T < 235K:
```math
\begin{equation}
  ln(p_{sol}) = 23.306 - 5.3465x + 12xw_h - 8.19xw_h^2 + \frac{1}{T}(-5814 + 928.9x - 1876.7xw_h)
\end{equation}
```
where ``x = \frac{wt%(H_2SO_4)}{100}``, ``w_h = 1.4408x``, and temperature is in Kelvins.

``p_{i,sat}`` in Pa follows [MurphyKoop2005](@cite) for temperatures in Kelvins above 110K:
```math
\begin{equation}
  p_{i,sat} = exp[9.550426 - \frac{1}{T}5723.265 + 3.53068 ln(T) - 0.00728332T]
\end{equation}
```

Once ``J_{het}`` is calculated, it can be used to determine the ice production rate, ``P_{ice}``, per minute via immersion freezing.
```math
\begin{equation}
  P_{ice} = J_{het}A(N_{tot}-N_{ice})
\end{equation}
```
where ``A`` is surface area of an individual ice nuclei, ``N_{tot}`` is total number of ice nuclei, and ``N_{ice}`` is number of ice crystals already in the system. 

## Homogeneous Freezing for Sulphuric Acid Containing Droplets
Based off [Koop2000](@cite), this parameterization determines a homoegneous nucleation rate coefficient, ``J_{hom}``, using water activity. First, ``a_w^i(T,P = 0)``, ``\int [v_w(c,T,P) - v^i(T,P)] dP``, and ``a_w(c,T,P = 0)`` are calculated on the assumption that the solution droplets contain an ideal solvent (water) and that ``v_w(c,T,P) \approx v_w^0(T,P)``. They are then used to calculate ``\Delta a_w(c,T,P)`` and ultimately, ``J_{hom}(\Delta a_w)`` with units of ``cm^{-3}s^{-1}``.

To find ``a_w^i(T,P = 0)``, the activity of water in solution in equilibrium with ice at ambient pressure,
```math
\begin{equation}
  a_w^i(T,P) = exp{[\mu_w^i(T,P) - \mu_w^0(T,P)]/RT}
\end{equation}
```
where:
  - ``R`` is the ideal gas constant,
  - ``T`` is temperature in Kelvins,
  - ``\mu_w^i(T,P)`` is chemical potential of pure ice,
  - ``\mu_w^0(T,P)`` is chemical potential of pure water.

Because we need water activity at ambient pressure, ``\mu_w^i(T,P) - \mu_w^0(T,P)`` can be determined empirically at ``P = 0`` through
```math
\begin{equation}
  \mu_w^i(T,0) - \mu_w^0(T,0) = 210368 + 131.438T - 3.32373 \times 10^6 T^{-1} - 41729.1ln(T)
\end{equation}
```

``\int [v_w(c,T,P) - v^i(T,P)] dP`` can be approximated as
```math
\begin{equation}
  \int [v_w(c,T,P) - v^i(T,P)] dP \approx v_w^0(T,0){P - \frac{1}{2}\kappa_T^0(T,P)P^2 - \frac{1}{6}[\frac{\delta \kappa_T^0(T,P)}{\delta P}]P^3} - v^i(T,0){P - \frac{1}{2}\kappa_T^i(T,P)P^2 - \frac{1}{6}[\frac{\delta \kappa_T^i(T,P)}{\delta P}]P^3}
\end{equation}
```
where:
  - ``v_w`` is partial molar volume of water in solution with units of ``cm^3mol^{-1}``,
  - ``v^i`` is molar volume of hexagonal ice,
  - ``v_w^0`` is molar volume of pure liquid water,
  - ``\kappa_T^0`` is isothermal compressibility of pure water
      (at ambient pressure (P=0), this is taken to be 1.6 GPa^{-1}),
  - ``\kappa_T^i`` is isothermal compressibility of pure ice
      (at ambient pressure (P=0), this is taken to be 0.22 GPa^{-1}),
  - ``\frac{\delta \kappa_T^0}{\delta P}`` is pressure dependence of ``\kappa_T^0``
      (this is taken to be -8.8 GPa^{-2}),
  - ``\frac{\delta \kappa_T^i}{\delta P}`` is pressure dependence of ``\kappa_T^i``
      (this is taken to be -0.17 GPa^{-2}).
Both ``v_w^0(T,0)`` and ``v^i(T,0)`` are computed empirically using
```math
\begin{equation}
  v_w^0(T,0) = -230.76 - 0.1478T + 4099.2T^{-1} + 48.8341ln(T)
\end{equation}
```
```math
\begin{equation}
  v^i(T,0) = 19.43 - 2.2\times 10^{-3}T + 1.08 \times 10^{-5}T^2
\end{equation}
```

``a_w(c,T,P = 0)`` can be found from the equation for chemical potential of solution.
```math
\begin{equation}
  RTln(a_w(c,T,P)) = \mu_w(c,T,0) + \int v_w(c,T,P)dP - \mu_w^0(T,0) - int v_w^0(T,P)dP
\end{equation}
```
!!! note

    Requires parameterization of \mu_w(c,T,0) - \mu_w^0(T,0)

``\Delta a_w(C,T,P)`` can now be computed using
```math
\begin{equation}
  \Delta a_w(C,T,P) = a_w(c,T,0)exp \{ \frac{1}{RT} \int [v_w(c,T,P) - v^i(T, P)]dP \} - a_w^i(T,0)
\end{equation}
```
And finally, the nucleation rate coefficient is determined with a cubic function
```math
\begin{equation}
  logJ_{hom} = -906.7 + 8502 \Delta a_w - 26924(\Delta a_w)^2 + 29180(\Delta a_w)^3
\end{equation}
```
This parameterization is valid for ``0.26 < \Delta a_w < 0.36``.