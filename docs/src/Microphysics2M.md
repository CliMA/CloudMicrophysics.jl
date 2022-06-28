# Microphysics 2M

The `Microphysics2M.jl` module defines a 2-moment bulk parameterization of
  autoconversion and accretion. The estimates for the rain autoconversion and accretion rate functions are a result of integrating stochastic collection equations (SCE).
The module is based on the ideas of the papers utilized in 
  [Wood2005](@cite).

The cloud microphysics variables are expressed as specific humidities:
  - `q_tot` - total water specific humidity,
  - `q_liq` - cloud water specific humidity,
  - `q_rai` - rain specific humidity.

## Rain Autoconversion Functions

### Khairoutdinov and Kogan (2000)

Rain autoconversion defines the rate of conversion form cloud liquid water to 
  rain water due to collisions between cloud droplets. It is parameterized 
  following 
  [KhairoutdinovKogan2000](@cite):

``` math
\begin{equation}
  A q_{liq}^a N_d^b \rho^c
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific humidity
  - ``N_d`` is the cloud droplet concentration
  - ``\rho`` is the rain water density
  - ``A = 7.42 \times 10^{13}``
  - ``a = 2.47``
  - ``b = -1.79``
  - ``c = -1.47``

### Beheng (1994)

Rain autoconversion defines the rate of conversion form cloud liquid water to 
  rain water due to collisions between cloud droplets. It is parameterized 
  following 
  [Beheng1994](@cite):

```math
\begin{equation}
  C d^a q_{liq}^b N_d^c
\end{equation}
```

where:
  - ``q_{liq}`` is the cloud liquid water specific humidity
  - ``N_d`` is the cloud droplet concentration
  - ``C = 3.0 \times 10^{34}``
  - ``d = 9.9`` for ``N_d < 200 cm^{-3}``
  - ``d = 3.9`` for ``N_d > 200 cm ^{-3}``
  - ``a = -1.7``
  - ``b = 4.7``
  - ``c = -3.3``

### Tripoli and Cotton (1980)

Rain autoconversion defines the rate of conversion form cloud liquid water to 
  rain water due to collisions between cloud droplets. It is parameterized 
  following 
  [TripoliCotton1980](@cite):

```math
\begin{equation}
  D q_{liq}^a N_d^b H(q_{liq} - q_{liq_threshold})
\end{equation}
```

where:
  - ``q_{liq}`` is cloud liquid water specific humidity
  - ``N_d`` is the cloud droplet concentration
  - ``H(x)`` is the Heaviside step function
  - ``q_{liq_threshold}`` is cloud liquid to rain water autoconversion threshold (``\frac{4}{3} \pi \rho_w N_d r_{cm}^3``)
  - ``r_{cm}^3 = 7 \times 10^{-6} m``
  - ``D = 3268``
  - ``a = \frac{7}{3}``
  - ``b = \frac{-1}{3}``

### Liu and Daum (2004)

Rain autoconversion defines the rate of conversion form cloud liquid water to 
  rain water due to collisions between cloud droplets. It is parameterized 
  following 
  [LiuDaum2004](@cite):

```math
\begin{equation}
  E q_{liq}^a N_d^b H(R_6 - R_{6C})
\end{equation}
```

where:
  - ``q_{liq}`` is cloud liquid water specific humiditiy
  - ``N_d`` is the cloud droplet concentration
  - ``E = 1.08 \times 10^{10} \beta_6^6``
  - ``R_6 = \beta_6 r_\nu``
  - ``\beta_6 = \left( \frac{r_\nu + 3}{r_\nu} \right)^{1/3}``
  - ``r_\nu`` is mean volume radius
  - ``R_{6C} = \frac{2.25}{q_{liq}^{1/6} R_6^{1/2}}``

According to [Wood2005](@cite), the sixth moment of the cloud droplet size
distribution is given by

```math
R_6 = \left{ \left[ \int_0^{r_0} r^6 n(r) \, dr \right]/N_d \right}^{1/6}
```

In this package, we defined the mean volume radius as defined in the 1-moment
parameterization module, 

```math
m(r) = \chi_m \, m_0 \left(\frac{r}{r_0}\right)^{m_e + \Delta_m}
```

where the ``\lambda`` parameter is defined as
```math
\lambda =
  \left(
    \frac{ \Gamma(m_e + \Delta_m + 1) \, \chi_m \, m_0 \, n_0}
         {q \, \rho \, (r_0)^{m_e + \Delta_m}}
  \right)^{\frac{1}{m_e + \Delta_m + 1}}
```
where:
 - ``q`` is rain, snow or ice specific humidity
 - ``\chi_m``, ``m_0``, ``m_e``, ``\Delta_m``, ``r_0``, and ``n_0``
   are the corresponding mass(radius) and size distribution parameters
 - ``\Gamma()`` is the gamma function.

## Rain Accretion Functions

### Khairoutdinov and Kogan (2000)
Accretion defines the rates of conversion between different categories due to
  collisions between particles. It is parameterized following 
  [KhairoutdinovKogan2000](@cite):

```math
\begin{equation}
  A (q_{liq} q_{rai})^a \rho^b
\end{equation}
```

where:
  - ``q_{liq}`` is cloud liquid water specific humiditiy
  - ``q_{rai}`` is rain specific humidity
  - ``A = 67``
  - ``a = 1.15``
  - ``b = -1.3``

### Beheng (1994)

Accretion defines the rates of conversion between different categories due to
  collisions between particles. It is parameterized following 
  [Beheng1994](@cite):

```math
\begin{equation}
  B q_{liq} q_{rai}
\end{equation}
```

where:
  - ``q_{liq}`` is cloud liquid water specific humiditiy
  - ``q_{rai}`` is rain specific humidity
  - ``B = 6.0``

### Tripoli and Cotton (1980)

Accretion defines the rates of conversion between different categories due to
  collisions between particles. It is parameterized following 
  [TripoliCotton1980](@cite):

```math
\begin{equation}
  C q_{liq} q_{rai}
\end{equation}
```

where:
  - ``q_{liq}`` is cloud liquid water specific humiditiy
  - ``q_{rai}`` is rain specific humidity
  - ``C = 4.7``

## Example Figures

```@example example_figures
import Plots

import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

const PL = Plots
const CMT = CloudMicrophysics.CommonTypes
const CM1 = CloudMicrophysics.Microphysics1M
const CM2 = CloudMicrophysics.Microphysics2M
const TD = Thermodynamics
const CP = CLIMAParameters
const CP_planet = CLIMAParameters.Planet

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

const liquid = CMT.LiquidType()
const ice = CMT.IceType()
const rain = CMT.RainType()
const snow = CMT.SnowType()

# example values
q_liq_range = range(1e-8, stop=1e-3, length=100)
q_rai_range = range(1e-8, stop=1e-3, length=100)
N_d_range = range(1e7, stop=1e9, length=100)
q_liq = 5e-4
q_rai = 5e-4
ρ_air = 1.2 # kg m^-3

PL.plot( q_liq_range * 1e3, [CM2.conv_q_liq_to_q_rai_KK2000(param_set, q_liq, 1.0, 1e8) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="KK2000")
PL.plot!( q_liq_range * 1e3, [CM2.conv_q_liq_to_q_rai_B1994(param_set, q_liq, 1e8) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="B1994")
PL.plot!( q_liq_range * 1e3, [CM2.conv_q_liq_to_q_rai_TC1980(param_set, q_liq, 1e8) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="TC1980")
PL.plot!( q_liq_range * 1e3, [CM2.conv_q_liq_to_q_rai_LD2004(param_set, q_liq, 1e8) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="LD2004")
PL.plot!( q_liq_range * 1e3, [CM1.conv_q_liq_to_q_rai(param_set, q_liq) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="K1969")
PL.savefig("acnv_rate_q_liq.svg") # hide

PL.plot( N_d_range * 1e3, [CM2.conv_q_liq_to_q_rai_KK2000(param_set, q_liq, 1.0, 1e8) for N_d in N_d_range],     linewidth=3, xlabel="N_d [cm^-3]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="KK2000")
PL.plot!( N_d_range * 1e3, [CM2.conv_q_liq_to_q_rai_B1994(param_set, q_liq, 1.0, 1e8) for N_d in N_d_range],     linewidth=3, xlabel="N_d [cm^-3]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="B1994")
PL.plot!( N_d_range * 1e3, [CM2.conv_q_liq_to_q_rai_TC1980(param_set, q_liq, 1.0, 1e8) for N_d in N_d_range],     linewidth=3, xlabel="N_d [cm^-3]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="TC1980")
PL.plot!( N_d_range * 1e3, [CM2.conv_q_liq_to_q_rai_LD2004(param_set, q_liq, 1.0, 1e8) for N_d in N_d_range],     linewidth=3, xlabel="N_d [cm^-3]", ylabel="autoconversion rate [kg m^-3 s^-1]", label="LD2004")
PL.savefig("acnv_rate_N_d.svg") # hide

PL.plot( q_liq_range * 1e3, [CM2.accretion_KK2000(param_set, q_liq, q_rai, 1.0) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg] (q_rai = 0.5e-4)", ylabel="accretion rate [?]", label="KK2000")
PL.plot!( q_liq_range * 1e3, [CM2.accretion_B1994(param_set, q_liq, q_rai, 1.0) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg] (q_rai = 0.5e-4)", ylabel="accretion rate [?]", label="B1994")
PL.plot!( q_liq_range * 1e3, [CM2.accretion_TC1980(param_set, q_liq, q_rai, 1.0) for q_liq in q_liq_range],     linewidth=3, xlabel="q_liq [g/kg] (q_rai = 0.5e-4)", ylabel="accretion rate [?]", label="TC1980")
PL.plot!( q_liq_range * 1e3, [CM1.accretion(param_set, liquid, rain, q_liq, q_rai, ρ_air) for q_rai in q_rai_range],     linewidth=3, xlabel="q_liq [g/kg] (q_rai = 0.5e-4)", ylabel="accretion rate [?]", label="K1969")
PL.savefig("acc_rate_q_liq.svg") # hide

PL.plot( q_rai_range * 1e3, [CM2.accretion_KK2000(param_set, q_liq, q_rai, 1.0) for q_rai in q_rai_range],     linewidth=3, xlabel="q_rai [g/kg] (q_liq = 0.5e-4)", ylabel="accretion rate [?]", label="KK2000")
PL.plot!( q_rai_range * 1e3, [CM2.accretion_B1994(param_set, q_liq, q_rai, 1.0) for q_rai in q_rai_range],     linewidth=3, xlabel="q_rai [g/kg] (q_liq = 0.5e-4)", ylabel="accretion rate [?]", label="B1994")
PL.plot!( q_rai_range * 1e3, [CM2.accretion_TC1980(param_set, q_liq, q_rai, 1.0) for q_rai in q_rai_range],     linewidth=3, xlabel="q_rai [g/kg] (q_liq = 0.5e-4)", ylabel="accretion rate [?]", label="TC1980")
PL.plot!( q_rai_range * 1e3, [CM1.accretion(param_set, liquid, rain, q_liq, q_rai, ρ_air) for q_rai in q_rai_range],     linewidth=3, xlabel="q_rai [g/kg] (q_liq = 0.5e-4)", ylabel="accretion rate [?]", label="K1969")
PL.savefig("acc_rate_q_rai.svg") # hide

```

![](acnv_rate_q_liq.svg)
![](acnv_rate_N_d.svg)
![](acc_rate_q_liq.svg)
![](acc_rate_q_rai.svg)
