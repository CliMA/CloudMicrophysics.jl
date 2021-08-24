# Aerosol Activation

The `AerosolActivation.jl` module contains parameterization
  of activation of aerosol particles into cloud droplets
  via deposition of water vapor.
Accompanying it, is the `AerosolDistribution.jl` module, which contains
  information about the aerosol size distribution and chemical properties.
The parameterization computes the total activated number and mass
  from a given aerosol size distribution.
It is based on [Köhler theory](https://en.wikipedia.org/wiki/K%C3%B6hler_theory)
  and assumes equilibrium thermodynamics.
The modules are an implementation of the parameterization
  from [Abdul-Razzaketal1998](@cite) and
  [Abdul-RazzakandGhan2000](@cite).


## Assumed aerosol size distribution and properties

Aerosol particles are assumed to follow a multi-mode lognormal
  size distribution.
Particles in each mode are assumed to be internally mixed.
The following table lists the parameters defining the aerosol
  physical and chemical properties.
The ``r_{dry}``, ``\sigma``, and ``N_{tot}`` are given for each mode.
Other parameters are defined for each component in each mode.

|    variable name     |         definition                                | units               |
|----------------------|---------------------------------------------------|---------------------|
|``r_{dry}``           | geometric mean dry radius                         | ``m``               |
|``\sigma``            | geometric standard deviation                      | ``-``               |
|``N_{tot}``           | total number concentration                        | ``m^{-3}``          |
|``r``                 | component mass mixing ratio                       | ``-``               |
|``\epsilon``          | mass fraction of water soluble material           | ``-``               |
|``\phi``              | osmotic coefficient                               | ``-``               |
|``M_a``               | molar mass                                        | ``kg \, mol^{-1}``  |
|``\nu``               | number of ions the salt dissociates into in water | ``-``               |
|``\rho_a``            | aerosol density                                   | ``kg \, m^{-3}``    |

!!! note

    The parameterization assumes that the solute is sufficiently soluble
    so that its concentration does not increase as the droplet grows.
    The effects of surfactants on surface tension are also not considered.


## Mean Hygroscopicity

Hygroscopicity describes the impact of solute on aerosol efficiency in taking up
  water vapor from the environment (i.e. the [Raoult's law](https://en.wikipedia.org/wiki/Raoult%27s_law)).
For a given aerosol, it is defined as in eq. (3) in [Abdul-RazzakandGhan2000](@cite):

```math
B = \frac{\nu \, \phi \, M_w \, \rho_a}{M_a \, \rho_w}
```
where:
  - ``\nu, \, \phi, \, M_a, \, \rho_a`` are the aerosol properties defined in the above table,
  - ``M_w`` is the molar mass of water,
  - ``\rho_w`` is the density of water.

The mean hygroscopicity for internally mixed mode ``i``
  made up of ``j`` aerosol species is computed
  according to eq. (4) in [Abdul-RazzakandGhan2000](@cite):

```math
\begin{equation}
\bar{B_i} = \frac{M_w \sum_{j = 1}^{J} \frac{r_{ij} \, \nu_{ij} \, \phi_{ij} \, \epsilon_{ij}}{M_{aij}}}{\rho_w \sum_{j = 1}^{J} \frac{r_{ij}}{\rho_{aij}}}
\end{equation}
```
where:
  - ``i = 1, 2, ..., I`` iterates over aerosol size distribution modes,
  - ``j = 1, 2, ..., J`` iterates over aerosol components within a given mode,
  - ``r_{ij}`` is the mass ratio of component ``j`` in mode ``i``.


## Critical supersaturation

Supersaturation ``S`` is the ratio of water vapor pressure to saturation vapor pressure.
Köhler theory defines ``S`` at which the growing aerosol particle
  is in equilibrium with the environment over a range of its sizes.
It takes into account the curvature effects and the solute effects.
Aerosol activation occurs when the threshold supersaturation,
  named critical supersaturation ``S_c``, is reached.
After reaching ``S_c``, even when ``S`` decreases,
  as long as the conditions remain saturated,
  the particle will continue to grow.
The critical supersaturation is defined by the maximum of the
  [Köhler curve](https://en.wikipedia.org/wiki/K%C3%B6hler_theory).
For example, eq.(9) in [Abdul-RazzakandGhan2000](@cite) defines the
  the critical supersaturation for a particle with
  dry radius equal to the mean mode radius:

```math
\begin{equation}
S_{ci} = \frac{2}{\sqrt{\bar{B_{i}}}} \bigg( \frac{A}{3 \, r_{dry, \, i}} \bigg)^{3/2}
\label{eq:Scriti}
\end{equation}
```
where:
  - ``A`` is the coefficient describing the curvature effects (i.e. [Kelvin effect](https://en.wikipedia.org/wiki/Kelvin_equation)),
  - ``r_{dry \, i}`` is the mean dry radius for mode ``i``.

Coefficient ``A`` is defined as in equation (5) in [Abdul-Razzaketal1998](@cite).

```math
\begin{equation}
A = \frac{2 \tau M_w}{\rho_w R T}
\end{equation}
```
where:
  - ``\tau`` is the surface tension of water,
  - ``R`` is the universal gas constant,
  - ``T`` is the temperature.

## Maximum Supersaturation

Maximum supersaturation reached by the system ``S_{max}`` governs
  what aerosol sizes are activated and what aren't.
We estimate ``S_{max}`` by considering a parcel of air raising adiabatically
  with a constant velocity, see for example [Rogers1975](@cite).
The time rate of change of ``S`` is given by eq (10) in [Rogers1975](@cite)

```math
\begin{equation}
  \frac{dS}{dt} = \alpha w - \gamma \frac{d\chi}{dt}
  \label{eq:Sevolution}
\end{equation}
```
where:
  - ``w`` is the vertical velocity,
  - ``d\chi / dt`` is the water condensation rate during aerosol activation and growth,
  - ``\alpha`` and ``\gamma`` are coefficients that do not depend on aerosol properties.

The parameters ``\alpha`` and ``\gamma`` are defined by eq. (11) and (12)
  in [Abdul-RazzakandGhan2000](@cite):
```math
\begin{equation}
\alpha = \frac{g \, M_w \, L_v}{c_p \, R \, T^2} - \frac{g \, M_{air}}{R  T}
\end{equation}
```
```math
\begin{equation}
\gamma = \frac{R T}{p_{vap}^{sat} \, M_w} + \frac{M_w \, L_v^2}{c_p \,p \, M_{air} \, T}
\end{equation}
```
where:
  - ``g`` is gravitational acceleration,
  - ``L_v`` is the latent heat of vaporization,
  - ``c_p`` is the specific heat of air,
  - ``p_{vap}^{sat}`` is the saturation vapor pressure,
  - ``p`` is the air pressure.

The maximum supersaturation is estimated from eq. (\ref{eq:Sevolution})
  assuming steady a state solution ``dS/dt = 0``.
[Abdul-Razzaketal1998](@cite) and [Abdul-RazzakandGhan2000](@cite)
  show how to derive an approximate solution for ``S_{max}``,
  since analytical solution is in general not possible.
They consider approximate solutions for very small and very large
  critical supersaturation values relative to maximum supersaturation,
  and combine them into a final expression for ``S_{max}``.
The final formula is presented in eq (6) in [Abdul-RazzakandGhan2000](@cite)
```math
\begin{equation}
S_{max} = \frac{1}{{\sum_{i=1}^{I} \frac{1}{S_{ci}^{2}} \bigg[ f_i \bigg( \frac{\zeta}{\eta_{i}} \bigg)^{\frac{3}{2}} + g_{i} \bigg(\frac{S_{ci}^{2}}{\eta_{i} + 3\zeta} \bigg)^{\frac{3}{4}}} \bigg]^{\frac{1}{2}}}
\end{equation}
```
where
  - ``S_{ci}`` is the critical supersaturation for mode ``i`` defined in eq. (\ref{eq:Scriti}),
  - ``f_i``, ``g_i``, ``\zeta``, ``\eta_i`` are the coefficients defined in eqs. (7, 8, 10 and 11) in  [Abdul-RazzakandGhan2000](@cite).

```math
\begin{equation}
f_i  = 0.5 \, \mathrm{exp} (2.5 \, \mathrm{ln}^{2} \sigma_{i})
\end{equation}
```
```math
\begin{equation}
g_i  = 1 + 0.25 \, \mathrm{ln} \sigma_i
\end{equation}
```
```math
\begin{equation}
\zeta = \frac{2A}{3} \bigg(\frac{\alpha w}{G}\bigg)^{\frac{1}{2}}
\end{equation}
```

```math
\begin{equation}
\eta_i = \bigg(\frac{\alpha w}{G}\bigg)^{\frac{3}{2}} \frac{1}{2 \pi \rho_w \gamma N_i}
\end{equation}
```
where:
 - ``G(T) = \frac{1}{\rho_w} \, \left(\frac{L_v}{KT} \left(\frac{L_v}{R_v T} - 1 \right) + \frac{R_v T}{p_{vap}^{sat} D} \right)^{-1}``
     combines the effects of thermal conductivity ``K`` and water diffusivity ``D``.


## Number and mass of activated particles

The total number ``N_{act}`` and mass ``M_{act}`` of activated aerosol particles
  can be computed by integrating their size distribution
  starting from the smallest activated size.
Following the derivations of
  [Abdul-Razzaketal1998](@cite) and [Abdul-RazzakandGhan2000](@cite)
  this can be expressed in terms of critical supersaturations of each
  size distribution mode ``S_{ci}`` and the maximum supersaturation ``S_{max}``.

```math
\begin{equation}
N_{act} = \sum_{i = 1}^{I} N_{i}\frac{1}{2}\bigg[1 - \mathrm{erf}(u_{i})\bigg]
\end{equation}
```
```math
\begin{equation}
M_{act} = \sum_{i = 1}^{I} M_{i}\frac{1}{2}\bigg[1 - \mathrm{erf}\bigg(u_{i} - \frac{3 \sqrt2}{2} ln(\sigma_i)\bigg)\bigg]
\end{equation}
```
where:

  - ``M_i`` is the average molar mass of aerosol particles in mode ``i``,
  - ``u_i`` is given in equation (15) in [Abdul-RazzakandGhan2000](@cite).

```math
\begin{equation}
u_i = \frac{2}{3\sqrt2 \, ln(\sigma_i)} ln\bigg( \frac{S_{ci}}{S_{max}} \bigg)
\end{equation}
```
where:
  - ``S_{ci}`` is the mode critical supersaturation,
  - ``S_{max}`` is the maximum supersaturation.

## Example figures

```@example example_figures
import Plots

import CloudMicrophysics
import CLIMAParameters

const PL = Plots
const AM = CloudMicrophysics.AerosolModel
const AA = CloudMicrophysics.AerosolActivation
const CP =  CLIMAParameters

struct EarthParameterSet <: CP.AbstractEarthParameterSet end
const param_set = EarthParameterSet()

# Atmospheric conditions
T = 294.0         # air temperature
p = 1000.0 *1e2   # air pressure
w = 0.5           # vertical velocity

# Abdul-Razzak and Ghan 2000 Figure 1 mode 1
# https://doi.org/10.1029/1999JD901161
r_dry_paper = 0.05 * 1e-6 # um
stdev_paper = 2.0         # -
N_1_paper = 100.0 * 1e6   # 1/m3

# Sulfate - universal parameters
M_sulfate = 0.132
ρ_sulfate = 1770.0
ϕ_sulfate = 1.0
ν_sulfate = 3.0
ϵ_sulfate = 1.0

paper_mode_1 = AM.Mode(
    r_dry_paper,
    stdev_paper,
    N_1_paper,
    (1.0,),
    (ϵ_sulfate,),
    (ϕ_sulfate,),
    (M_sulfate,),
    (ν_sulfate,),
    (ρ_sulfate,),
    1,
)

N_2_range = range(0, stop=5000 * 1e6, length=100)
N_act_frac = Vector{Float64}(undef, 100)

it = 1
for N_2_paper in N_2_range

        paper_mode_2 = AM.Mode(
          r_dry_paper,
          stdev_paper,
          N_2_paper,
          (1.0,),
          (ϵ_sulfate,),
          (ϕ_sulfate,),
          (M_sulfate,),
          (ν_sulfate,),
          (ρ_sulfate,),
          1,
        )

        AD =  AM.AerosolDistribution((paper_mode_1, paper_mode_2))
        N_act_frac[it] = AA.N_activated_per_mode(param_set, AD, T, p, w)[1] / N_1_paper

        global it += 1
end

# data read from Fig 1 in Abdul-Razzak and Ghan 2000
# using https://automeris.io/WebPlotDigitizer/
N_2_obs = [18.74716810149539, 110.41572270049846, 416.00589034889026, 918.1014952424102, 1914.816492976891, 4919.913910285455]
N_act_obs = [0.7926937018577255, 0.7161078386950611, 0.5953670140462167, 0.4850589034888989, 0.34446080652469424, 0.162630267331219]
N_2_paper_param = [54.6839601268689, 72.69483461712753, 109.48119619392855, 127.7469415496148, 155.44290892614436, 183.2238332578163,
                   220.2650657000454, 275.8269143633893, 312.86814680561884, 368.5999093792484, 433.50702310829183, 554.4007702763934,
                   628.9080199365658, 693.9850475758949, 740.5414589941097, 787.0129134571821, 852.3448119619393, 917.4218396012684,
                   973.3235160851837, 1047.915722700498, 1103.8173991844137, 1178.2396918894428, 1234.1413683733576,1280.7827367467153,
                   1364.5502945174449,1467.1782963298592,1513.734707748075, 1560.2911191662893, 1709.5604893520617, 1765.547122791119,
                   1849.4845944721342, 1989.408699592207, 2092.0367014046215, 2157.3685999093796, 2222.7004984141367, 2297.292705029451,
                   2353.364295423651, 2456.07725419121, 2549.2750339827826, 2605.26166742184, 2689.199139102854, 2810.6026280018123,
                   2997.3380154055285, 3099.8810602628, 3221.2845491617572, 3314.6522428636154, 3389.2444494789315, 3463.9216130493887,
                   3566.634571816946, 3650.657000453103, 3762.5453103760756, 3827.9621658359765, 3911.814680561849, 3995.922066153149,
                   4126.585863162664, 4247.904395106479, 4322.496601721793, 4397.258722247395, 4490.54145899411, 4574.563887630267,
                   4667.9315813321255, 4835.891481649297]
N_act_paper_param = [0.7307884005437245, 0.7016538287267784, 0.676110104213865, 0.657884005437245, 0.643271409152696, 0.6322949705482556,
                     0.6176597190756684, 0.5957068418667876, 0.5810715903942004, 0.5663910285455369, 0.5444154961486181, 0.5186678749433621,
                     0.5075781603987314, 0.49287494336202997, 0.48548935206162214, 0.4744676030811056, 0.47067285908473044, 0.45596964204802903,
                     0.4485613955595833, 0.4411078386950612, 0.43369959220661536, 0.41897371998187594, 0.4115654734934301, 0.4078160398731311,
                     0.39306751246035343, 0.3855459900317173, 0.37816039873130947, 0.37077480743090174, 0.3595038513819666, 0.3557317625736294,
                     0.3482555505210695, 0.33700724966017237, 0.32948572723153613, 0.325690983235161, 0.32189623923878585, 0.31444268237426376,
                     0.31430675124603535, 0.31042138649750806, 0.2992863615768011, 0.2955142727684641, 0.2880380607159041, 0.28410738559130055,
                     0.27638196647032176, 0.265224286361577, 0.2612936112369735, 0.25743090167648397, 0.249977344811962, 0.24615994562754873,
                     0.24227458087902143, 0.23843452650657015, 0.22725419120978696, 0.22709560489352065, 0.21598323516085194, 0.21577933846850939,
                     0.2081898504757591, 0.20062301767104684, 0.19316946080652464, 0.19298821930222032, 0.1854893520616221, 0.18164929768917082,
                     0.17778658812868153, 0.16647032170367027]


PL.plot(N_2_range * 1e-6, N_act_frac, label="CliMA", xlabel="Mode 2 aerosol number concentration [1/cm3]", ylabel="Mode 1 number fraction activated")
PL.scatter!(N_2_obs, N_act_obs, label="paper observations")
PL.plot!(N_2_paper_param, N_act_paper_param, label="paper parameterization")

PL.savefig("Abdul-Razzak_and_Ghan_fig_1.svg")
```
![](Abdul-Razzak_and_Ghan_fig_1.svg)
