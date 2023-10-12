# 0-dimensional Box Model

This model solves for number concentration of ice crystals
    produced via water-activity based immersion freezing
    in a 0-dimensional setup with a prescribed constant cooling rate.
It is used to validate the immersion freezing parameterization
  and is based on Figure 4 from [Alpert2016](@cite).

Note that because the cooling rate is prescribed, the temperature change does not
   depend on the freezing rate and the associated latent heat release.
The production rate of ice is given by
```math
\begin{equation}
\frac{dN_{i}}{dt} = J \; N_{l} \; A_{d}
\end{equation}
```
where
- ``J`` - is the immersion freezing rate coefficient,
- ``N_{l}`` - is number of liquid droplets,
- ``A_{d}`` - is the surface area available for freezing for a single droplet.

``A_{d}`` is currently assumed to be either constant among all droplets
or sampled from a lognormal distribution with user defined mean and width.

In the plots below, we first reproduce Fig. 4 from [Alpert2016](@cite)
  by computing the apparent and actual nucleation rates based on the paper defined
  frozen fraction.
As as second step we run the box model for constant and variable ``A_d``
  and compare the results.
Note that the variable ``A_d`` scenario is closer to what is observed, while
  constant ``A_d`` is closer to what we can currently implement in our
  bulk microphysics scheme for ESM.
```@example
include("../../box/Alpert_Knopf_2016_backward.jl")
include("../../box/Alpert_Knopf_2016_forward.jl")
```
![](Alpert_Knopf_2016_backward.svg)
![](Alpert_Knopf_2016_forward.svg)
