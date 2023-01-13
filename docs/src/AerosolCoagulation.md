# Aerosol Coagulation

Aerosol coagulation describes the formation of a new particle via collision of two particles from the same or different modes.
In MAM3 (and our implementation), coagulation rates and effects are computed by moment and mode. 
According to [Liu2012](@cite), coagulation involving modes with sizes larger than the accumulation mode is very slow
and is therefore neglected.
Intramodal coagulation reduces number and leaves mass unchanged.
For intermodal coagulation, mass is transferred to the mode with the larger mean size. The number concentration for the smaller mode decreases, while the larger mode stays constant.
The change to a given mode and moment is given by the following integral over the distribution sizes of the two given modes:
``` math
\frac{\partial }{\partial t}(M_{k_{i}}) = \Gamma_T \int_{0}^{\infty}\int_{0}^{\infty}
f_k(d_{p_1}, d_{p_2})\gamma(d_{p_1}, d_{p_2})n_i(d_{p_1})n_j(d_{p_2})
d(\text{ln}d_{p_1})d(\text{ln}d_{p_2})
```
where 
 - ``\Gamma_T`` is the non size-dependent term
 - ``f_k`` is the moment factor (``d^p``, where ``p`` is the moment)
 - ``\gamma`` is the size-dependent coagulation rate, equivalent to ``\beta`` below
 - ``n_i`` and ``n_j`` are the number concentrations of the given modes (``i``,``j``) at the given diameters (``d_{p_1}``, ``d_{p_2}``)

#### Intramodal Coagulation
The change to the 0-th moment due to intramodal coagulation is as follows:
``` math
\begin{equation}
\frac{\partial }{\partial t}(M_{0_{i}}) = -\frac{1}{2}\int_{0}^{\infty}\int_{0}^{\infty}
d^{0}_{p_{1}}\beta_(d_{p_{1}},d_{p_{2}})
n_{i}(d_{p_{1}})n_{i}(d_{p_{2}})\text{d}d_{p_{1}}\text{d}d_{p_{2}}
\end{equation}
```
The change to the 2nd moment due to intramodal coagulation is given by:
``` math
\begin{equation}
\frac{\partial }{\partial t}(M_{2_{i}}) = \int_{0}^{\infty}\int_{0}^{\infty}
d^{0}_{p_{1}}\beta_(d_{p_{1}},d_{p_{2}})
n_{i}(d_{p_{1}})n_{i}(d_{p_{2}})\text{d}d_{p_{1}}\text{d}d_{p_{2}}
\end{equation}
```
Since intramodal coagulation conserves mass, the 3rd moment change is zero.
#### Intermodal coagulation
Intermodal coagulation integrals take the following form for the ``k``-th moment:
``` math
\begin{equation}
\frac{\partial }{\partial t}(M_{k_{i}}) = 
-\int_{0}^{\infty}\int_{0}^{\infty}
d^{k}_{p_{1}}\beta(d_{p_{1}},d_{p_{2}})
n_{i}(d_{p_{1}})n_{j}(d_{p_{2}})\text{d}d_{p_{1}}\text{d}d_{p_{2}}
\end{equation}
```

The `Coagulation.jl` module contains analytical expressions for the coagulation integrals of the modal dynamics equations.
The module has two methods for approximating these integrals: quadrature, and approximation of the size-dependent coagulation rate terms.

## Quadrature
Quadrature is computed using the [Cubature.jl](https://github.com/JuliaMath/Cubature.jl) package.
Following Binkowski and Shankar, 1995, ``\beta(d_{p_{1}},d_{p_{2}})`` has two forms, dependent on the Knudsen regime. 
For the free-molecule regime,
``` math
\begin{equation}
\beta_{fm}(d_{p_{1}},d_{p_{2}})=
\sqrt{\frac{3k_BT}{\rho_p}}
\sqrt{\frac{1}{d_{p_1}^3}+\frac{1}{d_{p_2}^3}}
(d_{p_1}+d_{p_2})^2
\end{equation}
```
where
 - ``k_B`` is the Boltzmann constant
 - ``T`` is temperature (K)
 - ``\rho_p`` is the particle density for the given mode.

For the near-continuum regime, 
``` math
\begin{equation}
\beta_{nc}(d_{p_{1}},d_{p_{2}})=
\frac{2k_BT}{3\mu}(d_{p_1} + d_{p_2})

[\frac{1}{d_{p_1}}+\frac{1}{d_{p_2}}+
2.492\lambda(\frac{1}{d_{p_1}^2}+\frac{1}{d_{p_2}^2})
]
\end{equation}
```

where
 - ``\mu`` is the gas viscosity
 - ``\lambda`` is the mean free path
## Substitution-based Approximation
For the substitution-based approximation, a full derivation of these expressions is found in Appendix H of [Whitby1991](@cite), 
but an example of 0-th moment change due to intramodal coagulation for the ``i``-th mode is included below:
``` math
\begin{equation}
\frac{\partial }{\partial t}(M_{0_{i}}) = -\frac{1}{2}\int_{0}^{\infty}\int_{0}^{\infty}
d^{0}_{p_{1}}\beta(d_{p_{1}},d_{p_{2}})
n_{i}(d_{p_{1}})n_{i}(d_{p_{2}})\text{d}d_{p_{1}}\text{d}d_{p_{2}}
\end{equation}
```
where ``d_p`` is the particle diameter (m). ``beta_{fm}(d_{p_{1}},d_{p_{2}})`` is the same as in the Quadrature section, and
``\beta_{nc}(d_{p_{1}},d_{p_{2}})``, the coagulation coefficient ``(m^3/s)`` is given by:
``` math
\begin{equation}
\beta_{nc}(d_{p_{1}},d_{p_{2}}) = 2\pi (D_{p_1} + D_{p_2})(d_{p_1} + d_{p_2})
\end{equation}
```
where  ``D_p`` is the particle diffusion coefficient ``(m^2/s)``:
``` math
D_p = [k_B\cdot T\cdot C_C /(3\pi\mu d_p)]
\\
C_C \approx 1 + A \cdot Kn
\\
A = 1.392Kn_{g}^{0.0783}
\\
Kn_g = 2\lambda/D_{gn}
```
where: 
 - ``Kn`` is the Knudsen number for air at the given temperature and pressure (m)
 - ``\lambda`` is the mean free path ``(m)``
 - ``D_{gn}`` is the average diameter of the ``i``-th mode ``(m)``
Substituting ``(2)`` into ``(1)`` yields:
``` math
\frac{1}{2}\int_{0}^{\infty}\int_{0}^{\infty}
d^{0}_{p_{1}}\beta_{nc}(d_{p_{1}},d_{p_{2}})
n_{i}(d_{p_{1}})n_{i}(d_{p_{2}})\text{d}d_{p_{1}}\text{d}d_{p_{2}}
=
N_{i}^{2}K_{nc} [1 + \text{ESG}_i^8 + A_i Kn_{g_i}(\text{ESG}_i^{20}+\text{ESG}_i^4)]
```
where:
 - ``N_i`` is the number concentration of the ``i``-th mode ``(1/m^3)``
 - ``K_{nc} = \sqrt{2k_B T / 3\mu}`` is the non size-dependent term for the continuum/near-continuum coagulation coefficient ``(m^3/s)``, and ``\mu`` is the gas viscosity ``[kg/(m\cdot s)]``
 - ``ESG_i = exp(1/8 \cdot ln^2 \sigma_g)``, where ``\sigma_g`` is the geometric standard deviation of the ``i``-th mode