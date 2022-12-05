# Aerosol Coagulation

Aerosol coagulation describes the formation of a new particle via collision of two particles from the same or different modes.
According to [Liu2012](@cite), coagulation involving modes with sizes larger than the accumulation mode is very slow
and is therefore neglected.
Intramodal coagulation reduces number and leaves mass unchanged.
For intermodal coagulation, mass is transferred to the mode with the larger mean size.

The `Coagulation.jl` module contains analytical expressions for the coagulation integrals of the modal dynamics equations.
A full derivation of these expressions is found in Appendix H of [Whitby1991](@cite), 
but an example of 0-th moment change due to intramodal coagulation under the continuum/near-continuum regime for the ``i``-th mode is included below:
``` math
\begin{equation}
\frac{\partial }{\partial t}(M_{0_{i}}) = -\frac{1}{2}\int_{0}^{\infty}\int_{0}^{\infty}
d^{0}_{p_{1}}\beta_{nc}(d_{p_{1}},d_{p_{2}})
n_{i}(d_{p_{1}})n_{i}(d_{p_{2}})\text{d}d_{p_{1}}\text{d}d_{p_{2}}
\end{equation}
```
where ``d_p`` is the particle diameter (m), and 
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