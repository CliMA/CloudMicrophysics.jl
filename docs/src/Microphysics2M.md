# Microphysics 2M

The `Microphysics2M.jl` module provides 2-moment warm rain bulk parameterization of cloud microphysical processes including autoconversion, accretion, cloud droplets and raindrops self-collection, raindrops breakup, mean terminal velocity and rain evaporation. Autoconversion defines the rate of transfer from cloud liquid water to rain water due to collisions between cloud droplets. Accretion defines the rate of transfer from cloud liquid water to rain water due to collisions between cloud droplets and rain drops. Cloud self-collection defines the rate of change of cloud droplets number density due to collisions between cloud droplets, and similarly, rain self-collection defines the rate of change of raindrops number density due to collsions between raindrops. Rain drops breakup defines the rate of new raindrops production due to the disintegration of larger raindrops. Mean terminal velocity represents the mean fall speed of raindrops, and rain evaporation describes the rate of the transformation of rainwater to water vapor. Specifically, `Microphysics2M.jl` implements:
 - the double-moment [SeifertBeheng2006](@cite) parametrization, which includes autoconversion, accretion, cloud and rain self-collection rates, breakup, terminal velocity and evaporation;
 - and other double-moment autoconversion and accretions schemes from [Wood2005](@cite) based on the works of [KhairoutdinovKogan2000](@cite), [Beheng1994](@cite), [TripoliCotton1980](@cite) and [LiuDaum2004](@cite).

The microphysics variables are expressed as specific contents [kg substance / kg air] and number densities [1 / m$^3$ air]:
  - `q_liq` - cloud water specific content,
  - `q_rai` - rain specific content,
  - `N_liq` - cloud droplets number density,
  - `N_rai` - raindrops number density.
The default values of free parameters are defined in
  [ClimaParams](https://github.com/CliMA/ClimaParams.jl)
  and can be overwritten using the `toml` files.

## The Seifert and Beheng (2006) parametrization

The [SeifertBeheng2006](@cite) parametrization provides process rates for autoconversion, accretion, self-collection of cloud droplets and raindrops, raindrops breakup, raindrops mean fall speed, and rain evaporation. This parametrization is directly derived from the stochastic collection equation (SCE) with a piecewise polynomial collsion kernel. 
It assumes
- a generalized gamma distribution for cloud droplets (as a function of mass)
- an exponential distribution for raindrops (as a function of diameter)

The piece-wise polynomial collection kernel, used for the derivation of the parametrization, is given by:
```math
\begin{align}
    K(x,y) =
    \begin{cases}
    k_{cc}(x^2+y^2), \quad & x\wedge y < x^*\\
    k_{cr}(x+y), \quad & x\oplus y \geq x^*,\\
    k_{rr}(x+y)\exp[-\kappa_{rr} (x^{1/3} + y^{1/3})], \quad & x\wedge y \geq x^*,
    \end{cases}
\end{align}
```
where ``x`` [kg] and ``y`` [kg] are drop masses and ``x^*`` [kg] is the mass threshold chosen to separate the cloud and rain portions of the mass distribution. ``K`` has units of m$^3$ s$^{-1}$, and the constants are:

|   symbol      | default value                                  |
|---------------|------------------------------------------------|
|``k_{cc}``     | ``4.44 × 10^9`` m$^3$ kg$^{-2}$ s$^{-1}$       |
|``k_{cr}``     | ``5.25`` m$^3$ kg$^{-1}$ s$^{-1}$              |
|``k_{rr}``     | ``7.12`` m$^3$ kg$^{-1}$ s$^{-1}$              |
|``κ_{rr}``     | ``60.7`` kg$^{-1/3}$                           |
|``x^*``        | ``6.54 × 10^{-11}`` kg                         |

Assuming spherical raindrops, the mass ``x`` is related to the drop radius ``r`` by ``x = \frac{4π}{3}ρ_w r^3``, then the default value of ``x^*=6.54\times 10^{-11}`` kg corresponds to the drop radius ``r^* ≈ 25`` μm.

### Assumed size distributions

In this section, we describe the assumed size distributions for cloud droplets and raindrops. We use the symbol ``f`` to denote the number distribution as a function of mass, ``x`` [kg], and the symbol ``n`` to denote the number distribution as a function of diameter, ``D`` [m]. We will show how to convert between the two representations. Symbols related to cloud droplets and raindrops are denoted with the subscript ``c`` and ``r``, respectively (e.g. ``f_c(x)`` and ``f_r(D)``).

#### Cloud droplets

The cloud droplet number distribution, as a function of mass ``x`` [kg], is assumed to follow a generalized Gamma distribution[^1]

[^1]: c.f. Eq. (79) in [SeifertBeheng2006](@cite), but using the symbol `B` instead of `λ`

!!! todo "Form of size distribution parameters"
    It would simplify math (and intuition?) if we instead defined the size distribution parameters in the form:

    ```math
    f(x) = \frac{N}{x̄} ⋅ \left(\frac{x}{x̄}\right)^ν ⋅ A ⋅ \exp\left(- \left(B\frac{x}{x̄}\right)^μ\right)
    ```
    where
    ```math
    A = \frac{μ ⋅ B^{ν+1}}{Γ(\frac{ν + 1}{μ})}, \qquad
    B = \frac{Γ(\frac{ν + 2}{μ})}{Γ(\frac{ν + 1}{μ})}
    ```

    At least, it would imply that ``x`` is always normalized by ``x̄``, and that ``A`` and ``B`` are non-dimensional.

```math
\begin{align}
    f_c(x) = A_c x^{ν_c} \exp\left(-B_c x^{μ_c} \right), 
\end{align}
```
We assume that ``ν_c`` and ``μ_c`` are fixed constants. Our default choice for these parameters are ``ν_c = 2`` and ``μ_c = 1``.

The free parameters for cloud droplets (``A_c``, ``B_c``) can be found analytically by integrating over the assumed 
mass distribution to find the prognostic variables
```math
\begin{align}
    N_c = ∫_0^∞ f_c(x) dx, \quad  L_c = ρ_a q_c = ∫_0^∞ x \; f_c(x) dx
\end{align}
```
where ``L_c = ρ_a q_c`` [kg / m$^3$] is the cloud liquid water content, ``ρ_a`` [kg / m$^3$] is the air density, and ``\overline{x}_c = \tfrac{L_c}{N_c}`` [kg] is the mean droplet mass.

!!! details "Derivation of the expressions for ``A_c`` and ``B_c``"

    TODO: Add derivation

Then[^2], 
```math
\begin{align}
    B_c = \left(\overline{x}_c \; \frac{Γ\left(\frac{ν_c+1}{μ_c}\right)}{Γ\left(\frac{ν_c+2}{μ_c}\right)}\right)^{-μ_c}, \qquad
    A_c = \frac{μ_c \, N_c \, }{Γ\left(\frac{ν_c+1}{μ_c}\right)} \; B_c^{\frac{ν_c+1}{μ_c}}.
\end{align}
```

[^2]: c.f. Eq. (80) in [SeifertBeheng2006](@cite)

##### Interactions with the P3 scheme
For interactions (e.g. collisions) with frozen particles, as parameterized in the P3 scheme, 
bulk rates are computed by integrating over diameter instead of mass.

!!! details "Change of variables between mass (x) and diameter (D)"
    Assuming spherical raindrops, the raindrop mass ``x`` is related to the raindrop diameter ``D`` by
    ```math
    x = k_m D^3= \frac{π}{6} ρ_w D^3,
    ```
    where ``k_m = \tfrac{ρ_w π}{6}``, and ``ρ_w`` is the density of water. In terms of ``x``, we have ``D = k_m^{-1/3} x^{1/3}``. To do a change of variable
    from ``x`` to ``D``, we need to consider the transformation of ``dx`` to ``dD``, which would appear in an integral against ``f``.
    We can calculate ``\tfrac{dx}{dD} = 3k_m D^2`` and thus ``dx = 3k_m D^2 dD``.
    Thus, an integral expression with ``f(x)dx`` becomes
    ```math
    f(x)dx \rightarrow f(x(D)) 3k_m D^2 dD
    ```
    Similarly, starting from an integral expression ``n(D)dD``, in terms of ``x`` we have
    ``dD = \tfrac{1}{3} k_m^{-1/3} x^{-2/3} dx``, so
    ```math
    n(D)dD \rightarrow n(D(x)) \frac{1}{3} k_m^{-1/3} x^{-2/3} dx
    ```

    We remark that if the integration bounds are at ``0`` and ``∞``,  then they remain unchanged by the 
    change of variable. For any other bounds, we would need to consider the transformation of the bounds.

The resulting cloud droplet distribution in terms of diameter is
```math
\begin{align}
    n_c(D) 
    = f_c(x(D)) \frac{ρ_w π}{2} D^2 
    = A_c x(D)^{\nu_c} \exp\left(-B_c x(D)^{\mu_c} \right) 3k_m D^2
\end{align}
```
where ``x(D) = k_m D^3 = \tfrac{π}{6} ρ_w D^3``, and ``ρ_w`` is the density of water. Explicitly written out, we have
```math
\begin{align}
    n_c(D) = 3 A_c k_m^{ν_c + 1} D^{3ν_c + 2} \exp\left(-B_c k_m^{μ_c} D^{3μ_c} \right)
\end{align}
```

!!! note "Undeveloped cloud droplet spectrum"
    In the derivation of the parametrization, it is assumed that the cloud droplet distribution
    ``f_c(x)`` does not contain a significant number of droplets with masses almost equal
    or larger than ``x^*``. This is reffered to as the undeveloped cloud droplet spectrum assumption.
    Similarly the raindrop distribution does not contain a significant number of rain drops
    with masses almost equal or smaller than ``x^*``. These assumptions allow us
    to simplify the calculation of moments of the distributions by integrating from zero to infinity.

#### Raindrops

!!! todo "Notation fixes"
    This section describes the psd in terms of `λ_r`. In the code, we use `D_mean ≡ 1/λ_r`.
    Once this change has percolated through the code, we should use `D_mean` consistently, including in the equations below.

The raindrop number distribution, as a function of diameter ``D`` [m], is assumed to follow an exponential distribution
```math
\begin{align}
    n_r(D) = N_0 \exp\left(- \frac{D}{\overline{D}_r} \right).
\end{align}
```

To find the free parameters ``N_0`` and ``\overline{D}_r``, we write the expressions for the raindrop number density and the raindrop liquid water content:
```math
\begin{align}
    N_r = ∫_0^∞ n_r(D) dD, \quad  
    L_r = ρ_a q_r = k_m ∫_0^∞ D^3 \; n_r(D) dD
\end{align}
```
where ``L_r`` [kg / m$^3$] is the raindrop content. The second expression comes from integrating the raindrop distribution against the raindrop mass ``x(D) = k_m D^3 = \tfrac{ρ_w π}{6} D^3``.

!!! details "Derivation of the expressions for N₀ and D̄ᵣ"

    The moments of the exponential distribution are given by
    ```math
    M^k =∫_0^∞ D^k n_r(D) dD = N_0 ∫_0^∞ D^k \exp\left(- \frac{D}{\overline{D}_r} \right) dD = N_0 Γ(k+1) \overline{D}_r^{k+1}.
    ```
    for any ``k>-1``, ``\overline{D}_r > 0``.
    This allows us to solve the integrals above for ``N_0`` and ``\overline{D}_r``:
    ```math
    N_r = N_0 \overline{D}_r, \quad  
    L_r = k_m N_0 Γ(4) \overline{D}_r^4
        = 6 k_m N_r \overline{D}_r^3
    ```
    where the last equality uses ``N_0 = N_r / \overline{D}_r``.
    This implies that
    ```math
        N_0            = \frac{N_r}{\overline{D}_r}, \quad
        \overline{D}_r = \left( \frac{L_r}{6 k_m N_r} \right)^{1/3}
                       = \left( \frac{\overline{x}_r}{6 k_m} \right)^{1/3}
                       = \left( \frac{\overline{x}_r}{ρ_w π} \right)^{1/3}
    ```
    
we find
```math
\begin{align}
    \overline{D}_r 
      = \left( \frac{\overline{x}_r}{6 k_m} \right)^{1/3}
      = \left( \frac{\overline{x}_r}{π ρ_w} \right)^\frac{1}{3}, \qquad
    N_0 = \frac{N_r}{\overline{D}_r}.
\end{align}
```
where ``\overline{x}_r = L_r / N_r`` is the mean raindrop mass. This implies a relation between the mean mass and mean diameter 
```math
\begin{align}
    \overline{x}_r = 6 k_m \overline{D}_r^3 = ρ_w π \overline{D}_r^3.
\end{align}
```

!!! details "Raindrop distribution in terms of mass"
    To express the raindrop distribution in terms of mass, we use the change of variable machinery from the previous section, noting in particular ``D = k_m^{-1/3} x^{1/3}`` where ``k_m = \tfrac{ρ_w π}{6}``, and ``dD = \tfrac{1}{3} k_m^{-1/3} x^{-2/3} dx``.
    In terms of mass, the distribution is given by
    ```math
    \begin{align*}
        f_r(x) 
        &= n(D(x)) \frac{1}{3} k_m^{-1/3} x^{-2/3} \\
        &= \frac{N_0 k_m^{-1/3}}{3} x^{-2/3} \exp\left(- \frac{k_m^{-1/3} x^{1/3}}{\overline{D}_r}\right) \\
        &= \frac{N_r k_m^{-1/3}}{3\overline{D}_r} x^{-2/3} \exp\left(- \frac{k_m^{-1/3} x^{1/3}}{\overline{D}_r}\right) \\
    \end{align*}
    ```
    We identify ``B_r`` by ``B_r = \tfrac{k_m^{-1/3}}{\overline{D}_r} = \left(\tfrac{\overline{x}_r}{6}\right)^{-1/3}``. This implies that ``A_r = \frac{N_r}{3}B_r``.
    
In the model code, we provide two options for calculating the raindrop number distribution parameters ``N_0`` and ``\overline{D}_r``. One option is to use the expressions above. The other option, described below, limits the range of ``\overline{x}_r``, ``N_0`` and ``\overline{D}_r`` to avoid numerical artifacts.

##### Limiting the range of D̄ᵣ

!!! todo "Comment on the limiting procedure"
    Given that each step limits the range of a "problematic" quantity, which is then used in the next limiting step, it is odd that we need to consider additional limits in subsequent steps. We should check whether the subsequent limiting steps are reached, and thus whether those limits are needed at all.

[WackerSeifert2001](@cite) showed that every one- or two-moment scheme may experience numerical artifacts, especially as ``N→0`` and ``L→0``, which results in ``\overline{x}_r = L_r / N_r`` being ill-defined. This manifests, for example, in the bulk rates for sedimentation and evaporation. To avoid this, [SeifertBeheng2006](@cite) proposed limiting the range of ``λ_r ≡ \tfrac{1}{\overline{D}_r}``, ``N_0``, and ``\overline{x}_r`` by a sequence of steps.

First, compute a limited mean mass by
```math
\begin{align}
    \tilde{x}_r → \overline{x}_{r, \text{min}} ≤ \frac{L_r}{N_r} ≤ \overline{x}_{r, \text{max}}.
\end{align}
```
Then, limit ``N_0`` to be consistent with the limited mean mass 

!!! details "Derivation of the limit on N₀"
    Substitute the expression ``\overline{D}_r = \left(\tfrac{\overline{x}_r}{ρ_w π}\right)^{1/3}`` into the expression for ``N_0``:
    ```math
    N_0 = \frac{N_r}{\overline{D}_r} = N_r \left(\frac{π ρ_w}{\tilde{x}_r}\right)^\frac{1}{3}.
    ```

The resulting expression is
```math
\begin{align}
    N_0 → N_{0, \text{min}} ≤ N_r \left(\frac{π ρ_w}{\tilde{x}_r}\right)^\frac{1}{3} ≤ N_{0, \text{max}}.
\end{align}
```
With the limited ``N_0``, we then limit ``λ_r ≡ \tfrac{1}{\overline{D}_r}``.

!!! details "Derivation of the limit on λᵣ"
    !!! todo "TODO: Change from limiting λᵣ to limiting D̄ᵣ"
        The limit on ``λ_r`` is what is proposed in [SeifertBeheng2006](@cite).
        However, we use ``\overline{D}_r`` in the code, so we should limit ``\overline{D}_r`` instead.
        This requires a change to the `ClimaParams` parameters.
      
    Next, limit ``λ_r`` by considering the integral expression for ``L_r``:
    ```math
    \begin{align*}
        L_r &= k_m ∫_0^∞ D^3 n_r(D) dD \\
            &= k_m N_0 ∫_0^∞ D^3 \exp\left(-λ_r D\right) dD \\
            &= k_m N_0 \frac{Γ(3+1)}{λ_r^{3+1}} 
             = \frac{6 k_m N_0}{λ_r^4}.
    \end{align*}
    ```
    which implies
    ```math
    λ_r = \left( \frac{6 k_m N_0}{L_r} \right)^\frac{1}{4}
        = \left( \frac{ρ_w π N_0}{L_r} \right)^\frac{1}{4}.
    ```

The resulting expression is
```math
\begin{align}
    λ_r → λ_{r, \text{min}} ≤ \left( \frac{ρ_w π N_0}{L_r} \right)^\frac{1}{4} ≤ λ_{r, \text{max}}.
\end{align}
```
Finally, we obtain the limited mean mass ``\overline{x}_r`` from the relation ``N_0 = λ_r N_r = λ_r \tfrac{L_r}{\overline{x}_r}``, which is consistent with the limited ``N_0`` and ``λ_r``:
```math
\begin{align}
    \overline{x}_r → \overline{x}_{r, \text{min}} ≤ \frac{λ_r L_r}{N_0} ≤ \overline{x}_{r, \text{max}}.
\end{align}
```
We obtain the limited ``\overline{D}_r`` by ``\overline{D}_r = \tfrac{1}{λ_r}``.

### Autoconversion

The autoconversion rate can be estimated by looking at variations in the second moment of the particle mass spectrum. Specifically, the sum of variations in the second moment of the cloud droplets and raindrops spectrum equals the variations in the second moment of the particle mass spectrum:
```math
\begin{align}
\frac{\partial Z}{\partial t} = \frac{\partial Z_c}{\partial t} + \frac{\partial Z_r}{\partial t},
\end{align}
```
where ``Z`` represents the second moment, and ``c`` and ``r`` subscripts denote cloud and rain categories respectively. In the early stages of rain evolution, an estimate of the variations in the second moment of the particle mass spectrum is obtained from the stochastic collection equation: ``\tfrac{∂Z}{∂t} ≈ 2k_c L_c M_c^{(3)}``, where ``M_c^{(3)}`` is the third moment of the cloud droplets spectrum. Using these equations, along with computing ``Z_c``, ``M_c^{(3)}``, ``Z_r`` directly by integrating the distribution functions, allows us to derive an equation for the autoconversion rate. To simplify the derivation, we assume that in the initial stage of the rain evolution raindrops have sizes of the order of ``x^*`` and the mean radius of cloud droplets is much less than ``x^*``. This approach yields an approximation of the autoconversion rate in the early stages of rain evolution. The early stage rain evolution assumption is then relaxed by means of a universal function that depends on a process time scale.

The rate of change of rain specific content by autoconversion is finally expressed as
``` math
\begin{equation}
  \left. \frac{∂q_{rai}}{∂t} \right|_\text{acnv} 
  = \frac{k_{cc}}{20 \; x^* \; ρ} \frac{(ν+2)(ν+4)}{(ν+1)^2} (q_{liq} ρ)^2 \overline{x}_c^2 \left(1+\frac{ϕ_\text{acnv}(τ)}{1-τ^2}\right)\frac{ρ_0}{ρ},
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``ρ`` is the moist air density,
  - ``ρ_0 = 1.225`` [kg / m$^3$] is the air density at surface conditions,
  - ``k_{cc}`` is the cloud-cloud collection kernel constant,
  - ``ν`` is the cloud droplet gamma distribution parameter,
  - ``x^*`` is the drop mass separating the cloud and rain categories
  - ``\overline{x}_c = (q_{liq} ρ) / N_{liq}`` is the cloud droplet mean mass with ``N_{liq}`` denoting the cloud droplet number density. Here, to ensure numerical stability, we limit ``\overline{x}_c`` by the upper bound of ``x^*``.

The function ``ϕ_\text{acnv}(τ)`` is used to correct the autoconversion rate for the undeveloped cloud droplet spectrum and the early stage rain evolution assumptions. This is a universal function which is obtained by fitting to numerical results of the SCE:
```math
\begin{equation}
  ϕ_\text{acnv}(τ) = A τ^a(1-τ^a)^b,
\end{equation}
```
where
  - ``τ = 1 - \tfrac{q_{liq}}{q_{liq} + q_{rai}}`` is a dimensionless internal time scale with ``q_{rai}`` being the cloud liquid water specific content.

The default free parameter values are:

|   symbol   | default value |
|------------|---------------|
|``ν``       | ``2``         |
|``A``       | ``400``       |
|``a``       | ``0.7``       |
|``b``       | ``3``         |

The rate of change of raindrops number density is
``` math
\begin{equation}
  \left. \frac{∂N_{rai}}{∂t} \right|_\text{acnv} 
  = \frac{ρ}{x^*} \left. \frac{d \, q_{rai}}{dt} \right|_\text{acnv},
\end{equation}
```
and the rate of change of liquid water specific content and cloud droplets number density are
``` math
\begin{align}
  \left. \frac{∂q_{liq}}{∂t} \right|_{acnv} 
    = - \left. \frac{∂q_{rai}}{∂t} \right|_{acnv},\\
  \left. \frac{∂N_{liq}}{∂t} \right|_{acnv} 
    = -2 \left. \frac{∂N_{rai}}{∂t} \right|_{acnv}.
\end{align}
```
!!! note
    The Seifert and Beheng parametrization is formulated for the rate of change of liquid water content ``L = ρ q``. Here, we assume constant ``ρ`` and divide the rates by ``ρ`` to derive the equations for the rate of change of specific contents.

### Accretion
An approximation for the accretion rate is obtained by directly evaluating the integral:
```math
\begin{align}
    \left. \frac{∂q_{rai}}{∂t} \right|_\text{accr} = \frac{1}{ρ} ∫_{x=0}^∞ ∫_{y=0}^∞ f_c(x) f_r(y) K(x,y) x dy dx.
\end{align}
```
Similar to the autoconversion rate, the accretion rate is modified by a universal function. Thus, the rate of change of rain specific content by accretion becomes
```math
\begin{align}
  \left. \frac{∂ q_{rai}}{∂t} \right|_\text{accr} = & \frac{k_{cr}}{ρ} (q_{liq} ρ) (q_{rai} ρ) ϕ_\text{accr}(τ),\nonumber\\
   = & k_r ρ q_{liq} q_{rai} ϕ_\text{accr}(τ) \left(\frac{ρ_0}{ρ}\right)^{1/2},
\end{align}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``q_{rai}`` is the rain liquid water specific content,
  - ``ρ`` is the moist air density,
  - ``ρ_0`` is the air density at surface conditions,
  - ``k_{cr}`` is the cloud-rain collection kernel constant.

The universal function ``ϕ_\text{accr}(τ)`` is used to correct the accretion rate for the assumption of collsion efficiency being one. Fitting to numerical solutions of the SCE obtains:
```math
\begin{equation}
  ϕ_\text{accr}(τ) = \left(\frac{τ}{τ+τ_0}\right)^c,
\end{equation}
```
where
  - ``τ = 1 - \tfrac{q_{liq}}{q_{liq} + q_{rai}}`` is a dimensionless internal time scale.

The default free parameter values are:

|   symbol   | default value   |
|------------|-----------------|
|``τ_0``     | ``5 × 10^{-5}`` |
|``c``       | ``4``           |

The rate of change of raindrops number density by accretion is zero, and the rate of change of liquid water specific content and cloud droplets number density are
``` math
\begin{align}
  \left. \frac{∂q_{liq}}{∂t} \right|_\text{accr} = - \left. \frac{∂q_{rai}}{∂t} \right|_\text{accr},\\
  \left. \frac{∂N_{liq}}{∂t} \right|_\text{accr} = \frac{ρ}{\overline{x}_c} \left. \frac{∂q_{liq}}{∂t} \right|_\text{accr},
\end{align}
```
where ``\overline{x}_c = \tfrac{q_{liq}}{N_{liq}}`` is the cloud droplet mean mass.

### Cloud droplet self-collection

An approximation for the self-collection rate of cloud droplets is obtained by the following equation:
```math
\begin{align}
   \left. \frac{∂N_{liq}}{∂t} \right|_\text{sc} 
     &= \left. \frac{∂N_{liq}}{∂t} \right|_\text{acnv, sc} 
      - \left. \frac{∂q_{rai}}{∂t} \right|_\text{acnv}\nonumber\\
     &= - \frac{1}{2}∫_{x=0}^{∞}∫_{y=0}^{∞} f_c(x) f_c(y) K(x,y) dy dx 
        - \left. \frac{d \, q_{rai}}{dt} \right|_\text{acnv}.
\end{align}
```
Direct evaluation of the integral results in the following approximation of the rate of change of cloud droplets number density due to self-collection
``` math
\begin{equation}
  \left. \frac{∂N_{liq}}{∂t} \right|_\text{sc} 
    = -k_{cc} \frac{ν + 2}{ν + 1} \frac{ρ_0}{ρ} (q_{liq} ρ)^2 
    - \left. \frac{∂N_{liq}}{∂t} \right|_\text{acnv},
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``ρ`` is the moist air density,
  - ``ρ_0`` is the air density at surface conditions,
  - ``k_{cc}`` is the Long's collection kernel constant,
  - ``ν`` is the cloud droplet gamma distribution parameter,
  - ``\left. \frac{d \, N_{liq}}{dt} \right|_\text{acnv}`` is the rate of change of cloud droplets number density by autoconversion.

### Raindrop self-collection

An approximation for rate of change of raindrops number density due to self-collection is obtained by directly evaluating the integral:
```math
\begin{align}
    \left. \frac{∂N_{rai}}{∂t} \right|_\text{sc}= -\frac{1}{2} ∫_{x=0}^∞ ∫_{y=0}^∞ f_r(x) f_r(y) K(x,y) dy dx.
\end{align}
```
This yields,
```math
\begin{equation}
  \left. \frac{∂N_{rai}}{∂t} \right|_\text{sc} = -k_{rr} N_{rai} (q_{rai} ρ) \left(1+\frac{κ_{rr}}{B_r} \right)^d \left(\frac{ρ_0}{ρ}\right)^{1/2},
\end{equation}
```
where:
  - ``q_{rai}`` is the rain water specific content,
  - ``ρ`` is the moist air density,
  - ``ρ_0`` is the air density at surface conditions,
  - ``N_{rai}`` is the raindrops number density,
  - ``k_{rr}`` and ``κ_{rr}`` are the rain-rain collection kernel constants.
  - ``B_r`` is the raindrops mass distribution parameter ``B_r = \left(\frac{6}{\overline{x}_r}\right)^{1/3}``.

The default constant value is:

|   symbol   | default value                       |
|------------|-------------------------------------|
|``d``       | ``-5``                              |

!!! note
    In the paper ``d=-9`` which seems to be a mistake! Evaluating the integral for derving the self-collection rate results in ``d=-5``.

!!! note
    For the same numerical instabilities which in the paper are mentioned for terminal velocity and evaporation, here for rain self-collection, the value of ``\lambda_r`` is bounded within a range. In fact we first compute the bounded ``\lambda_r`` based on drop diameter by the algorithm given in the paper and then convert it to ``\lambda_r`` based on mass (the conversion can be done by multiplying to a constant value).

### Raindrops breakup

Raindrops breakup is modeled by assuming that in a precipitation event coalescence and breakup ultimately reach an equilibrium with a self-similar equilibrium size distribution. As a result, the breakup process can be coupled to raindrops self-collection by the following parameterization

```math
\begin{equation}
  \left. \frac{\partial N_{rai}}{\partial t} \right|_{br} = -[\Phi_{br}(\Delta \overline{D}_r) + 1] \left. \frac{\partial N_{rai}}{\partial t} \right|_{sc},
\end{equation}
```
where ``\Delta \overline{D}_r = \overline{D}_r - \overline{D}_{eq}`` with ``\overline{D}_r`` denoting the mean volume raindrop diameter and ``\overline{D}_{eq}`` being the equilibrium mean diameter. The function ``\Phi_{br}(\Delta \overline{D}_r)`` is given by
```math
  \begin{align}
    \Phi_{br}(\Delta \overline{D}_r) =
    \begin{cases}
    -1, \quad & \overline{D}_r < \overline{D}_{threshold},\\
    k_{br} \Delta \overline{D}_r, \quad & \overline{D}_{threshold} < \overline{D}_r < \overline{D}_{eq},\\
    2 (exp(\kappa_{br} \Delta \overline{D}_r) -1), \quad & \overline{D}_{eq} < \overline{D}_r.
    \end{cases}
  \end{align}
```

The default free parameter values are:

|   symbol                   | default value                       |
|----------------------------|-------------------------------------|
|``k_{br}``                  | ``1000 \, m^{-1}``                  |
|``\kappa_{br}``             | ``2300 \, m^{-1}``                  |
|``\overline{D}_{threshold}``| ``0.35 \times 10^{-3}  \, m``       |
|``\overline{D}_{eq}``       | ``0.9 \times 10^{-3}  \, m``        |

!!! note
    In the paper for ``\overline{D}_{eq} < \overline{D}_r`` the equation ``\Phi_{br}(\Delta \overline{D}_r) = 2 exp(\kappa_{br} \Delta \overline{D}_r) -1`` is given. This equations seems to be missing parentheses as the equation must be continuous at ``\Delta \overline{D}_r = 0`` as shown in Fig. 2 of the paper.

### Rain evaporation

The parametrization of rain evaporation is obtained by considering the time scale of evaporation of individual raindrops:
```math
\begin{equation}
  \tau_{eva} = \frac{x_r}{\frac{dx_r}{dt}\bigg|_{eva}} = \frac{x_r}{2 \pi G_{lv}(T, p) S D_r(x_r) F_v(x_r)},
\end{equation}
```
where
```math
\begin{equation}
  G_{lv}(T, p) = \left[\frac{R_v T}{p_{lv}(T) D_v} + \frac{L_{lv}}{K_T T} \left(\frac{L_{lv}}{R_v T}-1\right)\right]^{-1}
\end{equation}
```
with temperature ``T``, thermal conductivity ``K_T``, diffucivity of water vapor ``D_v``, specific gas constant for water vapor ``R_v``, latent heat of evaporation ``L_{lv}`` and liquid-vapor saturation pressure ``p_{lv}``. The ventilation factor is given by ``F_v(x_r) = a_v + b_v N_{Sc}^{1/3} N_{Re}(x_r)^{1/2}`` where ``N_{Sc} = \nu_{air} / D_v`` is the Schmidt number and ``N_{Re}(x_r) = \frac{v_r (x_r) D_r (x_r)}{\nu_{air}}`` is the Reynolds number with kinematic viscosity of air ``\nu_{air}``. The average evaporation rates are obtained from the following integral:
```math
\begin{equation}
  \frac{\partial M_r^k}{\partial t}\bigg|_{eva} = \int_0^\infty \frac{x^k f_r(x)}{\tau_{eva}} dx = 2 \pi G_{lv}(T, p) S \int_0^\infty D_r(x) F_v(x) f_r(x) x^{k-1} dx.
\end{equation}
```
where the superscript ``k`` indicates the moment number, ``k=0`` for number-weighted and ``k=1`` for mass-weighted average. Here in the computation of the evaporation rate, a power-law fall speed is assumed: ``v_r(x) \cong \alpha_r x^{\beta_r} \left(\frac{\rho_0}{\rho}\right)^{\frac{1}{2}}``. Evaluating the integral results in:
```math
\begin{equation}
  \frac{\partial M_r^k}{\partial t}\bigg|_{eva} = 2 \pi G_{lv}(T, p) S N_{rai} D_r(\overline{x}_r) \overline{F}_{v,\, k}(\overline{x}_r) \overline{x}_r^{k-1},
\end{equation}
```
where ``\overline{F}_{v,\, k}`` is an average ventilation factor for the ``k``-th moment:
```math
\begin{equation}
\overline{F}_{v,\, k}(\overline{x}_r) = a_{vent,\, k} + b_{vent,\, k} N_{Sc}^{1/3} N_{Re}(\overline{x}_r)^{1/2},
\end{equation}
```
with
```math
\begin{align}
  a_{vent,\, k} &= a_v 6^{2/3-k} \Gamma(3k-1),\\
  b_{vent,\, k} &= b_v 6^{1/2-\beta_r/2-k} \Gamma(3k-1/2+3\beta_r/2).
\end{align}
```

!!! note
    For ``k = 0`` the integral for computing the mean evaporation rate does not converge. In this case it is reasonable to change the lower bound of the integral to ``x=x^*``. The results remain the same except that the Gamma functions in the equations for ``a_{vent,\, 0}`` and ``b_{vent,\, 0}``, which are ``\Gamma(-1)`` and ``\Gamma(-1/2+3\beta_r/2)``, are replaced by the upper incomplete gamma function ``\Gamma(-1, (6 x^* / \overline{x}_r)^{1/3})`` and ``\Gamma(-1/2+3\beta_r/2, (6 x^* / \overline{x}_r)^{1/3})``, respectively. This issue and the suggested workaround are not mentioned in the paper.

The two-moment parametrization of evaporation suffers from the similar numerical instability issues as the sedimentation scheme. Thus, the same limiting algorithm as the sedimentation scheme is applied here to bound size distribution parameters. These limited parameters are then used to compute the mean raindrop mass by the following equation:
```math
\begin{equation}
  \overline{x}_r = max \left(\overline{x}_{r,\, min} , min \left(\overline{x}_{r,\, max} , \frac{\rho q_{rai} \lambda_r}{N_0}\right)\right).
\end{equation}
```
This mean mass is used for computing the evaporation rate.

The default free parameter values are:

|   symbol                   | default value                                   |
|----------------------------|-------------------------------------------------|
|``a_v``                     | ``0.78``                                        |
|``b_v``                     | ``0.308``                                       |
|``\alpha_r``                | ``159 \, m \cdot s^{-1} \cdot kg^{-\beta_r}``   |
|``\alpha_r``                | ``0.266``                                       |

!!! note
    In our implementation we approximate the incomplete gamma function in order to get good
    performance on the GPU. Below we show the evaporation rates for the number concnetration
    and mass using both the exact and approximated gamma function

```@example
include("plots/RainEvapoartionSB2006.jl")
```
![](SB2006_rain_evaporation.svg)

### Terminal velocity

The number- and mass-weighted terminal velocities for rain water are obtained by calculating the following integral:
```math
\begin{equation}
  \overline{v}_{r,\, k} = \frac{1}{M_r^k} \int_0^\infty x^k f_r(x) v(x) dx,
\end{equation}
```
where the superscript ``k`` indicates the moment number,  with``k=0`` for number density and ``k=1`` for mass.

The individual terminal velocity of particles is approximated by
```math
v(D) = \left(\rho_0/\rho\right)^{\frac{1}{2}} [a_R - b_R exp(-c_R D)]
```
where ``a_R``, ``b_R`` and ``c_R`` are free parameters and ``D`` is the particle diameter.
Evaluating the integral results in
```math
\begin{equation}
  \overline{v}_{r,\, k} = \left(\frac{\rho_0}{\rho}\right)^{\frac{1}{2}}\left[a_R - b_R \left(1+\frac{c_R}{\lambda_r}\right)^{-(3k+1)}\right],
\label{eq:SBTerminalVelocity}
\end{equation}
```
where ``\lambda_r`` is the raindrop size distribution parameter (based on diameter): ``\lambda_r = (\phi \rho_w/\overline{x}_r)^{1/3}``. To avoid numerical instabilities, especially when ``N_{rai} \rightarrow 0`` and ``q_{rai} \rightarrow 0``, ``\lambda_r`` is bounded. The limiting algorithm is as follows:
```math
\begin{align}
\widetilde{x}_r &= \max \left[ \overline{x}_{r,\, \min} , \min \left(\overline{x}_{r,\, \max} , \frac{\rho q_{rai}}{N_{rai}}\right)\right] ,\\
N_0 &= \max \left[ N_{0,\, \min} , \min \left(N_{0,\, \max} , N_{rai}\left(\frac{\pi \rho_w}{\widetilde{x}_r}\right)^{\frac{1}{3}}\right)\right],\\
\lambda_r &= \max \left[ \lambda_{\min} , \min \left(\lambda_{\max} , \left(\frac{\pi \rho_w N_0}{\rho q_{rai}}\right)^{\frac{1}{4}}\right)\right].
\end{align}
```

When the limiting algorithm is not applied, the terminal velocity given by eq. (\ref{eq:SBTerminalVelocity}) can become negative for small mean radius values. This occurs because the equation for the individual particle's terminal velocity may predict negative values for small particles. To avoid negative terminal velocities (or a sudden transition to zero velocity if we return zero instead of negative values), we need to adjust the integration bounds. Specifically, the integrals should be evaluated from the radius at which the individual terminal velocity is zero to infinity. This adjustment leads to the following equation for the terminal velocity:
```math
\begin{align}
  \overline{v}_{r,\, k} &= \frac{1}{M_r^k} \int_{r_c}^\infty x^k f_r(x) v(x) dx \nonumber\\
  &= \left(\frac{\rho_0}{\rho}\right)^{\frac{1}{2}}\left[a_R \frac{\Gamma(3k+1, 2r_c\lambda_r)}{\Gamma(3k+1)} - b_R \frac{\Gamma(3k+1, 2r_c(\lambda_r + c_R))}{\Gamma(3k+1)} \left(1+\frac{c_R}{\lambda_r}\right)^{-(3k+1)}\right],
\label{eq:SBModifiedTerminalVelocity}
\end{align}
```
where ``r_c = -\ln(a_R / b_R)/(2c_R)``.

The default parameter values are:

|   symbol                   | default value                       |
|----------------------------|-------------------------------------|
|``a_R``                     | ``9.65 \, m \cdot s^{-1}``          |
|``b_R``                     | ``10.3 \, m \cdot s^{-1}``          |
|``c_R``                     | ``600 \, m^{-1}``                   |
|``\overline{x}_{r,\, \min}`` | ``6.54 \times 10^{-11} \, m``        |
|``\overline{x}_{r,\, \max}`` | ``5 \times 10^{-6}  \, m``          |
|``N_{0,\, \min}``            | ``3.5 \times 10^{5}  \, m^{-4}``    |
|``N_{0,\, \max}``            | ``2 \times 10^{10}  \, m^{-4}``      |
|``\lambda_{\min}``           | ``1 \times 10^{3}  \, m^{-1}``      |
|``\lambda_{\max}``           | ``4 \times 10^{4}  \, m^{-1}``      |

Below we compare number-weighted (left) and mass-weighted (right) terminal velocities
  for the original parameterization of SB2006 [SeifertBeheng2006](@cite)
  and the modifed parameterzation given by eq. (\ref{eq:SBModifiedTerminalVelocity}),
  without the limiting distribution shape factor $\lambda_r$.
```@example
include("plots/TerminalVelocity2M.jl")
```
![](2M_terminal_velocity_comparisons.svg)

For the Chen et al. [Chen2022](@cite) terminal velocity parameterization,
  the number- and mass-weighted group terminal velocities are:
```math
\begin{equation}
  \overline{v_k} = \frac{\int_0^\infty v(D) \, D^k \, n(D) \, dD}
             {\int_0^\infty D^k \, n(D) \, dD}
             = (\phi)^{\kappa} \Sigma_{i} \frac{a_i \lambda^\delta \Gamma(b_i + \delta)}
             {(\lambda + c_i)^{b_i + \delta} \; \Gamma(\delta)},
\end{equation}
```
where $\Gamma$ is the gamma function, $\delta = k + 1$, and
$\lambda$ is the size distribution parameter.
The velocity $\overline{v_k}$ is the number-weighted mean terminal velocity when k = 0,
and the mass-weighted mean terminal velocity when k = 3.

## Additional 2-moment microphysics options

The other autoconversion and accretion rates in the `Microphysics2M.jl` module
  are implemented after Table 1 from [Wood2005](@cite) and are based on the works of
  [KhairoutdinovKogan2000](@cite),
  [Beheng1994](@cite),
  [TripoliCotton1980](@cite) and
  [LiuDaum2004](@cite) respectively.
From the above works:
  (i) the [KhairoutdinovKogan2000](@cite) parameterisation is based
  on a fit to drop spectrum resolving scheme and designed to work
  for stratocumulus topped boundary layers,
  (ii) the [Beheng1994](@cite) parameterisation is based on a fit
  to stochastic collection equation,
  (iii) the [TripoliCotton1980](@cite) parameterisation is developed
  for a deep convective case, and
  (iv) the [LiuDaum2004](@cite) parameterisation is derived to
  include the effects of relative dispersion
  of the cloud droplet size distribution on precipitation formation rates
  and assumes a modified gamma distribution.

### Autoconversion

#### Khairoutdinov and Kogan (2000)

``` math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{acnv} = A \; q_{liq}^a \; N_d^b \; \rho^c
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``N_d`` is the cloud droplet concentration,
  - ``\rho`` is the air density,

and the default free parameter values are:

|   symbol   | default value             |
|------------|---------------------------|
|``A``       | ``7.42 \times 10^{13} ``  |
|``a``       | ``2.47``                  |
|``b``       | ``-1.79``                 |
|``c``       | ``-1.47``                 |


#### Beheng (1994)

``` math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{acnv} = \frac{C \; d^a \; (q_{liq} \rho)^b \; N_d^c}{\rho}
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``N_d`` is the cloud droplet number concentration,

and the default free parameter values are:

|   symbol   | default value             |
|------------|---------------------------|
|``C``       | ``3 \times 10^{34} ``     |
|``a``       | ``-1.7``                  |
|``b``       | ``4.7``                   |
|``c``       | ``-3.3``                  |
|``d``       | ``9.9`` for ``N_d < 200 cm^{-3}``,  ``3.9`` for ``N_d > 200 cm ^{-3}`` |


#### Tripoli and Cotton (1980)

```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{acnv} =
    D \; q_{liq}^a \; N_d^b \; \mathrm{H}(q_{liq} - q_{liq\_threshold})
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``q_{liq_threshold}`` is the cloud liquid to rain water threshold,
  - ``N_d`` is the cloud droplet number concentration,
  - ``\mathrm{H}(x)`` is the Heaviside step function.

The cloud liquid to rain water autoconversion threshold is defined
  assuming spherical liquid water drops of radius equal to ``7 \mu m``:
```math
\begin{equation}
  q_{liq\_threshold} = \frac{4}{3} \pi \rho_w N_d r_{cm}^3
\end{equation}
```
where:
  - ``\rho_w`` is the liquid water density,

and the default free parameter values are:

|   symbol   | default value             |
|------------|---------------------------|
|``D``       | ``3268 ``                 |
|``a``       | ``\frac{7}{3}``           |
|``b``       | ``\frac{-1}{3}``          |
|``r_{cm}``  | ``7 \times 10^{-6} m``    |

#### Liu and Daum (2004)

```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{acnv} =
    \frac{E \; (q_{liq} \; \rho)^3 \; \mathrm{H}(R_6 - R_{6C})}{N_d \; \rho}
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``N_d`` is the cloud droplet number concentration,
  - ``\rho`` is the air density.

The parameterisation is formulated using mean volume radius
  ``r_{vol}`` expressed in ``\mu m`` which we compute as
```math
\begin{equation}
    r_{vol} = \left(\frac{\rho q_{liq}}{4/3 \pi \; \rho_w \; N_d}\right)^{1/3} 10^6
\end{equation}
```
where:
  - ``\rho_w`` is the liquid water density.

Then the ``R_6`` and ``R_{6C}`` are defined as
  - ``R_6 = \beta_6 \; r_{vol}``
  - ``R_{6C} = \frac{R_{C0}}{(q_{liq} \rho)^{1/6} R_6^{1/2}}``
  - ``\beta_6 = \left( \frac{r_{vol} + 3}{r_{vol}} \right)^{1/3}``
  - ``E = E_0 \beta_6^6``

|   symbol   | default value             |
|------------|---------------------------|
|``R_{C0}``  | ``7.5 ``                  |
|``E_0``     | ``1.08 \times 10^{10}``   |

#### Auto-conversion with time scale depending on number density

```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{acnv} =
    \frac{q_{liq}}{\tau_{acnv,\, 0} \left(\frac{N_d}{100\, cm^{-3}}\right)^{\alpha_{acnv}}}
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``N_d`` is the cloud droplet number concentration,
  - ``\tau_{acnv,\, 0}`` is the auto-conversion time scale at ``N_d = 100 cm^{-3}``.

The default free parameter values are:

|   symbol             | default value             |
|----------------------|---------------------------|
|``\tau_{acnv,\, 0}``  | ``1000\ s``               |
|``\alpha_{acnv}``     | ``1``                     |

### Accretion

#### Khairoutdinov and Kogan (2000)

```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{accr} = A \; (q_{liq} q_{rai})^a \; \rho^b
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``q_{rai}`` is the rain water specific content,
  - ``\rho``    is the air density,

and the default free parameter values are:

|   symbol   | default value             |
|------------|---------------------------|
|``A``       | ``67 ``                   |
|``a``       | ``1.15``                  |
|``b``       | ``-1.3``                  |


#### Beheng (1994)

```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{accr} = A \; q_{liq} \; q_{rai} \; \rho
\end{equation}
```
where:
  - ``q_{liq}`` is the cloud liquid water specific content,
  - ``q_{rai}`` is the rain specific content,
  - ``\rho``    is the air density,

and the default free parameter values are:

|   symbol   | default value             |
|------------|---------------------------|
|``A``       | ``6 ``                    |


#### Tripoli and Cotton (1980)

```math
\begin{equation}
  \left. \frac{d \, q_{rai}}{dt} \right|_{accr} = A \; q_{liq} \; q_{rai}
\end{equation}
```
where:
  - ``q_{liq}`` is cloud liquid water specific content
  - ``q_{rai}`` is rain specific content

and the default free parameter values are:

|   symbol   | default value             |
|------------|---------------------------|
|``A``       | ``4.7``                   |

## Number concentration adjustment (Horn 2012)
To ensure that the number concentration ``N`` remains physically consistent with the mass mixing ratio ``q``, the `Microphysics2M.jl` module provides a pair of relaxation tendencies that bound the implied mean particle mass ``x = \rho q / N`` within user-defined limits ``[x_{min}, x_{max}]``. This is useful in two-moment schemes to prevent unphysical droplet masses and improve numerical stability.

The tendencies are applied over a relaxation timescale ``\tau``, following the method described in [Horn2012](@cite). The total tendency can be written as:
```math
\frac{\partial N}{\partial t} = \frac{1}{\tau}\left[max\left(0, \frac{\rho q}{x_{max}}-N\right) + min\left(0, \frac{\rho q}{x_{min}}-N\right)\right].
```
This correction increases ``N`` when droplets are too heavy (i.e., ``x > x_{max}``) and decreases ``N`` when they are too light (``x < x_{min}``), while applying no correction when ``x`` lies within bounds. 

This formulation is implemented as two separate functions:
- `number_increase_for_mass_limit`: returns the rate needed to increase `N` when particles are too heavy (`x > x_{max}`).
- `number_decrease_for_mass_limit`: returns the rate needed to decrease `N` when particles are too light (`x < x_{min}`).
No correction is applied when x lies within bounds. The default parameter value is

|   symbol   | default value                       |
|------------|-------------------------------------|
|``\tau``    | ``100 \, s``                        |

The adjustment is assumed to exchange number concentration with a background reservoir of cloud condensation nuclei (CCN), so that any increase in droplet number is accompanied by a corresponding decrease in CCN number, and vice versa:
```math
\frac{\partial N_{CCN}}{\partial t} = -\frac{\partial N}{\partial t}
```

!!! note
    In the reference paper, this approach is given for the number concentration of cloud droplets; however, the same formulation can also be used for raindrops.

## Example figures

```@example
include("plots/Microphysics2M_plots.jl")
```
![](Autoconversion_accretion.svg)
