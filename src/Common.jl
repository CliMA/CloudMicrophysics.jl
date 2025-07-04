"""
    Module for functions shared by different parameterizations.
"""
module Common

import SpecialFunctions as SF

import ..Parameters as CMP
import ..ThermodynamicsInterface as TDI
const HPS = CMP.H2SO4SolutionParameters

export G_func_liquid
export G_func_ice
export H2SO4_soln_saturation_vapor_pressure
export a_w_xT
export a_w_eT
export a_w_ice
export Chen2022_vel_coeffs
export Chen2022_monodisperse_pdf
export Chen2022_exponential_pdf
export ventilation_factor

"""
    G_func_liquid(air_props, tps, T)

 - `air_props` - struct with air parameters
 - `tps` - struct with thermodynamics parameters
 - `T` - air temperature

Utility function combining thermal conductivity and vapor diffusivity effects
for vapor to liquid phase change.
"""
function G_func_liquid(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TDI.PS,
    T::FT,
) where {FT}
    R_v = TDI.Rᵥ(tps)
    L = TDI.Lᵥ(tps, T)
    p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)

    return FT(1) /
           (L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs)
end

"""
    G_func_ice(air_props, tps, T)

 - `air_props` - struct with air parameters
 - `tps` - struct with thermodynamics parameters
 - `T` - air temperature

Utility function combining thermal conductivity and vapor diffusivity effects.
for vapor to ice phase change.
"""
function G_func_ice(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TDI.PS,
    T::FT,
) where {FT}
    R_v = TDI.Rᵥ(tps)
    L = TDI.Lₛ(tps, T)
    p_vs = TDI.saturation_vapor_pressure_over_ice(tps, T)

    return FT(1) /
           (L / K_therm / T * (L / R_v / T - FT(1)) + R_v * T / D_vapor / p_vs)
end

"""
    A Heaviside step function
"""
heaviside(x::FT) where {FT} = FT(x > 0)

"""
    logistic_function(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the curve, characterizing steepness of the transition

Returns the value of the logistic function for smooth transitioning at thresholds. This is
a normalized curve changing from 0 to 1 while x varies from 0 to Inf (for positive k). For
x < 0 the value at x = 0 (zero) is returned. For x_0 = 0 H(x) is returned.
"""
function logistic_function(x::FT, x_0::FT, k::FT) where {FT}

    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return FT(1)
    end

    return FT(1) / (FT(1) + exp(-k * (x / x_0 - x_0 / x)))
end

"""
    logistic_function_integral(x, x_0, k)

 - `x` - independent variable
 - `x_0` - threshold value for x
 - `k` - growth rate of the logistic function, characterizing steepness of the transition

Returns the value of the indefinite integral of the logistic function, for smooth transitioning
of piecewise linear profiles at thresholds. This curve smoothly transition from y = 0
for 0 < x < x_0 to y = x - x_0 for x_0 < x.
"""
function logistic_function_integral(x::FT, x_0::FT, k::FT) where {FT}
    @assert k > 0
    @assert x_0 >= 0
    x = max(0, x)

    if abs(x) < eps(FT)
        return FT(0)
    elseif abs(x_0) < eps(FT)
        return x
    end

    # translation of the curve in x and y to enforce zero at x = 0
    trnslt = -log(FT(1) - exp(-k)) / k

    kt = k * (x / x_0 - FT(1) + trnslt)
    return (kt > FT(40)) ? x - x_0 : (log(FT(1) + exp(kt)) / k - trnslt) * x_0
end

"""
    H2SO4_soln_saturation_vapor_pressure(prs, x, T)

 - `prs` - a struct with H2SO4 solution free parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns the saturation vapor pressure above a sulphuric acid solution droplet in Pa.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight
"""
function H2SO4_soln_saturation_vapor_pressure(
    (;
        T_max,
        T_min,
        w_2,
        c1,
        c2,
        c3,
        c4,
        c5,
        c6,
        c7,
    )::CMP.H2SO4SolutionParameters{FT},
    x::FT,
    T::FT,
) where {FT}

    @assert T < T_max
    @assert T > T_min

    w_h = w_2 * x
    p_sol =
        exp(
            c1 - c2 * x + c3 * x * w_h - c4 * x * w_h^2 +
            (c5 + c6 * x - c7 * x * w_h) / T,
        ) * 100 # * 100 converts mbar --> Pa
    return p_sol
end

"""
    a_w_xT(H2SO4_prs, tps, x, T)

 - `H2SO4_prs` - a struct with H2SO4 solution free parameters
 - `tps` - a struct with thermodynamics parameters
 - `x` - wt percent sulphuric acid [unitless]
 - `T` - air temperature [K].

Returns water activity of H2SO4 containing droplet.
`x` is, for example, 0.1 if droplets are 10 percent sulphuric acid by weight.
"""
function a_w_xT(H2SO4_prs::HPS, tps::TDI.PS, x::FT, T::FT) where {FT}

    p_sol = H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x, T)
    p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)

    return p_sol / p_sat
end

"""
    a_w_eT(tps, e, T)

 - `tps` - struct with thermodynamics parameters
 - `e` - partial pressure of water [Pa]
 - `T` - air temperature [K].

Returns water activity of pure water droplet.
Valid when droplet is in equilibrium with surroundings.
"""
function a_w_eT(tps::TDI.PS, e::FT, T::FT) where {FT}
    # RH
    return e / TDI.saturation_vapor_pressure_over_liquid(tps, T)
end

"""
    a_w_ice(tps, T)

- `tps` - struct with thermodynamics parameters
 - `T` - air temperature [K].

Returns water activity of ice.
"""
function a_w_ice(tps::TDI.PS, T::FT) where {FT}

    return TDI.saturation_vapor_pressure_over_ice(tps, T) /
           TDI.saturation_vapor_pressure_over_liquid(tps, T)
end

"""
    Chen2022_vel_coeffs(coeffs, ρₐ)
    Chen2022_vel_coeffs(coeffs, ρₐ, ρᵢ)

Compute the coefficients for the Chen 2022 terminal velocity parametrization.

# Arguments
 - `coeffs`: a struct with terminal velocity free parameters
    - [`CMP.Chen2022VelTypeRain`](@ref): Fetch from Table B1
    - [`CMP.Chen2022VelTypeSmallIce`](@ref): Fetch from Table B2
    - [`CMP.Chen2022VelTypeLargeIce`](@ref): Fetch from Table B4
 - `ρₐ`: air density [kg/m³]
 - `ρᵢ`: apparent density of ice particles [kg/m³],
    only used for [`CMP.Chen2022VelTypeSmallIce`](@ref) and [`CMP.Chen2022VelTypeLargeIce`](@ref)

See [Chen2022](@cite) for more details.
"""
function Chen2022_vel_coeffs(coeffs::CMP.Chen2022VelTypeRain, ρₐ)
    (; ρ0, a, a3_pow, b, b_ρ, c) = coeffs
    # Table B1
    q = exp(ρ0 * ρₐ)
    ai = (a[1] * q, a[2] * q, a[3] * q * ρₐ^a3_pow)
    bi = (b[1] - b_ρ * ρₐ, b[2] - b_ρ * ρₐ, b[3] - b_ρ * ρₐ)
    ci = (c[1], c[2], c[3])
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return (aiu, bi, ciu)
end
function Chen2022_vel_coeffs(coeffs::CMP.Chen2022VelTypeSmallIce, ρₐ, ρᵢ)
    FT = eltype(coeffs)
    (; A, B, C, E, F, G) = coeffs
    # Table B3
    log_ρᵢ = log(ρᵢ)
    As = A[2] * log_ρᵢ^2 − A[3] * log_ρᵢ + A[1]
    Bs = FT(1) / (B[1] + B[2] * log_ρᵢ + B[3] / sqrt(ρᵢ))
    Cs = C[1] + C[2] * exp(C[3] * ρᵢ) + C[4] * sqrt(ρᵢ)
    Es = E[1] - E[2] * log_ρᵢ^2 + E[3] * sqrt(ρᵢ)
    Fs = -exp(F[1] - F[2] * log_ρᵢ^2 + F[3] * log_ρᵢ)
    Gs = FT(1) / (G[1] + G[2] / log_ρᵢ - G[3] * log_ρᵢ / ρᵢ)
    # Table B2
    ai = (Es * ρₐ^As, Fs * ρₐ^As)
    bi = (Bs + ρₐ * Cs, Bs + ρₐ * Cs)
    ci = (FT(0), Gs)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return (aiu, bi, ciu)
end
function Chen2022_vel_coeffs(coeffs::CMP.Chen2022VelTypeLargeIce, ρₐ, ρᵢ)
    FT = eltype(coeffs)
    (; A, B, C, E, F, G, H) = coeffs
    # Table B5
    log_ρᵢ = log(ρᵢ)
    Al = A[1] + A[2] * log_ρᵢ + A[3] * ρᵢ^FT(-3 / 2)
    Bl = exp(B[1] + B[2] * log_ρᵢ^2 + B[3] * log_ρᵢ)
    Cl = exp(C[1] + C[2] / log_ρᵢ + C[3] / ρᵢ)
    El = E[1] + E[2] * log_ρᵢ * sqrt(ρᵢ) + E[3] * sqrt(ρᵢ)
    Fl = F[1] + F[2] * log_ρᵢ - exp(log(-F[3]) - ρᵢ)
    Gl = (G[1] + G[2] * log_ρᵢ * sqrt(ρᵢ) + G[3] / sqrt(ρᵢ))^(-1)
    Hl = H[1] + H[2] * ρᵢ^FT(5 / 2) + exp(log(-H[3]) - ρᵢ)
    # Table B4
    ai = (Bl * ρₐ^Al, El * ρₐ^Al * exp(Hl * ρₐ))
    bi = (Cl, Fl)
    ci = (FT(0), Gl)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return (aiu, bi, ciu)
end

"""
    Chen2022_monodisperse_pdf(a, b, c)

# Arguments
 - `a`, `b`, `c`: free parameters defined in [Chen2022](@cite)

# Returns
 - `pdf(D)`: The monodisperse particle distribution function as a function of diameter, `D`, in [m/s].
"""
function Chen2022_monodisperse_pdf(a, b, c)
    return pdf(D) = a * D^b * exp(-c * D)
end

"""
    Chen2022_exponential_pdf(a, b, c, λ_inv, k)

 - a, b, c, - free parameters defined in Chen etl al 2022
 - λ_inv - inverse of the size distribution parameter
 - k - size distribution moment for which we compute the bulk fall speed

Returns the addends of the bulk fall speed of rain or ice particles
following Chen et al 2022 DOI: 10.1016/j.atmosres.2022.106171 in [m/s].
Assuming exponential size distribution and hence μ=0.
"""
function Chen2022_exponential_pdf(a::FT, b::FT, c::FT, λ_inv::FT, k::Int) where {FT}
    μ = 0 # Exponential instead of gamma distribution
    δ = FT(μ + k + 1)
    return a * exp(δ * log(1 / λ_inv) - (b + δ) * log(1 / λ_inv + c)) * SF.gamma(b + δ) / SF.gamma(δ)
end

"""
    particle_terminal_velocity(velocity_params, ρₐ)
    particle_terminal_velocity(velocity_params, ρₐ, ρᵢ)

Return a function `v_term(D)` that computes the particle terminal velocity

# Arguments
- `velocity_params`: a struct with terminal velocity parameters from Chen 2022
- `ρₐ`: air density [kg/m³]
- `ρᵢ`: apparent density of ice particles [kg/m³],
    only used for [`CMP.Chen2022VelTypeSmallIce`](@ref) and [`CMP.Chen2022VelTypeLargeIce`](@ref)

# Returns
- `v_term(D)`: The terminal velocity of a particle as a function of its size (diameter, `D`)

Needed for numerical integrals in the P3 scheme.

!!! note
    We use the same terminal velocity parametrization for cloud and rain water.
"""
function particle_terminal_velocity(velocity_params::CMP.TerminalVelocityType, ρs...)
    (ai, bi, ci) = Chen2022_vel_coeffs(velocity_params, ρs...)
    v_terms = Chen2022_monodisperse_pdf.(ai, bi, ci)  # tuple of functions
    v_term(D) = sum(@. D |> v_terms)
    return v_term
end


"""
    volume_sphere_D(D)

Calculate the volume of a sphere with diameter D.

```math
V = D^3 * π / 6
```

See also [`volume_sphere_R`](@ref).
"""
volume_sphere_D(D) = D^3 * π / 6

"""
    volume_sphere_R(R)

Calculate the volume of a sphere with radius R.

```math
V = (2R)^3 * π / 6
```

See also [`volume_sphere_D`](@ref).
"""
volume_sphere_R(R) = volume_sphere_D(2R)

"""
    ventilation_factor(vent, aps, v_term)

Returns a function that computes the ventilation factor for a particle as a function of its diameter, `D`.

The ventilation factor parameterizes the increase in the mass and heat exchange for falling particles.

# Arguments
- `vent`: Ventilation parameterization constants, [`CMP.VentilationFactor`](@ref)
- `aps`: Parameters with air properties, [`CMP.AirProperties`](@ref)
- `v_term`: A function `v_term(D)` that returns the terminal velocity of a particle with diameter `D`

# Returns
- `F_v(D)`: The ventilation factor as a function of diameter, `D`

See e.g. [SeifertBeheng2006](@cite) Eq. (24), or [`CMP.VentilationFactor`](@ref),
for the definition of the ventilation factor.
"""
function ventilation_factor(vent, aps, v_term)
    (; aᵥ, bᵥ) = vent
    (; ν_air, D_vapor) = aps
    N_sc = ν_air / D_vapor           # Schmidt number
    N_Re(D) = D * v_term(D) / ν_air  # Reynolds number
    F_v(D) = aᵥ + bᵥ * ∛(N_sc) * √(N_Re(D))  # Ventilation factor
    return F_v
end

end # module end
