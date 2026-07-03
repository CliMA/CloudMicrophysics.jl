"""
    Module for functions shared by different parameterizations.
"""
module Common

import SpecialFunctions as SF
import LogExpFunctions as LEF
using UnrolledUtilities

import ..Parameters as CMP
import ..ThermodynamicsInterface as TDI
import ..Utilities as UT
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

Utility function combining thermal conductivity and vapor diffusivity effects
for vapor to liquid phase change.

Includes guards against division by zero using ϵ_numerics thresholds for
numerical robustness.

# Arguments
- `air_props`: air parameters struct (contains `K_therm`, `D_vapor`)
- `tps`: thermodynamics parameters struct
- `T`: air temperature [K]

# Returns
- G function for liquid [kg/m²/s/Pa]

# Notes
- Division by `K_therm`, `D_vapor`, and `p_vs` are guarded with `max(value, UT.ϵ_numerics(FT))`
- This prevents NaN/Inf propagation from near-zero denominators
"""
@inline function G_func_liquid(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TDI.PS,
    T,
) where {FT}
    R_v = TDI.Rᵥ(tps)
    L = TDI.Lᵥ(tps, T)
    p_vs = TDI.saturation_vapor_pressure_over_liquid(tps, T)

    # Guard against division by zero
    p_vs_safe = max(p_vs, UT.ϵ_numerics(FT))
    D_vapor_safe = max(D_vapor, UT.ϵ_numerics(FT))
    K_therm_safe = max(K_therm, UT.ϵ_numerics(FT))

    return 1 /
           (L / K_therm_safe / T * (L / R_v / T - 1) + R_v * T / D_vapor_safe / p_vs_safe)
end

"""
    G_func_ice(air_props, tps, T)

Utility function combining thermal conductivity and vapor diffusivity effects
for vapor to ice phase change.

Includes guards against division by zero using ϵ_numerics thresholds for
numerical robustness.

# Arguments
- `air_props`: air parameters struct (contains `K_therm`, `D_vapor`)
- `tps`: thermodynamics parameters struct
- `T`: air temperature [K]

# Returns
- G function for ice [kg/m²/s/Pa]

# Notes
- Division by `K_therm`, `D_vapor`, and `p_vs` are guarded with `max(value, UT.ϵ_numerics(FT))`
- This prevents NaN/Inf propagation from near-zero denominators
"""
@inline function G_func_ice(
    (; K_therm, D_vapor)::CMP.AirProperties{FT},
    tps::TDI.PS,
    T,
) where {FT}
    R_v = TDI.Rᵥ(tps)
    L = TDI.Lₛ(tps, T)
    p_vs = TDI.saturation_vapor_pressure_over_ice(tps, T)

    # Guard against division by zero
    p_vs_safe = max(p_vs, UT.ϵ_numerics(FT))
    D_vapor_safe = max(D_vapor, UT.ϵ_numerics(FT))
    K_therm_safe = max(K_therm, UT.ϵ_numerics(FT))

    return 1 /
           (L / K_therm_safe / T * (L / R_v / T - 1) + R_v * T / D_vapor_safe / p_vs_safe)
end

"""
    A Heaviside step function
"""
@inline heaviside(x::FT) where {FT} = FT(x > 0)

"""
    logistic_function(x, x_0, k)

Returns the value of the logistic function for smooth transitioning at thresholds.
This is a normalized curve changing from 0 to 1 while x varies from 0 to Inf (for positive k).

For x < 0 the value at x = 0 is returned. For x_0 = 0, H(x) is returned.

# Arguments
- `x`: independent variable
- `x_0`: threshold value for x
- `k`: growth rate of the curve, characterizing steepness of the transition

# Returns
- Logistic function value (dimensionless, range [0,1])
"""
@inline function logistic_function(x::FT, x_0::FT, k::FT) where {FT}
    # GPU-compatible implementation with edge case handling
    x = max(FT(0), x)

    # Edge cases: x ≈ 0 → return 0, x_0 ≈ 0 → return 1 (if x > 0)
    x_safe = max(x, UT.ϵ_numerics(FT))
    x_0_safe = max(x_0, UT.ϵ_numerics(FT))

    # σ(z) = 1 / (1 + exp(-z)) = exp(-log1pexp(-z))
    z = k * (x_safe / x_0_safe - x_0_safe / x_safe)
    result = exp(-LEF.log1pexp(-z))

    # Handle edge cases with ifelse (branchless)
    return ifelse(x < UT.ϵ_numerics(FT), FT(0), ifelse(x_0 < UT.ϵ_numerics(FT), FT(1), result))
end

"""
    logistic_function_integral(x, x_0, k)

Returns the value of the indefinite integral of the logistic function, for smooth transitioning
of piecewise linear profiles at thresholds.

This curve smoothly transitions from y = 0 for 0 < x < x_0 to y = x - x_0 for x_0 < x.

# Arguments
- `x`: independent variable
- `x_0`: threshold value for x
- `k`: growth rate of the logistic function, characterizing steepness of the transition

# Returns
- Integral of logistic function (same units as x)
"""
@inline function logistic_function_integral(x, x_0, k)
    FT = UT.promote_typeof(x, x_0, k)
    # Branchless GPU-compatible implementation
    x = max(FT(0), x)
    x_safe = max(x, UT.ϵ_numerics(FT))
    x_0_safe = max(x_0, UT.ϵ_numerics(FT))

    # translation of the curve in x and y to enforce zero at x = 0
    # Using log1mexp for numerical stability: log1mexp(x) = log(1 - exp(x))
    trnslt = -LEF.log1mexp(-k) / k
    kt = k * (x_safe / x_0_safe - 1 + trnslt)
    # log1pexp handles all kt values correctly
    result = (LEF.log1pexp(kt) / k - trnslt) * x_0_safe

    # Handle edge cases: if x ≈ 0, return 0; if x_0 ≈ 0, return x
    return ifelse(x < UT.ϵ_numerics(FT), FT(0), ifelse(x_0 < UT.ϵ_numerics(FT), x, result))
end

"""
    H2SO4_soln_saturation_vapor_pressure(prs, x, T)

Returns the saturation vapor pressure above a sulphuric acid solution droplet.

# Arguments
- `prs`: H2SO4 solution free parameters struct
- `x`: weight percent sulphuric acid (dimensionless, e.g., 0.1 for 10%)
- `T`: air temperature [K]

# Returns
- Saturation vapor pressure [Pa]
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

Returns water activity of H2SO4 containing droplet.

# Arguments
- `H2SO4_prs`: H2SO4 solution free parameters struct
- `tps`: thermodynamics parameters struct
- `x`: weight percent sulphuric acid (dimensionless, e.g., 0.1 for 10%)
- `T`: air temperature [K]

# Returns
- Water activity (dimensionless)
"""
function a_w_xT(H2SO4_prs::HPS, tps::TDI.PS, x::FT, T::FT) where {FT}

    p_sol = H2SO4_soln_saturation_vapor_pressure(H2SO4_prs, x, T)
    p_sat = TDI.saturation_vapor_pressure_over_liquid(tps, T)

    return p_sol / p_sat
end

"""
    a_w_eT(tps, e, T)

Returns water activity of pure water droplet.
Valid when droplet is in equilibrium with surroundings.

# Arguments
- `tps`: thermodynamics parameters struct
- `e`: partial pressure of water [Pa]
- `T`: air temperature [K]

# Returns
- Water activity (dimensionless, equivalent to relative humidity)
"""
function a_w_eT(tps::TDI.PS, e::FT, T::FT) where {FT}
    # RH
    return e / TDI.saturation_vapor_pressure_over_liquid(tps, T)
end

"""
    a_w_ice(tps, T)

Returns water activity of ice.

# Arguments
- `tps`: thermodynamics parameters struct
- `T`: air temperature [K]

# Returns
- Water activity of ice (dimensionless)
"""
function a_w_ice(tps::TDI.PS, T::FT) where {FT}

    return TDI.saturation_vapor_pressure_over_ice(tps, T) /
           TDI.saturation_vapor_pressure_over_liquid(tps, T)
end

include("TerminalVelocityCurves.jl")

"""
    volume_sphere_D(D)

Calculate the volume of a sphere with diameter D.

```math
V = D^3 * π / 6
```

See also [`volume_sphere_R`](@ref).
"""
@inline volume_sphere_D(D) = D^3 * π / 6

"""
    volume_sphere_R(R)

Calculate the volume of a sphere with radius R.

```math
V = (2R)^3 * π / 6
```

See also [`volume_sphere_D`](@ref).
"""
@inline volume_sphere_R(R) = volume_sphere_D(2R)

"""
    ventilation_factor(vent, aps, v_term)

Returns a function that computes the ventilation factor for a particle as a function of its diameter.

The ventilation factor parameterizes the increase in the mass and heat exchange for falling particles.
See e.g., Seifert and Beheng (2006), https://doi.org/10.1007/s00703-005-0112-4, Eq. (24).

# Arguments
- `vent`: ventilation parameterization constants (contains `aᵥ`, `bᵥ`)
- `aps`: air properties parameters (contains `ν_air`, `D_vapor`)
- `v_term`: function `v_term(D)` that returns terminal velocity [m/s] for diameter D [m]

# Returns
- `F_v(D)`: ventilation factor function (dimensionless)
"""
@inline function ventilation_factor(vent, aps, v_term)
    (; aᵥ, bᵥ) = vent
    (; ν_air, D_vapor) = aps
    N_sc = ν_air / D_vapor           # Schmidt number
    cbrt_N_sc = cbrt(N_sc)           # loop-invariant over D; hoist out of F_v
    N_Re(D) = D * v_term(D) / ν_air  # Reynolds number
    F_v(D) = aᵥ + bᵥ * cbrt_N_sc * sqrt(N_Re(D))  # Ventilation factor
    return F_v
end

end # module end
