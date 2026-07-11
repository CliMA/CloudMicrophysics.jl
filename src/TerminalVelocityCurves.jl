# Chen 2022 terminal-velocity curves, part of the `Common` module (see Common.jl).

"""
    Chen2022_vel_coeffs(coeffs, ρₐ)
    Chen2022_vel_coeffs(coeffs, ρₐ, ρᵢ)
    Chen2022_vel_coeffs(v::Chen2022VelocityCurve)

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

# Returns
- `terms::NTuple{N, NTuple{3, FT}}`: the per-term coefficient triples
  `(aₖ, bₖ, cₖ)` of the velocity curve. The method on a
  [`Chen2022VelocityCurve`](@ref) returns the stored `terms` directly.
"""
@inline function Chen2022_vel_coeffs(coeffs::CMP.Chen2022VelTypeRain, ρₐ)
    (; ρ0, a, a3_pow, b, b_ρ, c) = coeffs
    ρₐ = max(ρₐ, zero(ρₐ))
    # Table B1
    q = exp(ρ0 * ρₐ)
    ai = (a[1] * q, a[2] * q, a[3] * q * ρₐ^a3_pow)
    bi = (b[1] - b_ρ * ρₐ, b[2] - b_ρ * ρₐ, b[3] - b_ρ * ρₐ)
    ci = (c[1], c[2], c[3])
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return map(tuple, aiu, bi, ciu)
end

@inline function Chen2022_vel_coeffs(coeffs::CMP.Chen2022VelTypeSmallIce, ρₐ, ρᵢ)
    FT = eltype(coeffs)
    (; A, B, C, E, F, G) = coeffs
    ρₐ = max(ρₐ, zero(ρₐ))
    # Table B3 - cache sqrt for reuse
    log_ρᵢ = log(ρᵢ)
    sqrt_ρᵢ = sqrt(ρᵢ)
    As = A[2] * log_ρᵢ^2−A[3] * log_ρᵢ + A[1]
    Bs = 1 / (B[1] + B[2] * log_ρᵢ + B[3] / sqrt_ρᵢ)
    Cs = C[1] + C[2] * exp(C[3] * ρᵢ) + C[4] * sqrt_ρᵢ
    Es = E[1] - E[2] * log_ρᵢ^2 + E[3] * sqrt_ρᵢ
    Fs = -exp(F[1] - F[2] * log_ρᵢ^2 + F[3] * log_ρᵢ)
    Gs = 1 / (G[1] + G[2] / log_ρᵢ - G[3] * log_ρᵢ / ρᵢ)
    # Table B2
    ai = (Es * ρₐ^As, Fs * ρₐ^As)
    bi = (Bs + ρₐ * Cs, Bs + ρₐ * Cs)
    ci = (FT(0), Gs)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return map(tuple, aiu, bi, ciu)
end

@inline function Chen2022_vel_coeffs(coeffs::CMP.Chen2022VelTypeLargeIce, ρₐ, ρᵢ)
    FT = eltype(coeffs)
    (; A, B, C, E, F, G, H) = coeffs
    ρₐ = max(ρₐ, zero(ρₐ))
    # Table B5 - cache sqrt for reuse
    log_ρᵢ = log(ρᵢ)
    sqrt_ρᵢ = sqrt(ρᵢ)
    Al = A[1] + A[2] * log_ρᵢ + A[3] / (ρᵢ * sqrt_ρᵢ)  # ρᵢ^(-3/2) = 1/(ρᵢ * sqrt(ρᵢ))
    Bl = exp(B[1] + B[2] * log_ρᵢ^2 + B[3] * log_ρᵢ)
    Cl = exp(C[1] + C[2] / log_ρᵢ + C[3] / ρᵢ)
    El = E[1] + E[2] * log_ρᵢ * sqrt_ρᵢ + E[3] * sqrt_ρᵢ
    Fl = F[1] + F[2] * log_ρᵢ - exp(log(-F[3]) - ρᵢ)
    Gl = 1 / (G[1] + G[2] * log_ρᵢ * sqrt_ρᵢ + G[3] / sqrt_ρᵢ)
    Hl = H[1] + H[2] * ρᵢ^2 * sqrt_ρᵢ + exp(log(-H[3]) - ρᵢ)  # ρᵢ^(5/2) = ρᵢ^2 * sqrt(ρᵢ)
    # Table B4
    ai = (Bl * ρₐ^Al, El * ρₐ^Al * exp(Hl * ρₐ))
    bi = (Cl, Fl)
    ci = (FT(0), Gl)
    # unit conversions
    aiu = ai .* 1000 .^ bi
    ciu = ci .* 1000
    return map(tuple, aiu, bi, ciu)
end

"""
    Chen2022VelocityCurve(velocity_params, ρₐ)
    Chen2022VelocityCurve(velocity_params, ρₐ, ρᵢ)
    Chen2022VelocityCurve(terms)

Callable holding the Chen 2022 terminal-velocity coefficients as
`terms::NTuple{N, NTuple{3, FT}}` of `(aₖ, bₖ, cₖ)`, from
[`Chen2022_vel_coeffs`](@ref) evaluated at air density `ρₐ` (and apparent ice
density `ρᵢ` for the small/large-ice tables). Evaluating it at a diameter `D`
returns `∑ₖ aₖ D^bₖ exp(-cₖ D)` [m/s].
"""
struct Chen2022VelocityCurve{N, FT}
    terms::NTuple{N, NTuple{3, FT}}
    # The implicit constructor leaves `FT` unbound for empty coefficient tuples
    Chen2022VelocityCurve{N, FT}(terms) where {N, FT} = new{N, FT}(terms)
end
function Chen2022VelocityCurve(velocity_params::CMP.TerminalVelocityType, ρₐ)
    return Chen2022VelocityCurve(Chen2022_vel_coeffs(velocity_params, ρₐ))
end
function Chen2022VelocityCurve(velocity_params::CMP.TerminalVelocityType, ρₐ, ρᵢ)
    return Chen2022VelocityCurve(Chen2022_vel_coeffs(velocity_params, ρₐ, ρᵢ))
end
function Chen2022VelocityCurve(terms::NTuple{N, NTuple{3, Any}}) where {N}
    FT = promote_type(map(t -> promote_type(map(typeof, t)...), terms)...)
    return Chen2022VelocityCurve{N, FT}(map(t -> map(FT, t), terms))
end
@inline (v::Chen2022VelocityCurve)(D) =
    unrolled_sum(t -> t[1] * D^t[2] * exp(-t[3] * D), v.terms)

@inline Chen2022_vel_coeffs(v::Chen2022VelocityCurve) = v.terms

"""
    Chen2022_monodisperse_pdf(a, b, c)

# Arguments
 - `a`, `b`, `c`: free parameters defined in [Chen2022](@cite)

# Returns
 - `pdf(D)`: The monodisperse particle distribution function as a function of diameter, `D`, in [m/s].
"""
@inline function Chen2022_monodisperse_pdf(a, b, c)
    # Fuse D^b * exp(-c*D) = exp(b*log(D) - c*D) into a single exp (D^b alone is a
    # pow = log+exp for runtime-float b, so this drops one exp per term). Evaluated
    # per quadrature node in the P3 velocity integrands, the hottest GPU path.
    return pdf(D) = a * exp(muladd(b, log(D), -c * D))
end

"""
    Chen2022_exponential_pdf(a, b, c, λ_inv, k)

Returns the addends of the bulk fall speed of rain or ice particles
following Chen et al. (2022), https://doi.org/10.1016/j.atmosres.2022.106171.
Assumes exponential size distribution (μ=0).

# Arguments
- `a`, `b`, `c`: free parameters defined in Chen et al. (2022)
 - `λ_inv`: inverse of the size distribution parameter [m]
 - `k`: size distribution moment for which we compute the bulk fall speed

# Returns
- Bulk fall speed component [m/s]
"""
@inline function Chen2022_exponential_pdf(a, b, c, λ_inv, k::Int)
    FT = UT.promote_typeof(a, b, c, λ_inv)
    # μ = 0 for exponential distribution, δ = k + 1
    δ = FT(k + 1)
    # Γ(δ) = k! for integer δ = k + 1. `UT.fac` computes the product directly;
    # `Base.factorial` reads a host-memory table, which is not GPU-compatible.
    gamma_delta = FT(UT.fac(k))
    return a * exp(-δ * log(λ_inv) - (b + δ) * log(1 / λ_inv + c)) * SF.gamma(b + δ) / gamma_delta
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
particle_terminal_velocity(velocity_params::CMP.TerminalVelocityType, ρs...) =
    Chen2022VelocityCurve(velocity_params, ρs...)

"""
    particle_terminal_velocity(velocity_params::CMP.StokesRegimeVelType{FT}, ρ::FT)

 - `velocity_params` - set with free parameters
 - `ρ` - air density

Returns a function `v_term(D)` that computes the analytical fall speed of a cloud droplet as a function of
its size (diameter, `D`) in the Stokes regime (Re < 1)
"""
function particle_terminal_velocity(velocity_params::CMP.StokesRegimeVelType, ρ)
    (; ρw, grav, ν_air) = velocity_params
    FT = eltype(ρ)
    terminal_velocity_prefactor = FT(1 / 18) * (ρw / ρ - 1) * grav / ν_air
    v_term(D) = terminal_velocity_prefactor * D^2
    return v_term
end
