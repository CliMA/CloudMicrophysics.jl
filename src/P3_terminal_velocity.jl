###
### P3 terminal velocity functions
###

# `ParametersP3` stores the aspect-ratio choice, but the functor methods need `ϕᵢ`,
# which is defined in this module, so the methods are defined here.
@inline (::CMP.Oblate)(particles, D) = cbrt(ϕᵢ(particles, D))
@inline (::CMP.NoAspectRatio)(particles, D) = one(D)

# Callable returned by `ice_particle_terminal_velocity`: piecewise small/large-ice
# Chen 2022 velocity scaled by the aspect-ratio factor from `particles.params`.
struct P3IceParticleVelocityFunctor{FT, VS, VL, S} <: Function
    v_term_small::VS
    v_term_large::VL
    D_cutoff::FT
    particles::S
end
@inline function (f::P3IceParticleVelocityFunctor)(D)
    vₜ = ifelse(D <= f.D_cutoff, f.v_term_small(D), f.v_term_large(D))
    return vₜ * f.particles.params.aspect_ratio(f.particles, D)
end

"""
    velocity_breakpoints(v_term)

Return the diameters where the terminal-velocity closure `v_term` changes functional form.

Integrals with `v_term` in the integrand place these on subinterval boundaries,
see [`velocity_integral_bounds`](@ref).
"""
velocity_breakpoints(f::P3IceParticleVelocityFunctor) = (f.D_cutoff,)
velocity_breakpoints(::CO.Chen2022VelocityCurve) = ()

"""
    ice_particle_terminal_velocity(velocity_params, particles, ρₐ)

Return the Chen 2022 terminal velocity `v_term(D)` of an ice particle of maximum dimension `D`.

The velocity is scaled by the aspect-ratio factor selected by `particles.params.aspect_ratio`.

# Arguments
 - `velocity_params`: A [`CMP.Chen2022VelType`](@ref)
 - `particles`: A [`P3Particles`](@ref)
 - `ρₐ`: Air density [kg/m³]
"""
@inline function ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    (; small_ice, large_ice) = velocity_params
    D_cutoff = small_ice.cutoff
    (; ρ_i) = particles.params
    v_term_small = CO.particle_terminal_velocity(small_ice, ρₐ, ρ_i)
    v_term_large = CO.particle_terminal_velocity(large_ice, ρₐ, ρ_i)
    return P3IceParticleVelocityFunctor(v_term_small, v_term_large, D_cutoff, particles)
end

# Integrand of `ice_terminal_velocities`: `(n, n v, n m, n m v)` at diameter `D`
struct P3VelocityIntegrand{N, V, P} <: Function
    n::N
    v_term::V
    particles::P
end
@inline function (f::P3VelocityIntegrand)(D)
    n = f.n(D)
    v = f.v_term(D)
    m = ice_mass(f.particles, D)
    return SA.SVector(n, n * v, n * m, n * m * v)
end

"""
    ice_terminal_velocities(velocity_params, state, ρₐ; [p])

Return the number- and mass-weighted mean terminal velocities of the ice size distribution.

The velocities are

- `v_n = ∫ n(D) v(D) dD / ∫ n(D) dD`, the number-weighted velocity [m/s],
- `v_m = ∫ n(D) m(D) v(D) dD / ∫ n(D) m(D) dD`, the mass-weighted velocity [m/s].

All four integrals come from one quadrature over the size distribution, so each velocity is a
weighted mean of `v(D)` over the quadrature nodes. Both are zero without ice number.

# Arguments
- `velocity_params`: [`CMP.Chen2022VelType`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]

# Keyword arguments
- `p`: tolerance parameter for the integral bounds, see [`integral_bounds`](@ref). Default is 1e-6.

# Returns
- `(; v_n, v_m)`: the number- and mass-weighted velocities [m/s]
"""
function ice_terminal_velocities(velocity_params, state, ρₐ; p = 1e-6)
    (; particles, integration, dist) = state
    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    n = DT.size_distribution(dist)
    bnds = velocity_integral_bounds(particles, dist, v_term; p)
    (∫n, ∫nv, ∫nm, ∫nmv) = integrate(P3VelocityIntegrand(n, v_term, particles), bnds, integration)
    return (; v_n = UT.guarded_quotient(∫nv, ∫n), v_m = UT.guarded_quotient(∫nmv, ∫nm))
end
