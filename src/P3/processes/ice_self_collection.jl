"""
    ice_self_collection(vel, state, ρₐ)

Compute the ice self-collection (aggregation) rate.

Self-collection decreases the ice number concentration and leaves mass, rime mass, and rime volume
unchanged.

The method of `ice_self_collection(vel, integration, state, ρₐ)` for the type of
`state.integration` evaluates the collision integral.

# Arguments
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection `[1/m³/s]` (always positive or zero, represents a loss rate), the fractional rate of the size distribution times the prognostic ice number
"""
@inline ice_self_collection(vel, state, ρₐ) = ice_self_collection(vel, state.integration, state, ρₐ)

@inline function ice_self_collection(vel, integration::QuadratureRule, state, ρₐ)
    (; particles, dist) = state
    n_i = DT.size_distribution(dist)
    v_ice = ice_particle_terminal_velocity(vel, particles, ρₐ)

    p = eps(one(ρₐ))
    ice_bounds = velocity_integral_bounds(particles, dist, v_ice; p)

    function inner_integral(D_1)
        v1 = v_ice(D_1)
        r1 = r_eff_ice(particles, D_1)
        # Volumetric collision rate integrand: E · K · |v₁ − v₂| · n(D₂)
        # where E = 1 (collision efficiency),
        #       K = π (r₁ + r₂)² is the geometric collision cross-section,
        #       |v₁ − v₂| is the differential sedimentation speed,
        #       r = √(A/π) is the effective radius from projected ice area A.
        integrand = D_2 -> begin
            v2 = v_ice(D_2)
            r2 = r_eff_ice(particles, D_2)
            K = π * (r1 + r2)^2                    # collision cross section
            return K * abs(v1 - v2) * n_i(D_2)     # E = 1 implied
        end
        # Integrate over D_2 ≥ D_1, so the triangle counts each unordered particle pair once.
        # The |v1 - v2| cusp at D_2 = D_1 is then an end point of the inner integral.
        rate_at_D1 = integrate(integrand, map(D -> max(D, D_1), ice_bounds), integration)
        return rate_at_D1 * n_i(D_1)
    end

    total_rate = integrate(inner_integral, ice_bounds, integration)

    dNdt = total_rate * UT.guarded_quotient(state.ρn_ice, dist.ρn)
    return (; dNdt)
end
