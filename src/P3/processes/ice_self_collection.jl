"""
    ice_self_collection(vel, state, ρₐ)

Compute the ice self-collection (aggregation) rate.

Self-collection decreases the ice number concentration and leaves mass, rime mass, and rime volume
unchanged.

# Arguments
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection `[1/m³/s]` (always positive or zero, represents a loss rate), the fractional rate of the size distribution times the prognostic ice number
"""
@inline function ice_self_collection(vel, state, ρₐ)
    (; particles, integration, dist) = state
    n_i = DT.size_distribution(dist)
    v_ice = ice_particle_terminal_velocity(vel, particles, ρₐ)

    p = eps(one(ρₐ))
    ice_bounds = integral_bounds(particles, dist; p)

    function inner_integral(D_1)
        v1 = v_ice(D_1)
        r1 = sqrt(ice_area(particles, D_1) / π)
        # Volumetric collision rate integrand: E · K · |v₁ − v₂| · n(D₂)
        # where E = 1 (collision efficiency),
        #       K = π (r₁ + r₂)² is the geometric collision cross-section,
        #       |v₁ − v₂| is the differential sedimentation speed,
        #       r = √(A/π) is the effective radius from projected ice area A.
        integrand = D_2 -> begin
            v2 = v_ice(D_2)
            r2 = sqrt(ice_area(particles, D_2) / π)
            K = π * (r1 + r2)^2                    # collision cross section
            return K * abs(v1 - v2) * n_i(D_2)     # E = 1 implied
        end
        # Split the inner integral at the |v1 - v2| cusp (D_2 = D_1, where the relative
        # fall speed vanishes and the integrand has a kink). Each half is then smooth, so
        # the quadrature converges like the other (cusp-free) P3 integrals rather than
        # being cusp-limited, letting a much lower node count reach target accuracy.
        D_lo, D_hi = first(ice_bounds), last(ice_bounds)
        rate_at_D1 = integrate(integrand, (D_lo, D_1), integration) + integrate(integrand, (D_1, D_hi), integration)
        return rate_at_D1 * n_i(D_1)
    end

    total_rate = integrate(inner_integral, ice_bounds, integration)

    # The 0.5 factor accounts for double-counting in self-collection
    FT = eltype(state)
    dNdt = FT(0.5) * total_rate * UT.guarded_quotient(state.ρn_ice, dist.ρn)
    return (; dNdt)
end
