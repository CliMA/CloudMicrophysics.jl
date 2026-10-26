"""
    ice_point(params, vel, generator, (logλ, w, ρ_rim, logρₐ))

Return the integrals of the ice table at a point, per unit number, or per unit number squared for the self-collection rate.

The value is an `SVector` of the ventilation moment `M_v` of [`ventilation_moment`](@ref)
[m² s^(-1/2)], the self-collection rate of [`node_self_collection`](@ref) [m³/s], and the
number- and mass-weighted velocities of [`ice_terminal_velocities`](@ref) [m/s], for the state of
[`table_state`](@ref) and the rule `generator`. The table stores their logarithms.
"""
function ice_point(params, vel, generator, (logλ, w, ρ_rim, logρₐ))
    state = table_state(params, generator, logλ, w, ρ_rim)
    ρₐ = exp(logρₐ)
    (; v_n, v_m) = ice_terminal_velocities(vel, state, ρₐ)
    return SA.SVector(ventilation_moment(vel, state, ρₐ), node_self_collection(vel, generator, state, ρₐ), v_n, v_m)
end

"""
    ventilation_moment(vel, state, ρₐ)

Return the ventilation moment `M_v = ∫ D^(3/2) v(D)^(1/2) N′(D) dD` [m⁻¹ s^(-1/2)] by quadrature.

With the coefficients `(a, b)` of [`CO.ventilation_coefficients`](@ref), the ventilation integral of
[`ice_ventilation_integral`](@ref) is `a M_1 + b M_v`, with the first moment `M_1 = ∫ D N′(D) dD`.
"""
function ventilation_moment(vel, state, ρₐ)
    (; particles, integration, dist) = state
    v_term = ice_particle_terminal_velocity(vel, particles, ρₐ)
    N′ = size_distribution(dist)
    bnds = velocity_integral_bounds(particles, dist, v_term; p = eltype(state)(1e-6))
    return integrate(D -> (DN′ = D * N′(D); DN′ * sqrt(D * v_term(D))), bnds, integration)
end

"""
    node_self_collection(vel, generator, state, ρₐ)

Return the self-collection rate of [`ice_self_collection`](@ref) per unit number squared, with the node values of each subinterval evaluated once.

For a `GaussLegendre` rule, the quadrature is that of [`ice_self_collection`](@ref), with the
values at the nodes of the subintervals after the one of the outer diameter reused.

# Arguments
- `vel`: [`CMP.Chen2022VelType`](@ref)
- `generator`: `QuadratureRule`
- `state`: [`P3State`](@ref) of unit ice number
- `ρₐ`: air density [kg/m³]
"""
node_self_collection(vel, generator, state, ρₐ) = ice_self_collection(vel, state, ρₐ).dNdt
function node_self_collection(vel, generator::GaussLegendre{<:Any, N}, state, ρₐ) where {N}
    (; particles, dist) = state
    n_i = DT.size_distribution(dist)
    v_ice = ice_particle_terminal_velocity(vel, particles, ρₐ)
    bounds = velocity_integral_bounds(particles, dist, v_ice; p = eps(one(ρₐ)))
    M = length(bounds) - 1
    # The diameters of the nodes of each subinterval, and the weights times the size distribution,
    # the effective radius and the terminal velocity at them
    nodes = ntuple(Val(M)) do i
        h = (bounds[i + 1] - bounds[i]) / 2
        c = (bounds[i] + bounds[i + 1]) / 2
        D = map(y -> h * y + c, generator.nodes)
        return (;
            D,
            wn = map((w, D) -> h * w * n_i(D), generator.weights, D),
            r = map(D -> r_eff_ice(particles, D), D),
            v = map(v_ice, D),
        )
    end
    pair(r₁, v₁, r₂, v₂) = π * (r₁ + r₂)^2 * abs(v₁ - v₂)
    total = zero(eltype(first(nodes).wn))
    for i in 1:M, a in 1:N
        (; D, r, v) = nodes[i]
        D₁, r₁, v₁ = D[a], r[a], v[a]
        inner =
            integrate(D₂ -> pair(r₁, v₁, r_eff_ice(particles, D₂), v_ice(D₂)) * n_i(D₂), D₁, bounds[i + 1], generator)
        for j in (i + 1):M, b in 1:N
            node = nodes[j]
            inner += node.wn[b] * pair(r₁, v₁, node.r[b], node.v[b])
        end
        total += nodes[i].wn[a] * inner
    end
    return total * UT.guarded_quotient(state.ρn_ice, dist.ρn)
end

"""
    table_state(params, generator, logλ, w, ρ_rim)

Return the [`P3State`](@ref) of unit ice number with the slope `logλ` [log(1/m)], the rime coordinate `w` and the rime density `ρ_rim` [kg/m³], integrated with the rule `generator`.
"""
function table_state(params, generator, logλ, w, ρ_rim)
    F_rim = rime_mass_fraction_at(w)
    particles = P3Particles(params, F_rim, ρ_rim)
    dist = P3Distribution(params.slope, one(logλ), logλ)
    x̄ = exp(logLdivN(particles, logλ))
    return P3State(particles, generator, x̄, one(x̄), F_rim * x̄ / ρ_rim, dist)
end

"""
    ice_table_key(state, ρₐ)

Return the coordinates `(logλ, w, ρ_rim, logρₐ)` of `state` in the ice tables.

`logλ` is a primal value, as for the quadrature, and the other coordinates keep their partials.
"""
@inline function ice_table_key(state, ρₐ)
    (; particles, dist) = state
    return (
        primal_value(dist.logλ),
        rime_coordinate(particles.F_rim),
        particles.ρ_rim,
        log(max(ρₐ, floatmin(primal_value(ρₐ)))),
    )
end

"""
    tabulated_ice_integrals(tables, state, ρₐ)

Return the `NamedTuple` `(; M_v, self_collection, v_n, v_m)` of the ice integral table of `tables` at the coordinates of `state`, see [`ice_point`](@ref).

The integrals are per unit number, or per unit number squared for `self_collection`.
"""
@inline function tabulated_ice_integrals(tables, state, ρₐ)
    (M_v, self_collection, v_n, v_m) = exp.(lookup(tables.groups.ice, ice_table_key(state, ρₐ)...))
    return (; M_v, self_collection, v_n, v_m)
end

@inline function ice_ventilation_integral(velocity_params, aps, state::TabulatedP3State, ρₐ)
    (; ρn, μ, logλ) = state.dist
    (a, b) = CO.ventilation_coefficients(state.particles.params.vent, aps)
    M_1 = ρn * (μ + 1) * exp(-logλ)
    M_v = tabulated_ice_integrals(state.integration, state, ρₐ).M_v * ρn
    return a * M_1 + b * M_v
end

@inline function ice_self_collection(vel, state::TabulatedP3State, ρₐ)
    K = tabulated_ice_integrals(state.integration, state, ρₐ).self_collection
    return (; dNdt = K * state.dist.ρn * state.ρn_ice)
end

function ice_terminal_velocities(velocity_params, state::TabulatedP3State, ρₐ; p = eltype(state)(1e-6))
    (; v_n, v_m) = tabulated_ice_integrals(state.integration, state, ρₐ)
    present = primal_value(state.dist.ρn) > 0
    return (; v_n = ifelse(present, v_n, zero(v_n)), v_m = ifelse(present, v_m, zero(v_m)))
end

"""
    ∫liquid_ice_collisions(
        psd_c, psd_r, aps, tps, vel, m_liq, state::TabulatedP3State,
        L_c, N_c, L_r, N_r, ρₐ, T,
    )

Return the liquid-ice collision rates of [`∫liquid_ice_collisions`](@ref) for a tabulated state, by quadrature with the rule of its tables.
"""
@inline function ∫liquid_ice_collisions(
    psd_c, psd_r, aps, tps, vel, m_liq, state::TabulatedP3State,
    L_c, N_c, L_r, N_r, ρₐ, T,
)
    (; particles, integration, ρq_ice, ρn_ice, ρb_rim, dist) = state
    quadrature = P3State(particles, integration.generator, ρq_ice, ρn_ice, ρb_rim, dist)
    return ∫liquid_ice_collisions(psd_c, psd_r, aps, tps, vel, m_liq, quadrature, L_c, N_c, L_r, N_r, ρₐ, T)
end
