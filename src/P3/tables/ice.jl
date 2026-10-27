"""
    IcePoint{P, V, G}

The point function of the ice integral table, the integrals at a point `(logλ, w, ρ_rim, logρₐ)`.

Each integral is a size-distribution integral of the ice per unit number, or per unit number
squared for a collision rate. Per unit number, an integral depends only on the shape of the
distribution, which `logλ` sets, on the particle properties, which the rime coordinate `w` (see
[`rime_coordinate`](@ref)) and the rime density `ρ_rim` set, and on the terminal velocity, which
the log air density `logρₐ` sets. The value is an `SVector` of

1. `M_v`, the ventilation moment of [`ice_ventilation_moments`](@ref) [m² s^(-1/2)],
2. `self_collection`, the self-collection rate of [`ice_self_collection`](@ref) [m³/s], see
   [`node_self_collection`](@ref),
3. `v_n`, the number-weighted velocity of [`ice_terminal_velocities`](@ref) [m/s],
4. `v_m`, the mass-weighted velocity of [`ice_terminal_velocities`](@ref) [m/s],

each integrated with the rule `generator` for the state of [`table_state`](@ref).

# Construction

    IcePoint(params, terminal_velocity, generator)

# Fields
$(FIELDS)
"""
struct IcePoint{P, V, G}
    "The P3 parameters, [`CMP.ParametersP3`](@ref)"
    params::P
    "The terminal velocity, [`CMP.Chen2022VelType`](@ref)"
    terminal_velocity::V
    "The `QuadratureRule` of the integrals"
    generator::G
end
function (f::IcePoint)((logλ, w, ρ_rim, logρₐ))
    vel = f.terminal_velocity
    state = table_state(f.params, f.generator, logλ, w, ρ_rim)
    ρₐ = exp(logρₐ)
    (; v_n, v_m) = ice_terminal_velocities(vel, state, ρₐ)
    M_v = ice_ventilation_moments(vel, state, ρₐ).M_v
    return SA.SVector(M_v, node_self_collection(vel, f.generator, state, ρₐ), v_n, v_m)
end

"""
    group_table(::Val{:ice}, values, grid, FT)

Return the ice integral table of a [`P3Tables`](@ref), the logarithms of the four elements of the values of [`IcePoint`](@ref) on the axes `(logλ, w, ρ_rim, logρₐ)` of `grid.ice`.

The other ventilation moment, `M_1`, has a closed form and no table. [`ice_table_key`](@ref)
gives the coordinates of a state, and [`tabulated_ice_integrals`](@ref) reads the table.
"""
group_table(::Val{:ice}, values, grid, ::Type{FT}) where {FT} =
    UniformTable(values, table_axes(grid).ice, FT; logs = (true, true, true, true))

"""
    node_self_collection(vel, generator, state, ρₐ)

Return the self-collection rate of [`ice_self_collection`](@ref) per unit number squared, with the node values of each subinterval evaluated once.

For a `GaussLegendre` rule, the rate is the same quadrature over the same subintervals as the
`QuadratureRule` method of [`ice_self_collection`](@ref). For an outer diameter `D₁` in a
subinterval, the inner integral over each later subinterval takes the nodes of that subinterval,
which do not depend on `D₁`, so the terminal velocity, the effective radius and the size
distribution at these nodes are evaluated once for every `D₁`. Only the inner integral over the
rest of the subinterval of `D₁` evaluates them at nodes that depend on `D₁`. For any other rule,
the rate is that of [`ice_self_collection`](@ref).

# Arguments
- `vel`: [`CMP.Chen2022VelType`](@ref)
- `generator`: `QuadratureRule`
- `state`: [`P3State`](@ref) of unit ice number
- `ρₐ`: air density [kg/m³]
"""
node_self_collection(vel, generator, state, ρₐ) = ice_self_collection(vel, generator, state, ρₐ).dNdt
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

Return the [`P3State`](@ref) of unit ice number with the slope `logλ` at a point of the ice tables.

The state has the rime mass fraction at the rime coordinate `w` (see
[`rime_mass_fraction_at`](@ref)), the rime density `ρ_rim`, the mean particle mass of the
distribution, and the integration `generator`.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `generator`: `QuadratureRule`
- `logλ`: log of the slope parameter [log(1/m)]
- `w`: rime coordinate [-]
- `ρ_rim`: rime density [kg/m³]
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

`logλ` is the slope of `state.dist`, a primal value as for every rule, so the shape of the
distribution is held under differentiation. `w` is the [`rime_coordinate`](@ref) of the rime mass
fraction, and `ρ_rim` the rime density. They and `logρₐ` keep their partials, so a lookup gives
the derivatives with respect to the rime and the air density along their axes, as the quadrature
does through the particle properties and the terminal velocity.

# Arguments
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]
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

Return the `NamedTuple` `(; M_v, self_collection, v_n, v_m)` of the ice integral table of `tables` at the coordinates of `state`, see [`IcePoint`](@ref).

The integrals are per unit number, or per unit number squared for `self_collection`.

# Arguments
- `tables`: [`P3Tables`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]
"""
@inline function tabulated_ice_integrals(tables, state, ρₐ)
    (M_v, self_collection, v_n, v_m) = exp.(lookup(tables.groups.ice, ice_table_key(state, ρₐ)...))
    return (; M_v, self_collection, v_n, v_m)
end

@inline function ice_ventilation_moments(velocity_params, tables::P3Tables, state, ρₐ)
    (; ρn, μ, logλ) = state.dist
    # The closed form of the first moment, ∫ D N′ dD = ρn (μ + 1) / λ
    M_1 = ρn * (μ + 1) * exp(-logλ)
    M_v = tabulated_ice_integrals(tables, state, ρₐ).M_v * ρn
    return (; M_1, M_v)
end

@inline function ice_ventilation_integral(velocity_params, aps, tables::P3Tables, state, ρₐ)
    (a, b) = CO.ventilation_coefficients(state.particles.params.vent, aps)
    (; M_1, M_v) = ice_ventilation_moments(velocity_params, tables, state, ρₐ)
    return a * M_1 + b * M_v
end

@inline function ice_self_collection(vel, tables::P3Tables, state, ρₐ)
    K = tabulated_ice_integrals(tables, state, ρₐ).self_collection
    return (; dNdt = K * state.dist.ρn * state.ρn_ice)
end

function ice_terminal_velocities(velocity_params, tables::P3Tables, state, ρₐ; p = 1e-6)
    (; v_n, v_m) = tabulated_ice_integrals(tables, state, ρₐ)
    present = primal_value(state.dist.ρn) > 0
    return (; v_n = ifelse(present, v_n, zero(v_n)), v_m = ifelse(present, v_m, zero(v_m)))
end
