"""
    P3Tables{G, R}

Lookup tables that evaluate the P3 size-distribution integrals and the slope solve.

A [`P3State`](@ref) whose `integration` is a `P3Tables` reads every ice integral group from the
tables, and its [`P3Distribution`](@ref) takes the slope from the inverse slope table.

The tables are built from a complete set of P3 parameters, and they depend on the mass, area and
slope laws, the ice density, the nucleated crystal size, the largest mean particle mass, the
terminal velocities, and the cloud and rain size distributions. `fingerprint` identifies these
parameters, so [`tables_are_current`](@ref) detects tables built from other parameters.
[`with_tables`](@ref) builds the tables and returns parameters that use them.

# Construction

    P3Tables(ice_params, grid, values)

assembles the tables from the arrays `values` of the values of the point functions of
[`table_point_functions`](@ref) at the points of the axes of [`table_axes`](@ref), a `NamedTuple`
with the keys of [`table_axes`](@ref), see [`tabulate`](@ref). The tables store the values in the
float type of `ice_params`, in arrays of the type of `values`, with the logarithm of the values
that scale with a number, see [`UniformTable`](@ref). It throws an `ArgumentError` for a cloud or
rain size distribution the tables do not support (see [`liquid_mass_bounds`](@ref)) or for
missing table groups, and a `DimensionMismatch` for values that do not match the grid.

# Arguments
- `ice_params`: [`CMP.P3IceParams`](@ref), the parameters of the point functions
- `grid`: [`P3TableGrid`](@ref)
- `values`: `NamedTuple` of arrays

# Fields
$(FIELDS)
"""
struct P3Tables{G <: NamedTuple, R}
    "Table of each table group, a [`UniformTable`](@ref), with the keys of [`table_axes`](@ref)"
    groups::G
    "Gauss-Legendre rule of the freeze and shed partition over the wet set, see [`tabulated_wet_set`](@ref)"
    partition_rule::R
    "Hash of the parameters the tables were built from"
    fingerprint::UInt64
end

Adapt.@adapt_structure P3Tables

function P3Tables(ice_params::CMP.P3IceParams, grid::P3TableGrid, values::NamedTuple)
    check_liquid_pdfs(ice_params)
    names = keys(table_axes(grid))
    issetequal(keys(values), names) ||
        throw(ArgumentError("the table values have the groups $(keys(values)), not $names"))
    return assemble_tables(ice_params, grid, name -> values[name])
end

"""
    assemble_tables(ice_params, grid, values_of)

Return the [`P3Tables`](@ref) of `ice_params` on `grid` with the values `values_of(name)` of each table group `name`, see [`table_axes`](@ref).

Each group is stored by [`group_table`](@ref) before the values of the next group are requested.
"""
function assemble_tables(ice_params, grid, values_of::F) where {F}
    FT = eltype(ice_params.scheme)
    names = keys(table_axes(grid))
    groups = ntuple(i -> group_table(Val(names[i]), values_of(names[i]), grid, FT), Val(length(names)))
    return P3Tables(NamedTuple{names}(groups), partition_rule(FT), table_fingerprint(ice_params))
end

"""
    group_table(::Val{name}, values, grid, FT)

Return the [`UniformTable`](@ref) of the table group `name` from the array `values` of the values of its point function at the points of `grid`, see [`table_axes`](@ref).
"""
function group_table end

"""
    table_point_functions(ice_params, grid; generator = GaussLegendre(Float64, 16))

Return the point functions of the table groups of the [`P3Tables`](@ref), a `NamedTuple` with the keys of [`table_axes`](@ref).

A point function maps the coordinates of one point of its group, an `NTuple` of `Float64`, to
the value of its tables at that point by quadrature with the rule `generator`. The value is a
number or an `SVector` with one element per table or per point of an axis that the function
evaluates in one pass. The point functions are

- `slope`: [`SlopePoint`](@ref),
- `ice`: [`IcePoint`](@ref),
- `cloud_collision`, `rain_collision`: [`CollisionPoint`](@ref),
- `cloud_rime_density`, `rain_rime_density`: [`RimeDensityPoint`](@ref) on the points of the
  `logR` axis,
- `collision_bounds`: [`CollisionBoundsPoint`](@ref),
- `cloud_collected`, `rain_collected`: [`CollectedPoint`](@ref),
- `capacity`: [`CapacityPoint`](@ref),
- `cloud_fraction`, `rain_fraction`: [`FractionPoint`](@ref) on the points of the
  `number_logit` axis,
- `capacity_fraction`: [`CapacityFractionPoint`](@ref) on the points of the `number_logit` axis,
- `cloud_particle_rime_density`, `rain_particle_rime_density`: [`ParticleRimeDensityPoint`](@ref)
  on the points of the `logR` axis.

Each point function is an `isbits` callable object, so it can be mapped over the points of its
axes in any array type, see [`tabulate`](@ref), and the values assemble into [`P3Tables`](@ref).

# Arguments
- `ice_params`: [`CMP.P3IceParams`](@ref)
- `grid`: [`P3TableGrid`](@ref)

# Keyword arguments
- `generator`: the `QuadratureRule` that evaluates the integrals
"""
function table_point_functions(ice_params, grid; generator = GaussLegendre(Float64, 16))
    check_liquid_pdfs(ice_params)
    (; scheme, terminal_velocity, cloud_pdf, rain_pdf) = ice_params
    logR = Tuple(axis_points(grid.rime_density.logR))
    particle_logR = Tuple(axis_points(grid.particle_rime_density.logR))
    logits = Tuple(axis_points(grid.wet_set.number_logit))
    return (;
        slope = SlopePoint(scheme),
        ice = IcePoint(scheme, terminal_velocity, generator),
        cloud_collision = CollisionPoint(scheme, terminal_velocity, cloud_pdf, generator),
        rain_collision = CollisionPoint(scheme, terminal_velocity, rain_pdf, generator),
        cloud_rime_density = RimeDensityPoint(scheme, terminal_velocity, cloud_pdf, generator, logR),
        rain_rime_density = RimeDensityPoint(scheme, terminal_velocity, rain_pdf, generator, logR),
        collision_bounds = CollisionBoundsPoint(),
        cloud_collected = CollectedPoint(scheme, terminal_velocity, cloud_pdf, generator),
        rain_collected = CollectedPoint(scheme, terminal_velocity, rain_pdf, generator),
        capacity = CapacityPoint(scheme, terminal_velocity),
        cloud_fraction = FractionPoint(scheme, terminal_velocity, cloud_pdf, generator, logits),
        rain_fraction = FractionPoint(scheme, terminal_velocity, rain_pdf, generator, logits),
        capacity_fraction = CapacityFractionPoint(scheme, terminal_velocity, generator, logits),
        cloud_particle_rime_density =
        ParticleRimeDensityPoint(scheme, terminal_velocity, cloud_pdf, generator, particle_logR),
        rain_particle_rime_density =
        ParticleRimeDensityPoint(scheme, terminal_velocity, rain_pdf, generator, particle_logR),
    )
end

"""
    with_tables(mp, ArrayType = Array; grid = P3TableGrid(mp), generator = GaussLegendre(Float64, 16))

Return a copy of the parameters `mp` whose size-distribution integrals evaluate from [`P3Tables`](@ref) built for them.

The tables are filled at the points of `grid`: the point functions of
[`table_point_functions`](@ref) are mapped over the points of the axes of [`table_axes`](@ref)
with [`tabulate`](@ref) in arrays of type `ArrayType`, and each group is stored as a table before
the next one is filled. The tables are therefore built where `ArrayType` places its elements, and
their arrays have that type. For an `Array`, the points are evaluated on all threads. The values
are computed in `Float64` with the rule `generator` and stored in the float type of `mp`.

    with_tables(mp, grid, values)

returns the copy with the tables that [`P3Tables`](@ref) assembles from the arrays `values` of
the point functions at the points of `grid`, for a host that maps the point functions itself.

These are the only ways to obtain parameters with tables: the constructors of
[`CMP.P3IceParams`](@ref) and [`CMP.Microphysics2MParams`](@ref) from a parameter dictionary
throw an `ArgumentError` for tables built from other parameters.

# Arguments
- `mp`: [`CMP.Microphysics2MParams`](@ref) with P3 ice, or [`CMP.P3IceParams`](@ref)
- `ArrayType`: the array type of the tables, such as `Array`

# Keyword arguments
- `grid`: [`P3TableGrid`](@ref)
- `generator`: the `QuadratureRule` that evaluates the integrals
"""
function with_tables(
    ice_params::CMP.P3IceParams,
    ArrayType = Array;
    grid = P3TableGrid(ice_params),
    generator = GaussLegendre(Float64, 16),
)
    points = table_point_functions(ice_params, grid; generator)
    axes = table_axes(grid)
    tables = assemble_tables(ice_params, grid, name -> tabulate(points[name], axes[name], ArrayType))
    return with_integration(ice_params, tables)
end
with_tables(ice_params::CMP.P3IceParams, grid::P3TableGrid, values::NamedTuple) =
    with_integration(ice_params, P3Tables(ice_params, grid, values))
with_tables(mp::CMP.Microphysics2MParams, args...; kwargs...) =
    CMP.Microphysics2MParams(mp.warm_rain, with_tables(ice_parameters(mp), args...; kwargs...))

"""
    with_integration(ice_params, integration)

Return a copy of the [`CMP.P3IceParams`](@ref) `ice_params` with the field `integration` replaced.
"""
function with_integration(ice_params, integration)
    fields = map(name -> name === :integration ? integration : getfield(ice_params, name), fieldnames(CMP.P3IceParams))
    return CMP.P3IceParams(fields...)
end

"""
    partition_rule(FT)

Return the Gauss-Legendre rule of the freeze and shed partition over the wet set, see [`tabulated_wet_set`](@ref).
"""
partition_rule(::Type{FT}) where {FT} = GaussLegendre(FT, 8)

"""
    table_fingerprint(ice_params)

Return a hash of the parameters of [`CMP.P3IceParams`](@ref) that the [`P3Tables`](@ref) depend on.
"""
table_fingerprint((; scheme, terminal_velocity, cloud_pdf, rain_pdf)) =
    hash((scheme, terminal_velocity, cloud_pdf, rain_pdf))

"""
    tables_are_current(ice_params)

Return whether the size-distribution integrals of `ice_params` evaluate with its own parameters.

This is `false` only if `ice_params.integration` is a [`P3Tables`](@ref) built from other
parameters, for example after a host changed the parameters of a copy without rebuilding the
tables.

# Arguments
- `ice_params`: [`CMP.P3IceParams`](@ref), or [`CMP.Microphysics2MParams`](@ref), which is
  current without ice
"""
tables_are_current(mp::CMP.Microphysics2MParams) = isnothing(mp.ice) || tables_are_current(mp.ice)
tables_are_current(ice_params::CMP.P3IceParams) = tables_are_current(ice_params.integration, ice_params)
tables_are_current(::QuadratureRule, ice_params) = true
tables_are_current(tables::P3Tables, ice_params) = tables.fingerprint == table_fingerprint(ice_params)

CMP.check_integration(tables::P3Tables, ice_params) =
    tables_are_current(tables, ice_params) ? nothing :
    throw(
        ArgumentError("the P3 tables were built from other parameters, build the parameters with P3Scheme.with_tables"),
    )

Base.broadcastable(tables::P3Tables) = tuple(tables)
Base.eltype(tables::P3Tables) = eltype(tables.groups.slope)
Base.show(io::IO, tables::P3Tables) = print(io, "P3Tables{", eltype(tables), "}")
function Base.show(io::IO, ::MIME"text/plain", tables::P3Tables)
    print(io, "P3Tables{", eltype(tables), "}")
    for (name, table) in pairs(tables.groups)
        print(io, "\n  ", name, ": ", table)
    end
end
