"""
    P3Tables{G}

Lookup tables that evaluate the P3 size-distribution integrals and the slope solve.

A [`P3State`](@ref) whose `integration` is a `P3Tables` reads every ice integral group from the
tables, and its [`P3Distribution`](@ref) takes the slope from the inverse slope table.

[`with_tables`](@ref) builds the tables and returns parameters that use them, and
[`tables_are_current`](@ref) detects tables built from other parameters.

The tables hold the integrals for the terminal velocities and the liquid size distributions they
were built with, so the table methods ignore the velocity and liquid mass arguments of the process
functions.

# Fields
$(FIELDS)
"""
struct P3Tables{G <: NamedTuple}
    "Table of each table group, a [`UniformTable`](@ref), with the keys of the groups of [`P3TableGrid`](@ref)"
    groups::G
    "Hash of the parameters the tables were built from"
    fingerprint::UInt64
end

Adapt.@adapt_structure P3Tables

"""
    TabulatedP3State

A [`P3State`](@ref) whose size-distribution integrals evaluate from [`P3Tables`](@ref).
"""
const TabulatedP3State = P3State{<:Any, <:Any, <:P3Tables}

"""
    table_groups(ice_params, grid, generator)

Return the point function and the storage of each table group, a `NamedTuple` with the keys of `grid.groups`.

Each entry is a `NamedTuple` `(; point, logs, inner)` of the point function, which maps the
coordinates of a grid point to the values of the group there by quadrature with the rule
`generator`, and the keyword arguments `logs` and `inner` of [`UniformTable`](@ref).
"""
function table_groups(ice_params, grid, generator)
    (; scheme, terminal_velocity, cloud_pdf, rain_pdf) = ice_params
    vel = terminal_velocity
    logR = Tuple(axis_points(first(grid.groups.cloud_rime_density)))
    logits = Tuple(axis_points(first(grid.groups.cloud_fraction)))
    group(point; logs = nothing, inner = false) = (; point, logs, inner)
    return (;
        slope = group(x -> slope_point(scheme, x)),
        ice = group(x -> ice_point(scheme, vel, generator, x); logs = (true, true, true, true)),
        cloud_collision = group(x -> collision_point(scheme, vel, cloud_pdf, generator, x); logs = (true, true, false)),
        rain_collision = group(x -> collision_point(scheme, vel, rain_pdf, generator, x); logs = (true, true, false)),
        cloud_rime_density = group(x -> rime_density_point(scheme, vel, cloud_pdf, generator, logR, x); inner = true),
        rain_rime_density = group(x -> rime_density_point(scheme, vel, rain_pdf, generator, logR, x); inner = true),
        collision_bounds = group(collision_bounds_point),
        cloud_collected = group(x -> collected_point(scheme, vel, cloud_pdf, generator, x); logs = (true, false)),
        rain_collected = group(x -> collected_point(scheme, vel, rain_pdf, generator, x); logs = (true, false)),
        capacity = group(x -> capacity_point(scheme, vel, x); logs = (true,)),
        cloud_fraction = group(x -> fraction_point(scheme, vel, cloud_pdf, generator, logits, x); inner = true),
        rain_fraction = group(x -> fraction_point(scheme, vel, rain_pdf, generator, logits, x); inner = true),
        capacity_fraction = group(x -> capacity_fraction_point(scheme, vel, generator, logits, x); inner = true),
    )
end

"""
    with_tables(mp, ArrayType = Array; grid = P3TableGrid(mp), generator = GaussLegendre(FT, 16))

Return a copy of the parameters `mp` whose size-distribution integrals evaluate from [`P3Tables`](@ref) built for them.

The point function of each group of [`table_groups`](@ref) is evaluated at the points of `grid`
by [`tabulate`](@ref) in arrays of type `ArrayType`, so the tables are built where `ArrayType`
places its elements, and each group is stored before the next one is evaluated. The values are
computed and stored in the float type `FT` of `mp`, with the rule `generator`. The constructors of
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
    generator = GaussLegendre(eltype(ice_params.scheme), 16),
)
    check_liquid_pdfs(ice_params)
    FT = eltype(ice_params.scheme)
    groups = map(table_groups(ice_params, grid, generator), grid.groups) do (; point, logs, inner), axes
        UniformTable(tabulate(point, inner ? Base.tail(axes) : axes, ArrayType), axes, FT; logs, inner)
    end
    return with_integration(ice_params, P3Tables(groups, table_fingerprint(ice_params)))
end
with_tables(mp::CMP.Microphysics2MParams, args...; kwargs...) =
    CMP.Microphysics2MParams(mp.warm_rain, with_tables(ice_parameters(mp), args...; kwargs...))

function with_integration(ice_params, integration)
    fields = map(name -> name === :integration ? integration : getfield(ice_params, name), fieldnames(CMP.P3IceParams))
    return CMP.P3IceParams(fields...)
end

table_fingerprint((; scheme, terminal_velocity, cloud_pdf, rain_pdf)) =
    hash((scheme, terminal_velocity, cloud_pdf, rain_pdf))

"""
    tables_are_current(ice_params)

Return whether the size-distribution integrals of `ice_params` evaluate with its own parameters, `false` only for [`P3Tables`](@ref) built from other parameters.

`ice_params` is a [`CMP.P3IceParams`](@ref), or a [`CMP.Microphysics2MParams`](@ref), which is
current without ice.
"""
tables_are_current(mp::CMP.Microphysics2MParams) = isnothing(mp.ice) || tables_are_current(mp.ice)
tables_are_current(ice_params::CMP.P3IceParams) =
    !(ice_params.integration isa P3Tables) || ice_params.integration.fingerprint == table_fingerprint(ice_params)

CMP.check_integration(tables::P3Tables, ice_params) =
    tables.fingerprint == table_fingerprint(ice_params) ? nothing :
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
