"""
    P3Tables{G, R}

Lookup tables that evaluate the P3 size-distribution integrals and the slope solve.

A [`P3State`](@ref) whose `integration` is a `P3Tables` takes the slope of its
[`P3Distribution`](@ref) from the inverse slope table, and each integral group has a method for
the type of the integration.

The tables are built from a complete set of P3 parameters, and they depend on the mass, area and
slope laws, the ice density, the nucleated crystal size, the largest mean particle mass, the
terminal velocities, and the cloud and rain size distributions. `fingerprint` identifies these
parameters, see [`table_fingerprint`](@ref).

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

"""
    group_table(::Val{name}, values, grid, FT)

Return the [`UniformTable`](@ref) of the table group `name` from the array `values` of the values of its point function at the points of `grid`, see [`table_axes`](@ref).
"""
function group_table end

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

Base.broadcastable(tables::P3Tables) = tuple(tables)
Base.eltype(tables::P3Tables) = eltype(tables.groups.slope)
Base.show(io::IO, tables::P3Tables) = print(io, "P3Tables{", eltype(tables), "}")
function Base.show(io::IO, ::MIME"text/plain", tables::P3Tables)
    print(io, "P3Tables{", eltype(tables), "}")
    for (name, table) in pairs(tables.groups)
        print(io, "\n  ", name, ": ", table)
    end
end
