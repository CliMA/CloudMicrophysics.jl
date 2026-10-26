"""
    rime_coordinate(F_rim)

Return the rime coordinate `w = 1 - √(1 - F_rim)` of the tables.

The coordinate spreads the rime mass fractions near one, where the ice properties change fastest
with `F_rim`. Its inverse is [`rime_mass_fraction_at`](@ref).
"""
@inline rime_coordinate(F_rim) = 1 - sqrt(1 - F_rim)

"""
    rime_mass_fraction_at(w)

Return the rime mass fraction `F_rim = w (2 - w)` at the rime coordinate `w`, below one.

It inverts [`rime_coordinate`](@ref), and returns the largest float below one at `w = 1`.
"""
@inline rime_mass_fraction_at(w) = min(w * (2 - w), prevfloat(one(w)))

"""
    P3TableGrid

The axes of each table group of the [`P3Tables`](@ref), with bounds derived from the parameters.

`groups` is a `NamedTuple` with a tuple of [`TableAxis`](@ref) objects for each table group, in
the float type of the parameters. The coordinates are

- `logλ` [log(1/m)]: within [`integral_logλ_bounds`](@ref), and within [`logλ_bounds`](@ref) for
  the inverse slope table,
- `w`: the rime coordinate of [`rime_coordinate`](@ref), from 0 to 1,
- `ρ_rim` [kg/m³]: within [`rime_density_bounds`](@ref),
- `logρₐ` [log(kg/m³)]: within the log of `ρₐ_bounds`.

A coordinate outside an axis reads the table at the nearest end of the axis.

# Construction

    P3TableGrid(ice_params; ρₐ_bounds = (0.005, 1.7), points...)

The keyword arguments `points` replace the numbers of points of the families of
[`TABLE_POINTS`](@ref). `P3TableGrid(mp; kwargs...)` takes the ice parameters of a
[`CMP.Microphysics2MParams`](@ref).

# Arguments
- `ice_params`: [`CMP.P3IceParams`](@ref)

# Fields
$(FIELDS)
"""
struct P3TableGrid{G <: NamedTuple}
    "Axes of each table group"
    groups::G
end

"""
    TABLE_POINTS

The default numbers of points of the axes of the [`P3TableGrid`](@ref), for each family of table groups.
"""
const TABLE_POINTS =
    (; slope = (; logλ = 128, w = 16, ρ_rim = 16), ice = (; logλ = 128, w = 33, ρ_rim = 12, logρₐ = 18))

function P3TableGrid(ice_params::CMP.P3IceParams; ρₐ_bounds = (0.005, 1.7), points...)
    issubset(keys(points), keys(TABLE_POINTS)) ||
        throw(ArgumentError("the table families are $(keys(TABLE_POINTS)), got $(keys(points))"))
    (; slope, ice) = merge(TABLE_POINTS, values(points))
    params = ice_params.scheme
    FT = eltype(params)
    bounds = (;
        logλ = integral_logλ_bounds(params),
        w = (zero(FT), one(FT)),
        ρ_rim = rime_density_bounds(params),
        logρₐ = log.(FT.(ρₐ_bounds)),
    )
    axes(family, names...) = map(name -> TableAxis(bounds[name]..., family[name]), names)
    return P3TableGrid((;
        slope = (TableAxis(logλ_bounds(params)..., slope.logλ), axes(slope, :w, :ρ_rim)...),
        ice = axes(ice, :logλ, :w, :ρ_rim, :logρₐ),
    ))
end

P3TableGrid(mp::CMP.Microphysics2MParams; kwargs...) = P3TableGrid(ice_parameters(mp); kwargs...)

"""
    integral_logλ_bounds(params)

Return the smallest and largest `logλ` [log(1/m)] of a size distribution with a mean particle mass within the limits of [`bounded_ice_number`](@ref).

The largest is that of [`logλ_bounds`](@ref), and the smallest is the slope of unrimed ice of the
mean mass [`ice_mean_particle_mass_max`](@ref).
"""
function integral_logλ_bounds(params)
    unrimed = P3Particles(params, zero(params.ρ_i), params.ρ_i)
    largest = get_distribution_logλ(unrimed, ice_mean_particle_mass_max(params), one(params.ρ_i))
    return (largest, last(logλ_bounds(params)))
end

ice_parameters(mp::CMP.Microphysics2MParams) = ice_parameters(mp.ice)
ice_parameters(ice_params::CMP.P3IceParams) = ice_params
ice_parameters(::Nothing) = throw(ArgumentError("the P3 tables need parameters with P3 ice"))
