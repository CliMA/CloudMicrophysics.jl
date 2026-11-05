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
the float type of the parameters. A group whose point function evaluates every point of an axis
in one pass has that axis first. The coordinates are

- `logλ` [log(1/m)]: within [`integral_logλ_bounds`](@ref), and within [`logλ_bounds`](@ref) for
  the inverse slope table,
- `w`: the rime coordinate of [`rime_coordinate`](@ref), from 0 to 1,
- `ρ_rim` [kg/m³]: within [`rime_density_bounds`](@ref),
- `logρₐ` [log(kg/m³)]: within the log of `ρₐ_bounds`,
- `cloud_logx̄`, `rain_logx̄` [log(kg)]: within [`liquid_mass_bounds`](@ref),
- `logR`: the log of the ratio of the impact scale to the temperature below freezing, see
  [`rime_density_point`](@ref), from 6 below to 6 above the range from 1 to 12 of
  [`rime_density_at`](@ref),
- `logD` [log(m)] and `μ`: the ice diameters and shape parameters of the collision integrals, see
  [`collision_diameter_bounds`](@ref),
- `number_logit`: the logit of the number fraction of the ice below a diameter, see
  [`number_logit`](@ref), within the quantiles of [`liquid_ice_collision_quantile`](@ref).

A coordinate outside an axis reads the table at the nearest end of the axis.

# Construction

    P3TableGrid(ice_params; ρₐ_bounds = (0.005, 1.7), points...)

The keyword arguments `points` replace the numbers of points of the families of
[`TABLE_POINTS`](@ref). `P3TableGrid(mp; kwargs...)` takes the ice parameters of a
[`CMP.Microphysics2MParams`](@ref).

# Arguments
- `ice_params`: [`CMP.P3IceParams`](@ref) with a [`CMP.CloudParticlePDF_SB2006`](@ref) and a
  [`CMP.RainParticlePDF_SB2006_limited`](@ref)

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
const TABLE_POINTS = (;
    slope = (; logλ = 128, w = 16, ρ_rim = 16),
    ice = (; logλ = 128, w = 33, ρ_rim = 12, logρₐ = 18),
    cloud_collision = (; logλ = 96, w = 17, ρ_rim = 8, logρₐ = 10, cloud_logx̄ = 24),
    rain_collision = (; logλ = 64, w = 33, ρ_rim = 12, logρₐ = 10, rain_logx̄ = 12),
    rime_density = (; logλ = 16, w = 5, ρ_rim = 4, logρₐ = 5, cloud_logx̄ = 8, rain_logx̄ = 6, logR = 32),
    wet_set = (;
        logλ = 32,
        w = 9,
        ρ_rim = 4,
        logρₐ = 5,
        cloud_logx̄ = 8,
        rain_logx̄ = 8,
        logD = 96,
        number_logit = 25,
        μ = 64,
    ),
)

function P3TableGrid(ice_params::CMP.P3IceParams; ρₐ_bounds = (0.005, 1.7), points...)
    check_liquid_pdfs(ice_params)
    issubset(keys(points), keys(TABLE_POINTS)) ||
        throw(ArgumentError("the table families are $(keys(TABLE_POINTS)), got $(keys(points))"))
    (; slope, ice, cloud_collision, rain_collision, rime_density, wet_set) = merge(TABLE_POINTS, values(points))
    params = ice_params.scheme
    FT = eltype(params)
    λ_bounds = integral_logλ_bounds(params)
    (; logD_bounds, μ_bounds) = collision_diameter_bounds(params, λ_bounds)
    p = liquid_ice_collision_quantile(FT)
    bounds = (;
        logλ = λ_bounds,
        w = (zero(FT), one(FT)),
        ρ_rim = rime_density_bounds(params),
        logρₐ = log.(FT.(ρₐ_bounds)),
        cloud_logx̄ = log.(liquid_mass_bounds(ice_params.cloud_pdf)),
        rain_logx̄ = log.(liquid_mass_bounds(ice_params.rain_pdf)),
        logR = (FT(-6), log(FT(12)) + 6),
        logD = logD_bounds,
        number_logit = (log(p) - log(1 - p), log(1 - p) - log(p)),
        μ = μ_bounds,
    )
    axes(family, names...) = map(name -> TableAxis(bounds[name]..., family[name]), names)
    ice_axes = (:logλ, :w, :ρ_rim, :logρₐ)
    particle_axes = (:logD, :w, :ρ_rim, :logρₐ)
    return P3TableGrid((;
        slope = (TableAxis(logλ_bounds(params)..., slope.logλ), axes(slope, :w, :ρ_rim)...),
        ice = axes(ice, ice_axes...),
        cloud_collision = axes(cloud_collision, ice_axes..., :cloud_logx̄),
        rain_collision = axes(rain_collision, ice_axes..., :rain_logx̄),
        cloud_rime_density = axes(rime_density, :logR, ice_axes..., :cloud_logx̄),
        rain_rime_density = axes(rime_density, :logR, ice_axes..., :rain_logx̄),
        collision_bounds = axes(wet_set, :μ),
        cloud_collected = axes(wet_set, particle_axes..., :cloud_logx̄),
        rain_collected = axes(wet_set, particle_axes..., :rain_logx̄),
        capacity = axes(wet_set, particle_axes...),
        cloud_fraction = axes(wet_set, :number_logit, ice_axes..., :cloud_logx̄),
        rain_fraction = axes(wet_set, :number_logit, ice_axes..., :rain_logx̄),
        capacity_fraction = axes(wet_set, :number_logit, ice_axes...),
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

function check_liquid_pdfs((; cloud_pdf, rain_pdf))
    foreach(liquid_mass_bounds, (cloud_pdf, rain_pdf))
    return nothing
end

"""
    collision_diameter_bounds(params, (logλ_lo, logλ_hi); n = 257)

Return the bounds of the ice diameters of the liquid-ice collisions over the slopes from `logλ_lo` to `logλ_hi`.

The collision integrals extend between the quantiles of [`liquid_ice_collision_quantile`](@ref)
of the size distribution, see [`integral_bounds`](@ref). The result is a `NamedTuple` of the
bounds of the log diameter `logD_bounds` [log(m)] and of the shape parameter `μ_bounds` over `n`
uniform values of `logλ`. A constant shape parameter gives the bounds `(μ, μ + 1)`.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `(logλ_lo, logλ_hi)`: bounds of `logλ` [log(1/m)]
"""
function collision_diameter_bounds(params, (logλ_lo, logλ_hi); n = 257)
    logλs = range(logλ_lo, logλ_hi; length = n)
    μs = map(logλ -> get_μ(params.slope, logλ), logλs)
    logξ_lo = map(collision_logξ_lo, μs)
    logξ_hi = map(collision_logξ_hi, μs)
    μ_lo, μ_hi = extrema(μs)
    return (;
        logD_bounds = (minimum(logξ_lo .- logλs), maximum(logξ_hi .- logλs)),
        μ_bounds = μ_lo < μ_hi ? (μ_lo, μ_hi) : (μ_lo, μ_lo + 1),
    )
end

"""
    collision_logξ_lo(μ)
    collision_logξ_hi(μ)

Return the log of the lower or upper bound of the scaled diameter `ξ = D λ` of the liquid-ice collisions.

The bounds are the quantiles `p` and `1 - p` of the gamma distribution of shape `μ + 1`, with
`p` from [`liquid_ice_collision_quantile`](@ref), as in [`integral_bounds`](@ref).
"""
function collision_logξ_lo(μ)
    p = liquid_ice_collision_quantile(typeof(μ))
    return log(UT.gamma_inc_inv(μ + 1, p, 1 - p))
end
function collision_logξ_hi(μ)
    p = liquid_ice_collision_quantile(typeof(μ))
    return log(UT.gamma_inc_inv(μ + 1, 1 - p, p))
end

"""
    liquid_mass_bounds(pdf)

Return the smallest and largest mean particle mass [kg] of a liquid size distribution the tables support.

For [`CMP.CloudParticlePDF_SB2006`](@ref), these are the mass limits of the distribution. The
slope of [`CMP.RainParticlePDF_SB2006_limited`](@ref) lies within its limits `λ_min` and `λ_max`,
so its mean mass `π ρw / λ³` lies within the bounds. Any other distribution throws an
`ArgumentError`.
"""
liquid_mass_bounds((; xc_min, xc_max)::CMP.CloudParticlePDF_SB2006) = (xc_min, xc_max)
liquid_mass_bounds((; ρw, λ_min, λ_max)::CMP.RainParticlePDF_SB2006_limited) = (π * ρw / λ_max^3, π * ρw / λ_min^3)
liquid_mass_bounds(pdf) =
    throw(ArgumentError("the P3 tables do not support the liquid size distribution $(nameof(typeof(pdf)))"))
