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

It inverts [`rime_coordinate`](@ref), and returns the largest `Float64` below one at `w = 1`.
"""
@inline rime_mass_fraction_at(w) = min(w * (2 - w), prevfloat(one(w)))

"""
    P3TableGrid

The axes of the [`P3Tables`](@ref), with bounds derived from the parameters.

Each field is a `NamedTuple` of [`TableAxis`](@ref) objects:

- `slope`: the inverse slope table, on `logλ` [log(1/m)] within [`logλ_bounds`](@ref), the
  bounds of the slope solve, the rime coordinate `w` (see [`rime_coordinate`](@ref)) from 0 to 1,
  and the rime density `ρ_rim` [kg/m³] within [`rime_density_bounds`](@ref).
- `ice`: the ice integral tables, on `logλ` within [`integral_logλ_bounds`](@ref), the slopes of
  the mean particle masses within the limits of the number, the rime coordinate `w`, the rime
  density `ρ_rim`, and the log air density `logρₐ` [log(kg/m³)].
- `liquid_ice`: the collected number and mass and the impact scale of the liquid-ice collisions,
  on the four axes of `ice` and the log mean mass of cloud droplets `cloud_logx̄` and of raindrops
  `rain_logx̄` [log(kg)] within [`liquid_mass_bounds`](@ref).
- `rime_density`: the mean local rime density of the collisions, on the axes of `liquid_ice` and
  `logR`, the log of the ratio of the impact scale to the temperature below freezing, see
  [`group_table`](@ref). `logR` extends 6 beyond the range from 1 to 12 where the local
  rime density depends on the ratio, see [`rime_density_at`](@ref).
- `wet_set`: the tables of the freeze and shed partition, see [`tabulated_wet_set`](@ref), on the
  axes of `liquid_ice`, the log ice diameter `logD` [log(m)], the logit of the number fraction
  of the ice below a diameter `number_logit` (see [`number_logit`](@ref)), and the shape
  parameter `μ`. The bounds of `logD` and `μ` span the integration bounds of the collisions, and
  those of `number_logit` are the quantiles of [`liquid_ice_collision_quantile`](@ref).
- `particle_rime_density`: the mean local rime density of the mass that one particle collects,
  on the axes `logD`, `w`, `ρ_rim`, `logρₐ`, `cloud_logx̄`, `rain_logx̄` and `logR`, with the
  bounds of `wet_set` and `rime_density`.

A coordinate outside the bounds of an axis reads the table at the nearest end of the axis.

# Construction

    P3TableGrid(ice_params;
        slope = (; logλ = 128, w = 16, ρ_rim = 16),
        ice = (; logλ = 128, w = 33, ρ_rim = 12, logρₐ = 18),
        liquid_ice = (; logλ = 96, w = 17, ρ_rim = 8, logρₐ = 10, cloud_logx̄ = 24, rain_logx̄ = 12),
        rime_density = (; logλ = 16, w = 5, ρ_rim = 4, logρₐ = 5, cloud_logx̄ = 8, rain_logx̄ = 6, logR = 32),
        wet_set = (;
            logλ = 32, w = 9, ρ_rim = 4, logρₐ = 5, cloud_logx̄ = 8, rain_logx̄ = 8,
            logD = 96, number_logit = 25, μ = 64,
        ),
        particle_rime_density = (; logD = 32, w = 5, ρ_rim = 4, logρₐ = 4, cloud_logx̄ = 6, rain_logx̄ = 6, logR = 16),
        ρₐ_bounds = (0.005, 1.7),
    )

The keyword arguments other than `ρₐ_bounds` give the number of points of each axis.

    P3TableGrid(mp; kwargs...)

takes the ice parameters of a [`CMP.Microphysics2MParams`](@ref).

# Arguments
- `ice_params`: [`CMP.P3IceParams`](@ref). Its cloud size distribution must be a
  [`CMP.CloudParticlePDF_SB2006`](@ref) and its rain size distribution a
  [`CMP.RainParticlePDF_SB2006_limited`](@ref).

# Keyword arguments
- `ρₐ_bounds`: the smallest and largest air density [kg/m³]

# Fields
$(FIELDS)
"""
struct P3TableGrid{S, I, L, R, W, P}
    "Axes of the inverse slope table"
    slope::S
    "Axes of the ice integral tables"
    ice::I
    "Axes of the collected number and mass and the impact scale tables"
    liquid_ice::L
    "Axes of the rime density tables"
    rime_density::R
    "Axes of the wet-set tables"
    wet_set::W
    "Axes of the rime density tables of one particle"
    particle_rime_density::P
end

function P3TableGrid(
    ice_params::CMP.P3IceParams;
    slope = (; logλ = 128, w = 16, ρ_rim = 16),
    ice = (; logλ = 128, w = 33, ρ_rim = 12, logρₐ = 18),
    liquid_ice = (; logλ = 96, w = 17, ρ_rim = 8, logρₐ = 10, cloud_logx̄ = 24, rain_logx̄ = 12),
    rime_density = (; logλ = 16, w = 5, ρ_rim = 4, logρₐ = 5, cloud_logx̄ = 8, rain_logx̄ = 6, logR = 32),
    wet_set = (;
        logλ = 32, w = 9, ρ_rim = 4, logρₐ = 5, cloud_logx̄ = 8, rain_logx̄ = 8,
        logD = 96, number_logit = 25, μ = 64,
    ),
    particle_rime_density = (; logD = 32, w = 5, ρ_rim = 4, logρₐ = 4, cloud_logx̄ = 6, rain_logx̄ = 6, logR = 16),
    ρₐ_bounds = (0.005, 1.7),
)
    check_liquid_pdfs(ice_params)
    params = ice_params.scheme
    ρ_rim_bounds = rime_density_bounds(params)
    λ_bounds = integral_logλ_bounds(params)
    logρₐ_bounds = log.(ρₐ_bounds)
    cloud_bounds = log.(liquid_mass_bounds(ice_params.cloud_pdf))
    rain_bounds = log.(liquid_mass_bounds(ice_params.rain_pdf))
    logR_bounds = (-6, log(12) + 6)
    (; logD_bounds, μ_bounds) = collision_diameter_bounds(params, λ_bounds)
    p = liquid_ice_collision_quantile(Float64)
    # The axes of the particle properties and the air density, from the numbers of points `n`
    particle_axes(n) = (;
        w = TableAxis(0, 1, n.w),
        ρ_rim = TableAxis(ρ_rim_bounds..., n.ρ_rim),
        logρₐ = TableAxis(logρₐ_bounds..., n.logρₐ),
        cloud_logx̄ = TableAxis(cloud_bounds..., n.cloud_logx̄),
        rain_logx̄ = TableAxis(rain_bounds..., n.rain_logx̄),
    )
    ice_axes(n) = (;
        logλ = TableAxis(λ_bounds..., n.logλ),
        w = TableAxis(0, 1, n.w),
        ρ_rim = TableAxis(ρ_rim_bounds..., n.ρ_rim),
        logρₐ = TableAxis(logρₐ_bounds..., n.logρₐ),
    )
    liquid_axes(n) = (; logλ = TableAxis(λ_bounds..., n.logλ), particle_axes(n)...)
    return P3TableGrid(
        (;
            logλ = TableAxis(logλ_bounds(params)..., slope.logλ),
            w = TableAxis(0, 1, slope.w),
            ρ_rim = TableAxis(ρ_rim_bounds..., slope.ρ_rim),
        ),
        ice_axes(ice),
        liquid_axes(liquid_ice),
        (; liquid_axes(rime_density)..., logR = TableAxis(logR_bounds..., rime_density.logR)),
        (;
            liquid_axes(wet_set)...,
            logD = TableAxis(logD_bounds..., wet_set.logD),
            number_logit = TableAxis(log(p) - log(1 - p), log(1 - p) - log(p), wet_set.number_logit),
            μ = TableAxis(μ_bounds..., wet_set.μ),
        ),
        (;
            logD = TableAxis(logD_bounds..., particle_rime_density.logD),
            particle_axes(particle_rime_density)...,
            logR = TableAxis(logR_bounds..., particle_rime_density.logR),
        ),
    )
end

P3TableGrid(mp::CMP.Microphysics2MParams; kwargs...) = P3TableGrid(ice_parameters(mp); kwargs...)

"""
    integral_logλ_bounds(params)

Return the smallest and largest `logλ` [log(1/m)] of a size distribution with a mean particle mass within the limits of [`bounded_ice_number`](@ref).

The largest is that of [`logλ_bounds`](@ref), and the smallest is the slope of unrimed ice of the
mean mass [`ice_mean_particle_mass_max`](@ref), the largest particles, see
[`get_distribution_logλ`](@ref). Rime adds mass at a given diameter, so a rimed distribution of
that mean mass has a larger slope.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
"""
function integral_logλ_bounds(params)
    unrimed = P3Particles(params, zero(params.ρ_i), params.ρ_i)
    largest = get_distribution_logλ(unrimed, ice_mean_particle_mass_max(params), one(params.ρ_i))
    return (Float64(largest), Float64(last(logλ_bounds(params))))
end

"""
    ice_parameters(mp)

Return the [`CMP.P3IceParams`](@ref) of the [`CMP.Microphysics2MParams`](@ref) `mp`, or throw an `ArgumentError` for parameters without ice.
"""
ice_parameters(mp::CMP.Microphysics2MParams) = ice_parameters(mp.ice)
ice_parameters(ice_params::CMP.P3IceParams) = ice_params
ice_parameters(::Nothing) = throw(ArgumentError("the P3 tables need parameters with P3 ice"))

"""
    check_liquid_pdfs(ice_params)

Throw an `ArgumentError` unless the tables support the cloud and rain size distributions of `ice_params`, see [`liquid_mass_bounds`](@ref).
"""
function check_liquid_pdfs((; cloud_pdf, rain_pdf))
    foreach(liquid_mass_bounds, (cloud_pdf, rain_pdf))
    return nothing
end

"""
    table_axes(grid)

Return the axes of the point function of each table group of the [`P3TableGrid`](@ref) `grid`, a `NamedTuple` of tuples of [`TableAxis`](@ref) objects.

A table group is the set of tables filled by one point function of
[`table_point_functions`](@ref), and its axes are the coordinates of that function, see
[`tabulate`](@ref). The groups are

- `slope`: `(logλ, w, ρ_rim)` of `grid.slope`,
- `ice`: `(logλ, w, ρ_rim, logρₐ)` of `grid.ice`,
- `cloud_collision`, `rain_collision`: `(logλ, w, ρ_rim, logρₐ, logx̄)` of `grid.liquid_ice`,
  with `logx̄` the `cloud_logx̄` or `rain_logx̄` axis,
- `cloud_rime_density`, `rain_rime_density`: `(logλ, w, ρ_rim, logρₐ, logx̄)` of
  `grid.rime_density`, whose point functions return the values at every point of `logR`,
- `collision_bounds`: `(μ,)` of `grid.wet_set`,
- `cloud_collected`, `rain_collected`: `(logD, w, ρ_rim, logρₐ, logx̄)` of `grid.wet_set`,
- `capacity`: `(logD, w, ρ_rim, logρₐ)` of `grid.wet_set`,
- `cloud_fraction`, `rain_fraction`: `(logλ, w, ρ_rim, logρₐ, logx̄)` of `grid.wet_set`, whose
  point functions return the values at every point of `number_logit`,
- `capacity_fraction`: `(logλ, w, ρ_rim, logρₐ)` of `grid.wet_set`, whose point function
  returns the values at every point of `number_logit`,
- `cloud_particle_rime_density`, `rain_particle_rime_density`: `(logD, w, ρ_rim, logρₐ, logx̄)`
  of `grid.particle_rime_density`, whose point functions return the values at every point of
  `logR`.
"""
function table_axes((; slope, ice, liquid_ice, rime_density, wet_set, particle_rime_density)::P3TableGrid)
    ice_axes(g) = (g.logλ, g.w, g.ρ_rim, g.logρₐ)
    particle_axes(g) = (g.logD, g.w, g.ρ_rim, g.logρₐ)
    return (;
        slope = (slope.logλ, slope.w, slope.ρ_rim),
        ice = ice_axes(ice),
        cloud_collision = (ice_axes(liquid_ice)..., liquid_ice.cloud_logx̄),
        rain_collision = (ice_axes(liquid_ice)..., liquid_ice.rain_logx̄),
        cloud_rime_density = (ice_axes(rime_density)..., rime_density.cloud_logx̄),
        rain_rime_density = (ice_axes(rime_density)..., rime_density.rain_logx̄),
        collision_bounds = (wet_set.μ,),
        cloud_collected = (particle_axes(wet_set)..., wet_set.cloud_logx̄),
        rain_collected = (particle_axes(wet_set)..., wet_set.rain_logx̄),
        capacity = particle_axes(wet_set),
        cloud_fraction = (ice_axes(wet_set)..., wet_set.cloud_logx̄),
        rain_fraction = (ice_axes(wet_set)..., wet_set.rain_logx̄),
        capacity_fraction = ice_axes(wet_set),
        cloud_particle_rime_density = (particle_axes(particle_rime_density)..., particle_rime_density.cloud_logx̄),
        rain_particle_rime_density = (particle_axes(particle_rime_density)..., particle_rime_density.rain_logx̄),
    )
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
    logλs = range(Float64(logλ_lo), Float64(logλ_hi); length = n)
    μs = map(logλ -> Float64(get_μ(params.slope, logλ)), logλs)
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
