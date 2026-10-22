# TODO: Add the liquid fraction `F_liq` to `P3Particles`
"""
    P3Particles{FT}

Properties of a P3 ice particle, which depend on its maximum dimension `D` but not on the amount of ice.

The struct holds:

- the P3 parameterizations `params`,
- the rime mass fraction `F_rim` and the rime density `ρ_rim`,
- the derived graupel density `ρ_g`,
- the size thresholds `D_th`, `D_gr`, and `D_cr` of the mass regimes.

Together they define the mass [`ice_mass`](@ref), the area [`ice_area`](@ref), and the aspect
ratio [`ϕᵢ`](@ref) of a particle.

# Construction

    P3Particles(params, F_rim, ρ_rim)

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `F_rim`: rime mass fraction [-]
- `ρ_rim`: rime density [kg/m³]

# Fields
$(FIELDS)
"""
struct P3Particles{FT, PARAMS <: CMP.ParametersP3}
    "[`CMP.ParametersP3`](@ref) object"
    params::PARAMS
    "Rime mass fraction"
    F_rim::FT
    "Rime density [kg/m³]"
    ρ_rim::FT
    "Graupel density [kg/m³] — `NaN` when `F_rim = 0` (no graupel regime)"
    ρ_g::FT
    "Critical size separating spherical and nonspherical ice [m]"
    D_th::FT
    "Size of equal mass for graupel and unrimed ice [m] — `Inf` when `F_rim = 0`"
    D_gr::FT
    "Size of equal mass for graupel and partially rimed ice [m] — `Inf` when `F_rim = 0`"
    D_cr::FT
end

function P3Particles(params::CMP.ParametersP3, F_rim, ρ_rim)
    FT = UT.promote_typeof(F_rim, ρ_rim)
    (; mass, ρ_i) = params
    ρ_g = FT(get_ρ_g(mass, F_rim, ρ_rim))
    D_th = FT(get_D_th(mass, ρ_i))
    D_gr = FT(ifelse(iszero(F_rim), FT(Inf), get_D_gr(mass, ρ_g)))
    D_cr = FT(ifelse(iszero(F_rim), FT(Inf), get_D_cr(mass, F_rim, ρ_g)))
    return P3Particles(params, FT(F_rim), FT(ρ_rim), ρ_g, D_th, D_gr, D_cr)
end

@inline _primal_particles(particles) =
    P3Particles(map(FD.value, ntuple(i -> getfield(particles, i), Val(fieldcount(typeof(particles)))))...)

"""
    P3Distribution{FT, PT}

Size distribution of the P3 ice, `N′(D) = N₀ D^μ e^{-λD}`.

# Construction

    P3Distribution(particles, integration, ρq_ice, ρn_ice)

The construction has three steps:

1. The ice number is bounded by [`bounded_ice_number`](@ref), so the mean particle mass lies
   within the scheme limits.
2. The slope `logλ` is solved by [`get_distribution_logλ`](@ref) from the ice mass, the bounded
   number, and the mass relation of the particles.
3. The shape parameter `μ` follows from the slope law of `particles.params`, and the intercept
   `logN₀` from the bounded number.

The method dispatches on the type of `integration`; for a `QuadratureRule` it returns the
distribution above.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `integration`: how the size-distribution integrals are evaluated, a `QuadratureRule`
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]

# Automatic differentiation

`logλ` and `μ` are solved from the primal values and stored as primal values, so they are
constant under differentiation with respect to the state; `ρn` and `logN₀` are not.

# Fields
$(FIELDS)
"""
struct P3Distribution{FT <: Real, PT <: Real}
    "Volumetric ice number concentration [1/m³], with the mean particle mass bounded"
    ρn::FT
    "Shape parameter [-], a primal value"
    μ::PT
    "Log of the slope parameter [log(1/m)], a primal value"
    logλ::PT
    "Log of the intercept parameter [log(1/m⁴)]"
    logN₀::FT
end

function P3Distribution(particles::P3Particles, ::QuadratureRule, ρq_ice, ρn_ice)
    ρn = bounded_ice_number(particles.params, ρq_ice, ρn_ice)
    logλ = get_distribution_logλ(particles, ρq_ice, ρn)
    μ, logλ = promote(get_μ(particles.params.slope, logλ), logλ)
    ρn, logN₀ = promote(ρn, get_logN₀(ρn, μ, logλ))
    return P3Distribution(ρn, μ, logλ, logN₀)
end

"""
    P3State{FT}

State of the P3 ice in one cell.

The struct holds:

- the [`P3Particles`](@ref) `particles`,
- the evaluation of the size-distribution integrals, `integration`,
- the ice mass `ρq_ice`, the prognostic ice number `ρn_ice`, and the prognostic rime volume `ρb_rim`,
- the [`P3Distribution`](@ref) `dist`.

# Construction

    P3State(params, integration, ρq_ice, ρn_ice, ρq_rim, ρb_rim)

The constructor builds each part once:

1. The prognostic variables are clamped to be nonnegative.
2. The particles are built from the rime mass fraction `F_rim = ρq_rim / ρq_ice` and the rime
   density `ρ_rim`. For `ρ_rim`, the rime volume is moved to the nearest value for which
   `ρq_rim / ρb_rim` lies within [`rime_density_bounds`](@ref), see
   [`UT.nearest_admissible_b`](@ref).
3. The distribution is [`P3Distribution`](@ref)`(particles, integration, ρq_ice, ρn_ice)`.

The field `ρb_rim` keeps the prognostic rime volume, not the moved one.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `integration`: how the size-distribution integrals are evaluated, a `QuadratureRule`
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]
- `ρq_rim`: rime mass concentration [kg/m³]
- `ρb_rim`: rime volume concentration [m³/m³]

# Fields
$(FIELDS)
"""
struct P3State{FT, P <: P3Particles, R, D}
    "[`P3Particles`](@ref) of the ice"
    particles::P
    "How the size-distribution integrals are evaluated"
    integration::R
    "Volumetric ice mass concentration [kg/m³]"
    ρq_ice::FT
    "Volumetric prognostic ice number concentration [1/m³]"
    ρn_ice::FT
    "Volumetric prognostic rime volume [m³/m³]"
    ρb_rim::FT
    "Size distribution, a [`P3Distribution`](@ref) for a `QuadratureRule` integration"
    dist::D
end

function P3State(params::CMP.ParametersP3, integration, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    FT = UT.promote_typeof(ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    # Clamp to the physical domain
    ρq_ice = UT.clamp_to_nonneg(FT(ρq_ice))
    ρn_ice = UT.clamp_to_nonneg(FT(ρn_ice))
    ρq_rim = UT.clamp_to_nonneg(FT(ρq_rim))
    ρb_rim = UT.clamp_to_nonneg(FT(ρb_rim))
    ρb_rim_admissible = UT.nearest_admissible_b(ρq_rim, ρb_rim, rime_density_bounds(params))
    F_rim = min(UT.rime_mass_fraction(ρq_rim, ρq_ice), FT(1) - eps(FT))
    # ρ_rim ≤ ρ_i gives ρ_g ≤ ρ_i, which the threshold ordering D_th ≤ D_gr requires
    ρ_rim = min(UT.rime_density(ρq_rim, ρb_rim_admissible), params.ρ_i)
    particles = P3Particles(params, F_rim, ρ_rim)
    dist = P3Distribution(particles, integration, ρq_ice, ρn_ice)
    return P3State(particles, integration, ρq_ice, ρn_ice, ρb_rim, dist)
end

"""
    ice_mean_particle_mass_min(params)

Return the smallest mean ice particle mass [kg].

It is the mass of a newly nucleated crystal, a solid-ice sphere of diameter `D_nuc`.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
"""
@inline ice_mean_particle_mass_min((; ρ_i, D_nuc)) = ρ_i * CO.volume_sphere_D(D_nuc)

"""
    ice_mean_particle_mass_max(params)

Return the largest mean ice particle mass [kg], `params.m_max`.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
"""
@inline ice_mean_particle_mass_max((; m_max)) = m_max

"""
    ice_moments_are_admissible(params, ρq_ice, ρn_ice)

Return whether the ice mass [kg/m³] and number [1/m³] describe a population the scheme can represent.

Both must be positive, with the mean particle mass between [`ice_mean_particle_mass_min`](@ref)
and [`ice_mean_particle_mass_max`](@ref), the limits that [`bounded_ice_number`](@ref) applies.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]
"""
@inline function ice_moments_are_admissible(params, ρq_ice, ρn_ice)
    x_min = ice_mean_particle_mass_min(params)
    x_max = ice_mean_particle_mass_max(params)
    # Products instead of the mean mass, so no division by a vanishing number
    return (ρq_ice > 0) & (ρn_ice > 0) & (ρn_ice * x_min <= ρq_ice) & (ρq_ice <= ρn_ice * x_max)
end

"""
    rime_density_bounds(params)

Return the bounds `(ρ_rim_min, ρ_rim_max)` [kg/m³] of the bulk rime density.

The lower bound is `ρ′_rim(1)`, the smallest Cober and List (1993) local rime density, and the
upper bound is the solid-ice density `ρ_i`. The [`P3State`](@ref) construction moves the rime
volume so that the rime density lies between them.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
"""
@inline rime_density_bounds((; ρ_rim_local, ρ_i)) = (ρ_rim_local(one(ρ_i)), ρ_i)

"""
    bounded_ice_number(params, ρq_ice, ρn_ice)

Return the ice number concentration [1/m³] with the mean particle mass bounded.

The bounds are [`ice_mean_particle_mass_min`](@ref) and [`ice_mean_particle_mass_max`](@ref).
Without ice mass or without ice number, the result is zero.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]
"""
@inline function bounded_ice_number(params, ρq_ice, ρn_ice)
    FT = typeof(ρn_ice)
    x_min = FT(ice_mean_particle_mass_min(params))
    x_max = FT(ice_mean_particle_mass_max(params))
    n_bounded = clamp(ρn_ice, ρq_ice / x_max, ρq_ice / x_min)
    # Compare primal values, since a dual number with zero value and nonzero partials exceeds zero
    present = (FD.value(ρq_ice) > 0) & (FD.value(ρn_ice) > 0)
    return ifelse(present, n_bounded, zero(FT))
end

Base.show(io::IO, mime::MIME"text/plain", x::P3State) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
ShowMethods.field_units(::P3State) = (; ρq_ice = "kg/m³", ρn_ice = "1/m³", ρb_rim = "m³/m³")
Base.show(io::IO, mime::MIME"text/plain", x::P3Particles) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
ShowMethods.field_units(::P3Particles) = (; ρ_rim = "kg/m³", ρ_g = "kg/m³", D_th = "m", D_gr = "m", D_cr = "m")
Base.show(io::IO, mime::MIME"text/plain", x::P3Distribution) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
ShowMethods.field_units(::P3Distribution) = (; ρn = "1/m³", logλ = "log(1/m)", logN₀ = "log(1/m⁴)")

Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

Base.eltype(::P3Particles{FT}) where {FT} = FT
Base.broadcastable(particles::P3Particles) = tuple(particles)

Base.eltype(::P3Distribution{FT}) where {FT} = FT
Base.broadcastable(dist::P3Distribution) = tuple(dist)

"""
    isunrimed(particles)

Return `true` if the particle is unrimed, i.e. `F_rim = 0`.

# Arguments
- `particles`: [`P3Particles`](@ref)
"""
isunrimed(particles) = iszero(particles.F_rim)

@inline exprel1(x) = expm1(x) / x            # exprel₁ = (exp(x)-1)/x
@inline _exprel2(x) = (expm1(x) - x) / (x * x)
@inline _exprel2_small(x) =
    evalpoly(x, ntuple(i -> inv(oftype(x, UT.fac(i + 1))), Val(8)))
@inline function exprel2(x)                  # exprel₂ = (exp(x)-1-x)/x²
    abs(x) < oftype(x, 1 / 5) && return _exprel2_small(x)
    return _exprel2(x)
end

"""
    exprel(x, ::Val{k})

Compute the relative exponential `exprelₖ(x) = Σₙ xⁿ/(n+k)!`.

# Arguments
- `x`: real argument.
- `k`: order of the function, passed as `Val(k)`.

# Details

`exprelₖ` is one of the `φ`-functions `φₖ(x)` that appear in exponential integrators.
It is implemented for `k = 1` (`(eˣ-1)/x`) and `k = 2` (`(eˣ-1-x)/x²`); other values of `k`
throw an `ArgumentError`. For `k = 2` at small `|x|`, a Taylor series is used to avoid
catastrophic loss of precision near `x = 0`.

The value `k` is passed as `Val(k)`, so the order resolves at compile time.

# References
- [Exponential integrators](https://en.wikipedia.org/wiki/Exponential_integrator)
- [Niesen & Wright (2009), A Krylov subspace algorithm for evaluating the φ-functions appearing in exponential integrators](https://arxiv.org/abs/0907.4631)
"""
@inline function exprel(x, ::Val{k}) where {k}
    k == 1 && return exprel1(x)
    k == 2 && return exprel2(x)
    throw(ArgumentError("exprel is only implemented for k = 1 and k = 2"))
end

"""
    get_ρ_d(mass, F_rim, ρ_rim)

Return the density of the unrimed portion of the particle [kg/m³].

The result is the exact solution as a function of the rime mass fraction `F_rim`, the mass power
law parameters `mass`, and the rime density `ρ_rim`.

For the derivation of the numerically stable form used here, see the
([P3 scheme documentation](@ref P3-assumed-particle-size-relationships)).

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref) parameters
- `F_rim`: rime mass fraction
- `ρ_rim`: rime density

# Returns
- `ρ_d`: density of the unrimed portion of the particle [kg/m³]

# Examples

```jldoctest
julia> import CloudMicrophysics.Parameters as CMP,
              ClimaParams as CP,
              CloudMicrophysics.P3Scheme as P3

julia> FT = Float64;

julia> mass = CMP.MassPowerLaw(CP.create_toml_dict(FT));

julia> F_rim, ρ_rim = FT(0.5), FT(916.7);

julia> ρ_d = P3.get_ρ_d(mass, F_rim, ρ_rim)
488.9120789986414
```
"""
function get_ρ_d((; β_va), F_rim, ρ_rim)
    p = 1 / (3 - β_va)
    logFᵤ = log1p(-F_rim)            # = log(1 - F_rim)
    φ₁ = exprel(logFᵤ, Val(1))
    φ₁₋ₚ = exprel((1 - p) * logFᵤ, Val(1))
    H = -p * exprel(-p * logFᵤ, Val(2)) - (1 - p) * exprel((1 - p) * logFᵤ, Val(2))
    G = H - φ₁₋ₚ * φ₁
    return -(ρ_rim * φ₁ * φ₁₋ₚ) / G
end

"""
    get_ρ_g(mass, F_rim, ρ_rim)

Return the density of total (deposition + rime) ice mass for graupel [kg/m³]

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref) parameters
- `F_rim`: rime mass fraction (`L_rim / L_ice`) [-]
- `ρ_rim`: rime density (`L_rim / B_rim`) [kg/m³]

# Returns
- `ρ_g`: density of total (deposition + rime) ice mass for graupel [kg/m³]

# Notes:
The density of the unrimed portion of the particle is [`get_ρ_d`](@ref).
See Eq. 16 in [MorrisonMilbrandt2015](@cite).
"""
get_ρ_g(mass, F_rim, ρ_rim) = weighted_average(F_rim, ρ_rim, get_ρ_d(mass, F_rim, ρ_rim))

"""
    _get_threshold(params, ρ)

All thresholds are on the form

```math
\\left( \\frac{6α_{va}}{π ρ} \\right)^\\frac{1}{3 - β_{va}}
```

where for the different thresholds, `ρ` is:
- `D_th`: `ρ = ρ_i` (see [`get_D_th`](@ref))
- `D_gr`: `ρ = ρ_g` (see [`get_D_gr`](@ref))
- `D_cr`: `ρ = ρ_g * (1 - F_rim)` (see [`get_D_cr`](@ref))

# Arguments
- `params`: [`CMP.MassPowerLaw`](@ref) parameters
- `ρ`: (ice/graupel) density [kg/m³]
"""
_get_threshold((; α_va, β_va), ρ) =
# A floor on ρ far below any physical density keeps the fractional power finite
    (6α_va / (π * max(ρ, oftype(ρ, 1e-4))))^(1 / (3 - β_va))

"""
    get_D_th(mass, ρ_i)

Return the critical size separating spherical and nonspherical ice [meters]

See Eq. 8 in [MorrisonMilbrandt2015](@cite).

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref)
- `ρ_i`: density of solid ice [kg/m³]
"""
get_D_th(mass, ρ_i) = _get_threshold(mass, ρ_i)

"""
    get_D_gr(mass, ρ_g)

Return the size of equal mass for graupel and unrimed ice [meters]

See Eq. 15 in [MorrisonMilbrandt2015](@cite).

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref)
- `ρ_g`: graupel density [kg/m³]
"""
get_D_gr(mass, ρ_g) = _get_threshold(mass, ρ_g)

"""
    get_D_cr(mass, F_rim, ρ_g)

Return the size of equal mass for graupel and partially rimed ice [meters]

See Eq. 14 in [MorrisonMilbrandt2015](@cite).

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref)
- `F_rim`: rime mass fraction [-]
- `ρ_g`: graupel density [kg/m³]
"""
get_D_cr(mass, F_rim, ρ_g) = _get_threshold(mass, ρ_g * (1 - F_rim))

"""
    segment_boundaries(particles, D_min = 0, D_max = Inf)

Return the P3 mass-regime boundaries `(D_min, D_th, D_gr, D_cr, D_max)`.

The thresholds are clamped into the integration window `[D_min, D_max]`. The result is suitable as
the `bnds` argument to [`integrate`](@ref) / [`subintervals`](@ref).

If `F_rim = 0`, `particles.D_gr` and `particles.D_cr` are `Inf`; the clamp
collapses them to `D_max`, producing zero-width upper segments —
correct for the unrimed regime where only `(D_min, D_th)` and
`(D_th, D_max)` carry mass.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `D_min`: lower end of the integration window [m]
- `D_max`: upper end of the integration window [m]
"""
function segment_boundaries(particles, D_min = zero(eltype(particles)), D_max = eltype(particles)(Inf))
    D_th = clamp(particles.D_th, D_min, D_max)
    D_gr = clamp(particles.D_gr, D_min, D_max)
    D_cr = clamp(particles.D_cr, D_min, D_max)
    return (D_min, D_th, D_gr, D_cr, D_max)
end

"""
    weighted_average(f_a, a, b)

Return the weighted average of `a` and `b` with fraction `f_a`,

```math
f_a ⋅ a + (1 - f_a) ⋅ b
```
"""
function weighted_average(f_a, a, b)
    return f_a * a + (1 - f_a) * b
end

"""
    regime_value(particles, D, small, unrimed, dense_rimed, graupel, partially_rimed)

Select the value for the P3 mass/area regime that the maximum dimension `D` falls in.

# Arguments
 - `particles`:         [`P3Particles`](@ref)
 - `D`:                 maximum particle dimension [m]
 - `small`:             small spherical ice     (`D < D_th`),
 - `unrimed`:           large unrimed ice       (`F_rim = 0` ∧ `D_th ≤ D`),
 - `dense_rimed`:       dense rimed ice         (`D_th ≤ D < D_gr`),
 - `graupel`:           graupel (rimed)         (`D_gr ≤ D < D_cr`),
 - `partially_rimed`:   partially rimed ice     (`D_cr ≤ D`).

The five values are positional, in the order above. They are promoted to a
common type so the selection is concretely typed.
"""
@inline function regime_value(particles, D, small, unrimed, dense_rimed, graupel, partially_rimed)
    (; F_rim, D_th, D_gr, D_cr) = particles
    small, unrimed, dense_rimed, graupel, partially_rimed =
        promote(small, unrimed, dense_rimed, graupel, partially_rimed)
    #! format: off
    return ifelse(D < D_th,      small,             # small spherical ice
           ifelse(iszero(F_rim), unrimed,           # large nonspherical unrimed ice
           ifelse(D < D_gr,      dense_rimed,       # dense nonspherical rimed ice
           ifelse(D < D_cr,      graupel,           # graupel (rimed)
                                 partially_rimed,   # partially rimed ice
    ))))
    #! format: on
end

"""
    ice_mass_coeffs(particles, D)

Return the coefficients for the ice mass power law at diameter `D`.

# Arguments
 - `particles`: The [`P3Particles`](@ref)
 - `D`: maximum particle dimension [m]

# Returns
 - `(a, b)`: coefficients for the ice mass power law, `a D^b`
"""
function ice_mass_coeffs(particles, D)
    FT = eltype(particles)
    (; params, F_rim, ρ_g) = particles
    (; ρ_i) = params
    (; α_va, β_va) = params.mass
    # unrimed fraction; eps is the float resolution of 1 - F_rim near 1,
    # so the floor acts only at F_rim = 1
    Fu = max(1 - F_rim, eps(FT))
    #! format: off
    #                          small        unrimed  rimed  graupel      partially-rimed
    a = regime_value(particles, D, ρ_i * π / 6, α_va,    α_va,  ρ_g * π / 6, α_va / Fu)
    b = regime_value(particles, D, FT(3),       β_va,    β_va,  FT(3),       β_va)
    #! format: on
    return (a, b)
end

"""
    ice_mass(particles, D)

Return the mass of a particle with diameter `D`

# Arguments
 - `particles`: The [`P3Particles`](@ref)
 - `D`: maximum particle dimension [m]
"""
function ice_mass(particles, D)
    (a, b) = ice_mass_coeffs(particles, D)
    return a * D^b
end

"""
    ice_density(particles, D)

Return the density of a particle with diameter `D`

# Arguments
 - `particles`: The [`P3Particles`](@ref)
 - `D`: maximum particle dimension [m]

# Notes:
 The density of nonspherical particles is assumed to be the particle mass divided
 by the volume of a sphere with the same D [MorrisonMilbrandt2015](@cite).
"""
ice_density(particles, D) = ice_mass(particles, D) / CO.volume_sphere_D(D)

"""
    ice_area(particles, D)

Return the cross-sectional area of a particle of maximum dimension `D` [m²].

The area follows the mass regime that `D` falls in.

# Arguments
 - `particles`: [`P3Particles`](@ref)
 - `D`: maximum particle dimension [m]
"""
function ice_area(particles, D)
    (; params, F_rim) = particles
    (; γ, σ) = params.area
    spherical = D^2 * π / 4
    nonspherical = γ * D^σ
    return regime_value(
        particles, D, spherical, nonspherical, nonspherical, spherical,
        weighted_average(F_rim, spherical, nonspherical),
    )
end

"""
    ϕ_material_density(particles, D)

Return the material density of the solid of an ice particle at diameter `D` [kg/m³].

It is used in the aspect-ratio closure [`ϕᵢ`](@ref). This is `ρ_i` in every mass regime except
graupel (`D_gr ≤ D < D_cr`), where it is the graupel density `ρ_g`.

This is the density of the actual solid material, not the size-dependent effective
density `mᵢ / (π D³ / 6)` returned by [`ice_density`](@ref); see the
[P3 documentation](@ref "Aspect ratio") for the distinction.

# Arguments
- `particles`: [`P3Particles`](@ref)
- `D`: maximum particle dimension [m]
"""
function ϕ_material_density(particles, D)
    (; params, ρ_g) = particles
    (; ρ_i) = params
    # small, unrimed, dense_rimed → ρ_i; graupel → ρ_g; partially_rimed → ρ_i
    return regime_value(particles, D, ρ_i, ρ_i, ρ_i, ρ_g, ρ_i)
end

"""
    ϕᵢ(particles, D)

Return the oblate aspect ratio of an ice particle of maximum dimension `D`.

The aspect ratio is `ϕ = 3√π mᵢ / (4 ρ aᵢ^{3/2})` (`κ = 1/3`), with mass `mᵢ = ice_mass(particles, D)`,
projected area `aᵢ = ice_area(particles, D)`, and material density
`ρ = ϕ_material_density(particles, D)`. Assumes zero liquid fraction.

# Arguments
 - `particles`: The [`P3Particles`](@ref)
 - `D`: maximum dimension of ice particle [m]

See also [`ϕ_material_density`](@ref) and the
[aspect-ratio section of the P3 documentation](@ref "Aspect ratio") for the
spheroid derivation and the residual `ϕ > 1` band above `D_th`.
"""
@inline function ϕᵢ(particles, D)
    FT = eltype(D)
    mᵢ = ice_mass(particles, D)
    aᵢ = ice_area(particles, D)
    ρ = ϕ_material_density(particles, D)

    # Oblate aspect ratio (κ = 1/3)
    ϕ_ob = 3 * sqrt(FT(π)) * mᵢ / (4 * ρ * aᵢ * sqrt(aᵢ))
    #ϕ_pr = 16 * ρ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2)  # prolate, κ = -1/6

    return ifelse(iszero(D), zero(ϕ_ob), ϕ_ob)
end
