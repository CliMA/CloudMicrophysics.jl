# TODO: Implement `F_liq` as another P3State-like struct
"""
    P3State{FT}

State of the P3 scheme.

This struct bundles the P3 parameterizations `params`, the provided rime state
(`F_rim`, `ρ_rim`), and the cached derived threshold variables
`D_th`, `D_gr`, `D_cr`, and `ρ_g` — computed once at construction.

# Construction

  - [`state_from_prognostic`](@ref): Main entry point.
    Accepts the volumetric prognostic variables `(ρq_ice, ρn_ice, ρq_rim, ρb_rim)`,
    regularises them into `(F_rim, ρ_rim)`, and returns the constructed state.

# Fields
$(FIELDS)
"""
struct P3State{FT, PARAMS <: CMP.ParametersP3}
    "[`CMP.ParametersP3`](@ref) object"
    params::PARAMS

    "Volumetric ice mass concentration [kg/m³]"
    ρq_ice::FT
    "Volumetric ice number concentration [1/m³]"
    ρn_ice::FT
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

function P3State(params::CMP.ParametersP3, ρq_ice, ρn_ice, F_rim, ρ_rim)
    FT = UT.promote_typeof(ρq_ice, ρn_ice, F_rim, ρ_rim)
    (; mass, ρ_i) = params
    # Clamp to the physical domain so the threshold formulas never evaluate a
    # power of a non-positive density. Inert on physical inputs.
    ρq_ice = UT.clamp_to_nonneg(FT(ρq_ice))
    ρn_ice = UT.clamp_to_nonneg(FT(ρn_ice))
    F_rim = clamp(FT(F_rim), FT(0), FT(1) - eps(FT))
    # Bound the rime density by the solid-ice density. Since `ρ_g ≤ ρ_rim`, this
    # keeps `ρ_g ≤ ρ_i`, preserving the `D_th ≤ D_gr` threshold ordering.
    ρ_rim = clamp(FT(ρ_rim), FT(0), ρ_i)
    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    # Floored to a numerically small positive value, far below any physical ice or
    # graupel density (order 100-917 kg/m³): a nonpositive `ρ_g` would give a DomainError
    # in `log(ρ_g·π/6)`, the graupel mass coefficient in the size-distribution moment
    # integral, and a negative aspect-ratio material density. Inert on physical inputs;
    # only engages for degenerate near-zero-mass states.
    ρ_g = max(get_ρ_g(F_rim, ρ_rim, ρ_d), oftype(F_rim, 1e-4))
    D_th = get_D_th(mass, ρ_i)
    D_gr = ifelse(iszero(F_rim), FT(Inf), get_D_gr(mass, ρ_g))
    D_cr = ifelse(iszero(F_rim), FT(Inf), get_D_cr(mass, F_rim, ρ_g))
    return P3State(
        params,
        ρq_ice, ρn_ice, F_rim, ρ_rim,
        FT(ρ_g), FT(D_th), FT(D_gr), FT(D_cr),
    )
end

Base.show(io::IO, mime::MIME"text/plain", x::P3State) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
ShowMethods.field_units(::P3State) = (;
    ρq_ice = "kg/m³", ρn_ice = "1/m³", ρ_rim = "kg/m³",
    ρ_g = "kg/m³", D_th = "m", D_gr = "m", D_cr = "m",
)

"""
    state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)

Construct a [`P3State`](@ref) from the volumetric prognostic ice variables directly,
computing the (clamped, regularised) rime mass fraction and rime density.

`(ρq_rim, ρb_rim)` is projected onto the admissible density cone
`ρ_rim_min·b ≤ q ≤ ρ_rim_max·b` via [`UT.nearest_admissible_b`](@ref) before
either ratio is computed, rather than clamped component-wise: an independent
clamp can leave positive mass paired with zero volume or the reverse, an
orphan pair for which the mass-to-volume ratio computed below does not
describe the state the mass and number moments still carry.

The regularised ratios come from [`UT.rime_mass_fraction`](@ref) and
[`UT.rime_density`](@ref), which smoothly go to zero when their
denominators are near machine precision, avoiding the discontinuity
at `q_ice = ϵ` / `b_rim = ϵ`. The upper clamps `F_rim < 1 - ε` and
`ρ_rim ≤ ρ_i` (the solid-ice density) keep the result inside the domain of the
threshold formulas evaluated by the [`P3State`](@ref) constructor. The rime
density bound is ordering-critical: since `ρ_g = F_rim·ρ_rim + (1-F_rim)·ρ_d`
satisfies `ρ_g ≤ ρ_rim`, bounding `ρ_rim ≤ ρ_i` gives `ρ_g ≤ ρ_i`, which the
`D_th ≤ D_gr` regime ordering requires. Only an unphysical `ρ_rim > ρ_i` can
breach it.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]
- `ρq_rim`: rime mass concentration [kg/m³]
- `ρb_rim`: rime volume concentration [m³/m³]
"""
function state_from_prognostic(params::CMP.ParametersP3, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    # Floor the prognostic moments so the regularised ratios stay non-negative;
    # F_rim and ρ_rim are bounded in the `P3State` constructor. `(ρq_rim, ρb_rim)`
    # is projected onto the admissible density cone jointly, not clamped
    # component-wise, so a violating pair cannot reach `rime_density` as an
    # orphan (positive mass, zero volume, or the reverse).
    ρq_ice = UT.clamp_to_nonneg(ρq_ice)
    ρn_ice = UT.clamp_to_nonneg(ρn_ice)
    ρq_rim = UT.clamp_to_nonneg(ρq_rim)
    ρ_rim_min, ρ_rim_max = rime_density_bounds(params)
    ρb_rim = UT.nearest_admissible_b(ρq_rim, ρb_rim, ρ_rim_min, ρ_rim_max)
    F_rim = UT.rime_mass_fraction(ρq_rim, ρq_ice)
    ρ_rim = UT.rime_density(ρq_rim, ρb_rim)
    return P3State(params, ρq_ice, ρn_ice, F_rim, ρ_rim)
end

"""
    rime_density_bounds(params::CMP.ParametersP3)

The interval `[ρ_rim_min, ρ_rim_max]` that the bulk rime density `ρq_rim/ρb_rim` is confined to
whenever rime is present, returned as a tuple.

The pair `(ρq_rim, ρb_rim)` describes one population twice, once by mass and once by volume, so
the physicality of their quotient is a consistency property of every increment pair: an addition
at an implied density inside the interval moves the quotient toward that density and never past
it (the mediant property), and a removal along the ray through the origin of the pair plane
preserves it exactly. With every source depositing inside the interval and every sink ray-form,
the interval is forward-invariant under the microphysics and no downstream clamp does physical
work. Transport is not covered: non-monotone horizontal operators can still leave it at
composition gradients.

Both endpoints are DERIVED from parameters the scheme already carries, so no new tuned constant
enters and they move with a recalibration.

  - `ρ_rim_max = ρ_i`, the solid-ice density: no parcel of rime can be denser than the lattice it
    freezes into. Bigg immersion, rain freezing and wet-growth densification all deposit at
    exactly this endpoint.
  - `ρ_rim_min = ρ′_rim(1)`, the Cober and List (1993) local rime density at the lower end of its
    own CLAMPED domain `Rᵢ ∈ [1, 12]`. Dry riming is the only source below solid ice, so with
    ray-form sinks nothing softer is reachable.

# What the lower endpoint is a property of

It is the infimum of the implied densities of the ADMITTED SOURCES under the clamped CL93 closure:
a property of THIS SCHEME'S CLOSURE SET, not a statement about rime. Rime softer than 159.5 kg/m³
occurs in nature, and the P3 Fortran's `ρ_rimeMin = 50 kg/m³` exists for that reason. This scheme
cannot produce it, because its only sub-solid-ice source is dry riming through a law whose argument
is clamped at `Rᵢ = 1`. The interval states what is REACHABLE, which is exactly what a
forward-invariance argument needs: a bound the sources cannot cross.

The two ways to move it are therefore not equivalent. Extending the CL93 closure past
`Rᵢ ∈ [1, 12]` and calibrating it against observations is the legitimate route, because it changes
what the scheme can PRODUCE and the endpoint then follows from the changed law without being
touched. Adopting the Fortran's 50 kg/m³ instead would take a lookup-table AXIS bound as a physical
one, widening the admissible set without widening the producing set - declaring states physical
that no process here can make, which is what the invariant exists to exclude. That is a physics
decision rather than a refactor and is deliberately not taken.

The endpoints follow the derivation automatically, being computed from [`LocalRimeDensity`](@ref
CloudMicrophysics.Parameters.LocalRimeDensity) rather than written down. If the law and this
derivation ever fall out of step the falsifier fails loudly rather than drifting, because it
asserts that `ρ_rim_min` IS the law's infimum over its whole clamped domain, not merely a number
that happens to sit below it.

Calibration hook, low priority: the clamp range `Rᵢ ∈ [1, 12]` and the CL93 coefficients are
calibration candidates and should eventually be named parameters rather than hardcoded, so that a
calibrated closure moves the interval without a code edit.

The value is needed by falsifier bounds, documentation and any restorer, not by any process rate.
"""
@inline function rime_density_bounds(params::CMP.ParametersP3)
    (; ρ_rim_local, ρ_i) = params
    return (ρ_rim_local(one(ρ_i)), ρ_i)
end

Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

"""
    isunrimed(state::P3State)

Return `true` if the particle is unrimed, i.e. `F_rim = 0`.

Reads the VALUE lane of `F_rim`: under ForwardDiff a differentiated zero (a `Dual` whose value is
zero at `q_rim = 0` but whose partials are not) fails `iszero` and selects the rimed branch, where
the `ρ_g = NaN` threshold sentinel then enters the thresholds. On plain floats this is identical.
"""
isunrimed(state::P3State) = iszero(FD.value(state.F_rim))

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
    get_ρ_d(mass::MassPowerLaw, F_rim, ρ_rim)
    get_ρ_d(state::P3State)

Exact solution for the density of the unrimed portion of the particle as
    function of the rime mass fraction `F_rim`, mass power law parameters `mass`,
    and rime density `ρ_rim`.

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
function get_ρ_d((; β_va)::CMP.MassPowerLaw, F_rim, ρ_rim)
    p = 1 / (3 - β_va)
    logFᵤ = log1p(-F_rim)            # = log(1 - F_rim)
    φ₁ = exprel(logFᵤ, Val(1))
    φ₁₋ₚ = exprel((1 - p) * logFᵤ, Val(1))
    H = -p * exprel(-p * logFᵤ, Val(2)) - (1 - p) * exprel((1 - p) * logFᵤ, Val(2))
    G = H - φ₁₋ₚ * φ₁
    return -(ρ_rim * φ₁ * φ₁₋ₚ) / G
end
get_ρ_d((; params, F_rim, ρ_rim)::P3State) = get_ρ_d(params.mass, F_rim, ρ_rim)

"""
    get_ρ_g(F_rim, ρ_rim, ρ_d)
    get_ρ_g(mass::MassPowerLaw, F_rim, ρ_rim)
    get_ρ_g(state::P3State)

Return the density of total (deposition + rime) ice mass for graupel [kg/m³]

# Arguments
- `F_rim`: rime mass fraction (`L_rim / L_ice`) [-]
- `ρ_rim`: rime density (`L_rim / B_rim`) [kg/m³]
- `ρ_d`: density of the unrimed portion of the particle [kg/m³], see [`get_ρ_d`](@ref)

# Returns
- `ρ_g`: density of total (deposition + rime) ice mass for graupel [kg/m³]

# Notes:
See Eq. 16 in [MorrisonMilbrandt2015](@cite).
"""
get_ρ_g(F_rim, ρ_rim, ρ_d) = weighted_average(F_rim, ρ_rim, ρ_d)
function get_ρ_g(mass::CMP.MassPowerLaw, F_rim, ρ_rim)
    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    return get_ρ_g(F_rim, ρ_rim, ρ_d)
end
get_ρ_g((; params, F_rim, ρ_rim)::P3State) = get_ρ_g(params.mass, F_rim, ρ_rim)

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
# The density is floored to a numerically small positive value, far below any physical ice
# or graupel density (order 100-917 kg/m³): the threshold is a fractional power of `1/ρ`, so
# a nonpositive `ρ` would give a DomainError. Inert on physical inputs; only engages for
# degenerate near-zero-mass states.
_get_threshold((; α_va, β_va)::CMP.MassPowerLaw, ρ) =
    (6α_va / (π * max(ρ, oftype(ρ, 1e-4))))^(1 / (3 - β_va))

"""
    get_D_th(mass::MassPowerLaw, ρ_i)
    get_D_th(params::ParametersP3)

Return the critical size separating spherical and nonspherical ice [meters]

See Eq. 8 in [MorrisonMilbrandt2015](@cite).
"""
get_D_th(mass::CMP.MassPowerLaw, ρ_i) = _get_threshold(mass, ρ_i)
get_D_th((; mass, ρ_i)::CMP.ParametersP3) = get_D_th(mass, ρ_i)

"""
    get_D_gr(mass::MassPowerLaw, ρ_g)

Return the size of equal mass for graupel and unrimed ice [meters]

See Eq. 15 in [MorrisonMilbrandt2015](@cite).
"""
get_D_gr(mass::CMP.MassPowerLaw, ρ_g) = _get_threshold(mass, ρ_g)

"""
    get_D_cr(mass::MassPowerLaw, F_rim, ρ_g)

Return the size of equal mass for graupel and partially rimed ice [meters]

See Eq. 14 in [MorrisonMilbrandt2015](@cite).
"""
get_D_cr(mass::CMP.MassPowerLaw, F_rim, ρ_g) = _get_threshold(mass, ρ_g * (1 - F_rim))

"""
    segment_boundaries(state::P3State, D_min = 0, D_max = Inf)

Return the 5-tuple `(D_min, D_th, D_gr, D_cr, D_max)` of P3 mass-regime
boundaries clamped into the requested integration window
`[D_min, D_max]`. Suitable as the `bnds` argument to
[`integrate`](@ref) / [`subintervals`](@ref).

If `F_rim = 0`, `state.D_gr` and `state.D_cr` are `Inf`; the clamp
collapses them to `D_max`, producing zero-width upper segments —
correct for the unrimed regime where only `(D_min, D_th)` and
`(D_th, D_max)` carry mass.
"""
function segment_boundaries(state::P3State{FT}, D_min = FT(0), D_max = FT(Inf)) where {FT}
    D_th = clamp(state.D_th, D_min, D_max)
    D_gr = clamp(state.D_gr, D_min, D_max)
    D_cr = clamp(state.D_cr, D_min, D_max)
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
    regime_value(state::P3State, D, small, unrimed, dense_rimed, graupel, partially_rimed)

Select the value for the P3 mass/area regime that the maximum dimension `D` falls in:
 - `small`:             small spherical ice     (`D < D_th`),
 - `unrimed`:           large unrimed ice       (`F_rim = 0` ∧ `D_th ≤ D`),
 - `dense_rimed`:       dense rimed ice         (`D_th ≤ D < D_gr`),
 - `graupel`:           graupel (rimed)         (`D_gr ≤ D < D_cr`),
 - `partially_rimed`:   partially rimed ice     (`D_cr ≤ D`).

The five values are positional, in the order above. They are promoted to a
common type so the selection is concretely typed.
"""
@inline function regime_value(state::P3State, D, small, unrimed, dense_rimed, graupel, partially_rimed)
    (; F_rim, D_th, D_gr, D_cr) = state
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
    ice_mass_coeffs(state::P3State, D)

Return the coefficients for the ice mass power law at diameter `D`.

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum particle dimension [m]

# Returns
 - `(a, b)`: coefficients for the ice mass power law, `a D^b`
"""
function ice_mass_coeffs(state::P3State, D)
    FT = eltype(state)
    (; params, F_rim, ρ_g) = state
    (; ρ_i) = params
    (; α_va, β_va) = params.mass
    ϵB = UT.ϵ_numerics_P3_B(FT)
    Fu = max(1 - F_rim, ϵB)  # unrimed fraction
    #! format: off
    #                          small        unrimed  rimed  graupel      partially-rimed
    a = regime_value(state, D, ρ_i * π / 6, α_va,    α_va,  ρ_g * π / 6, α_va / Fu)
    b = regime_value(state, D, FT(3),       β_va,    β_va,  FT(3),       β_va)
    #! format: on
    return (a, b)
end

"""
    ice_mass(state::P3State, D)

Return the mass of a particle with diameter `D`

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum particle dimension [m]
"""
function ice_mass(state::P3State, D)
    (a, b) = ice_mass_coeffs(state, D)
    return a * D^b
end

"""
    ice_density(state::P3State, D)

Return the density of a particle with diameter `D`

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum particle dimension [m]

# Notes:
 The density of nonspherical particles is assumed to be the particle mass divided
 by the volume of a sphere with the same D [MorrisonMilbrandt2015](@cite).
"""
ice_density(state::P3State, D) = ice_mass(state, D) / CO.volume_sphere_D(D)

function get_∂mass_∂D_coeffs(state::P3State, D)
    (a, b) = ice_mass_coeffs(state, D)
    return a * b, b - 1
end

"""
    ∂ice_mass_∂D(state::P3State, D)

Return the derivative of the ice mass with respect to the particle diameter `D`.

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum particle dimension [m]
"""
function ∂ice_mass_∂D(state::P3State, D)
    (a, b) = get_∂mass_∂D_coeffs(state, D)
    return a * D^b
end

"""
    ice_area(state::P3State, D)

Return the cross-sectional area of a particle based on where it falls in the
    particle-size-based properties regime.

# Arguments
 - `state`: [`P3State`](@ref) object
 - `D`: maximum particle dimension [m]
"""
function ice_area(state::P3State, D)
    (; params, F_rim) = state
    (; γ, σ) = params.area
    spherical = D^2 * π / 4
    nonspherical = γ * D^σ
    return regime_value(
        state, D, spherical, nonspherical, nonspherical, spherical,
        weighted_average(F_rim, spherical, nonspherical),
    )
end

"""
    ϕ_material_density(state::P3State, D)

Return the material density of the ice particle's solid at diameter `D`, used in
the aspect-ratio closure [`ϕᵢ`](@ref). This is `ρ_i` in every mass regime except
graupel (`D_gr ≤ D < D_cr`), where it is the graupel density `ρ_g`.

This is the density of the actual solid material, not the size-dependent effective
density `mᵢ / (π D³ / 6)` returned by [`ice_density`](@ref); see the
[P3 documentation](@ref "Aspect ratio") for the distinction.
"""
function ϕ_material_density(state::P3State, D)
    (; params, ρ_g) = state
    (; ρ_i) = params
    # small, unrimed, dense_rimed → ρ_i; graupel → ρ_g; partially_rimed → ρ_i
    return regime_value(state, D, ρ_i, ρ_i, ρ_i, ρ_g, ρ_i)
end

"""
    ϕᵢ(state::P3State, D)

Return the oblate aspect ratio `ϕ = 3√π mᵢ / (4 ρ aᵢ^{3/2})` (`κ = 1/3`) for an
ice particle of maximum dimension `D`, with mass `mᵢ = ice_mass(state, D)`,
projected area `aᵢ = ice_area(state, D)`, and material density
`ρ = ϕ_material_density(state, D)`. Assumes zero liquid fraction.

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum dimension of ice particle [m]

See also [`ϕ_material_density`](@ref) and the
[aspect-ratio section of the P3 documentation](@ref "Aspect ratio") for the
spheroid derivation and the residual `ϕ > 1` band above `D_th`.
"""
@inline function ϕᵢ(state::P3State, D)
    FT = eltype(D)
    mᵢ = ice_mass(state, D)
    aᵢ = ice_area(state, D)
    ρ = ϕ_material_density(state, D)

    # Oblate aspect ratio (κ = 1/3)
    ϕ_ob = 3 * sqrt(FT(π)) * mᵢ / (4 * ρ * aᵢ * sqrt(aᵢ))
    #ϕ_pr = 16 * ρ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2)  # prolate, κ = -1/6

    return ifelse(iszero(D), zero(ϕ_ob), ϕ_ob)
end
