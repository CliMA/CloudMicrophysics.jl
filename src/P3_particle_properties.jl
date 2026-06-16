# TODO: Implement `F_liq` as another P3State-like struct
"""
    P3State{FT}

State of the P3 scheme.

This struct bundles the P3 parameterizations `params`, the provided rime state
(`F_rim`, `ρ_rim`), and the cached derived threshold variables
`thresholds = (; D_th, D_gr, D_cr, ρ_g)` — computed once at construction.

# Construction

Two positional, broadcast-friendly, GPU-clean entry points; both wrap the
default inner constructor and auto-compute `thresholds` so the cached
derivatives stay consistent with the rime state. The field type `FT` is the
promotion of the four state inputs' types — so e.g. `ForwardDiff.Dual`
prognostics yield a Dual-valued state over plain-float `params` — while
mixed-type direct construction of the full field set errors.

  - `P3State(params, L_ice, N_ice, F_rim, ρ_rim)` — direct construction
    from explicit `(F_rim, ρ_rim)`. Used by tests, audit prototypes, and
    documentation plots that want exact rime-state inputs. **Performs no
    validation.** Out-of-range inputs (`F_rim ∉ [0, 1)`, `ρ_rim ∉ [0,
    ρ_l]`) produce a state whose downstream evaluation typically yields
    NaN or non-physical values.
  - [`state_from_prognostic`](@ref) — production / upstream-model entry
    point. Accepts the volumetric prognostic variables `(ρq_ice,
    ρn_ice, ρq_rim, ρb_rim)`, regularises them into `(F_rim, ρ_rim)`,
    and returns the constructed state.

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
    "Rime density"
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

# Positional 5-arg constructor — derives the cached thresholds and
# graupel density from `(F_rim, ρ_rim)`. The production hot path
# (`state_from_prognostic`) wraps this with input regularisation;
# tests / docs call this directly with explicit `(F_rim, ρ_rim)`.
#
# When `F_rim = 0`, `ρ_g` collapses to `NaN` and `D_gr = D_cr = Inf`
# (the "no graupel regime" sentinel). Downstream code is expected to
# gate on the unrimed branch before reading those fields.
function P3State(params::CMP.ParametersP3, ρq_ice, ρn_ice, F_rim, ρ_rim)
    FT = UT.promote_typeof(ρq_ice, ρn_ice, F_rim, ρ_rim)
    (; mass, ρ_i) = params
    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    ρ_g = get_ρ_g(F_rim, ρ_rim, ρ_d)
    D_th = get_D_th(mass, ρ_i)
    D_gr = ifelse(iszero(F_rim), FT(Inf), get_D_gr(mass, ρ_g))
    D_cr = ifelse(iszero(F_rim), FT(Inf), get_D_cr(mass, F_rim, ρ_g))
    return P3State(params, FT.((ρq_ice, ρn_ice, F_rim, ρ_rim, ρ_g, D_th, D_gr, D_cr))...)
end

Base.show(io::IO, mime::MIME"text/plain", x::P3State) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
ShowMethods.field_units(::P3State) = (;
    ρq_ice = "kg/m³", ρn_ice = "1/m³", ρ_rim = "kg/m³",
    ρ_g = "kg/m³", D_th = "m", D_gr = "m", D_cr = "m",
)

"""
    state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)

Construct a [`P3State`](@ref) from the volumetric prognostic ice
variables directly, computing the (clamped, regularised) rime mass
fraction and rime density. **The production / upstream-model entry
point.** Positional and broadcast-friendly.

The regularised ratios come from `Utilities.rime_mass_fraction` and
`Utilities.rime_density`, which smoothly go to zero when their
denominators are near machine precision — avoiding the hard discontinuity
at `q_ice = ϵ` / `b_rim = ϵ`. The upper clamps `F_rim < 1 − ε` and
`ρ_rim ≤ 0.8·ρ_l ≈ 730 kg/m³` keep the result inside the domain of validity
of the threshold formulas evaluated by the [`P3State`](@ref) constructor.

!!! note "TODO — revisit the `ρ_rim ≤ 0.8·ρ_l` cap"
    The closed-form graupel density `ρ_g = F_rim·ρ_rim + (1−F_rim)·ρ_d`
    can mathematically exceed `ρ_l` — `ρ_d` (the unrimed portion's
    density, [`get_ρ_d`](@ref)) is linear in `ρ_rim` with no built-in
    upper clamp, so feeding `ρ_rim` near `ρ_l` can produce `ρ_g > ρ_l`.
    That breaks the threshold ordering `D_th < D_gr < D_cr` that the P3
    partitioning assumes (`D_gr ∝ (6α_va/(π·ρ_g))^{1/(3−β_va)}` shrinks
    as `ρ_g` grows; eventually `D_gr < D_th`). The 0.8 cushion keeps
    `ρ_g` comfortably below `ρ_l` for the realistic `(F_rim, ρ_rim)`
    regime — Macklin-type rime-density formulations rarely give
    `ρ_rim > 700 kg/m³` anyway, so the cushion costs us nothing
    physically. To lift the cap to `ρ_l` we'd need to (i) explicitly
    bound `ρ_g` (e.g. `min(ρ_g, ρ_l)`) or rederive `ρ_d` so it's
    monotone-bounded by `ρ_l`, and (ii) accept that bulk rime densities
    800–917 kg/m³ are off the calibration domain of the original P3 fit.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]
- `ρq_rim`: rime mass concentration [kg/m³]
- `ρb_rim`: rime volume concentration [m³/m³]
"""
function state_from_prognostic(params::CMP.ParametersP3, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    FT = eltype(ρq_ice)
    F_rim = min(UT.rime_mass_fraction(ρq_rim, ρq_ice), one(FT) - eps(FT))
    ρ_rim = min(UT.rime_density(ρq_rim, ρb_rim), FT(0.8) * params.ρ_l)
    return P3State(params, ρq_ice, ρn_ice, F_rim, ρ_rim)
end

Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

"""
    isunrimed(state::P3State)

Return `true` if the particle is unrimed, i.e. `F_rim = 0`.
"""
isunrimed(state::P3State) = iszero(state.F_rim)

exprel_coeffs(n, k) = 1 // factorial(n + k)  # exprelₖ series coefficient = 1/(n+k)!

@inline exprel1(x) = expm1(x) / x            # exprel₁ = (exp(x)-1)/x
@inline _exprel2(x) = (expm1(x) - x) / (x * x)
@inline _exprel2_small(x) = evalpoly(x, ntuple(np1 -> exprel_coeffs(np1 - 1, 2), Val(8)))
@inline function exprel2(x)                  # exprel₂ = (exp(x)-1-x)/x²
    abs(x) < 1 // 5 && return _exprel2_small(x)
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

For the derivation of the numerically stable form used here, see the P3 scheme
documentation ([Assumed particle size relationships](@ref)).

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
_get_threshold((; α_va, β_va)::CMP.MassPowerLaw, ρ) = (6α_va / (π * ρ))^(1 / (3 - β_va))

"""
    get_D_th(mass::MassPowerLaw, ρ_i)

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

Return the flat 5-tuple `(D_min, D_th, D_gr, D_cr, D_max)` of P3 mass-
regime boundaries clamped into the requested integration window
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
    ice_mass_coeffs(state::P3State, D)

Return the coefficients for the ice mass power law at diameter `D`.

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum particle dimension [m]

# Returns
 - `(a, b)`: coefficients for the ice mass power law, `a D^b`
"""
function ice_mass_coeffs(state::P3State, D)
    FT = promote_type(eltype(state), eltype(D))
    (; params, F_rim, ρ_g, D_th, D_gr, D_cr) = state
    (; ρ_i) = params
    (; α_va, β_va) = params.mass
    a, b = if D < D_th       # small spherical ice
        (ρ_i * π / 6, 3)
    elseif iszero(F_rim)     # large nonspherical unrimed ice
        (α_va, β_va)
    elseif D_th ≤ D < D_gr   # dense nonspherical rimed ice
        (α_va, β_va)
    elseif D_gr ≤ D < D_cr   # graupel (rimed)
        (ρ_g * π / 6, 3)
    else # D_cr ≤ D          # partially rimed ice
        (α_va / (1 - F_rim), β_va)
    end
    return FT(a), FT(b)
end

"""
    ice_mass(state, D)

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
 Needed for aspect ratio calculation, so we assume zero liquid fraction.
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
    (; params, F_rim, D_th, D_gr, D_cr) = state
    (; γ, σ) = params.area
    spherical_area(D) = D^2 * π / 4
    nonspherical_area(D) = γ * D^σ
    return if D < D_th       # small spherical ice
        spherical_area(D)
    elseif iszero(F_rim)     # large nonspherical unrimed ice
        nonspherical_area(D)
    elseif D_th ≤ D < D_gr   # dense nonspherical rimed ice
        nonspherical_area(D)
    elseif D_gr ≤ D < D_cr   # graupel (rimed)
        spherical_area(D)
    else # D_cr ≤ D          # partially rimed ice
        weighted_average(F_rim, spherical_area(D), nonspherical_area(D))
    end
end

"""
    ϕᵢ(state::P3State, D)

Returns the aspect ratio (ϕ) for an ice particle with diameter `D`

# Arguments
 - `state`: The [`P3State`](@ref)
 - `D`: maximum dimension of ice particle [m]

# Notes
 The density of nonspherical particles is assumed to be equal to the particle mass
 divided by the volume of a spherical particle with the same D_max [MorrisonMilbrandt2015](@cite).
 Assuming zero liquid fraction and oblate shape.
"""
function ϕᵢ(state::P3State, D)
    FT = eltype(D)
    mᵢ = ice_mass(state, D)
    aᵢ = ice_area(state, D)
    ρᵢ = mᵢ / CO.volume_sphere_D(D)

    # TODO - prolate or oblate?
    ϕ_ob = min(1, 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^FT(1.5))) # κ =  1/3
    #ϕ_pr = max(1, 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2))       # κ = -1/6

    return ifelse(D == 0, zero(ϕ_ob), ϕ_ob)
end
