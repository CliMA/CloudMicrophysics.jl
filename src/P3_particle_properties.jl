# TODO: Implement `F_liq` as another P3State-like struct
"""
    P3State{FT}

State of the P3 scheme.

This struct bundles the P3 parameterizations `params`, the provided rime state
(`F_rim`, `ρ_rim`), and the cached derived threshold variables
`thresholds = (; D_th, D_gr, D_cr, ρ_g)` — computed once at construction.

To obtain a `P3State` object, use the [`get_state`](@ref) function (which
validates rime inputs) or the constructor `P3State(params, L_ice, N_ice,
F_rim, ρ_rim)` (no validation).

# Fields
$(FIELDS)
"""
struct P3State{FT, PARAMS <: CMP.ParametersP3{FT}, THRESH}
    "[`CMP.ParametersP3`](@ref) object"
    params::PARAMS

    "The ice mass concentration [kg/m³]"
    L_ice::FT
    "The ice number concentration [1/m³]"
    N_ice::FT
    "Rime mass fraction"
    F_rim::FT
    "Rime density"
    ρ_rim::FT

    "Cached derived thresholds `(; D_th, D_gr, D_cr, ρ_g)` — see [`get_thresholds_ρ_g`](@ref)"
    thresholds::THRESH
end

# Positional convenience constructor — auto-computes thresholds from
# (params, F_rim, ρ_rim). Use this or `get_state` instead of calling the
# 6-field inner constructor directly, so cached thresholds stay consistent
# with the rime state.
function P3State(params::CMP.ParametersP3, L_ice, N_ice, F_rim, ρ_rim)
    thresholds = get_thresholds_ρ_g(params, F_rim, ρ_rim)
    return P3State(params, L_ice, N_ice, F_rim, ρ_rim, thresholds)
end

# Keyword convenience (preserves the @kwdef-style call site).
P3State(; params, L_ice, N_ice, F_rim, ρ_rim) = P3State(params, L_ice, N_ice, F_rim, ρ_rim)

Base.show(io::IO, mime::MIME"text/plain", x::P3State) =
    ShowMethods.verbose_show_type_and_fields(io, mime, x)
ShowMethods.field_units(::P3State) = (; L_ice = "kg/m³", N_ice = "1/m³", ρ_rim = "kg/m³")

"""
    get_state(params; L_ice, N_ice, F_rim, ρ_rim)

Create a [`P3State`](@ref) from [`CMP.ParametersP3`](@ref) and rime state parameters.

# Arguments
 - `params`: [`CMP.ParametersP3`](@ref) object

# Keyword arguments
 - `L_ice`: ice mass concentration [kg/m³]
 - `N_ice`: ice number concentration [1/m³]
 - `F_rim`: rime mass fraction [-], `F_rim = L_rim / L_ice`
 - `ρ_rim`: rime density [kg/m³],   `ρ_rim = L_rim / B_rim`

# Notes

The returned state caches `thresholds = (; D_th, D_gr, D_cr, ρ_g)` —
computed once via [`get_thresholds_ρ_g`](@ref)
"""
function get_state(params::CMP.ParametersP3; L_ice, N_ice, F_rim, ρ_rim)
    # Rime mass fraction must always be non-negative AND strictly < 1
    # (fully-rimed ice would leave no unrimed part).
    0 ≤ F_rim < 1 || throw(DomainError(F_rim,
        "Rime mass fraction, `F_rim`, must be between 0 and 1"))
    # Rime density must be non-negative; as a bulk ice density it cannot
    # exceed the density of liquid water `ρ_l`.
    0 ≤ ρ_rim ≤ params.ρ_l || throw(DomainError(ρ_rim,
        "Rime density, `ρ_rim`, must be between 0 and ρ_l"))
    return P3State(; params, L_ice, N_ice, F_rim, ρ_rim)
end

"""
    get_state_from_prognostic(params, ρq_ice, ρn_ice, ρq_rim, ρb_rim)

Construct a [`P3State`](@ref) from the prognostic ice variables directly,
computing the (clamped, regularised) rime mass fraction and rime density.


The regularised ratios come from [`Utilities.rime_mass_fraction`](@ref) and
[`Utilities.rime_density`](@ref), which smoothly go to zero when their
denominators are near machine precision — avoiding the hard discontinuity
at `q_ice = ϵ` / `b_rim = ϵ`. The upper clamps `F_rim < 1 − ε` and
`ρ_rim ≤ 0.8·ρ_l ≈ 730 kg/m³` keep the result inside the domain of validity of
[`get_thresholds_ρ_g`](@ref).

!!! note "TODO — revisit the `ρ_rim ≤ 0.8·ρ_l` cap"
    The validated constructor [`get_state`](@ref) allows `ρ_rim ≤ ρ_l`.
    We use the tighter `0.8·ρ_l` here because the closed-form graupel
    density `ρ_g = F_rim·ρ_rim + (1−F_rim)·ρ_d` can mathematically exceed
    `ρ_l` — `ρ_d` (the unrimed portion's density, [`get_ρ_d`](@ref)) is
    linear in `ρ_rim` with no built-in upper clamp, so feeding `ρ_rim`
    near `ρ_l` can produce `ρ_g > ρ_l`. That breaks the threshold
    ordering `D_th < D_gr < D_cr` that the P3 partitioning assumes
    (`D_gr ∝ (6α_va/(π·ρ_g))^{1/(3−β_va)}` shrinks as `ρ_g` grows;
    eventually `D_gr < D_th`). The 0.8 cushion keeps `ρ_g` comfortably
    below `ρ_l` for the realistic `(F_rim, ρ_rim)` regime — Macklin-type
    rime-density formulations rarely give `ρ_rim > 700 kg/m³` anyway,
    so the cushion costs us nothing physically. To lift the cap to
    `ρ_l` we'd need to (i) explicitly bound `ρ_g` (e.g. `min(ρ_g, ρ_l)`)
    or rederive `ρ_d` so it's monotone-bounded by `ρ_l`, and
    (ii) accept that bulk rime densities 800–917 kg/m³ are off the
    calibration domain of the original P3 fit.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `ρq_ice`: ice mass concentration [kg/m³]
- `ρn_ice`: ice number concentration [1/m³]
- `ρq_rim`: rime mass concentration [kg/m³]
- `ρb_rim`: rime volume concentration [m³/m³]
"""
function get_state_from_prognostic(params::CMP.ParametersP3, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    FT = typeof(ρq_ice)
    F_rim = min(UT.rime_mass_fraction(ρq_rim, ρq_ice), one(FT) - eps(FT))
    ρ_rim = min(UT.rime_density(ρq_rim, ρb_rim),       FT(0.8) * params.ρ_l)
    return P3State(params, ρq_ice, ρn_ice, F_rim, ρ_rim)
end

Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

"""
    isunrimed(state::P3State)

Return `true` if the particle is unrimed, i.e. `F_rim = 0`.
"""
isunrimed(state::P3State) = iszero(state.F_rim)

"""
    get_ρ_d(mass::MassPowerLaw, F_rim, ρ_rim)
    get_ρ_d(state::P3State)

Exact solution for the density of the unrimed portion of the particle as
    function of the rime mass fraction `F_rim`, mass power law parameters `mass`,
    and rime density `ρ_rim`.

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
488.9120789986412
```
"""
function get_ρ_d((; β_va)::CMP.MassPowerLaw, F_rim, ρ_rim)
    k = (1 - F_rim)^(-1 / (3 - β_va))
    num = ρ_rim * F_rim
    den = (β_va - 2) * (k - 1) / ((1 - F_rim) * k - 1) - (1 - F_rim)
    return num / den
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

!!! note
    In the case where `ρ` is negative or zero, it is regularised to 1 kg/m³
    before the threshold is computed. This avoids taking a fractional power
    of a negative base, while remaining physically meaningful.
"""
_get_threshold((; α_va, β_va)::CMP.MassPowerLaw, ρ) =
    (6α_va / (π * max(ρ, 1)))^(1 / (3 - β_va))

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
    get_thresholds_ρ_g(state::P3State)
    get_thresholds_ρ_g(params::CMP.ParametersP3, F_rim, ρ_rim)

# Returns
- `(; D_th, D_gr, D_cr, ρ_g)`: The thresholds for the size distribution,
    and the density of total (deposition + rime) ice mass for graupel [kg/m³]
    If `F_rim = 0`, then we set `D_gr = D_cr = Inf`, and `ρ_g = NaN` (should not be used).

See [`get_D_th`](@ref), [`get_D_gr`](@ref), [`get_D_cr`](@ref), and [`get_ρ_g`](@ref) for more details.
"""
get_thresholds_ρ_g(state::P3State) = state.thresholds
function get_thresholds_ρ_g(params::CMP.ParametersP3, F_rim, ρ_rim)
    FT = eltype(F_rim)
    (; mass, ρ_i) = params
    D_th = get_D_th(mass, ρ_i)
    
    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    ρ_g = get_ρ_g(F_rim, ρ_rim, ρ_d)
    D_gr = get_D_gr(mass, ρ_g)
    D_cr = get_D_cr(mass, F_rim, ρ_g)

    # If unrimed, set D_gr and D_cr to infinity
    D_gr = ifelse(iszero(F_rim), FT(Inf), D_gr)
    D_cr = ifelse(iszero(F_rim), FT(Inf), D_cr)

    return (; D_th, D_gr, D_cr, ρ_g)
end

function get_bounded_thresholds(
    state::P3State{FT}, D_min::FT = FT(0), D_max::FT = FT(Inf)
) where {FT}
    (; D_th, D_gr, D_cr) = get_thresholds_ρ_g(state)
    return clamp.((D_min, D_th, D_gr, D_cr, D_max), D_min, D_max)
end

"""
    get_segments(state::P3State)

Return the segments of the size distribution as a tuple of intervals.

# Arguments
- `state`: [`P3State`](@ref) object

# Returns
- `segments`: tuple of tuples, each containing the lower and upper bounds of a segment

For example, if the thresholds are `(D_th, D_gr, D_cr)`, then the segments are:
- `(0, D_th)`, `(D_th, D_gr)`, `(D_gr, D_cr)`, `(D_cr, Inf)`
"""
function get_segments(state::P3State{FT}, D_min::FT = FT(0), D_max::FT = FT(Inf)) where {FT}
    (D_min, D_th, D_gr, D_cr, D_max) = get_bounded_thresholds(state, D_min, D_max)
    segments = ((D_min, D_th), (D_th, D_gr), (D_gr, D_cr), (D_cr, D_max))
    return segments
end

"""
    weighted_average(f_a, a, b)

Return the weighted average of `a` and `b` with fraction `f_a`,

```math
f_a ⋅ a + (1 - f_a) ⋅ b
```
"""
weighted_average(f_a, a, b) = f_a * a + (1 - f_a) * b

"""
    ice_mass_coeffs(state::P3State, D)
    ice_mass_coeffs(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Return the coefficients for the ice mass power law at diameter `D`.

# Arguments
 - `state`: The [`P3State`](@ref), or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum particle dimension [m]

# Returns
 - `(a, b)`: coefficients for the ice mass power law, `a D^b`
"""
function ice_mass_coeffs(state::P3State, D)
    (; params, F_rim, ρ_rim) = state
    FT = eltype(D)
    (; D_th, D_gr, D_cr, ρ_g) = get_thresholds_ρ_g(state)
    (; ρ_i) = params
    (; α_va, β_va) = params.mass

    return if D < D_th       # small spherical ice
        (ρ_i * π / 6, FT(3))
    elseif iszero(F_rim)     # large nonspherical unrimed ice
        (α_va, β_va)
    elseif D_th ≤ D < D_gr   # dense nonspherical rimed ice
        (α_va, β_va)
    elseif D_gr ≤ D < D_cr   # graupel (rimed)
        (ρ_g * π / 6, FT(3))
    else # D_cr ≤ D          # partially rimed ice
        (α_va / (1 - F_rim), β_va)
    end
end
ice_mass_coeffs(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    ice_mass_coeffs(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)

"""
    ice_mass(state, D)
    ice_mass(params, F_rim, ρ_rim, D)

Return the mass of a particle based on where it falls in the particle-size-based properties regime.

# Arguments
 - `state`: The [`P3State`](@ref), or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum particle dimension [m]
"""
function ice_mass(state::P3State, D)
    (a, b) = ice_mass_coeffs(state, D)
    return a * D^b
end
ice_mass(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    ice_mass(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)

"""
    ice_density(state::P3State, D)
    ice_density(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Return the density of a particle at diameter D

# Arguments
 - `state`: The [`P3State`](@ref), or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum particle dimension [m]

# Notes:
 The density of nonspherical particles is assumed to be the particle mass divided
 by the volume of a sphere with the same D [MorrisonMilbrandt2015](@cite).
 Needed for aspect ratio calculation, so we assume zero liquid fraction.
"""
ice_density(state::P3State, D) = ice_mass(state, D) / CO.volume_sphere_D(D)
ice_density(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    ice_density(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)

function get_∂mass_∂D_coeffs(state::P3State, D)
    (a, b) = ice_mass_coeffs(state, D)
    return a * b, b - 1
end
get_∂mass_∂D_coeffs(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    get_∂mass_∂D_coeffs(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)

"""
    ∂ice_mass_∂D(state::P3State, D)
    ∂ice_mass_∂D(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Return the derivative of the ice mass with respect to the particle diameter.
"""
function ∂ice_mass_∂D(state::P3State, D)
    (a, b) = get_∂mass_∂D_coeffs(state, D)
    return a * D^b
end
∂ice_mass_∂D(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    ∂ice_mass_∂D(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)

"""
    ice_area(state::P3State, D)
    ice_area(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Return the cross-sectional area of a particle based on where it falls in the
    particle-size-based properties regime.

# Arguments
 - `state`: [`P3State`](@ref) object, or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum particle dimension [m]
"""
function ice_area(state::P3State, D)
    (; params, F_rim, ρ_rim) = state
    (; D_th, D_gr, D_cr) = get_thresholds_ρ_g(state)
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
ice_area(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    ice_area(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)

"""
    ϕᵢ(state::P3State, D)
    ϕᵢ(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Returns the aspect ratio (ϕ) for an ice particle

# Arguments
 - `state`: The [`P3State`](@ref), or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum dimension of ice particle [m]

# Notes
 The density of nonspherical particles is assumed to be equal to the particle mass
 divided by the volume of a spherical particle with the same D_max [MorrisonMilbrandt2015](@cite).
 Assuming zero liquid fraction and oblate shape.
"""
function ϕᵢ(state::P3State, D)
    FT = eltype(D)
    # Reuse the mass-regime `(a, b)` coefficients for both mᵢ and ρᵢ so we
    # only compute `a * D^b` once per call.
    (a, b) = ice_mass_coeffs(state, D)
    mᵢ = a * D^b
    aᵢ = ice_area(state, D)
    ρᵢ = mᵢ / CO.volume_sphere_D(D)

    # TODO - prolate or oblate?
    ϕ_ob = min(1, 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^FT(1.5))) # κ =  1/3
    #ϕ_pr = max(1, 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2))       # κ = -1/6

    return ifelse(D == 0, FT(0), ϕ_ob)
end
ϕᵢ(params::CMP.ParametersP3, F_rim, ρ_rim, D) =
    ϕᵢ(P3State(params, zero(F_rim), zero(F_rim), F_rim, ρ_rim), D)
