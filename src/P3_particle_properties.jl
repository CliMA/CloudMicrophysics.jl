# TODO: Implement `F_liq` as another P3State-like struct
"""
    P3State{FT}

State of the P3 scheme.

This struct bundles the P3 parameterizations `params`, the provided rime state (`F_rim`, `ρ_rim`), 
    and the derived threshold variables (`D_th`, `D_gr`, `D_cr`, `ρ_g`).

To obtain a `P3State` object, use the [`get_state`](@ref) function.

# Fields
$(FIELDS)
"""
@kwdef struct P3State{FT, PARAMS <: CMP.ParametersP3{FT}}
    "[`CMP.ParametersP3`](@ref) object"
    params::PARAMS

    "Rime mass fraction"
    F_rim::FT
    "Rime density"
    ρ_rim::FT
end

"""
    get_state(params; F_rim, ρ_rim)

Create a [`P3State`](@ref) from [`CMP.ParametersP3`](@ref) and rime state parameters.

# Arguments
 - `params`: [`CMP.ParametersP3`](@ref) object

# Keyword arguments
 - `F_rim`: rime mass fraction [-], `F_rim = L_rim / L_ice`
 - `ρ_rim`: rime density [kg/m³],   `ρ_rim = L_rim / B_rim`

# Examples

 ```jldoctest
 julia> import CloudMicrophysics.Parameters as CMP, 
               CloudMicrophysics.P3Scheme as P3
 
 julia> FT = Float32;
 
 julia> params = CMP.ParametersP3(FT);
 
 julia> state = P3.get_state(params; F_rim = FT(0.5), ρ_rim = FT(916.7))
 P3State{Float32}
 ├── params = {MassPowerLaw, AreaPowerLaw, SlopePowerLaw, VentilationSB2005}
 ├── F_rim = 0.5 [-]
 └── ρ_rim = 916.7 [kg/m^3]
 ```
"""
function get_state(params::CMP.ParametersP3; F_rim, ρ_rim)
    # rime mass fraction must always be non-negative ...
    # ... and there must always be some unrimed part
    @assert 0 ≤ F_rim < 1 "Rime mass fraction, `F_rim`, must be between 0 and 1"
    # rime density must be positive ...
    # ... and as a bulk ice density can't exceed the density of water
    @assert 0 < ρ_rim ≤ params.ρ_l "Rime density, `ρ_rim`, must be between 0 and ρ_l"
    return P3State(; params, F_rim, ρ_rim)
end

Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

get_parameters(state::P3State) = state.params

isunrimed(state::P3State) = iszero(state.F_rim)

"""
    threshold_tuple(state)

Return a tuple of the thresholds for the current state.

This function is useful for providing thresholds to quadgk, for example.

!!! note
    If the state is unrimed, there is only one threshold, `D_th`.
    Otherwise (rimed state), there are three thresholds, `D_th < D_gr < D_cr`.
"""
function threshold_tuple(state::P3State)
    (; D_th, D_gr, D_cr) = get_thresholds_ρ_g(state)
    if isunrimed(state)
        return (D_th,)
    else
        return (D_th, D_gr, D_cr)
    end
end

"""
    get_ρ_d(mass::MassPowerLaw, F_rim, ρ_rim)

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
    get_thresholds_ρ_g(state::P3State)
    get_thresholds_ρ_g(params::CMP.ParametersP3, F_rim, ρ_rim)

# Returns
- `(; D_th, D_gr, D_cr, ρ_g)`: The thresholds for the size distribution, 
    and the density of total (deposition + rime) ice mass for graupel [kg/m³]


See [`get_D_th`](@ref), [`get_D_gr`](@ref), [`get_D_cr`](@ref), and [`get_ρ_g`](@ref) for more details.
"""
function get_thresholds_ρ_g(state::P3State)
    (; params, F_rim, ρ_rim) = state
    return get_thresholds_ρ_g(params, F_rim, ρ_rim)
end
function get_thresholds_ρ_g(params::CMP.ParametersP3, F_rim, ρ_rim)
    (; mass, ρ_i) = params
    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    ρ_g = get_ρ_g(F_rim, ρ_rim, ρ_d)

    D_th = get_D_th(mass, ρ_i)
    D_gr = get_D_gr(mass, ρ_g)
    D_cr = get_D_cr(mass, F_rim, ρ_g)

    return (; D_th, D_gr, D_cr, ρ_g)
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
    ice_mass_coeffs(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Return the coefficients for the ice mass power law at diameter `D`.

# Arguments
 - `state`: The [`P3State`](@ref), or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum particle dimension [m]

# Returns
 - `(a, b)`: coefficients for the ice mass power law, `a D^b`
"""
ice_mass_coeffs((; params, F_rim, ρ_rim)::P3State, D) = ice_mass_coeffs(params, F_rim, ρ_rim, D)
function ice_mass_coeffs(params::CMP.ParametersP3, F_rim, ρ_rim, D)
    FT = eltype(params)
    (; D_th, D_gr, D_cr, ρ_g) = get_thresholds_ρ_g(params, F_rim, ρ_rim)
    (; ρ_i) = params
    (; α_va, β_va) = params.mass

    return if D < D_th       # small spherical ice
        (ρ_i * π / 6, FT(3))
    elseif iszero(F_rim)  # large nonspherical unrimed ice
        (α_va, β_va)
    elseif D_th ≤ D < D_gr   # dense nonspherical rimed ice
        (α_va, β_va)
    elseif D_gr ≤ D < D_cr   # graupel (rimed)
        (ρ_g * π / 6, FT(3))
    else # D_cr ≤ D          # partially rimed ice
        (α_va / (1 - F_rim), β_va)
    end
end

"""
    ice_mass(state, D)
    ice_mass(params, F_rim, ρ_rim, D)

Return the mass of a particle based on where it falls in the particle-size-based properties regime.

# Arguments
 - `state`: The [`P3State`](@ref), or
 - `params, F_rim, ρ_rim`: The [`CMP.ParametersP3`](@ref), rime mass fraction, and rime density,
 - `D`: maximum particle dimension [m]
"""
function ice_mass(args_D...)
    D = last(args_D)
    (a, b) = ice_mass_coeffs(args_D...)
    return a * D^b
end

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
function ice_density(args_D...)
    D = last(args_D)
    return ice_mass(args_D...) / CO.volume_sphere_D(D)
end

function get_∂mass_∂D_coeffs(args_D...)
    (a, b) = ice_mass_coeffs(args_D...)
    return a * b, b - 1
end

"""
    ∂ice_mass_∂D(state::P3State, D)
    ∂ice_mass_∂D(params::CMP.ParametersP3, F_rim, ρ_rim, D)

Return the derivative of the ice mass with respect to the particle diameter.
"""
function ∂ice_mass_∂D(args_D...)
    D = last(args_D)
    (a, b) = get_∂mass_∂D_coeffs(args_D...)
    return a * D^b
end

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
ice_area((; params, F_rim, ρ_rim)::P3State, D) = ice_area(params, F_rim, ρ_rim, D)
function ice_area(params::CMP.ParametersP3, F_rim, ρ_rim, D)
    (; D_th, D_gr, D_cr) = get_thresholds_ρ_g(params, F_rim, ρ_rim)
    (; γ, σ) = params.area
    spherical_area(D) = D^2 * π / 4
    nonspherical_area(D) = γ * D^σ
    return if D < D_th       # small spherical ice
        spherical_area(D)
    elseif iszero(F_rim)  # large nonspherical unrimed ice
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
function ϕᵢ(args_D...)
    D = last(args_D)
    FT = eltype(D)
    mᵢ = ice_mass(args_D...)
    aᵢ = ice_area(args_D...)
    ρᵢ = ice_density(args_D...)

    # TODO - prolate or oblate?
    ϕ_ob = min(1, 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^FT(1.5))) # κ =  1/3
    #ϕ_pr = max(1, 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2))       # κ = -1/6

    return ifelse(D == 0, 0, ϕ_ob)
end

### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

function Base.show(io::IO, state::P3State)
    FT = eltype(state)
    _name(state, field) = typeof(getfield(state.params, field)).name.name
    param_types = join(_name.(state, (:mass, :area, :slope, :vent)), ", ")
    println(io, "P3State{$FT}")
    println(io, "├── params = {$param_types}")
    println(io, "├── F_rim = $(state.F_rim) [-]")
    println(io, "└── ρ_rim = $(state.ρ_rim) [kg/m³]")
end
