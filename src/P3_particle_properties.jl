# TODO: Implement `F_liq` as another P3State-like struct
"""
    P3State{FT}

State of the P3 scheme.

This struct bundles the P3 parameterizations `params` and 
    the frozen state parameters (`L_ice`, `N_ice`, `F_rim`, `ρ_rim`).

To obtain a `P3State` object, use the [`get_state`](@ref) function.

# Fields
$(FIELDS)
"""
@kwdef struct P3State{FT, PARAMS <: CMP.ParametersP3{FT}}
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
end

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
 ├── ρ_rim = 916.7 [kg/m^3]
 ```
"""
function get_state(params::CMP.ParametersP3; L_ice, N_ice, F_rim, ρ_rim)
    # rime mass fraction must always be non-negative ...
    # ... and there must always be some unrimed part
    @assert 0 ≤ F_rim < 1 "Rime mass fraction, `F_rim`, must be between 0 and 1"
    # rime density must be positive ...
    # ... and as a bulk ice density can't exceed the density of water
    @assert 0 ≤ ρ_rim ≤ params.ρ_l
    return P3State(; params, L_ice, N_ice, F_rim, ρ_rim)
end


Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

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
get_ρ_g(F_rim, ρ_rim, ρ_d) = CO.weighted_average(F_rim, ρ_rim, ρ_d)
function get_ρ_g(mass::CMP.MassPowerLaw, F_rim, ρ_rim)
    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    return get_ρ_g(F_rim, ρ_rim, ρ_d)
end

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

function get_thresholds_ρ_g(state::P3State)
    (; params, F_rim, ρ_rim) = state
    (; mass, ρ_i) = params

    ρ_d = get_ρ_d(mass, F_rim, ρ_rim)
    ρ_g = get_ρ_g(F_rim, ρ_rim, ρ_d)

    D_th = get_D_th(mass, ρ_i)
    D_gr = get_D_gr(mass, ρ_g)
    D_cr = get_D_cr(mass, F_rim, ρ_g)

    return (; D_th, D_gr, D_cr, ρ_g)
end

"""
    get_mass_law(state::P3State)

Construct the piecewise power law for ice mass as a function of diameter.

Returns a [`CO.PiecewisePowerLaw`](@ref) representing the `mass(D)` relationship.

# Arguments
 - `state`: [`P3State`](@ref) object

# Note
 Currently implemented for [`CMP.MassPowerLaw`](@ref). If other mass parameterizations 
 are added in the future, this function would need to dispatch on the mass parameterization type.
"""
@inline function get_mass_law(state::P3State{FT}) where {FT}
    (; params, F_rim) = state
    (; mass, ρ_i) = params
    (; α_va, β_va) = mass
    (; D_th, D_gr, D_cr, ρ_g) = get_thresholds_ρ_g(state)
    
    if iszero(F_rim)
        # Unrimed: 2 pieces
        # spherical ice, nonspherical ice
        piece1 = CO.BoundedPowerLaw{FT}(; D_min = FT(0), D_max = D_th, a = ρ_i * π / 6, b = FT(3))
        piece2 = CO.BoundedPowerLaw{FT}(; D_min = D_th, D_max = FT(Inf), a = α_va, b = β_va)
        
        return CO.PiecewisePowerLaw{FT, 2}(; pieces = (piece1, piece2))
    else
        # Rimed: 4 pieces
        # spherical ice, nonspherical ice, graupel, partially rimed ice
        piece1 = CO.BoundedPowerLaw{FT}(; D_min = FT(0), D_max = D_th, a = ρ_i * π / 6, b = FT(3))
        piece2 = CO.BoundedPowerLaw{FT}(; D_min = D_th, D_max = D_gr, a = α_va, b = β_va)
        piece3 = CO.BoundedPowerLaw{FT}(; D_min = D_gr, D_max = D_cr, a = ρ_g * π / 6, b = FT(3))
        piece4 = CO.BoundedPowerLaw{FT}(; D_min = D_cr, D_max = FT(Inf), a = α_va / (1 - F_rim), b = β_va)
        
        return CO.PiecewisePowerLaw{FT, 4}(; pieces = (piece1, piece2, piece3, piece4))
    end
end
"""
    ice_mass(state, D)

Return the mass of a particle at diameter D

# Arguments
 - `state`: [`P3State`](@ref) object
 - `D`: maximum particle dimension [m]
"""
# function ice_mass(state::P3State, D)
#     mass_law = get_mass_law(state)
#     return mass_law(D)
# end

"""
    ice_density(state::P3State, D)

Return the density of a particle at diameter D

# Arguments
 - `state`: [`P3State`](@ref) object
 - `D`: maximum particle dimension [m]

# Notes:
 The density of nonspherical particles is assumed to be the particle mass divided 
 by the volume of a sphere with the same D [MorrisonMilbrandt2015](@cite). 
 Needed for aspect ratio calculation, so we assume zero liquid fraction.
"""
ice_density(state::P3State, D) = ice_mass(state, D) / CO.volume_sphere_D(D)

# """
#     p3_mass(state, D)

# Return mass(D) regime, used to create figures for the docs page.

# # Arguments
# - `state`: [`P3State`](@ref) object
# - `D`: maximum particle dimension [m]

# """
# function p3_mass(state::P3State, D)
#     # TODO: Refactor to be a parameterization for use with `F_liq`
#     (; ρ_l) = get_parameters(state)
#     m_liq = mass_spherical(ρ_l, D)
#     m_ice = ice_mass(state, D)
#     mass_avg = weighted_average(F_liq, m_liq, m_ice)
#     return mass_avg
# end


"""
    ∂ice_mass_∂D(state, D)

Calculate the derivative of the ice mass with respect to the maximum particle dimension.

Used in snow melting calculations.

# Arguments
- `state`: [`P3State`](@ref) object
- `D`: maximum dimension of the particle [m]
"""
function ∂ice_mass_∂D(state::P3State, D)
    mass_law = get_mass_law(state)
    piece = CO.get_power_law_at(mass_law, D)
    (a, b) = CO.get_coefficients(piece)
    return b * a * D^(b - 1)
end

"""
    ice_area(state, D)

Return the cross-sectional area of a particle based on where it falls in the 
    particle-size-based properties regime.

# Arguments
- `state`: [`P3State`](@ref) object
- `D`: maximum particle dimension [m]
"""
function ice_area(state::P3State, D)
    area_spherical(D) = D^2 * π / 4
    area_nonspherical((; γ, σ)::CMP.AreaPowerLaw, D) = γ * D^σ
    function area_rimed(area::CMP.AreaPowerLaw, F_rim, D)
        rimed_area = area_spherical(D)
        unrimed_area = area_nonspherical(area, D)
        return CO.weighted_average(F_rim, rimed_area, unrimed_area)
    end

    (; params, F_rim) = state
    (; area) = params
    (; D_th, D_gr, D_cr) = get_thresholds_ρ_g(state)
    
    return if D < D_th
        area_spherical(D)           # small spherical ice
    elseif isunrimed(state)
        area_nonspherical(area, D)  # large nonspherical unrimed ice
    elseif D_th ≤ D < D_gr
        area_nonspherical(area, D)  # dense nonspherical rimed ice
    elseif D_gr ≤ D < D_cr
        area_spherical(D)           # graupel (rimed)
    else # D_cr ≤ D
        area_rimed(area, F_rim, D)  # partially rimed ice
    end
end

# """
#     p3_area(state, D)

# Return area(D), used to create figures for the documentation.

# # Arguments
# - `state`: [`P3State`](@ref) object
# - `D`: maximum particle dimension [m]
# """
# function p3_area(state::P3State, D)
#     # TODO: Refactor to be a parameterization for use with `F_liq`
#     (; F_liq) = get_parameters(state)
#     area_liq = area_spherical(D)
#     area_ice = ice_area(state, D)
#     return weighted_average(F_liq, area_liq, area_ice)
# end

"""
    ϕᵢ(state, D)

Returns the aspect ratio (ϕ) for an ice particle

# Arguments
- `state`: [`P3State`](@ref) object
- `D`: maximum dimension of ice particle [m]

The density of nonspherical particles is assumed to be equal to the particle mass divided by the volume of a
spherical particle with the same D_max [MorrisonMilbrandt2015](@cite). Assuming zero liquid fraction and oblate shape.
"""
function ϕᵢ(state::P3State, D)
    FT = eltype(state)
    mᵢ = ice_mass(state, D)
    aᵢ = ice_area(state, D)
    ρᵢ = ice_density(state, D)

    # TODO - prolate or oblate?
    ϕ_ob = min(1, 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^FT(1.5))) # κ =  1/3
    #ϕ_pr = max(1, 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2))       # κ = -1/6

    return ifelse(D == 0, 0, ϕ_ob)
end

### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

function Base.show(io::IO, state::P3State{FT}) where {FT}
    _name(state, field) = typeof(getfield(state.params, field)).name.name
    param_types = join(_name.(state, (:mass, :area, :slope, :vent)), ", ")
    println(io, "P3State{$FT}")
    println(io, "├── params = {$param_types}")
    println(io, "├── F_rim = $(state.F_rim) [-]")
    println(io, "└── ρ_rim = $(state.ρ_rim) [kg/m³]")
end
