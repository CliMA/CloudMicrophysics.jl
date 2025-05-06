# TODO: Implement `F_liq` as another P3State-like struct
"""
    P3State{FT}

State of the P3 scheme.

This struct bundles the P3 parameterizations `params`, the provided rime state (`F_rim`, `ρ_r`), 
    and the derived threshold variables (`D_th`, `D_gr`, `D_cr`, `ρ_g`).

To obtain a `P3State` object, use the [`get_state`](@ref) function.

# Fields
$(FIELDS)
"""
@kwdef struct P3State{FT}
    "[`CMP.ParametersP3`](@ref) object"
    params::PSP3{FT}

    "Rime mass fraction"
    F_rim::FT
    "Rime density"
    ρ_r::FT
    "Graupel density"
    ρ_g::FT

    "Spherical ice threshold"
    D_th::FT
    "Graupel threshold"
    D_gr::FT
    "Partially rimed ice threshold"
    D_cr::FT
end

"""
    get_state(params; F_rim, ρ_r)

Create a [`P3State`](@ref) object from a [`CMP.ParametersP3`](@ref) object and rime state parameters.

# Arguments
- `params`: [`CMP.ParametersP3`](@ref) object
- `F_rim`: rime mass fraction
- `ρ_r`: rime density

# Examples

```jldoctest
julia> import CloudMicrophysics.Parameters as CMP, CloudMicrophysics.P3Scheme as P3

julia> FT = Float32;

julia> params = CMP.ParametersP3(FT);

julia> state = P3.get_state(params; F_rim = FT(0.5), ρ_r = FT(916.7))
P3State{Float32}
├── params = {MassPowerLaw, AreaPowerLaw, SlopePowerLaw, VentilationSB2005}
├── F_rim = 0.5 [-]
├── ρ_r = 916.7 [kg/m^3]
├── ρ_g = 702.8062 [kg/m^3]
├── D_th = 9.728088e-5 [m]
├── D_gr = 0.00012385941 [m]
└── D_cr = 0.00023259086 [m]
```
"""
function get_state(params::PSP3{FT}; F_rim, ρ_r) where {FT}
    # rime mass fraction must always be non-negative ...
    # ... and there must always be some unrimed part
    @assert 0 ≤ F_rim < 1 "Rime mass fraction, `F_rim`, must be between 0 and 1"

    (; mass, ρ_i, ρ_l) = params

    D_th = get_D_th(mass, ρ_i)

    if iszero(F_rim)
        ρ_g = D_gr = D_cr = FT(NaN)
    else
        # rime density must be positive ...
        # ... and as a bulk ice density can't exceed the density of water
        @assert 0 < ρ_r <= ρ_l

        ρ_d = get_ρ_d(mass, F_rim, ρ_r)
        ρ_g = get_ρ_g(ρ_r, F_rim, ρ_d)

        D_gr = get_D_gr(mass, ρ_g)
        D_cr = get_D_cr(mass, ρ_g, F_rim)
    end

    return P3State(; params, F_rim, ρ_r, ρ_g, D_th, D_gr, D_cr)
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
    if isunrimed(state)
        return (state.D_th,)
    else
        return (state.D_th, state.D_gr, state.D_cr)
    end
end

"""
    get_ρ_d(mass::MassPowerLaw, F_rim, ρ_r)

Exact solution for the density of the unrimed portion of the particle as 
    function of the rime mass fraction `F_rim`, mass power law parameters `mass`, 
    and rime density `ρ_r`.

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref) parameters
- `F_rim`: rime mass fraction
- `ρ_r`: rime density

# Returns
- `ρ_d`: density of the unrimed portion of the particle [kg/m³]

# Examples

```jldoctest
julia> import CloudMicrophysics.Parameters as CMP, 
              ClimaParams as CP, 
              CloudMicrophysics.P3Scheme as P3

julia> FT = Float64;

julia> mass = CMP.MassPowerLaw(CP.create_toml_dict(FT));

julia> F_rim, ρ_r = FT(0.5), FT(916.7);

julia> ρ_d = P3.get_ρ_d(mass, F_rim, ρ_r)
488.9120789986412
```
"""
function get_ρ_d((; β_va)::CMP.MassPowerLaw, F_rim, ρ_r)
    k = (1 - F_rim)^(-1 / (3 - β_va))
    num = ρ_r * F_rim
    den = (β_va - 2) * (k - 1) / ((1 - F_rim) * k - 1) - (1 - F_rim)
    return num / den
end

"""
    get_ρ_g(ρ_r, F_rim, ρ_d)

Return the density of total (deposition + rime) ice mass for graupel [kg/m³]

# Arguments
- `ρ_r`: rime density (`L_rim/B_rim`) [kg/m³]
- `F_rim`: rime mass fraction (`L_rim / L_ice`) [-]
- `ρ_d`: density of the unrimed portion of the particle [kg/m³], see [`get_ρ_d`](@ref)

# Returns
- `ρ_g`: density of total (deposition + rime) ice mass for graupel [kg/m³]

# Examples

```jldoctest
julia> import CloudMicrophysics.P3Scheme as P3

julia> FT = Float64; 

julia> ρ_r, F_rim, ρ_d = FT(916.7), FT(0.5), FT(916.7);

julia> ρ_g = P3.get_ρ_g(ρ_r, F_rim, ρ_d)
916.7
```

# Notes:
See Eq. 16 in [MorrisonMilbrandt2015](@cite).
"""
get_ρ_g(ρ_r, F_rim, ρ_d) = F_rim * ρ_r + (1 - F_rim) * ρ_d

"""
    get_threshold(params, ρ)

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
- `ρ`: density [kg/m³]
"""
get_threshold((; α_va, β_va)::CMP.MassPowerLaw{FT}, ρ::FT) where {FT} =
    (6α_va / (FT(π) * ρ))^(1 / (3 - β_va))

"""
    get_D_th(mass::MassPowerLaw, ρ_i)

Return the critical size separating spherical and nonspherical ice [meters]

See Eq. 8 in [MorrisonMilbrandt2015](@cite).
"""
get_D_th(mass::CMP.MassPowerLaw, ρ_i) = get_threshold(mass, ρ_i)

"""
    get_D_gr(mass::MassPowerLaw, ρ_g)

Return the size of equal mass for graupel and unrimed ice [meters]

See Eq. 15 in [MorrisonMilbrandt2015](@cite).
"""
get_D_gr(mass::CMP.MassPowerLaw, ρ_g) = get_threshold(mass, ρ_g)

"""
    get_D_cr(mass::MassPowerLaw, ρ_g, F_rim)

Return the size of equal mass for graupel and partially rimed ice [meters]

See Eq. 14 in [MorrisonMilbrandt2015](@cite).
"""
get_D_cr(mass::CMP.MassPowerLaw, ρ_g, F_rim) =
    get_threshold(mass, ρ_g * (1 - F_rim))

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
    mass_spherical(ρ, D)

Calculate the mass as a function of size for spherical particles, used for small ice, liquid, and graupel.

# Arguments
- `ρ`: bulk density [kg m⁻³]
- `D`: maximum particle dimension [m]

"""
mass_spherical(ρ, D) = ρ * CO.volume_sphere_D(D)

"""
    mass_nonspherical(mass, D)

Calculate the mass as a function of size for large nonspherical ice or dense nonspherical ice.

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref) parameters
- `D`: maximum particle dimension [m]
"""
mass_nonspherical((; α_va, β_va)::CMP.MassPowerLaw, D) = α_va * D^β_va

"""
    mass_rimed(mass, D, F_rim)

Calculate the mass as a function of size for partially rimed ice.

# Arguments
- `mass`: [`CMP.MassPowerLaw`](@ref) parameters
- `D`: maximum particle dimension [m]
- `F_rim`: rime mass fraction [`L_rim/L_ice`]
"""
mass_rimed((; α_va, β_va)::CMP.MassPowerLaw, D, F_rim) = α_va * D^β_va / (1 - F_rim)

"""
    ice_mass(state, D)

Return the mass of a particle based on where it falls in the particle-size-based properties regime.

# Arguments
- `state`: [`P3State`](@ref) object
- `D`: maximum particle dimension [m]
"""
function ice_mass(state::P3State, D)
    (; D_th, D_gr, D_cr, ρ_g, F_rim) = state
    (; mass, ρ_i) = get_parameters(state)
    return if D < D_th
        mass_spherical(ρ_i, D)      # small spherical ice
    elseif isunrimed(state)
        mass_nonspherical(mass, D)  # large nonspherical unrimed ice
    elseif D_th ≤ D < D_gr
        mass_nonspherical(mass, D)  # dense nonspherical rimed ice
    elseif D_gr ≤ D < D_cr
        mass_spherical(ρ_g, D)      # graupel (rimed)
    else # D_cr ≤ D
        mass_rimed(mass, D, F_rim)  # partially rimed ice
    end
end

"""
    ice_density(state::P3State, D)

Return the density of a particle based on where it falls in the particle-size-based properties regime. 

# Arguments
- `state`: [`P3State`](@ref) object
- `D`: maximum particle dimension [m]

The density of nonspherical particles is assumed to be the particle mass divided by the volume of a sphere with the 
    same D [MorrisonMilbrandt2015](@cite). Needed for aspect ratio calculation, so we assume zero liquid fraction.
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
    (; D_th, D_gr, D_cr, ρ_g, F_rim) = state
    (; mass, ρ_i) = get_parameters(state)

    ∂mass_spherical_∂D(ρ::FT, D::FT) where {FT} = ρ * D^2 * FT(π) / 2
    ∂mass_nonspherical_∂D((; α_va, β_va)::CMP.MassPowerLaw, D) =
        β_va * α_va * D^(β_va - 1)
    ∂mass_rimed_∂D((; α_va, β_va)::CMP.MassPowerLaw, D, F_rim) =
        β_va * α_va * D^(β_va - 1) / (1 - F_rim)

    return if D < D_th
        ∂mass_spherical_∂D(ρ_i, D)      # small spherical ice
    elseif isunrimed(state)
        ∂mass_nonspherical_∂D(mass, D)  # large nonspherical unrimed ice
    elseif D_th ≤ D < D_gr
        ∂mass_nonspherical_∂D(mass, D)  # dense nonspherical rimed ice
    elseif D_gr ≤ D < D_cr
        ∂mass_spherical_∂D(ρ_g, D)      # graupel (rimed)
    else # D_cr ≤ D
        ∂mass_rimed_∂D(mass, D, F_rim)  # partially rimed ice
    end
end


"""
    area_spherical(D)

Calculate the cross-sectional area of a spherical particle.

# Arguments
- `D`: maximum particle dimension [m]
"""
area_spherical(D::FT) where {FT} = D^2 * FT(π) / 4

"""
    area_nonspherical(area, D)

Calculate the cross-sectional area of a nonspherical particle.

# Arguments
- `area`: [`CMP.AreaPowerLaw`](@ref) parameters
- `D`: maximum particle dimension [m]
"""
area_nonspherical((; γ, σ)::CMP.AreaPowerLaw, D) = γ * D^σ

"""
    area_rimed(area, F_rim, D)

Calculate the cross-sectional area of a partially rimed particle.

# Arguments
- `area`: [`CMP.AreaPowerLaw`](@ref) parameters
- `F_rim`: rime mass fraction [`L_rim/L_ice`]
- `D`: maximum particle dimension [m]
"""
function area_rimed(area::CMP.AreaPowerLaw, F_rim, D)
    rimed_area = area_spherical(D)
    unrimed_area = area_nonspherical(area, D)
    return weighted_average(F_rim, rimed_area, unrimed_area)
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
    (; area) = get_parameters(state)
    (; D_th, D_gr, D_cr, F_rim) = state
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

"""
    ventilation_factor(state, vel, ρₐ, aps)

Returns a function that computes the ventilation factor for an ice particle 
    as a function of the maximum particle dimension, `D`.

# Arguments
- `state`: The [`P3State`](@ref)
- `vel`: The [`Chen2022VelType`](@ref)
- `ρₐ`: Air density [kg/m³]
- `aps`: [`CMP.AirProperties`](@ref)

# Returns
- `vent_factor(D)`: The ventilation factor as a function of `D`
"""
function ventilation_factor(state, vel, ρₐ, aps)
    (; ν_air, D_vapor) = aps
    (; vent_a, vent_b) = state.params.vent
    N_sc = ν_air / D_vapor
    v_term = ice_particle_terminal_velocity(state, vel, ρₐ)

    # Reynolds number
    N_Re(D) = D * v_term(D) / ν_air
    # Ventilation factor
    vent_factor(D) = vent_a + vent_b * cbrt(N_sc) * sqrt(N_Re(D))
    return vent_factor
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
    println(io, "├── ρ_r = $(state.ρ_r) [kg/m^3]")
    println(io, "├── ρ_g = $(state.ρ_g) [kg/m^3]")
    println(io, "├── D_th = $(state.D_th) [m]")
    println(io, "├── D_gr = $(state.D_gr) [m]")
    println(io, "└── D_cr = $(state.D_cr) [m]")
end
