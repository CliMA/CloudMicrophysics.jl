
# TODO: Consider order of properties, and also whether L, N should be in here
Base.@kwdef struct P3State{FT}
    # Core parameters
    params::PSP3{FT}

    # Rime state
    F_rim::FT
    F_liq::FT  # TODO: Make this a struct with liquid properties
    ρ_r::FT
    ρ_g::FT

    # Regime boundary thresholds
    D_th::FT
    D_gr::FT
    D_cr::FT
end

function get_state(params::PSP3{FT}; F_rim, F_liq, ρ_r) where {FT}
    # rime mass fraction must always be non-negative ...
    # ... and there must always be some unrimed part
    @assert 0 ≤ F_rim < 1

    (; mass, ρ_i, ρ_l) = params

    D_th = get_D_th(mass, ρ_i)

    if iszero(F_rim)
        ρ_g = D_gr = D_cr = FT(NaN)
    else
        # rime density must be positive ...
        # ... and as a bulk ice density can't exceed the density of water
        @assert 0 < ρ_r <= ρ_l

        ρ_d = get_ρ_d_exact(mass, F_rim, ρ_r)
        ρ_g = get_ρ_g(ρ_r, F_rim, ρ_d)

        D_gr = get_D_gr(mass, ρ_g)
        D_cr = get_D_cr(mass, ρ_g, F_rim)
    end

    return P3State(; params, F_rim, F_liq, ρ_r, ρ_g, D_th, D_gr, D_cr)
end

Base.eltype(::P3State{FT}) where {FT} = FT
Base.broadcastable(state::P3State) = tuple(state)

get_parameters(state::P3State) = state.params

isrimed(state::P3State) = state.F_rim > 0
isunrimed(state::P3State) = !isrimed(state)

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
    get_ρ_d_exact(F_rim, β_va, ρ_r)

Exact solution for the density of the unrimed portion of the particle

as function of the rime mass fraction `F_rim`, mass power law exponent `β_va`, and rime density `ρ_r`.
"""
function get_ρ_d_exact((; β_va)::MassPowerLaw, F_rim, ρ_r)
    k = (1 - F_rim) ^ (-1/(3 - β_va))
	num = ρ_r * F_rim
	den = (β_va - 2) * (k - 1) / ( (1 - F_rim) * k - 1 ) - (1 - F_rim)
	return num / den
end

"""
    get_ρ_g(ρ_r, F_rim, ρ_d)

- ρ_r   - rime density (L_rim/B_rim) [kg/m^3]
- F_rim - rime mass fraction (L_rim / L_ice) [-]
- ρ_d   - is the density of the unrimed portion of the particle [kg/m^3]

Return the density of total (deposition + rime) ice mass for graupel [kg/m3]

See Eq. 16 in Morrison and Milbrandt (2015).
"""
get_ρ_g(ρ_r, F_rim, ρ_d) = F_rim * ρ_r + (1 - F_rim) * ρ_d

"""
    get_threshold(params, ρ)

All thresholds are on the form
    ``(fac * 6α_va / (π ρ))^(1 / (3 - β_va))``

where for
    - D_th: ρ = ρ_i
    - D_gr: ρ = ρ_g
    - D_cr: ρ = ρ_g ⋅ (1 - F_rim)
"""
get_threshold((; α_va, β_va)::MassPowerLaw, ρ) = (6α_va / (π * ρ))^(1 / (3 - β_va))

"""
    get_D_th(mass::MassPowerLaw, ρ_i)

Return the critical size separating spherical and nonspherical ice [meters]

# Notes:
See Eq. 8 in Morrison and Milbrandt (2015).
"""
get_D_th(mass::MassPowerLaw, ρ_i) = get_threshold(mass, ρ_i)

"""
    get_D_gr(mass::MassPowerLaw, ρ_g)

Return the size of equal mass for graupel and unrimed ice [meters]

# Notes:
See Eq. 15 in Morrison and Milbrandt (2015).
"""
get_D_gr(mass::MassPowerLaw, ρ_g) = get_threshold(mass, ρ_g)

"""
    get_D_cr(mass::MassPowerLaw, ρ_g, F_rim)

Return the size of equal mass for graupel and partially rimed ice [meters]

# Notes:
See Eq. 14 in Morrison and Milbrandt (2015).
"""
get_D_cr(mass::MassPowerLaw, ρ_g, F_rim) = get_threshold(mass, ρ_g * (1 - F_rim))

struct P3Distribution{FT}
    # Particle state
    state::P3State{FT}

    # Distribution parameters
    log_λ::FT
    log_N₀::FT
end

get_state(dist::P3Distribution) = dist.state

Base.eltype(::P3Distribution{FT}) where {FT} = FT
Base.broadcastable(state::P3Distribution) = tuple(state)

get_parameters(dist::P3Distribution) = get_parameters(get_state(dist))

isrimed(dist::P3Distribution) = isrimed(get_state(dist))
isunrimed(dist::P3Distribution) = isunrimed(get_state(dist))

"""
    p3_F_liq_average(F_liq, X_ice, X_liq)

 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - X_ice - ice core parameterization (i.e. mass, etc)
 - X_liq - liquid part parameterization

Returns the liquid fraction weighted average of X_ice and X_liq.
"""
# function p3_F_liq_average(F_liq, X_ice, X_liq)
#     return (1 - F_liq) * X_ice + F_liq * X_liq
# end

"""
    weighted_average(f_a, a, b)

Return the weighted average of `a` and `b` with fraction `f_a`,
    
    ``f_a * a + (1 - f_a) * b``
"""
function weighted_average(f_a, a, b)
    return f_a * a + (1 - f_a) * b
end

"""
    volume_sphere(D)

Calculate the volume of a sphere with diameter D.
"""
volume_sphere(D) = D^3 * π / 6

"""
    mass_spherical(ρ, D)

Calculate the mass as a function of size for spherical particles, used for small ice, liquid, and graupel.

# Arguments
- `ρ`: bulk density [kg m⁻³]
- `D`: maximum particle dimension [m]

"""
mass_spherical(ρ, D) = ρ * volume_sphere(D)  # TODO: used to be: `mass_s` (update tests)

"""
    mass_nonspherical(mass, D)

Calculate the mass as a function of size for large nonspherical ice or dense nonspherical ice.

# Arguments
- `mass::MassPowerLaw`: mass power law parameters
- `D`: maximum particle dimension [m]
"""
mass_nonspherical((; α_va, β_va)::MassPowerLaw, D) = α_va * D^β_va  # TODO: used to be: `mass_nl`

"""
    mass_rimed(mass, D, F_rim)

Calculate the mass as a function of size for partially rimed ice.

# Arguments
- `mass::MassPowerLaw`: mass power law parameters
- `D`: maximum particle dimension [m]
- `F_rim`: rime mass fraction [L_rim/L_ice]
"""
mass_rimed((; α_va, β_va)::MassPowerLaw, D, F_rim) = α_va * D^β_va / (1 - F_rim)  # TODO: used to be: `mass_r`

"""
    ice_mass(state::P3State, D)

Return the mass of a particle based on where it falls in the particle-size-based properties regime.

# Arguments
- `state`: The P3 state
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

- `state`: The P3 state
- `D`: maximum particle dimension [m]

Return the density of a particle based on where it falls in the particle-size-based properties regime. 

Following Morrison and Milbrandt (2015), the density of nonspherical particles is assumed to be the particle mass divided 
by the volume of a sphere with the same D. Needed for aspect ratio calculation, so we assume zero liquid fraction.
"""
ice_density(state::P3State, D) = ice_mass(state, D) / volume_sphere(D)

"""
    p3_mass(p3, D, F_rim, F_liq, th)

Returns mass(D) regime, used to create figures for the docs page.

# Arguments
- `p3`: a struct with P3 scheme parameters
- `D`: maximum particle dimension [m]
- `F_rim`: rime mass fraction (`L_rim / L_ice`)
- `F_liq`: liquid fraction (`L_liq / L_ice`)
- `th`: P3 scheme thresholds() output tuple (D_cr, D_gr, ρ_g, ρ_d)

"""
function p3_mass(state::P3State, D)
    # TODO: Refactor to be a parameterization for use with `F_liq`
    (; ρ_l) = params = get_parameters(state)
    m_liq = mass_spherical(ρ_l, D)
    m_ice = ice_mass(state, D)
    mass_avg = p3_F_liq_average(F_liq, m_ice, m_liq)
    return mass_avg
end


"""
    ∂ice_mass_∂D(state, D)

Calculate the derivative of the ice mass with respect to the maximum particle dimension.

Used in snow melting calculations.

# Arguments
- `state::P3State`: P3 state
- `D`: maximum dimension of the particle [m]
"""
function ∂ice_mass_∂D(state::P3State, D)
    (; D_th, D_gr, D_cr, ρ_g, F_rim) = state
    (; mass, ρ_i) = get_parameters(state)

    ∂mass_spherical_∂D(ρ, D) = ρ * D^2 * π / 2
    ∂mass_nonspherical_∂D((; α_va, β_va)::MassPowerLaw, D) = β_va * α_va * D^(β_va - 1)
    ∂mass_rimed_∂D((; α_va, β_va)::MassPowerLaw, D, F_rim) = β_va * α_va * D^(β_va - 1) / (1 - F_rim)
    
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


# for spherical particles
area_spherical(D) =  D^2 * π / 4
# for nonspherical particles
area_nonspherical((; γ, σ)::AreaPowerLaw, D) = γ * D^σ
# partially rimed ice
area_rimed(area::AreaPowerLaw, F_rim, D) = weighted_average(F_rim, area_spherical(D), area_nonspherical(area, D))

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

"""
    p3_area(p3, D, F_rim, F_liq, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_rim - rime mass fraction (L_rim / L_ice)
 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - th - P3 scheme thresholds() output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns area(D), used to create figures for the documentation.
"""
function p3_area(state::P3State, D)
    # TODO: Refactor to be a parameterization for use with `F_liq`
    (; F_liq) = get_parameters(state)
    area_liq = area_spherical(D)
    area_ice = ice_area(state, D)
    return weighted_average(F_liq, area_liq, area_ice)
end

"""
    ϕᵢ(p3, D, F_rim, th)

 - p3 - a struct containing P3 parameters
 - D - maximum dimension of ice particle [m]
 - F_rim - rime mass fraction (L_rim/ L_ice) [-]
 - th - P3 particle properties thresholds

Returns the aspect ratio (ϕ) for an ice particle with mass, cross-sectional area,
and ice density determined using the size-dependent particle property regimes
following Morrison and Milbrandt (2015). The density of nonspherical
particles is assumed to be equal to the particle mass divided by the volume of a
spherical particle with the same D_max.
Assuming zero liquid fraction and oblate shape.
"""
# function ϕᵢ(p3::PSP3{FT}, D, F_rim, th) where {FT}
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
    println(io, "P3State{$FT}")
    println(io, "├── F_rim = $(state.F_rim) [-]")
    println(io, "├── F_liq = $(state.F_liq) [-]")
    println(io, "├── ρ_r = $(state.ρ_r) [kg/m^3]")
    println(io, "├── ρ_g = $(state.ρ_g) [kg/m^3]")
    println(io, "├── D_th = $(state.D_th) [m]")
    println(io, "├── D_gr = $(state.D_gr) [m]")
    println(io, "└── D_cr = $(state.D_cr) [m]")
end

function Base.show(io::IO, p3s::P3Distribution{FT}) where {FT}
    println(io, "P3Distribution{$FT}")
    println(io, "├── log_λ = $(p3s.log_λ) [log(1/m)]")
    println(io, "└── log_N₀ = $(p3s.log_N₀) [log(1/m^3)]")
end
