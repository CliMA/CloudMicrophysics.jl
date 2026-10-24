"""
    collision_cross_section_ice_liquid_coeffs(rᵢ)
    collision_cross_section_ice_liquid_coeffs(particles, Dᵢ)

Return the monomial coefficients `(k₀, k₁, k₂)` of the ice-liquid collision cross-section.

The cross-section is a polynomial in the liquid diameter `Dₗ`,

```math
σ(Dᵢ, Dₗ) = π (rᵢ + Dₗ/2)² = k₀ + k₁ Dₗ + k₂ Dₗ²,
```

with `k₀ = π rᵢ²`, `k₁ = π rᵢ`, `k₂ = π/4`, where the ice effective radius
is `rᵢ = √(ice_area(particles, Dᵢ)/π)`; see [`ice_area`](@ref).

Used in [`collision_cross_section_ice_liquid`](@ref)
"""
@inline collision_cross_section_ice_liquid_coeffs(rᵢ::FT) where {FT} =
    (π * rᵢ^2, π * rᵢ, FT(π / 4))
@inline collision_cross_section_ice_liquid_coeffs(particles, Dᵢ) =
    collision_cross_section_ice_liquid_coeffs(√(ice_area(particles, Dᵢ) / π))

"""
    collision_cross_section_ice_liquid(particles, Dᵢ, Dₗ)

Return the ice-liquid collision cross-section [m²], `π (rᵢ(Dᵢ) + Dₗ/2)²`.

It is evaluated by Horner from the shared [`collision_cross_section_ice_liquid_coeffs`](@ref).
"""
collision_cross_section_ice_liquid(particles, Dᵢ, Dₗ) =
    evalpoly(Dₗ, collision_cross_section_ice_liquid_coeffs(particles, Dᵢ))

"""
    volumetric_collision_rate_integrand(velocity_params, particles, ρₐ)

Return a function that computes the volumetric collision rate integrand for ice-liquid collisions [m³/s].

The returned function takes ice and liquid particle diameters as arguments.

# Arguments
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `particles`: [`P3Particles`](@ref)
- `ρₐ`: air density [kg/m³]

# Returns
A [`VolumetricCollisionRate`](@ref), a function `(D_ice, D_liq) -> E * K * |vᵢ - vₗ|` where:
- `D_ice` and `D_liq` are the (maximum) diameters of the ice and liquid particles
- `E` is the collision efficiency
- `K` is the collision cross section
- `vᵢ` and `vₗ` are the terminal velocities of ice and liquid particles

Note that `E`, `K`, `vᵢ` and `vₗ` are all, in general, functions of `D_ice` and `D_liq`.

This function is a component of integrals like

```math
∫ ∫ E * K * |vᵢ - vₗ| * N'_i * N'_l dD_i dD_l
```
"""
function volumetric_collision_rate_integrand(velocity_params, particles, ρₐ)
    v_ice = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    return VolumetricCollisionRate(particles, v_ice, v_liq)
end

"""
    VolumetricCollisionRate(particles, v_i, v_l)

The volumetric collision rate `(D_ice, D_liq) -> E * K * |vᵢ - vₗ|` [m³/s] of ice-liquid collisions.

See [`volumetric_collision_rate_integrand`](@ref).

# Fields
- `particles`: [`P3Particles`](@ref)
- `v_i`, `v_l`: ice and liquid terminal velocity functions [m/s]
"""
struct VolumetricCollisionRate{P, VI, VL} <: Function
    particles::P
    v_i::VI
    v_l::VL
end
function (∂ₜV::VolumetricCollisionRate)(D_ice::FT, D_liq::FT) where {FT}
    E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
    K = collision_cross_section_ice_liquid(∂ₜV.particles, D_ice, D_liq)
    return E * K * abs(∂ₜV.v_i(D_ice) - ∂ₜV.v_l(D_liq))
end

"""
    compute_max_freeze_rate(aps, tps, velocity_params, particles, ρₐ, Tₐ)

Return a function `max_freeze_rate(Dᵢ)` of the maximum possible freezing rate [kg/s].

The rate is that of an ice particle of diameter `Dᵢ` [m]. It evaluates to `0` if `T ≥ T_freeze`.

# Arguments
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `particles`: [`P3Particles`](@ref)
- `ρₐ`: air density [kg/m³]
- `Tₐ`: air temperature [K]

This rate represents the thermodynamic upper limit to collisional freezing,
which occurs when the heat transfer from the ice particle to the environment is
balanced by the latent heat of fusion.

From Eq (A7) in Musil (1970), [Musil1970](@cite).
"""
function compute_max_freeze_rate(aps, tps, velocity_params, particles, ρₐ, Tₐ)
    (; heat, denom, T_frz) = max_freeze_rate_scalars(aps, tps, ρₐ, Tₐ)
    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    F_v = CO.ventilation_factor(particles.params.vent, aps, v_term)
    function max_freeze_rate(Dᵢ)
        # fallback values typed by the promotion of the node and the captured particles
        # (mixed plain/Dual under differentiation)
        FT = UT.promote_typeof(Dᵢ, heat, denom)
        Tₐ ≥ T_frz && return zero(FT)  # No collisional freezing above the freezing temperature
        denom > 0 || return floatmax(FT)
        return 2 * (π * Dᵢ) * F_v(Dᵢ) * heat / denom
    end
    return max_freeze_rate
end

"""
    max_freeze_rate_scalars(aps, tps, ρₐ, Tₐ)

Return the factors of the maximum freezing rate of [`compute_max_freeze_rate`](@ref) that do not depend on the particle.

The maximum freezing rate of an ice particle of diameter `Dᵢ` is `2π Dᵢ F_v(Dᵢ) heat / denom`
below the freezing temperature `T_frz` if `denom > 0`, and unbounded if `denom ≤ 0`.

# Arguments
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `ρₐ`: air density [kg/m³]
- `Tₐ`: air temperature [K]

# Returns
A `NamedTuple` `(; heat, denom, T_frz)`, where:
- `heat`: heat flux from the particle surface per unit capacity, by conduction and by vapor
  diffusion [W/m]
- `denom`: heat to remove per unit mass of collected liquid, which warms from `Tₐ` to the
  surface temperature and freezes there [J/kg]; it is nonpositive for `T_frz - Tₐ ≥ L_f / cp_l
  ≈ 80 K`, where every collected droplet freezes
- `T_frz`: freezing temperature [K]
"""
@inline function max_freeze_rate_scalars(aps, tps, ρₐ, Tₐ)
    (; D_vapor, K_therm) = aps
    cp_l = TDI.cp_l(tps)
    T_frz = TDI.T_freeze(tps)
    Tₛ = T_frz  # the surface of the ice particle is assumed to be at the freezing temperature
    Lᵥ = TDI.Lᵥ(tps, Tₛ)
    L_f = TDI.Lf(tps, Tₛ)
    ΔT = Tₛ - Tₐ  # temperature difference between the surface of the ice particle and the air
    Δρᵥ_sat =
        ρₐ * (  # saturation vapor density difference between the surface of the ice particle and the air
            TDI.p2q(tps, Tₛ, ρₐ, TDI.saturation_vapor_pressure_over_ice(tps, Tₛ)) -
            TDI.p2q(tps, Tₐ, ρₐ, TDI.saturation_vapor_pressure_over_ice(tps, Tₐ))
        )
    heat = K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat
    denom = L_f - cp_l * ΔT
    return (; heat, denom, T_frz)
end

"""
    compute_local_rime_density(velocity_params, particles, ρₐ, T)

Return a function `ρ′_rim(Dᵢ, Dₗ)` that computes the local rime density [kg/m³].

The arguments are the ice particle diameter `Dᵢ` [m] and the liquid particle diameter `Dₗ` [m].

# Arguments
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `particles`: [`P3Particles`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Returns
A function that computes the local rime density [kg/m³] using the equation:

```math
ρ'_{rim} = a + b R_i + c R_i^2
```
where
```math
R_i = -\\frac{ 10^6 ⋅ D_{liq} ⋅ |v_{liq} - v_{ice}| }{ 2 T_{sfc} }
```
and ``T_{sfc} < 0`` is the sub-zero surface temperature [°C], ``D_{liq}`` is the liquid
particle diameter [m], ``v_{liq/ice}`` is the particle terminal velocity [m/s].
With the ``10^6`` factor converting ``D_{liq}`` from [m] to [μm], the units of
``R_i`` are [μm m s⁻¹ °C⁻¹]. The units of ``ρ'_{rim}`` are [kg/m³].

We assume for simplicity that ``T_{sfc}`` equals ``T``, the ambient air temperature.
For real graupel, ``T_{sfc}`` is slightly higher than ``T`` due to latent heat release
of freezing liquid particles onto the ice particle. Morrison & Milbrandt (2013)
found little sensitivity to "realistic" increases in ``T_{sfc}``.

See also [`LocalRimeDensity`](@ref CloudMicrophysics.Parameters.LocalRimeDensity).

# Extended help

 Implementation follows Cober and List (1993), Eq. 16 and 17.
 See also the P3 fortran code, `microphy_p3.f90`, Line 3315-3323,
 which extends the range of the calculation to ``R_i ≤ 12``, the upper limit of which
 then equals the solid bulk ice density, ``ρ_ice = 916.7 kg/m^3``.
 At and above ``T_{sfc} = 0``, the local rime density is ``ρ_ice``.

 Note that Morrison & Milbrandt (2015) [MorrisonMilbrandt2015](@cite) only uses this
 parameterization for collisions with cloud droplets.
 For rain drops, they use a value near the solid bulk ice density, ``ρ^* = 900 kg/m^3``.
 We do not consider this distinction, and use this parameterization for all liquid particles.
"""
function compute_local_rime_density(velocity_params, particles, ρₐ, T)
    (; T_freeze, ρ_rim_local) = particles.params
    T°C = T - T_freeze  # Convert to °C

    v_ice = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    return CollisionRimeDensity(ρ_rim_local, T°C, v_ice, v_liq)
end

"""
    CollisionRimeDensity(ρ_rim_local, T°C, v_i, v_l)

Local rime density `ρ′_rim(Dᵢ, Dₗ)` [kg/m³] of a collision, returned by [`compute_local_rime_density`](@ref).

# Fields
- `ρ_rim_local`: [`CMP.LocalRimeDensity`](@ref)
- `T°C`: temperature [°C]
- `v_i`, `v_l`: ice and liquid terminal velocity functions [m/s]
"""
struct CollisionRimeDensity{R, FT, VI, VL} <: Function
    ρ_rim_local::R
    T°C::FT
    v_i::VI
    v_l::VL
end
@inline (ρ′::CollisionRimeDensity)(Dᵢ, Dₗ) = rime_density_at(ρ′, abs(ρ′.v_i(Dᵢ) - ρ′.v_l(Dₗ)), Dₗ)

"""
    rime_density_at(ρ′, Δv, Dₗ)

Return the local rime density [kg/m³] of a collision at the fall-speed difference `Δv`.

# Arguments
- `ρ′`: [`CollisionRimeDensity`](@ref)
- `Δv`: fall-speed difference `|vᵢ - vₗ|` [m/s]
- `Dₗ`: drop diameter [m]
"""
@inline function rime_density_at(ρ′, Δv, Dₗ)
    μm = 1_000_000  # m to μm factor, c.f. units of rₘ in Eq. 16 in Cober and List (1993)
    Rᵢ = -(Dₗ * μm * Δv) / (2 * ρ′.T°C)  # Eq. 16 in Cober and List (1993), Rᵢ > 0 for T°C < 0
    # At and above the melting point, use the solid-ice limit Rᵢ = 12
    Rᵢ = ifelse(ρ′.T°C < 0, Rᵢ, oftype(Rᵢ, 12))
    return ρ′.ρ_rim_local(Rᵢ)
end


"""
    get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds, quad)

Return a function `liquid_integrals(Dᵢ)` of the liquid particle integrals at ice diameter `Dᵢ`.

For a [`VolumetricCollisionRate`](@ref) and a [`CollisionRimeDensity`](@ref), the integrals are split at the
fall-speed crossing, see [`crossing_integral_bounds`](@ref).

# Arguments
- `n`: liquid particle size distribution function `n(D)`
- `∂ₜV`: volumetric collision rate integrand function `∂ₜV(Dᵢ, D)`
- `m_liq`: liquid particle mass function `m_liq(D)`
- `ρ′_rim`: local rime density function `ρ′_rim(Dᵢ, D)`
- `liq_bounds`: integration bounds for liquid particles
- `quad`: quadrature rule, a `QuadratureRule`

# Notes
The function `liquid_integrals(Dᵢ)` returns a tuple `(∂ₜN_col, ∂ₜM_col, ∂ₜB_col)`
    of collision rates at `Dᵢ`, where:
- `∂ₜN_col`: number collision rate [1/s]
- `∂ₜM_col`: mass collision rate [kg/s]
- `∂ₜB_col`: rime volume collision rate [m³/s]
"""
@inline function get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds, quad)
    function liquid_integrals(Dᵢ)
        integrand = D -> begin
            ∂ₜN = ∂ₜV(Dᵢ, D) * n(D)
            ∂ₜM = ∂ₜN * m_liq(D)
            ∂ₜB = ∂ₜM / ρ′_rim(Dᵢ, D)
            return SA.SVector((∂ₜN, ∂ₜM, ∂ₜB))
        end
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, liq_bounds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end
@inline function get_liquid_integrals(
    n, ∂ₜV::VolumetricCollisionRate, m_liq, ρ′_rim::CollisionRimeDensity, liq_bounds, quad,
)
    function liquid_integrals(Dᵢ)
        v_i_at_Dᵢ = ∂ₜV.v_i(Dᵢ)
        coeffs = collision_cross_section_ice_liquid_coeffs(∂ₜV.particles, Dᵢ)
        integrand = liquid_collision_integrand(n, ∂ₜV.v_l, m_liq, ρ′_rim, v_i_at_Dᵢ, coeffs)
        bounds = crossing_integral_bounds(liq_bounds, ∂ₜV.v_l, v_i_at_Dᵢ)
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, bounds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end

"""
    liquid_collision_integrand(n, v_l, m_liq, ρ′_rim, v_i_at_Dᵢ, coeffs)

Return the integrand `D -> (∂ₜN, ∂ₜM, ∂ₜB)` of [`get_liquid_integrals`](@ref) at one ice diameter.

# Arguments
- `n`: liquid particle size distribution function `n(D)`
- `v_l`: liquid terminal velocity function [m/s]
- `m_liq`: liquid particle mass function `m_liq(D)`
- `ρ′_rim`: local rime density, a [`CollisionRimeDensity`](@ref)
- `v_i_at_Dᵢ`: terminal velocity of the ice particle [m/s]
- `coeffs`: cross-section coefficients from [`collision_cross_section_ice_liquid_coeffs`](@ref)
"""
@inline function liquid_collision_integrand(n, v_l, m_liq, ρ′_rim, v_i_at_Dᵢ, coeffs)
    return D -> begin
        Δv = abs(v_i_at_Dᵢ - v_l(D))
        ∂ₜV = evalpoly(D, coeffs) * Δv  # E = 1, as in `VolumetricCollisionRate`
        ∂ₜN = ∂ₜV * n(D)
        ∂ₜM = ∂ₜN * m_liq(D)
        ∂ₜB = ∂ₜM / rime_density_at(ρ′_rim, Δv, D)
        return SA.SVector((∂ₜN, ∂ₜM, ∂ₜB))
    end
end

"""
    crossover_diameter(v_target, v_l, D_min, D_max)

Find the diameter `D` in `[D_min, D_max]` where `v_l(D) = v_target`
"""
function crossover_diameter(v_target, v_l::F, D_min, D_max) where {F}
    FT = float(promote_type(typeof(v_target), typeof(D_min), typeof(D_max)))
    f(D) = v_l(D) - v_target
    maxiters = FT === Float32 ? 8 : 10
    sol = RS.find_zero(f,
        RS.BrentsMethod(FT(D_min), FT(D_max)), RS.CompactSolution(),
        FixedIterations{FT}(), maxiters,
    )
    return sol.root
end

"""
    crossing_integral_bounds((D_min, D_max), v_l, v_i_at_Dᵢ)

Return the liquid integration bounds `(D_min, D*, D_max)` split at the fall-speed crossing `D*`.

At `D*` a drop falls as fast as the ice particle, `v_l(D*) = v_i_at_Dᵢ`, and `|v_i_at_Dᵢ - v_l(D)|` is
not smooth.
"""
@inline function crossing_integral_bounds((D_min, D_max), v_l, v_i_at_Dᵢ)
    Dstar = crossover_diameter(v_i_at_Dᵢ, v_l, D_min, D_max)
    return (D_min, clamp(Dstar, D_min, D_max), D_max)
end

"""
    ∫liquid_ice_collisions(
        n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, quad,
    )

Computes the bulk collision rate integrands between ice and liquid particles.

# Arguments
- `n_i`: ice particle size distribution function n_i(D)
- `∂ₜM_max`: maximum freezing rate function ∂ₜM_max(Dᵢ)
- `cloud_integrals`: inner liquid integrals for cloud particles, e.g. from [`get_liquid_integrals`](@ref)
- `rain_integrals`: inner liquid integrals for rain particles, e.g. from [`get_liquid_integrals`](@ref)
- `ice_bounds`: integration bounds for ice particles, from [`integral_bounds`](@ref)
- `quad`: quadrature rule, a `QuadratureRule`

# Returns
A 9-element vector of integrated rates,
`(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL)`;
see the method below for the definition of each entry.
"""
@inline function ∫liquid_ice_collisions(
    n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, quad,
)
    function liquid_ice_collisions_integrands(Dᵢ)
        # Inner integrals over liquid particle diameters
        ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col = cloud_integrals(Dᵢ)
        ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col = rain_integrals(Dᵢ)

        # Partition the mass collisions between freezing and shedding
        ∂ₜM_col = ∂ₜM_c_col + ∂ₜM_r_col  # [kg / s]

        ∂ₜM_frz = min(∂ₜM_col, ∂ₜM_max(Dᵢ))
        f_frz = ∂ₜM_col ≥ floatmin(typeof(∂ₜM_col)) ? ∂ₜM_frz / ∂ₜM_col : zero(∂ₜM_frz)

        n = n_i(Dᵢ)
        # Integrating over `Dᵢ` gives another unit of `[m]`, so `[X / s / m]` --> `[X / s]`
        # ∂ₜX = ∫ ∂ₜX(Dᵢ) nᵢ(Dᵢ) dDᵢ
        return SA.SVector((
            n * ∂ₜM_c_col * f_frz,        # QCFRZ
            n * ∂ₜM_c_col * (1 - f_frz),  # QCSHD
            n * ∂ₜN_c_col,                # NCCOL
            n * ∂ₜM_r_col * f_frz,        # QRFRZ
            n * ∂ₜM_r_col * (1 - f_frz),  # QRSHD
            n * ∂ₜN_r_col,                # NRCOL
            n * ∂ₜM_col,                  # ∫M_col,      total collision rate
            n * ∂ₜB_c_col * f_frz,        # BCCOL,       ∂ₜB_rim source
            n * ∂ₜB_r_col * f_frz,        # BRCOL,       ∂ₜB_rim source
        ))
    end
    return integrate(liquid_ice_collisions_integrands, ice_bounds, quad)
end

"""
    ∫liquid_ice_collisions(
        psd_c, psd_r, aps, tps, vel, m_liq, state,
        L_c, N_c, L_r, N_r, ρₐ, T,
    )

Compute key liquid-ice collision rates and quantities. Used by [`bulk_liquid_ice_collision_sources`](@ref).

The method of

    ∫liquid_ice_collisions(
        psd_c, psd_r, aps, tps, vel, m_liq, integration, state,
        L_c, N_c, L_r, N_r, ρₐ, T,
    )

for the type of `state.integration` evaluates the integrals.

# Arguments
- `psd_c`: [`CMP.CloudParticlePDF_SB2006`](@ref)
- `psd_r`: [`CMP.RainParticlePDF_SB2006`](@ref)
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `vel`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `m_liq`: liquid particle mass function `m_liq(D)`
- `state`: [`P3State`](@ref)
- `L_c`: cloud liquid water content [kg/m³]
- `N_c`: cloud liquid water number concentration [1/m³]
- `L_r`: rain water content [kg/m³]
- `N_r`: rain number concentration [1/m³]
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Returns
A tuple `(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL)`, where:
1. `QCFRZ` - Cloud mass collision rate due to freezing [kg/s]
2. `QCSHD` - Cloud mass collision rate due to shedding [kg/s]
3. `NCCOL` - Cloud number collision rate [1/s]
4. `QRFRZ` - Rain mass collision rate due to freezing [kg/s]
5. `QRSHD` - Rain mass collision rate due to shedding [kg/s]
6. `NRCOL` - Rain number collision rate [1/s]
7. `∫M_col` - Total collision rate [kg/s]
8. `BCCOL` - Cloud rime volume source [m³/m³/s]
9. `BRCOL` - Rain rime volume source [m³/m³/s]
"""
@inline ∫liquid_ice_collisions(
    psd_c, psd_r, aps, tps, vel, m_liq, state,
    L_c, N_c, L_r, N_r, ρₐ, T,
) = ∫liquid_ice_collisions(
    psd_c, psd_r, aps, tps, vel, m_liq, state.integration, state,
    L_c, N_c, L_r, N_r, ρₐ, T,
)

@inline function ∫liquid_ice_collisions(
    psd_c, psd_r, aps, tps, vel, m_liq, integration::QuadratureRule, state,
    L_c, N_c, L_r, N_r, ρₐ, T,
)
    FT = eltype(state)
    (; particles, dist) = state

    # Particle size distributions
    n_c = DT.size_distribution(psd_c, L_c / ρₐ, ρₐ, N_c)  # n_c(Dₗ)
    n_r = DT.size_distribution(psd_r, L_r / ρₐ, ρₐ, N_r)  # n_r(Dₗ)
    n_i = DT.size_distribution(dist)                      # n_i(Dᵢ)

    p = liquid_ice_collision_quantile(FT)
    ice_bounds = velocity_integral_bounds(particles, dist, ice_particle_terminal_velocity(vel, particles, ρₐ); p)
    bounds_c = CM2.get_size_distribution_bounds(psd_c, L_c / ρₐ, ρₐ, N_c, p)
    bounds_r = CM2.get_size_distribution_bounds(psd_r, L_r / ρₐ, ρₐ, N_r, p)

    # Integrand components
    # NOTE: We assume collision efficiency, shape (spherical), and terminal velocity is the
    #   same for cloud and precipitating liquid particles ⟹ same volumetric collision rate, ∂ₜV
    ∂ₜV = volumetric_collision_rate_integrand(vel, particles, ρₐ)  # ∂ₜV(Dᵢ, Dₗ)
    ρ′_rim = compute_local_rime_density(vel, particles, ρₐ, T)  # ρ′_rim(Dᵢ, Dₗ)
    ∂ₜM_max = compute_max_freeze_rate(aps, tps, vel, particles, ρₐ, T)  # ∂ₜM_max(Dᵢ)

    cloud_integrals = get_liquid_integrals(n_c, ∂ₜV, m_liq, ρ′_rim, bounds_c, integration)  # (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col)
    rain_integrals = get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, integration)  # (∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)

    # Split the ice integral at the wet-growth onset diameters, where the freeze-shed partition is not smooth
    D_wet = wet_growth_onset_diameters(cloud_integrals, rain_integrals, ∂ₜM_max, first(ice_bounds), last(ice_bounds))
    ice_bounds = insert_breakpoints(ice_bounds, D_wet)

    return ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, integration)
end

"""
    liquid_ice_collision_quantile(FT)

Return the quantile `p` of the size distributions at which the liquid-ice collision integrals are cut.

The integrals extend from the quantile `p` to the quantile `1 - p` of the ice, cloud and rain size
distributions, see [`integral_bounds`](@ref).
"""
liquid_ice_collision_quantile(::Type{FT}) where {FT} = FT(1e-5)

"""
    wet_growth_onset_diameters(cloud_integrals, rain_integrals, ∂ₜM_max, D_lo, D_hi)

Return the two ice diameters in `[D_lo, D_hi]` at which the collected liquid mass rate equals `∂ₜM_max`.

The collected rate is the sum of the mass rates of `cloud_integrals` and `rain_integrals`, see
[`get_liquid_integrals`](@ref). The crossings are found by [`CloudMicrophysics.Utilities.find_zeros`](@ref) on a
log-spaced grid, and an absent crossing is returned as `D_lo`.
"""
function wet_growth_onset_diameters(cloud_integrals, rain_integrals, ∂ₜM_max, D_lo, D_hi)
    # The onset diameters enter as breakpoints only, so the search drops the dual partials
    excess_mass_rate(D) = FD.value(cloud_integrals(D)[2] + rain_integrals(D)[2] - ∂ₜM_max(D))
    D_lo, D_hi = FD.value(D_lo), FD.value(D_hi)
    FT = typeof(D_lo)
    # Scan 16 log-spaced diameters and refine each sign change by three regula falsi steps
    (; roots, n_found) = UT.find_zeros(
        l -> excess_mass_rate(exp(l)), RS.RegulaFalsiMethod, log(D_lo), log(D_hi),
        Val(16), Val(2), FixedIterations{FT}(), 3,
    )
    return ntuple(i -> i ≤ n_found ? exp(roots[i]) : D_lo, Val(2))
end

"""
    bulk_liquid_ice_collision_sources(
        psd_c, psd_r, aps, tps, vel, state,
        L_c, N_c, L_r, N_r, ρₐ, T,
    )

Computes the bulk rates for ice and liquid particle collisions.

# Arguments
- `psd_c`: a [`CMP.CloudParticlePDF_SB2006`](@ref)
- `psd_r`: a [`CMP.RainParticlePDF_SB2006`](@ref)
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: thermodynamics parameters
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `state`: [`P3State`](@ref)
- `L_c`: cloud liquid water content [kg/m³]
- `N_c`: cloud liquid water number concentration [1/m³]
- `L_r`: rain water content [kg/m³]
- `N_r`: rain number concentration [1/m³]
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Returns
A `NamedTuple` of `(; ∂ₜq_c, ∂ₜq_r, ∂ₜN_c, ∂ₜN_r, ∂ₜL_rim, ∂ₜL_ice, ∂ₜB_rim)`, where:
1. `∂ₜq_c`: cloud liquid water content tendency [kg/kg/s]
2. `∂ₜq_r`: rain water content tendency [kg/kg/s]
3. `∂ₜN_c`: cloud number concentration tendency [1/m³/s]
4. `∂ₜN_r`: rain number concentration tendency [1/m³/s]
5. `∂ₜL_rim`: riming mass tendency [kg/m³/s]
6. `∂ₜL_ice`: ice water content tendency [kg/m³/s]
7. `∂ₜB_rim`: rime volume tendency [m³/m³/s]
"""
@inline function bulk_liquid_ice_collision_sources(
    psd_c, psd_r, aps, tps, vel, state,
    L_c, N_c, L_r, N_r, ρₐ, T,
)
    FT = promote_type(eltype(state), UT.promote_typeof(L_c, N_c, L_r, N_r, ρₐ, T))
    (; τ_wet, ρ_i) = state.particles.params
    D_shd = FT(1e-3) # 1mm  # TODO: Externalize this parameter

    ρw = psd_c.ρw
    @assert ρw == psd_r.ρw "Cloud and rain should have the same liquid water density"
    m_liq(Dₗ) = ρw * CO.volume_sphere_D(Dₗ)

    rates = ∫liquid_ice_collisions(
        psd_c, psd_r, aps, tps, vel, m_liq, state,
        L_c, N_c, L_r, N_r, ρₐ, T,
    )
    (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫∂ₜM_col, BCCOL, BRCOL) = rates

    # The shed fraction drives the rime densification, which needs freezing
    f_shd = ifelse(T < TDI.T_freeze(tps), UT.guarded_quotient(QCSHD + QRSHD, ∫∂ₜM_col), zero(∫∂ₜM_col))

    # Shedding of rain
    # QRSHD = ∫∂ₜM_col - (QCFRZ + QRFRZ)
    NRSHD = QRSHD / m_liq(D_shd)
    # NCSHD = QCSHD / m_liq(D_shd)

    # Densification of rime: relax the prognostic `(L_rim, ρb_rim)` toward `(ρq_ice, ρq_ice / ρ_i)`
    (; ρq_ice, ρb_rim) = state
    (; F_rim) = state.particles
    QIWET = f_shd * ρq_ice * (1 - F_rim) / τ_wet   # densification of rime mass
    BIWET = f_shd * (ρq_ice / ρ_i - ρb_rim) / τ_wet  # densification of rime volume

    # Bulk rates
    ## Liquid phase
    ∂ₜq_c = (-QCFRZ - QCSHD) / ρₐ
    ∂ₜq_r = (-QRFRZ + QCSHD) / ρₐ
    ∂ₜN_c = -NCCOL
    ∂ₜN_r = -NRCOL + NRSHD
    ## Ice phase
    ∂ₜL_rim = QCFRZ + QRFRZ + QIWET
    ∂ₜL_ice = QCFRZ + QRFRZ
    # ∂ₜN_ice = 0
    ∂ₜB_rim = BCCOL + BRCOL + BIWET

    return @NamedTuple{∂ₜq_c::FT, ∂ₜq_r::FT, ∂ₜN_c::FT, ∂ₜN_r::FT, ∂ₜL_rim::FT, ∂ₜL_ice::FT, ∂ₜB_rim::FT}(
        (∂ₜq_c, ∂ₜq_r, ∂ₜN_c, ∂ₜN_r, ∂ₜL_rim, ∂ₜL_ice, ∂ₜB_rim)
    )
end
