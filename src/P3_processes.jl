"""
    het_ice_nucleation(aerosol, tps, q_lcl, N_lcl, RH, T, ρₐ)

Calculate the ice nucleation rate from heterogeneous freezing due to some `aerosol`

# Arguments
  - `aerosol`: aerosol parameters (supported types: desert dust, illite, kaolinite)
  - `tps`: thermodynamics parameters
  - `q_lcl`: cloud liquid water specific content
  - `N_lcl`: cloud droplet number concentration
  - `RH`: relative humidity
  - `T`: temperature
  - `ρₐ`: air density

# Returns
- A `NamedTuple` with the fields:
  - `dNdt`: ice number concentration change rate [m⁻³ s⁻¹]
  - `dLdt`: ice content change rate [kg m⁻³ s⁻¹]
"""
function het_ice_nucleation(
    aerosol::Union{CMP.DesertDust, CMP.Illite, CMP.Kaolinite},
    tps::TDI.PS,
    q_lcl, N_lcl, RH, T, ρₐ,
)
    FT = eltype(tps)
    #TODO - Also consider rain freezing

    # Immersion freezing nucleation rate coefficient
    J = CM_HetIce.ABIFM_J(aerosol, RH - CO.a_w_ice(tps, T))

    # Assumed erosol surface area
    # TODO - Make it a parameter of ABIFM scheme
    # We could consider making it a function of the droplet size distribution
    A_aer = FT(1e-10)

    # NaN guard: if ABIFM returns a non-finite J for a degenerate input,
    # treat non-finite J as "no nucleation" (0 rate), consistent with
    # the physical limit of an unresolvable state.
    JA_aer = ifelse(isfinite(J), J * A_aer, zero(J))

    dNdt = max(0, JA_aer * N_lcl)
    dLdt = max(0, JA_aer * q_lcl * ρₐ)

    return (; dNdt, dLdt)
end

"""
    ice_melt(velocity_params, aps, tps, state, ρₐ, Tₐ)

Compute the melting rates of ice content and number concentration.

The rates are QIMLT in [MorrisonMilbrandt2015](@cite).

# Arguments
 - `velocity_params`: [`CMP.Chen2022VelType`](@ref)
 - `aps`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `state`: a [`P3State`](@ref) object
 - `ρₐ`: air density
 - `Tₐ`: temperature (K)

# Returns
A `NamedTuple` `(; dNdt, dLdt)` with the melting rates of ice number
  concentration [1/m³/s] and ice content [kg/m³/s].

A particle melts by heat conduction at `dm/dt = 2π D K_therm (Tₐ - T_freeze) F_v(D) / L_f`,
for capacitance ``C = D/2`` and without the vapor diffusion contribution in subsaturated air.
Mass and number melt at the same fractional rate, that of the size distribution,

```math
f = \\frac{∫ \\frac{dm}{dt} N'(D) \\, dD}{∫ m(D) N'(D) \\, dD}, \\qquad
\\frac{dL}{dt} = L f, \\qquad \\frac{dN}{dt} = N f,
```
with `N` the prognostic ice number.
"""
@inline function ice_melt(velocity_params, aps, tps, state, ρₐ, Tₐ)
    FT = eltype(state)
    (; particles, integration, dist) = state
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; K_therm) = aps
    L_f = TDI.Lf(tps, Tₐ)

    (; ρq_ice) = state
    (; T_freeze, vent) = particles.params

    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    F_v = CO.ventilation_factor(vent, aps, v_term)
    N′ = size_distribution(dist)

    # Integrate
    bnds = velocity_integral_bounds(particles, dist, v_term; p = 1e-6)
    ∫DFvN = integrate(D -> D * F_v(D) * N′(D), bnds, integration)

    # only consider melting (not fusion)
    ΔT = max(Tₐ - T_freeze, 0)
    distribution_melt_rate = 2 * FT(π) * K_therm * ΔT / L_f * ∫DFvN  # ∫ dm/dt N′ dD
    distribution_mass = dist.ρn * exp(logLdivN(particles, dist.logλ))  # ∫ m N′ dD
    melt_frac = UT.guarded_quotient(distribution_melt_rate, distribution_mass)
    return (; dNdt = state.ρn_ice * melt_frac, dLdt = ρq_ice * melt_frac)
end

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
A function `(D_ice, D_liq) -> E * K * |vᵢ - vₗ|` where:
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
    function integrand(D_ice::FT, D_liq::FT) where {FT}
        E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
        K = collision_cross_section_ice_liquid(particles, D_ice, D_liq)
        return E * K * abs(v_ice(D_ice) - v_liq(D_liq))
    end

    return integrand
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
    (; D_vapor, K_therm) = aps
    cp_l = TDI.cp_l(tps)
    T_frz = TDI.T_freeze(tps)
    Lᵥ = TDI.Lᵥ(tps, Tₐ)
    L_f = TDI.Lf(tps, Tₐ)
    Tₛ = T_frz  # the surface of the ice particle is assumed to be at the freezing temperature
    ΔT = Tₛ - Tₐ  # temperature difference between the surface of the ice particle and the air
    Δρᵥ_sat =
        ρₐ * (  # saturation vapor density difference between the surface of the ice particle and the air
            TDI.p2q(tps, Tₛ, ρₐ, TDI.saturation_vapor_pressure_over_ice(tps, Tₛ)) -
            TDI.p2q(tps, Tₐ, ρₐ, TDI.saturation_vapor_pressure_over_ice(tps, Tₐ))
        )
    v_term = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    F_v = CO.ventilation_factor(particles.params.vent, aps, v_term)
    # Musil (1970) dry-growth formula: the denominator `(L_f - cp_l·ΔT)`
    # represents the *net* latent heat per unit mass available to freeze a
    # colliding droplet. At Tₐ ≲ 220 K (ΔT ≳ L_f/cp_l ≈ 53 K with
    # T-dependent L_f, see Eq. A7 in Musil 1970), the denominator flips
    # sign, making `max_freeze_rate < 0` — which is unphysical. Cold air
    # is *further from* the dry/wet-growth transition, not closer to it:
    # the physical answer is `f_frz → 1` (every colliding droplet
    # freezes). We enforce that by returning `floatmax(FT)` when the
    # denominator is non-positive, so `min(∂ₜM_col, ∂ₜM_max) = ∂ₜM_col` and
    # `f_frz = 1`.
    denom = L_f - cp_l * ΔT
    function max_freeze_rate(Dᵢ)
        # fallback values typed by the promotion of the node and the captured particles
        # (mixed plain/Dual under differentiation)
        FT = UT.promote_typeof(Dᵢ, ΔT, Δρᵥ_sat, denom)
        Tₐ ≥ T_frz && return zero(FT)  # No collisional freezing above the freezing temperature
        denom > 0 || return floatmax(FT)
        return 2 * (π * Dᵢ) * F_v(Dᵢ) * (K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat) / denom
    end
    return max_freeze_rate
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
    μm = 1_000_000  # m to μm factor, c.f. units of rₘ in Eq. 16 in Cober and List (1993)

    v_ice = ice_particle_terminal_velocity(velocity_params, particles, ρₐ)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    function ρ′_rim(Dᵢ, Dₗ)
        v_term = abs(v_ice(Dᵢ) - v_liq(Dₗ))
        Rᵢ = -(Dₗ * μm * v_term) / (2 * T°C)  # Eq. 16 in Cober and List (1993), Rᵢ > 0 for T°C < 0
        # At and above the melting point, select the solid-ice limit, Rᵢ = 12, directly
        Rᵢ = ifelse(T°C < 0, Rᵢ, oftype(Rᵢ, 12))
        return ρ_rim_local(Rᵢ)
    end
    return ρ′_rim
end

"""
    get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds, quad)

Return a function `liquid_integrals(Dᵢ)` of the liquid particle integrals at ice diameter `Dᵢ`.

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
            V_val = ∂ₜV(Dᵢ, D)
            n_val = n(D)
            m_val = m_liq(D)
            term1 = V_val * n_val
            term2 = term1 * m_val
            term3 = term2 / ρ′_rim(Dᵢ, D)
            return SA.SVector((term1, term2, term3))
        end
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, liq_bounds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
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
    closed_rain_inner_NM(
        Dᵢ, v_i_at_Dᵢ, v_l, rᵢ, ρw, ai, bi, ci, D_min, D_max, N₀r, Dr_mean,
    )

Closed-form `(∂ₜN_col, ∂ₜM_col)` for the rain inner integral at one outer `Dᵢ`.
"""
function closed_rain_inner_NM(Dᵢ, v_i_at_Dᵢ, v_l::F, rᵢ, ρw, ai, bi, ci, D_min, D_max, N₀r, Dr_mean) where {F}
    FT = float(eltype(ai))
    λ = inv(Dr_mean)  # rain PSD slope: n_r(D) ∝ e^{-λ D}
    Dstar = crossover_diameter(v_i_at_Dᵢ, v_l, D_min, D_max)

    # Compute rain PSD incomplete moments weighted by ice-liquid collision
    # cross-section `K`, and sedimentation velocity difference `|vᵢ - vₗ|`
    coeffs = SA.SVector(collision_cross_section_ice_liquid_coeffs(rᵢ))
    function Iᵖ(a, b, p, α)
        acc = @inbounds coeffs[1] * gamma_inc_moment(a, b, p, α)
        @inbounds for i in 2:lastindex(coeffs)
            acc += coeffs[i] * gamma_inc_moment(a, b, p + (i - 1), α)
        end
        return acc
    end
    function flux(a, b, p)  # ≡ ∫ₐᵇ K(Dᵢ, Dₗ) ⋅ (vᵢ(Dᵢ) - vₗ(Dₗ)) ⋅ n_r(Dₗ) dDₗ
        s = v_i_at_Dᵢ * Iᵖ(a, b, p, λ)  # vᵢ ⋅ ∫ₐᵇ K ⋅ n_r dDₗ
        @inbounds for j in eachindex(ai)  # - ∫ₐᵇ K ⋅ vₗ ⋅ n_r dDₗ
            s -= ai[j] * Iᵖ(a, b, p + bi[j], λ + ci[j])
        end
        return s
    end
    crossing(p) = flux(D_min, Dstar, p) - flux(Dstar, D_max, p)  # sign flip at Dstar
    mfac = ρw * CO.volume_sphere_D(one(FT))  # m_liq(D) = mfac Dₗ³
    return (N₀r * crossing(FT(0)), N₀r * mfac * crossing(FT(3)))  # number: D⁰, mass: D³
end

"""
    get_liquid_integrals_rain_closed(
        psd_r::RainParticlePDF_SB2006, vel::Chen2022VelType,
        n_r, ρₐ, L_r, N_r, particles, ∂ₜV, m_liq, ρ′_rim, bounds_r, quad,
    )

Return a function `liquid_integrals(Dᵢ) -> (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)` for rain.

N and M are the exact incomplete-gamma closed form, and B_rim is computed by quadrature.
"""
@inline function get_liquid_integrals_rain_closed(
    psd_r::CMP.RainParticlePDF_SB2006, vel::CMP.Chen2022VelType,
    n_r, ρₐ, L_r, N_r, particles, ∂ₜV, m_liq, ρ′_rim, bounds_r, quad,
)
    FT = promote_type(eltype(particles), UT.promote_typeof(ρₐ, L_r, N_r))
    ρw = psd_r.ρw
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    ai_t, bi_t, ci_t = CO.Chen2022_vel_coeffs(vel.rain, ρₐ)
    ai, bi, ci = SA.SVector(ai_t), SA.SVector(bi_t), SA.SVector(ci_t)
    v_l = CO.particle_terminal_velocity(vel.rain, ρₐ)
    v_i = ice_particle_terminal_velocity(vel, particles, ρₐ)
    D_min, D_max = bounds_r
    zero_rates = (zero(FT), zero(FT), zero(FT))
    function liquid_integrals(Dᵢ)
        if iszero(N₀r) || !(D_max > D_min)
            return zero_rates
        end
        v_i_at_Dᵢ = v_i(Dᵢ)
        rᵢ = sqrt(ice_area(particles, Dᵢ) / π)
        ∂ₜN_col, ∂ₜM_col = closed_rain_inner_NM(
            Dᵢ, v_i_at_Dᵢ, v_l, rᵢ, ρw, ai, bi, ci,
            D_min, D_max, N₀r, Dr_mean,
        )
        if !(isfinite(∂ₜN_col) && isfinite(∂ₜM_col))
            return zero_rates
        end
        ∂ₜB_col = integrate(
            D -> ∂ₜV(Dᵢ, D) * n_r(D) * m_liq(D) / ρ′_rim(Dᵢ, D),
            bounds_r,
            quad,
        )
        return (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)
    end
    return liquid_integrals
end

@inline _rain_inner_integrals(
    psd_r::CMP.RainParticlePDF_SB2006, vel::CMP.Chen2022VelType,
    n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, particles, quad,
) = get_liquid_integrals_rain_closed(
    psd_r, vel, n_r, ρₐ, L_r, N_r, particles, ∂ₜV, m_liq, ρ′_rim, bounds_r, quad,
)
@inline _rain_inner_integrals(
    ::Any, ::Any,
    n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, particles, quad,
) = get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, quad)

"""
    ∫liquid_ice_collisions(
        n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, quad,
    )

Computes the bulk collision rate integrands between ice and liquid particles.

# Arguments
- `n_i`: ice particle size distribution function n_i(D)
- `∂ₜM_max`: maximum freezing rate function ∂ₜM_max(Dᵢ)
- `cloud_integrals`: inner liquid integrals for cloud particles, e.g. from [`get_liquid_integrals`](@ref)
- `rain_integrals`: inner liquid integrals for rain particles, e.g. from [`get_liquid_integrals`](@ref) or `get_liquid_integrals_rain_closed`
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
        f_frz = iszero(∂ₜM_col) ? zero(∂ₜM_frz) : ∂ₜM_frz / ∂ₜM_col

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
@inline function ∫liquid_ice_collisions(
    psd_c, psd_r, aps, tps, vel, m_liq, state,
    L_c, N_c, L_r, N_r, ρₐ, T,
)
    FT = eltype(state)
    (; particles, integration, dist) = state

    # Particle size distributions
    n_c = DT.size_distribution(psd_c, L_c / ρₐ, ρₐ, N_c)  # n_c(Dₗ)
    n_r = DT.size_distribution(psd_r, L_r / ρₐ, ρₐ, N_r)  # n_r(Dₗ)
    n_i = DT.size_distribution(dist)                      # n_i(Dᵢ)

    # Initialize integration buffers by evaluating a representative integral
    p = FT(0.00001)
    ice_bounds = integral_bounds(particles, dist; p)
    bounds_c = CM2.get_size_distribution_bounds(psd_c, L_c / ρₐ, ρₐ, N_c, p)
    bounds_r = CM2.get_size_distribution_bounds(psd_r, L_r / ρₐ, ρₐ, N_r, p)

    # Integrand components
    # NOTE: We assume collision efficiency, shape (spherical), and terminal velocity is the
    #   same for cloud and precipitating liquid particles ⟹ same volumetric collision rate, ∂ₜV
    ∂ₜV = volumetric_collision_rate_integrand(vel, particles, ρₐ)  # ∂ₜV(Dᵢ, Dₗ)
    ρ′_rim = compute_local_rime_density(vel, particles, ρₐ, T)  # ρ′_rim(Dᵢ, Dₗ)
    ∂ₜM_max = compute_max_freeze_rate(aps, tps, vel, particles, ρₐ, T)  # ∂ₜM_max(Dᵢ)

    cloud_integrals = get_liquid_integrals(n_c, ∂ₜV, m_liq, ρ′_rim, bounds_c, integration)  # (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col)
    # Rain inner: exact closed form for the (SB2006-exp PSD, Chen-2022) pair
    # Numerical fallback for any other PSD/velocity type.
    rain_integrals = _rain_inner_integrals(
        psd_r, vel, n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r,
        ρₐ, L_r, N_r, particles, integration,
    )  # (∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)

    return ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds, integration)
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

"""
    ice_self_collection(vel, state, ρₐ)

Compute the ice self-collection (aggregation) rate.

Self-collection decreases the ice number concentration and leaves mass, rime mass, and rime volume
unchanged.

# Arguments
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `state`: [`P3State`](@ref)
- `ρₐ`: air density [kg/m³]

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection `[1/m³/s]` (always positive or zero, represents a loss rate), the fractional rate of the size distribution times the prognostic ice number
"""
@inline function ice_self_collection(vel, state, ρₐ)
    (; particles, integration, dist) = state
    n_i = DT.size_distribution(dist)
    v_ice = ice_particle_terminal_velocity(vel, particles, ρₐ)

    p = eps(one(ρₐ))
    ice_bounds = integral_bounds(particles, dist; p)

    function inner_integral(D_1)
        v1 = v_ice(D_1)
        r1 = sqrt(ice_area(particles, D_1) / π)
        # Volumetric collision rate integrand: E · K · |v₁ − v₂| · n(D₂)
        # where E = 1 (collision efficiency),
        #       K = π (r₁ + r₂)² is the geometric collision cross-section,
        #       |v₁ − v₂| is the differential sedimentation speed,
        #       r = √(A/π) is the effective radius from projected ice area A.
        integrand = D_2 -> begin
            v2 = v_ice(D_2)
            r2 = sqrt(ice_area(particles, D_2) / π)
            K = π * (r1 + r2)^2                    # collision cross section
            return K * abs(v1 - v2) * n_i(D_2)     # E = 1 implied
        end
        # Split the inner integral at the |v1 - v2| cusp (D_2 = D_1, where the relative
        # fall speed vanishes and the integrand has a kink). Each half is then smooth, so
        # the quadrature converges like the other (cusp-free) P3 integrals rather than
        # being cusp-limited, letting a much lower node count reach target accuracy.
        D_lo, D_hi = first(ice_bounds), last(ice_bounds)
        rate_at_D1 = integrate(integrand, (D_lo, D_1), integration) + integrate(integrand, (D_1, D_hi), integration)
        return rate_at_D1 * n_i(D_1)
    end

    total_rate = integrate(inner_integral, ice_bounds, integration)

    # The 0.5 factor accounts for double-counting in self-collection
    FT = eltype(state)
    dNdt = FT(0.5) * total_rate * UT.guarded_quotient(state.ρn_ice, dist.ρn)
    return (; dNdt)
end
