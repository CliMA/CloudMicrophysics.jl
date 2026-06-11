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
    ice_melt(velocity_params, aps, tps, Tₐ, ρₐ, state, logλ; ∫kwargs...)

# Arguments
 - `velocity_params`: [`CMP.Chen2022VelType`](@ref)
 - `aps`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `Tₐ`: temperature (K)
 - `ρₐ`: air density
 - `state`: a [`P3State`](@ref) object
 - `logλ`: the log of the slope parameter [log(1/m)]

# Keyword arguments
 - `quad`: quadrature rule, default is `ChebyshevGauss(100)`

Returns the melting rate of ice (QIMLT in Morrison and Mildbrandt (2015)).
"""
function ice_melt(
    velocity_params, aps::CMP.AirProperties, tps::TDI.PS,
    Tₐ, ρₐ, state::P3State, logλ;
    quad = ChebyshevGauss(100),
)
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; K_therm) = aps
    L_f = TDI.Lf(tps, Tₐ)

    (; ρq_ice, ρn_ice) = state
    (; T_freeze, vent) = state.params

    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    F_v = CO.ventilation_factor(vent, aps, v_term)
    N′ = size_distribution(state, logλ)

    # Integrate
    fac = 4 * K_therm / L_f * (Tₐ - T_freeze)
    bnds = integral_bounds(state, logλ; p = 1e-6)
    melt_integrand = D -> ∂ice_mass_∂D(state, D) * F_v(D) * N′(D) / D
    dLdt_unclamped = fac * integrate(melt_integrand, bnds, quad)

    # only consider melting (not fusion)
    dLdt = max(0, dLdt_unclamped)
    # compute change of N_ice proportional to change in mass
    dNdt = ρn_ice / ρq_ice * dLdt

    return (; dNdt, dLdt)
end

function collision_cross_section_ice_liquid(state, Dᵢ, Dₗ)
    rᵢ_eff(Dᵢ) = √(ice_area(state, Dᵢ) / π)
    return π * (rᵢ_eff(Dᵢ) + Dₗ / 2)^2  # collision cross section
end

"""
    volumetric_collision_rate_integrand(state, velocity_params, ρₐ)

Returns a function that computes the volumetric collision rate integrand for ice-liquid collisions [m³/s].
The returned function takes ice and liquid particle diameters as arguments.

# Arguments
- `state`: [`P3State`](@ref)
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density

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
function volumetric_collision_rate_integrand(velocity_params, ρₐ, state)
    v_ice = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    function integrand(D_ice::FT, D_liq::FT) where {FT}
        E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
        K = collision_cross_section_ice_liquid(state, D_ice, D_liq)
        return E * K * abs(v_ice(D_ice) - v_liq(D_liq))
    end

    return integrand
end

"""
    compute_max_freeze_rate(aps, tps, velocity_params, ρₐ, Tₐ, state)

Returns a function `max_freeze_rate(Dᵢ)` that returns the maximum possible freezing rate [kg/s]
    for an ice particle of diameter `Dᵢ` [m]. Evaluates to `0` if `T ≥ T_freeze`.

# Arguments
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `Tₐ`: air temperature [K]
- `state`: [`P3State`](@ref)

This rate represents the thermodynamic upper limit to collisional freezing,
which occurs when the heat transfer from the ice particle to the environment is
balanced by the latent heat of fusion.

From Eq (A7) in Musil (1970), [Musil1970](@cite).
"""
function compute_max_freeze_rate(aps, tps, velocity_params, ρₐ, Tₐ, state)
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
    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    F_v = CO.ventilation_factor(state.params.vent, aps, v_term)
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
        Tₐ ≥ T_frz && return zero(Dᵢ)  # No collisional freezing above the freezing temperature
        denom > 0 || return floatmax(typeof(Dᵢ))
        return 2 * (π * Dᵢ) * F_v(Dᵢ) * (K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat) / denom
    end
    return max_freeze_rate
end

"""
    compute_local_rime_density(velocity_params, ρₐ, T, state)

Provides a function `ρ′_rim(Dᵢ, Dₗ)` that computes the local rime density [kg/m³]
    for a given ice particle diameter `Dᵢ` [m] and liquid particle diameter `Dₗ` [m].

# Arguments
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]
- `state`: [`P3State`](@ref)

# Returns
A function that computes the local rime density [kg/m³] using the equation:

```math
ρ'_{rim} = a + b R_i + c R_i^2
```
where
```math
R_i = \\frac{ 10^6 ⋅ D_{liq} ⋅ |v_{liq} - v_{ice}| }{ 2 T_{sfc} }
```
and ``T_{sfc}`` is the surface temperature [°C], ``D_{liq}`` is the liquid particle
diameter [m], ``v_{liq/ice}`` is the particle terminal velocity [m/s].
So the units of ``R_i`` are [m² s⁻¹ °C⁻¹]. The units of ``ρ'_{rim}`` are [kg/m³].

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

 Note that Morrison & Milbrandt (2015) [MorrisonMilbrandt2015](@cite) only uses this
 parameterization for collisions with cloud droplets.
 For rain drops, they use a value near the solid bulk ice density, ``ρ^* = 900 kg/m^3``.
 We do not consider this distinction, and use this parameterization for all liquid particles.
"""
function compute_local_rime_density(velocity_params, ρₐ, T, state)
    (; T_freeze, ρ_rim_local) = state.params
    T°C = T - T_freeze  # Convert to °C
    μm = 1_000_000  # Note: m to μm factor, c.f. units of rₘ in Eq. 16 in Cober and List (1993)

    v_ice = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    function ρ′_rim(Dᵢ, Dₗ)
        v_term = abs(v_ice(Dᵢ) - v_liq(Dₗ))
        Rᵢ = (Dₗ * μm * v_term) / (2 * T°C)  # Eq. 16 in Cober and List (1993). Note: no `-` due to absolute value in v_term
        return ρ_rim_local(Rᵢ)
    end
    return ρ′_rim
end

"""
    get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; [quad])

Returns a function `liquid_integrals(Dᵢ)` that computes the liquid particle integrals
    for a given ice particle diameter `Dᵢ`.

# Arguments
- `n`: liquid particle size distribution function `n(D)`
- `∂ₜV`: volumetric collision rate integrand function `∂ₜV(Dᵢ, D)`
- `m_liq`: liquid particle mass function `m_liq(D)`
- `ρ′_rim`: local rime density function `ρ′_rim(Dᵢ, D)`
- `liq_bounds`: integration bounds for liquid particles

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Notes
The function `liquid_integrals(Dᵢ)` returns a tuple `(∂ₜN_col, ∂ₜM_col, ∂ₜB_col)`
    of collision rates at `Dᵢ`, where:
- `∂ₜN_col`: number collision rate [1/s]
- `∂ₜM_col`: mass collision rate [kg/s]
- `∂ₜB_col`: rime volume collision rate [m³/s]
"""
function get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; quad = ChebyshevGauss(100))
    function liquid_integrals(Dᵢ)
        integrand =
            D -> SA.SVector(
                # ∂ₜN_col = ∫ ∂ₜV ⋅ n ⋅ dD
                ∂ₜV(Dᵢ, D) * n(D),
                # ∂ₜM_col = ∫ ∂ₜV ⋅ n ⋅ m_liq ⋅ dD
                ∂ₜV(Dᵢ, D) * n(D) * m_liq(D),
                # ∂ₜB_col = ∫ ∂ₜV ⋅ n ⋅ m_liq / ρ′_rim ⋅ dD
                ∂ₜV(Dᵢ, D) * n(D) * m_liq(D) / ρ′_rim(Dᵢ, D),
            )
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(integrand, liq_bounds, quad)
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end

# ---------------------------------------------------------------------------
# Closed-form rain inner integral  (A3 / issue 003 / A4 Phase 3)
#
# The rain SB2006 PSD is pure-exponential `n_r(D)=N₀r e^{-D/Dr_mean}`
# and Chen-2022 rain `v_l(D)=Σⱼ Aⱼ D^{Bⱼ} e^{-Cⱼ D}`; the cross-section
# is a degree-2 poly in `Dₗ`. The only obstruction `|v_i(Dᵢ)−v_l(Dₗ)|`
# is removed EXACTLY by splitting the inner integral at the velocity
# crossover `D*(Dᵢ)` where `v_l(D*)=v_i(Dᵢ)`. Each piece is a finite sum
# of regularized-incomplete-gamma moments. This is exact (no physical
# approximation); the only numeric is the scalar `D*` root, found to
# ~1e-12 — far tighter than the ~1e-3 bulk-tendency budget. Validated
# vs ChebyshevGauss(1024): max relerr 5.4e-6 (N) / 2.0e-6 (M) over a
# 1944-cell grid; AD-clean under the frozen-logλ/Float64-param regime
# (`gamma_inc` 2nd-arg only). B-rim stays numerical (piecewise-rational
# `1/ρ′_rim`). See `audits/p3_audit/A3_closed_form_reductions.md` §4.

# The shared incomplete-gamma partial-moment kernel `gamma_inc_moment`
# (the linear-space twin of `loggamma_inc_moment`, used by the
# `closed_rain_inner_NM` sign-alternating split) is defined next to
# `loggamma_inc_moment` in `P3_size_distribution.jl` so the log/linear
# pair stays co-located and in sync (P1-2). See its docstring there.

"""
    crossover_diameter(v_target, v_l, D_min, D_max)

`D* ∈ [D_min,D_max]` with `v_l(D*) = v_target`, by bisection.
`v_l` is the rain terminal-velocity callable, passed in from the
**canonical** CM API `CO.particle_terminal_velocity(vel.rain, ρₐ)` —
the *exact same* Chen-2022 rain `v(D)` the numerical fallback uses in
`volumetric_collision_rate_integrand` / `compute_local_rime_density`
(`v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)`).
Reusing that callable here (rather than re-rolling the
`CO.Chen2022_monodisperse_pdf` summation) guarantees the closed-form
root solve cannot silently diverge from CM's Chen rain law if that
table/form ever changes.
Returns `D_min` if `v_target ≤ v_l(D_min)` / `D_max` if
`v_target ≥ v_l(D_max)` (the `|·|` never flips in the bracket).
Chen-2022 rain `v_l` is monotone except a tiny ρₐ-dependent dimple at
`D≈5–9 mm` with rel amplitude ≤1e-4 at the extreme PSD tail; measured
single-`D*` worst-case error ≤1.4e-6 (below the bulk budget), so a
single root is used. See A3 §3.4 for the general (deep-dimple /
non-Chen) multi-root fallback (defensive — not reached for Chen rain).
"""
function crossover_diameter(
    v_target, v_l::F, D_min, D_max;
    tol = nothing, max_iters::Int = 60,
) where {F}
    # AD-generic: `v_target` (∝ v̄ᵢ, Dual via state) and the bracket
    # `D_min,D_max` (Dual via the rain PSD) need not share a type — work
    # in their promotion (the over-tie was the latent AD blocker, P1-3).
    T = float(promote_type(typeof(v_target), typeof(D_min), typeof(D_max)))
    tol_ = tol === nothing ? T(1e-12) : T(tol)
    D_min, D_max = T(D_min), T(D_max)
    f(D) = v_l(D) - v_target
    flo, fhi = f(D_min), f(D_max)
    flo > 0 && return D_min          # v_target below the bracket — no flip
    fhi < 0 && return D_max          # v_target above the bracket — no flip
    a, b = D_min, D_max
    for _ in 1:max_iters
        m = (a + b) / 2
        fm = f(m)
        b - a < tol_ * max(b, one(T)) && return m
        if (flo < 0) == (fm < 0)
            a, flo = m, fm
        else
            b, fhi = m, fm
        end
    end
    return (a + b) / 2
end

"""
    closed_rain_inner_NM(Dᵢ, v_i_at_Dᵢ, v_l, rᵢ, ρ_w, ai, bi, ci,
                         D_min, D_max, N₀r, D̄r)

Closed-form `(∂ₜN_col, ∂ₜM_col)` for the rain inner integral at one
outer `Dᵢ`. `rᵢ = √(ice_area(state,Dᵢ)/π)`; `v_l` is the canonical
rain terminal-velocity callable `CO.particle_terminal_velocity(vel.rain,
ρₐ)` (reused for the `D*` crossover root — see [`crossover_diameter`]);
`(ai,bi,ci) = CO.Chen2022_vel_coeffs(vel.rain, ρₐ)` are the *same* Chen
coefficients that callable is built from, kept here because the closed
form's `αⱼ=α0+cⱼ` exponent merge and the `Aⱼ`/`Bⱼ+m` partial-moment
weights need the per-term `(aⱼ,bⱼ,cⱼ)` split, for which CM exposes no
existing API (`CO.Chen2022_exponential_pdf` is the *full-domain* `[0,∞)`
1M analogue — not a drop-in here, the `D*` split needs *partial* moments).
Rain PSD `n_r(D)=N₀r e^{-D/D̄r}`.
"""
function closed_rain_inner_NM(
    Dᵢ, v_i_at_Dᵢ, v_l::F, rᵢ, ρ_w, ai, bi, ci, D_min, D_max, N₀r, D̄r,
) where {F}
    # AD-generic working type. In the implicit-2M+P3 Jacobian (A2) the
    # `Dᵢ` quadrature node is Float64 while the state/PSD-derived inputs
    # (`v̄ᵢ, rᵢ, N₀r, D̄r`, the split bounds) are `Dual`; the inputs are a
    # *mix* of Float64 and Dual, so the working type is their promotion
    # (no `where {T}` over-tie — that was the latent AD blocker the
    # harness's central-FD evidence masked; P1-3 / P1-test gap 2).
    T = float(
        promote_type(
            typeof(v_i_at_Dᵢ), typeof(rᵢ), typeof(ρ_w),
            typeof(D_min), typeof(D_max), typeof(N₀r), typeof(D̄r),
            eltype(ai),
        ),
    )
    # `Tp` is the *non-Dual* exponent type (the Chen-coeff float type):
    # the `gamma_inc_moment` moment order `p` (hence `z = p+1`) MUST stay
    # a plain float, never a `Dual` — that is exactly the A2-safe
    # `SF.gamma_inc(z, α·D)` "2nd-arg-Dual only" invariant (a Dual
    # *first* arg `z` hits the missing `_gamma_inc` rule). The
    # cross-section / Chen / m_liq exponents are physical constants and
    # are never differentiated, so this is correct, not a workaround.
    Tp = float(eltype(ai))
    K0 = T(π) * rᵢ^2
    K1 = T(π) * rᵢ
    K2 = T(π) / 4
    α0 = inv(D̄r)
    Dstar = crossover_diameter(v_i_at_Dᵢ, v_l, D_min, D_max)
    mfac = ρ_w * T(π) / 6
    # `Tt`/`Tpt` are passed as `::Type` ARGUMENTS, not captured from the
    # enclosing scope. A computed *type value* (`T`, `Tp` above) captured by an
    # inner closure is boxed by Julia as `::DataType` (it loses its `Type{FT}`
    # precision), which makes `zero(Tt)`, `Tpt(0)`, … infer `::Any` and widens
    # this closure's return to `Tuple{Any,Any}` — the root of the 1.10 JET
    # failures and the GPU dynamic-dispatch / mpfr InvalidIRError on the
    # closed-form rain path. Taking them as type arguments keeps them concrete.
    function piece_NM(::Type{Tt}, ::Type{Tpt}, a, b, sgn) where {Tt, Tpt}
        b > a || return (zero(Tt), zero(Tt))
        vi_part_N =
            sgn * v_i_at_Dᵢ *
            (
                K0 * gamma_inc_moment(a, b, Tpt(0), α0) +
                K1 * gamma_inc_moment(a, b, Tpt(1), α0) +
                K2 * gamma_inc_moment(a, b, Tpt(2), α0)
            )
        vi_part_M =
            sgn * v_i_at_Dᵢ * mfac *
            (
                K0 * gamma_inc_moment(a, b, Tpt(3), α0) +
                K1 * gamma_inc_moment(a, b, Tpt(4), α0) +
                K2 * gamma_inc_moment(a, b, Tpt(5), α0)
            )
        v_part_N = zero(Tt)
        v_part_M = zero(Tt)
        @inbounds for j in eachindex(ai)
            αj = α0 + ci[j]
            Aj = ai[j]
            Bj = bi[j]
            v_part_N +=
                -sgn * Aj *
                (
                    K0 * gamma_inc_moment(a, b, Bj, αj) +
                    K1 * gamma_inc_moment(a, b, Bj + Tpt(1), αj) +
                    K2 * gamma_inc_moment(a, b, Bj + Tpt(2), αj)
                )
            v_part_M +=
                -sgn * Aj * mfac *
                (
                    K0 * gamma_inc_moment(a, b, Bj + Tpt(3), αj) +
                    K1 * gamma_inc_moment(a, b, Bj + Tpt(4), αj) +
                    K2 * gamma_inc_moment(a, b, Bj + Tpt(5), αj)
                )
        end
        return (vi_part_N + v_part_N, vi_part_M + v_part_M)
    end
    lower = piece_NM(T, Tp, D_min, Dstar, +one(T))
    upper = piece_NM(T, Tp, Dstar, D_max, -one(T))
    return (N₀r * (lower[1] + upper[1]), N₀r * (lower[2] + upper[2]))
end

"""
    get_liquid_integrals_rain_closed(psd_r, vel, n_r, ρₐ, L_r, N_r, state,
                                     ∂ₜV, m_liq, ρ′_rim, bounds_r; quad)

Closed-form rain inner integrals: a `liquid_integrals(Dᵢ)` returning
`(∂ₜN_col, ∂ₜM_col, ∂ₜB_col)` where N and M are the exact incomplete-
gamma closed form and B-rim falls back to the numerical 1-D `integrate`
(piecewise-rational `1/ρ′_rim`; see A3 §4 — deliberately not closed).
Dispatched only for `(RainParticlePDF_SB2006, Chen2022VelType)`.

The Chen-2022 rain `v_l` callable is sourced from the **canonical** CM
API `CO.particle_terminal_velocity(vel.rain, ρₐ)` — the *same*
construction the numerical fallback uses (`volumetric_collision_rate_
integrand`, `compute_local_rime_density`) — so the closed path tracks
any future change to CM's Chen law instead of forking from it. The
B-rim quadrature reuses the **passed-in `n_r` closure** (identical to
the numerical path's `∂ₜV·n_r·m_liq/ρ′_rim` integrand) rather than
re-deriving `N₀r e^{-D/D̄r}` by hand.

`α0 = inv(Dr_mean) > 0` is an invariant of the `(SB2006, Chen)`
bundle (`Dr_mean > 0` for any valid SB2006 PSD, Chen rain `cⱼ > 0` ⇒
`αⱼ = α0+cⱼ > 0`), so `gamma_inc_moment`'s `α≤0` NaN sentinel is
unreachable here. As a belt-and-suspenders matching the numerical
fallback's robustness (and the `iszero(N₀r)` guard above — one NaN
cell would otherwise abort a whole TRMM column), a degenerate
non-finite `(N,M)` is folded to `(0,0)` rather than NaN-poisoning the
outer integral.
"""
function get_liquid_integrals_rain_closed(
    psd_r::CMP.RainParticlePDF_SB2006, vel::CMP.Chen2022VelType,
    n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad,
)
    FT = eltype(state)
    ρ_w = psd_r.ρw
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(psd_r, L_r / ρₐ, ρₐ, N_r)
    ai, bi, ci = CO.Chen2022_vel_coeffs(vel.rain, ρₐ)
    # Canonical Chen-2022 rain v(D) — the exact same callable the
    # numerical fallback builds (`v_liq = CO.particle_terminal_velocity(
    # velocity_params.rain, ρₐ)`, see volumetric_collision_rate_integrand).
    v_l = CO.particle_terminal_velocity(vel.rain, ρₐ)
    v_i = ice_particle_terminal_velocity(vel, ρₐ, state)
    D_min, D_max = bounds_r
    function liquid_integrals(Dᵢ)
        if iszero(N₀r) || !(D_max > D_min)
            return (zero(FT), zero(FT), zero(FT))
        end
        vi_Dᵢ = v_i(Dᵢ)
        rᵢ = sqrt(ice_area(state, Dᵢ) / FT(π))
        ∂ₜN_col, ∂ₜM_col = closed_rain_inner_NM(
            FT(Dᵢ), vi_Dᵢ, v_l, rᵢ, ρ_w, ai, bi, ci,
            D_min, D_max, N₀r, Dr_mean,
        )
        # Degrade like the numerical path (→0) instead of propagating a
        # NaN through the outer 10-tuple if a degenerate state ever made
        # the unreachable `α≤0` sentinel fire (see docstring / P1-1).
        if !(isfinite(∂ₜN_col) && isfinite(∂ₜM_col))
            return (zero(FT), zero(FT), zero(FT))
        end
        # B-rim: piecewise-rational `1/ρ′_rim` — keep numerical (A3 §4).
        # Reuse the passed-in `n_r` closure (identical integrand to the
        # numerical `get_liquid_integrals` B-rim) — no hand re-derivation.
        ∂ₜB_col = integrate(
            D -> ∂ₜV(Dᵢ, D) * n_r(D) * m_liq(D) / ρ′_rim(Dᵢ, D),
            bounds_r,
            quad,
        )
        return (∂ₜN_col, ∂ₜM_col, ∂ₜB_col)
    end
    return liquid_integrals
end

# Dispatch: closed form for the SB2006-exp × Chen-2022 bundle, else the
# existing numerical path (behavior byte-unchanged for other bundles).
_rain_inner_integrals(
    psd_r::CMP.RainParticlePDF_SB2006, vel::CMP.Chen2022VelType,
    n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
) = get_liquid_integrals_rain_closed(
    psd_r, vel, n_r, ρₐ, L_r, N_r, state, ∂ₜV, m_liq, ρ′_rim, bounds_r;
    quad,
)
_rain_inner_integrals(
    ::Any, ::Any,
    n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r, ρₐ, L_r, N_r, state; quad,
) = get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r; quad)

"""
    ∫liquid_ice_collisions(
        n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; [quad]
    )

Computes the bulk collision rate integrands between ice and liquid particles.

# Arguments
- `n_i`: ice particle size distribution function n_i(D)
- `∂ₜM_max`: maximum freezing rate function ∂ₜM_max(Dᵢ)
- `cloud_integrals`: an instance of [`get_liquid_integrals`](@ref) for cloud particles
- `rain_integrals`: an instance of [`get_liquid_integrals`](@ref) for rain particles
- `ice_bounds`: integration bounds for ice particles, from [`integral_bounds`](@ref)

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Returns
A tuple of 8 integrands, see [`∫liquid_ice_collisions`](@ref) for details.
"""
function ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad = ChebyshevGauss(100))
    function liquid_ice_collisions_integrands(Dᵢ)
        # Inner integrals over liquid particle diameters
        ∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col = cloud_integrals(Dᵢ)
        ∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col = rain_integrals(Dᵢ)

        # Partition the mass collisions between freezing and shedding
        ∂ₜM_col = ∂ₜM_c_col + ∂ₜM_r_col  # [kg / s]

        ∂ₜM_frz = min(∂ₜM_col, ∂ₜM_max(Dᵢ))
        f_frz = iszero(∂ₜM_col) ? zero(∂ₜM_frz) : ∂ₜM_frz / ∂ₜM_col
        𝟙_wet = ∂ₜM_col > ∂ₜM_frz  # Used for wet densification

        n = n_i(Dᵢ)
        # Integrating over `Dᵢ` gives another unit of `[m]`, so `[X / s / m]` --> `[X / s]`
        # ∂ₜX = ∫ ∂ₜX(Dᵢ) nᵢ(Dᵢ) dDᵢ
        return SA.SVector(
            n * ∂ₜM_c_col * f_frz,        # QCFRZ
            n * ∂ₜM_c_col * (1 - f_frz),  # QCSHD
            n * ∂ₜN_c_col,                # NCCOL
            n * ∂ₜM_r_col * f_frz,        # QRFRZ
            n * ∂ₜM_r_col * (1 - f_frz),  # QRSHD
            n * ∂ₜN_r_col,                # NRCOL
            n * ∂ₜM_col,                  # ∫M_col,      total collision rate
            n * ∂ₜB_c_col * f_frz,        # BCCOL,       ∂ₜB_rim source
            n * ∂ₜB_r_col * f_frz,        # BRCOL,       ∂ₜB_rim source
            n * 𝟙_wet * ∂ₜM_col,          # ∫𝟙_wet_M_col, wet growth indicator
        )
    end
    return integrate(liquid_ice_collisions_integrands, ice_bounds, quad)
end

"""
    ∫liquid_ice_collisions(
        state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T, m_liq; [quad]
    )

Compute key liquid-ice collision rates and quantities. Used by [`bulk_liquid_ice_collision_sources`](@ref).

# Arguments
- `state`: [`P3State`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `psd_c`: [`CMP.CloudParticlePDF_SB2006`](@ref)
- `psd_r`: [`CMP.RainParticlePDF_SB2006`](@ref)
- `L_c`: cloud liquid water content [kg/m³]
- `N_c`: cloud liquid water number concentration [1/m³]
- `L_r`: rain water content [kg/m³]
- `N_r`: rain number concentration [1/m³]
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: `TDP.ThermodynamicsParameters`
- `vel`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]
- `m_liq`: liquid particle mass function `m_liq(D)`

# Keyword arguments
- `quad`: A `QuadratureRule` instance

# Returns
A tuple `(QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫M_col, BCCOL, BRCOL, ∫𝟙_wet_M_col)`, where:
1. `QCFRZ` - Cloud mass collision rate due to freezing [kg/s]
2. `QCSHD` - Cloud mass collision rate due to shedding [kg/s]
3. `NCCOL` - Cloud number collision rate [1/s]
4. `QRFRZ` - Rain mass collision rate due to freezing [kg/s]
5. `QRSHD` - Rain mass collision rate due to shedding [kg/s]
4. `NRCOL` - Rain number collision rate [1/s]
5. `∫M_col` - Total collision rate [kg/s]
6. `BCCOL` - Cloud rime volume source [m³/m³/s]
7. `BRCOL` - Rain rime volume source [m³/m³/s]
8. `∫𝟙_wet_M_col` - Wet growth indicator [kg/s]
"""
function ∫liquid_ice_collisions(
    state, logλ,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T, m_liq; quad,
)
    FT = eltype(state)

    # Particle size distributions
    n_c = DT.size_distribution(psd_c, L_c / ρₐ, ρₐ, N_c)  # n_c(Dₗ)
    n_r = DT.size_distribution(psd_r, L_r / ρₐ, ρₐ, N_r)  # n_r(Dₗ)
    n_i = DT.size_distribution(state, logλ)               # n_i(Dᵢ)

    # Initialize integration buffers by evaluating a representative integral
    p = FT(0.00001)
    ice_bounds = integral_bounds(state, logλ; p)
    bounds_c = CM2.get_size_distribution_bounds(psd_c, L_c / ρₐ, ρₐ, N_c, p)
    bounds_r = CM2.get_size_distribution_bounds(psd_r, L_r / ρₐ, ρₐ, N_r, p)

    # Integrand components
    # NOTE: We assume collision efficiency, shape (spherical), and terminal velocity is the
    #   same for cloud and precipitating liquid particles ⟹ same volumetric collision rate, ∂ₜV
    ∂ₜV = volumetric_collision_rate_integrand(vel, ρₐ, state)  # ∂ₜV(Dᵢ, Dₗ)
    ρ′_rim = compute_local_rime_density(vel, ρₐ, T, state)  # ρ′_rim(Dᵢ, Dₗ)
    ∂ₜM_max = compute_max_freeze_rate(aps, tps, vel, ρₐ, T, state)  # ∂ₜM_max(Dᵢ)

    cloud_integrals = get_liquid_integrals(n_c, ∂ₜV, m_liq, ρ′_rim, bounds_c; quad)  # (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col)
    # Rain inner: exact closed form for the (SB2006-exp PSD, Chen-2022)
    # bundle (N,M via incomplete gamma; B-rim numerical). Numerical
    # fallback for any other PSD/velocity type. Cloud inner + outer axis
    # unchanged (compose with GL40). See A3 §4 / A4 Phase 3 / issue 003.
    rain_integrals = _rain_inner_integrals(
        psd_r, vel, n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r,
        ρₐ, L_r, N_r, state; quad,
    )  # (∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)

    return ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; quad)
end

"""
    bulk_liquid_ice_collision_sources(
        params, logλ, L_ice, F_rim, ρ_rim,
        psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T,
    )

Computes the bulk rates for ice and liquid particle collisions.

# Arguments
- `params`: the [`CMP.ParametersP3`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `L_ice`: ice water content [kg/m³]
- `F_rim`: riming fraction
- `ρ_rim`: rime density [kg/m³]
- `psd_c`: a [`CMP.CloudParticlePDF_SB2006`](@ref)
- `psd_r`: a [`CMP.RainParticlePDF_SB2006`](@ref)
- `L_c`: cloud liquid water content [kg/m³]
- `N_c`: cloud liquid water number concentration [1/m³]
- `L_r`: rain water content [kg/m³]
- `N_r`: rain number concentration [1/m³]
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: thermodynamics parameters
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
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
function bulk_liquid_ice_collision_sources(
    state, logλ,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T; quad = ChebyshevGauss(100),
)
    FT = eltype(state)
    (; τ_wet, ρ_i) = state.params
    D_shd = FT(1e-3) # 1mm  # TODO: Externalize this parameter

    ρw = psd_c.ρw
    @assert ρw == psd_r.ρw "Cloud and rain should have the same liquid water density"
    m_liq(Dₗ) = ρw * CO.volume_sphere_D(Dₗ)

    rates = ∫liquid_ice_collisions(
        state, logλ,
        psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T, m_liq; quad,
    )
    (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫∂ₜM_col, BCCOL, BRCOL, ∫𝟙_wet_M_col) = rates

    # Bulk wet growth fraction
    f_wet = iszero(∫∂ₜM_col) ? zero(∫∂ₜM_col) : ∫𝟙_wet_M_col / ∫∂ₜM_col

    # Shedding of rain
    # QRSHD = ∫∂ₜM_col - (QCFRZ + QRFRZ)
    NRSHD = QRSHD / m_liq(D_shd)
    # NCSHD = QCSHD / m_liq(D_shd)

    # Densification of rime
    (; ρq_ice, F_rim, ρ_rim) = state
    B_rim = iszero(ρ_rim) ? zero(ρ_rim) : (ρq_ice * F_rim) / ρ_rim  # from: ρ_rim = L_rim / B_rim
    QIWET = f_wet * ρq_ice * (1 - F_rim) / τ_wet   # densification of rime mass
    BIWET = f_wet * (ρq_ice / ρ_i - B_rim) / τ_wet  # densification of rime volume

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


function collision_cross_section_ice_ice(state, D_1, D_2)
    r_eff(D) = √(ice_area(state, D) / π)
    return π * (r_eff(D_1) + r_eff(D_2))^2  # collision cross section
end

"""
    volumetric_ice_ice_collision_rate_integrand(state, velocity_params, ρₐ)

Returns a function that computes the volumetric collision rate integrand for ice-ice collisions [m³/s].

# Arguments
- `state`: [`P3State`](@ref)
- `velocity_params`: velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density

# Returns
A function `(D_1, D_2) -> E * K * |vᵢ(D_1) - vᵢ(D_2)|` where:
- `D_1` and `D_2` are the (maximum) diameters of the ice particles
- `E` is the collision efficiency
- `K` is the collision cross section
- `vᵢ` is the terminal velocity of ice particles
"""
function volumetric_ice_ice_collision_rate_integrand(velocity_params, ρₐ, state)
    v_ice = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    function integrand(D_1::FT, D_2::FT) where {FT}
        E = FT(1)  # Collision efficiency
        K = collision_cross_section_ice_ice(state, D_1, D_2)
        return E * K * abs(v_ice(D_1) - v_ice(D_2))
    end
    return integrand
end

"""
    ice_self_collection(state, logλ, vel, ρₐ; [quad])

Computes the ice self-collection (aggregation) rate, which decreases the ice number concentration
while leaving mass, rime mass, and rime volume unchanged.

# Arguments
- `state`: [`P3State`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]

# Keyword arguments
- `quad`: quadrature rule, default is `ChebyshevGauss(100)`

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection [1/m³/s] (always positive or zero, represents a loss rate)
"""
function ice_self_collection(state, logλ, vel, ρₐ; quad = ChebyshevGauss(100))
    n_i = DT.size_distribution(state, logλ)
    ∂ₜV = volumetric_ice_ice_collision_rate_integrand(vel, ρₐ, state)

    p = eps(one(ρₐ))
    ice_bounds = integral_bounds(state, logλ; p)

    function inner_integral(D_1)
        integrand = D_2 -> ∂ₜV(D_1, D_2) * n_i(D_2)
        rate_at_D1 = integrate(integrand, ice_bounds, quad)
        return rate_at_D1 * n_i(D_1)
    end

    total_rate = integrate(inner_integral, ice_bounds, quad)

    # The 0.5 factor accounts for double-counting in self-collection
    dNdt = (1 // 2) * total_rate
    return (; dNdt)
end
