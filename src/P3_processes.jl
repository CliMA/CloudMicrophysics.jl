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

    dNdt = J * A_aer * N_lcl
    dLdt = J * A_aer * q_lcl * ρₐ

    # nucleation rates are always positive definite...
    dNdt = max(0, dNdt)
    dLdt = max(0, dLdt)

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
 - `∫kwargs...`: Additional keyword arguments passed to the quadrature rule

Returns the melting rate of ice (QIMLT in Morrison and Mildbrandt (2015)).
"""
function ice_melt(
    velocity_params, aps::CMP.AirProperties, tps::TDI.PS,
    Tₐ, ρₐ, state::P3State, logλ;
    ∫kwargs...,
)
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; K_therm) = aps
    L_f = TDI.Lf(tps, Tₐ)

    (; L_ice, N_ice) = state
    (; T_freeze, vent) = state.params

    v_term = ice_particle_terminal_velocity(velocity_params, ρₐ, state)
    F_v = CO.ventilation_factor(vent, aps, v_term)
    N′ = size_distribution(state, logλ)

    # Integrate
    fac = 4 * K_therm / L_f * (Tₐ - T_freeze)
    bnds = integral_bounds(state, logλ; p = 1e-6)
    dLdt = fac * integrate(bnds...; ∫kwargs...) do D
        ∂ice_mass_∂D(state, D) * F_v(D) * N′(D) / D
    end

    # only consider melting (not fusion)
    dLdt = max(0, dLdt)
    # compute change of N_ice proportional to change in L
    dNdt = N_ice / L_ice * dLdt

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
    function max_freeze_rate(Dᵢ)
        Tₐ ≥ T_frz && return zero(Dᵢ)  # No collisional freezing above the freezing temperature
        return 2 * (π * Dᵢ) * F_v(Dᵢ) * (K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat) / (L_f - cp_l * ΔT)
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
    get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; ∫kwargs...)

Returns a function `liquid_integrals(Dᵢ)` that computes the liquid particle integrals
    for a given ice particle diameter `Dᵢ`.

# Arguments
- `n`: liquid particle size distribution function `n(D)`
- `∂ₜV`: volumetric collision rate integrand function `∂ₜV(Dᵢ, D)`
- `m_liq`: liquid particle mass function `m_liq(D)`
- `ρ′_rim`: local rime density function `ρ′_rim(Dᵢ, D)`
- `liq_bounds`: integration bounds for liquid particles

# Keyword arguments
- `∫kwargs...`: Additional keyword arguments passed to the quadrature rule

# Notes
The function `liquid_integrals(Dᵢ)` returns a tuple `(∂ₜN_col, ∂ₜM_col, ∂ₜB_col)`
    of collision rates at `Dᵢ`, where:
- `∂ₜN_col`: number collision rate [1/s]
- `∂ₜM_col`: mass collision rate [kg/s]
- `∂ₜB_col`: rime volume collision rate [m³/s]
"""
function get_liquid_integrals(n, ∂ₜV, m_liq, ρ′_rim, liq_bounds; ∫kwargs...)
    function liquid_integrals(Dᵢ)
        (∂ₜN_col, ∂ₜM_col, ∂ₜB_col) = integrate(liq_bounds...; ∫kwargs...) do D
            return SA.SVector(
                # ∂ₜN_col = ∫ ∂ₜV ⋅ n ⋅ dD
                ∂ₜV(Dᵢ, D) * n(D),
                # ∂ₜM_col = ∫ ∂ₜV ⋅ n ⋅ m_liq ⋅ dD
                ∂ₜV(Dᵢ, D) * n(D) * m_liq(D),
                # ∂ₜB_col = ∫ ∂ₜV ⋅ n ⋅ m_liq / ρ′_rim ⋅ dD
                ∂ₜV(Dᵢ, D) * n(D) * m_liq(D) / ρ′_rim(Dᵢ, D),
            )
        end
        return ∂ₜN_col, ∂ₜM_col, ∂ₜB_col
    end
    return liquid_integrals
end

"""
    ∫liquid_ice_collisions(
        n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; ∫kwargs...
    )

Computes the bulk collision rate integrands between ice and liquid particles.

# Arguments
- `n_i`: ice particle size distribution function n_i(D)
- `∂ₜM_max`: maximum freezing rate function ∂ₜM_max(Dᵢ)
- `cloud_integrals`: an instance of [`get_liquid_integrals`](@ref) for cloud particles
- `rain_integrals`: an instance of [`get_liquid_integrals`](@ref) for rain particles
- `ice_bounds`: integration bounds for ice particles, from [`integral_bounds`](@ref)

# Keyword arguments
- `∫kwargs...`: Additional keyword arguments passed to the quadrature rule

# Returns
A tuple of 8 integrands, see [`∫liquid_ice_collisions`](@ref) for details.
"""
function ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; ∫kwargs...)
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
    return integrate(liquid_ice_collisions_integrands, ice_bounds...; ∫kwargs...)
end

"""
    ∫liquid_ice_collisions(
        state, logλ, psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T, m_liq
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
    aps, tps, vel, ρₐ, T, m_liq; ∫kwargs...,
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

    cloud_integrals = get_liquid_integrals(n_c, ∂ₜV, m_liq, ρ′_rim, bounds_c; ∫kwargs...)  # (∂ₜN_c_col, ∂ₜM_c_col, ∂ₜB_c_col)
    rain_integrals = get_liquid_integrals(n_r, ∂ₜV, m_liq, ρ′_rim, bounds_r; ∫kwargs...)  # (∂ₜN_r_col, ∂ₜM_r_col, ∂ₜB_r_col)

    return ∫liquid_ice_collisions(n_i, ∂ₜM_max, cloud_integrals, rain_integrals, ice_bounds; ∫kwargs...)
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
    params, logλ, L_ice, N_ice, F_rim, ρ_rim,
    psd_c, psd_r, L_c, N_c, L_r, N_r,
    aps, tps, vel, ρₐ, T; ∫kwargs...,
)
    FT = eltype(params)
    (; τ_wet, ρ_i) = params
    D_shd = FT(1e-3) # 1mm  # TODO: Externalize this parameter

    ρw = psd_c.ρw
    @assert ρw == psd_r.ρw "Cloud and rain should have the same liquid water density"
    m_liq(Dₗ) = ρw * CO.volume_sphere_D(Dₗ)

    # state = get_state(params; L_ice, N_ice, F_rim, ρ_rim)
    state = P3State(params, L_ice, N_ice, F_rim, ρ_rim)

    rates = ∫liquid_ice_collisions(
        state, logλ,
        psd_c, psd_r, L_c, N_c, L_r, N_r,
        aps, tps, vel, ρₐ, T, m_liq; ∫kwargs...,
    )
    (QCFRZ, QCSHD, NCCOL, QRFRZ, QRSHD, NRCOL, ∫∂ₜM_col, BCCOL, BRCOL, ∫𝟙_wet_M_col) = rates

    # Bulk wet growth fraction
    f_wet = iszero(∫∂ₜM_col) ? zero(∫∂ₜM_col) : ∫𝟙_wet_M_col / ∫∂ₜM_col

    # Shedding of rain
    # QRSHD = ∫∂ₜM_col - (QCFRZ + QRFRZ)
    NRSHD = QRSHD / m_liq(D_shd)
    # NCSHD = QCSHD / m_liq(D_shd)

    # Densification of rime
    (; L_ice, F_rim, ρ_rim) = state
    B_rim = iszero(ρ_rim) ? zero(ρ_rim) : (L_ice * F_rim) / ρ_rim  # from: ρ_rim = L_rim / B_rim
    QIWET = f_wet * L_ice * (1 - F_rim) / τ_wet   # densification of rime mass
    BIWET = f_wet * (L_ice / ρ_i - B_rim) / τ_wet  # densification of rime volume

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
    ice_self_collection(state, logλ, aps, tps, vel, ρₐ, T; ∫kwargs...)

Computes the ice self-collection (aggregation) rate, which decreases the ice number concentration
while leaving mass, rime mass, and rime volume unchanged.

# Arguments
- `state`: [`P3State`](@ref)
- `logλ`: the log of the slope parameter [log(1/m)]
- `aps`: [`CMP.AirProperties`](@ref)
- `tps`: thermodynamics parameters
- `vel`: the velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Returns
A `NamedTuple` of `(; dNdt)`, where:
1. `dNdt`: ice number concentration tendency due to self-collection [1/m³/s] (always positive or zero, represents a loss rate)
"""
function ice_self_collection(state, logλ, aps, tps, vel, ρₐ, T; ∫kwargs...)
    n_i = DT.size_distribution(state, logλ)
    ∂ₜV = volumetric_ice_ice_collision_rate_integrand(vel, ρₐ, state)

    p = eps(one(ρₐ))
    ice_bounds = integral_bounds(state, logλ; p)

    function inner_integral(D_1)
        (rate_at_D1,) = integrate(ice_bounds...; ∫kwargs...) do D_2
            return SA.SVector(∂ₜV(D_1, D_2) * n_i(D_2))
        end
        return rate_at_D1 * n_i(D_1)
    end

    (total_rate,) = integrate(inner_integral, ice_bounds...; ∫kwargs...)

    # The 0.5 factor accounts for double-counting in self-collection
    dNdt = (1 // 2) * total_rate
    return (; dNdt)
end
