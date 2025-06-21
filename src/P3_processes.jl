"""
    het_ice_nucleation(pdf_c, p3, tps, q, N, T, ρₐ, p, aerosol)

 - aerosol - aerosol parameters (supported types: desert dust, illite, kaolinite)
 - tps - thermodynamics parameters
 - qₚ - phase partition
 - N_liq - cloud water number concentration
 - RH - relative humidity
 - T - temperature
 - ρₐ - air density
 - dt - model time step

Returns a named tuple with ice number concentration and ice content
hetergoeneous freezing rates from cloud droplets.
"""
function het_ice_nucleation(
    aerosol::Union{CMP.DesertDust, CMP.Illite, CMP.Kaolinite},
    tps::TDP.ThermodynamicsParameters{FT},
    qₚ::TD.PhasePartition{FT},
    N_liq::FT,
    RH::FT,
    T::FT,
    ρₐ::FT,
    dt::FT,
) where {FT}
    #TODO - Also consider rain freezing

    # Immersion freezing nucleation rate coefficient
    J = CM_HetIce.ABIFM_J(aerosol, RH - CO.a_w_ice(tps, T))

    # Assumed erosol surface area
    # TODO - Make it a parameter of ABIFM scheme
    # We could consider making it a function of the droplet size distribution
    A_aer = FT(1e-10)

    dNdt = J * A_aer * N_liq
    dLdt = J * A_aer * qₚ.liq * ρₐ

    # nucleation rates are always positive definite...
    dNdt = max(0, dNdt)
    dLdt = max(0, dLdt)
    # ... and dont exceed the available number and mass of water droplets
    dNdt = min(dNdt, N_liq / dt)
    dLdt = min(dLdt, qₚ.liq * ρₐ / dt)

    return (; dNdt, dLdt)
end

"""
    ice_melt(state, logλ, Chen2022, aps, tps, Tₐ, ρₐ, dt; ∫kwargs...)

# Arguments
 - `state`: a [`P3State`](@ref) object
 - `logλ`: the log of the slope parameter [log(1/m)]
 - `Chen2022`: struct containing Chen 2022 velocity parameters
 - `aps`: air properties
 - `tps`: thermodynamics parameters
 - `Tₐ`: temperature (K)
 - `ρₐ`: air density
 - `dt`: model time step (for limiting the tendnecy)

# Keyword arguments
 - `∫kwargs`: Named tuple of keyword arguments passed to [`∫fdD`](@ref)

Returns the melting rate of ice (QIMLT in Morrison and Mildbrandt (2015)).
"""
function ice_melt(
    state::P3State, logλ, Chen2022::CMP.Chen2022VelType,
    aps::CMP.AirProperties, tps::TDP.ThermodynamicsParameters,
    Tₐ, ρₐ, dt;
    ∫kwargs = (;),
)
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; K_therm) = aps
    L_f = TD.latent_heat_fusion(tps, Tₐ)

    (; L_ice, N_ice) = state
    (; T_freeze, vent) = state.params

    v_term = ice_particle_terminal_velocity(state, Chen2022, ρₐ)
    F_v = CO.ventilation_factor(vent, aps, v_term)
    N′ = size_distribution(state, logλ)

    # Integrate
    fac = 4 * K_therm / L_f * (Tₐ - T_freeze)
    dLdt = fac * ∫fdD(state, logλ; ∫kwargs...) do D
        ∂ice_mass_∂D(state, D) * F_v(D) * N′(D) / D
    end

    # only consider melting (not fusion)
    dLdt = max(0, dLdt)
    # compute change of N_ice proportional to change in L
    dNdt = N_ice / L_ice * dLdt

    # ... and don't exceed the available number and mass of water droplets
    dNdt = min(dNdt, N_ice / dt)  # TODO: Apply limiters in CA.jl
    dLdt = min(dLdt, L_ice / dt)
    return (; dNdt, dLdt)
end

"""
    compute_max_freeze_rate(state, aps, tps, velocity_params, ρₐ, Tₐ)

Returns a function `max_freeze_rate(Dᵢ)` that returns the maximum possible freezing rate [kg/s] 
    for an ice particle of diameter `Dᵢ` [m]. Evaluates to 0 if T ≥ T_freeze.

# Arguments
- `state`: [`P3State`](@ref)
- `aps`: [`CMP.AirProperties`](@ref) 
- `tps`: `TDP.ThermodynamicsParameters`
- `velocity_params`: The velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `Tₐ`: air temperature [K]

This rate represents the thermodynamic upper limit to collisional freezing, 
which occurs when the heat transfer from the ice particle to the environment is 
balanced by the latent heat of fusion.

From Eq (A7) in Musil (1970), [Musil1970](@cite).
"""
function compute_max_freeze_rate(state, aps, tps, velocity_params, ρₐ, Tₐ)
    (; D_vapor, K_therm) = aps
    cp_l = TDP.cp_l(tps)
    T_frz = TDP.T_freeze(tps)
    Lᵥ = TD.latent_heat_vapor(tps, Tₐ)
    L_f = TD.latent_heat_fusion(tps, Tₐ)
    Tₛ = T_frz  # the surface of the ice particle is assumed to be at the freezing temperature
    ΔT = Tₛ - Tₐ  # temperature difference between the surface of the ice particle and the air
    Δρᵥ_sat = ρₐ * (  # saturation vapor density difference between the surface of the ice particle and the air
        TD.q_vap_saturation(tps, Tₛ, ρₐ, TD.PhaseNonEquil) - 
        TD.q_vap_saturation(tps, Tₐ, ρₐ, TD.PhaseNonEquil)
    )
    v_term = ice_particle_terminal_velocity(state, velocity_params, ρₐ)
    F_v = CO.ventilation_factor(state.params.vent, aps, v_term)
    function max_freeze_rate(Dᵢ)
        Tₐ ≥ T_frz && return zero(Dᵢ)  # No collisional freezing above the freezing temperature
        return 2 * (π * Dᵢ) * F_v(Dᵢ) * (K_therm * ΔT + Lᵥ * D_vapor * Δρᵥ_sat) / (L_f - cp_l * ΔT)
    end
    return max_freeze_rate
end

"""
    compute_local_rime_density(state, velocity_params, ρₐ, T)

Provides a function `ρ′ᵣ(Dᵢ, Dₗ)` that computes the local rime density [kg/m³] 
    for a given ice particle diameter `Dᵢ` [m] and liquid particle diameter `Dₗ` [m].

From Cober and List (1993).

# Arguments
- `state`: a [`P3State`](@ref) object
- `velocity_params`: The velocity parameterization, e.g. [`CMP.Chen2022VelType`](@ref)
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Returns
A function that computes the local rime density [kg/m³] using the equation:

```math
ρ′_r = 51 + 114 R_i - 5.5 R_i^2
```
where
```math
R_i = ( D_{liq} ⋅ |v_{liq} - v_{ice}| ) / ( 2 T_{sfc} )
```
and ``T_{sfc}`` is the surface temperature [°C], ``D_{liq}`` is the liquid particle
diameter [m], ``v_{liq/ice}`` is the particle terminal velocity [m/s].
So the units of ``R_i`` are [m² s⁻¹ °C⁻¹]. The units of ``ρ′_r`` are [kg/m³].

They assume for simplicity that `T_sfc` equals `T`, the ambient air temperature.
For real graupel, `T_sfc` is slightly higher than `T` due to latent heat release 
of freezing liquid particles onto the ice particle. MM13 found little sensitivity 
to "realistic" increases in `T_sfc`.

Implementation follows Cober and List (1993), Eq. 16 and 17.
See also the P3 fortran code, microphy_p3.f90, Line 3315-3323,
which extends the range of the calculation to Rᵢ ≤ 12, the upper limit of which
then equals the solid bulk ice density, ``ρ^⭒ = 900 kg/m^3``.

Note that MM15 only uses this parameterization for collisions with cloud droplets.
For rain drops, they use the solid bulk ice density, ``ρ^* = 900 kg/m^3``.
We do not consider this distinction, and use this parameterization for all liquid particles.
"""
function compute_local_rime_density(state, velocity_params, ρₐ, T)
    (; T_freeze) = state.params
    T°C = T - T_freeze  # Convert to °C

    v_ice = ice_particle_terminal_velocity(state, velocity_params, ρₐ)
    v_liq = CO.particle_terminal_velocity(velocity_params.rain, ρₐ)
    function ρ′ᵣ(Dᵢ, Dₗ)
        v_term = abs(v_ice(Dᵢ) - v_liq(Dₗ))
        Rᵢ = (Dₗ * v_term) / (2 * T°C)  # Eq. 16 in Cober and List (1993). Note: no `-` due to absolute value in v_term
        return ρ′ᵣ_P3(Rᵢ)
    end
    return ρ′ᵣ
end

"""
    ρ′ᵣ_P3(Rᵢ)

Returns the local rime density [kg/m³] for a given Rᵢ [m² s⁻¹ °C⁻¹].

Based on Cober and List (1993), Eq. 16 and 17 (valid for 1 ≤ Rᵢ ≤ 8).
For 8 < Rᵢ ≤ 12, linearly interpolate between ρ′ᵣ(8) ≡ 611 kg/m³ and ρ⭒ = 900 kg/m³.
See also the P3 fortran code
"""
function ρ′ᵣ_P3(Rᵢ)
    # TODO: Externalize these parameters
    a, b, c = 51, 114, -11//2 # coeffs for Eq. 17 in Cober and List (1993), converted to [kg / m³]
    ρ⭒ = 900  # ρ^⭒: density of solid bulk ice
    
    Rᵢ = clamp(Rᵢ, 1, 12)  # P3 fortran code, microphy_p3.f90, Line 3315 clamps to 1 ≤ Rᵢ ≤ 12
    
    ρ′ᵣ_CL93(Rᵢ) = a + b * Rᵢ + c * Rᵢ^2  # Eq. 17 in Cober and List (1993), in [kg / m³], valid for 1 ≤ Rᵢ ≤ 8
    ρ′ᵣ = if Rᵢ ≤ 8
        ρ′ᵣ_CL93(Rᵢ)
    else
        # following P3 fortran code, microphy_p3.f90, Line 3323
        #   https://github.com/P3-microphysics/P3-microphysics/blob/main/src/microphy_p3.f90#L3323
        # for 8 < Rᵢ ≤ 12, linearly interpolate between ρ′ᵣ(8) ≡ 611 kg/m³ and ρ⭒ = 900 kg/m³
        ρ′ᵣ8 = ρ′ᵣ_CL93(8)
        f_ρ⭒ = (Rᵢ - 8) / (12 - 8)
        (1 - f_ρ⭒) * ρ′ᵣ8 + f_ρ⭒ * ρ⭒  # Linear interpolation beyond 8.
    end
    return float(ρ′ᵣ)
end
