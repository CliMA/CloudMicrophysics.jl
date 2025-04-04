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
    ice_melt(p3, Chen2022, aps, tps, Tₐ, ρₐ, dt; ∫kwargs...)

# Arguments
 - `dist`: a [`P3Distribution`](@ref) object
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
    dist::P3Distribution{FT},
    Chen2022::CMP.Chen2022VelType,
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    Tₐ::FT,
    ρₐ::FT,
    dt::FT;
    ∫kwargs = (;),
) where {FT}
    # Note: process not dependent on `F_liq`
    # (we want ice core shape params)
    # Get constants
    (; ν_air, D_vapor, K_therm) = aps
    L_f = TD.latent_heat_fusion(tps, Tₐ)
    N_sc = ν_air / D_vapor

    state = get_state(dist)
    params = get_parameters(dist)

    # Ice particle terminal velocity
    v(D) = ice_particle_terminal_velocity(state, D, Chen2022, ρₐ)
    # Reynolds number
    N_Re(D) = D * v(D) / ν_air
    # Ventillation factor
    (; vent_a, vent_b) = params.vent
    F_v(D) = vent_a + vent_b * N_sc^FT(1 / 3) * N_Re(D)^FT(1 / 2)

    # Integrate
    fac = 4 * K_therm / L_f * (Tₐ - params.T_freeze)
    dLdt = fac * ∫fdD(state; ∫kwargs...) do D
        ∂ice_mass_∂D(state, D) * F_v(D) * N′ice(dist, D) / D
    end

    # only consider melting (not fusion)
    dLdt = max(0, dLdt)
    # compute change of N_ice proportional to change in L
    (; N, L) = dist
    dNdt = N / L * dLdt

    # ... and don't exceed the available number and mass of water droplets
    dNdt = min(dNdt, N / dt)
    dLdt = min(dLdt, L / dt)
    return (; dNdt, dLdt)
end

function collision_cross_section_ice_liquid(state, Dᵢ, Dₗ)
    rᵢ_eff(Dᵢ) = √(ice_area(state, Dᵢ) / π)
    return π * (rᵢ_eff(Dᵢ) + Dₗ / 2)^2  # collision cross section  -- TODO: Check if this is correct
end

"""
    volumetric_collision_rate_integrand(D_ice, D_liq, ice_dist, vel, ρₐ)

Computes the integrand for the volumetric collision rate of ice with liquid [m³/s].

```math
E * K * |vᵢ - vₗ|
```

A component of integrals like

```math
∫ ∫ E * K * |vᵢ - vₗ| * N′_ice * N′_liq dDᵢ dDₗ
```

# Arguments
- `D_ice`: ice particle diameter
- `D_liq`: liquid particle diameter
- `ice_dist`: a [`P3Distribution`](@ref) object
- `vel`: a [`Chen2022VelType`](@ref) object
- `ρₐ`: air density
"""
function volumetric_collision_rate_integrand(D_ice, D_liq, ice_dist, vel, ρₐ)
    FT = eltype(ice_dist)
    (; state) = ice_dist
    E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
    K = collision_cross_section_ice_liquid(state, D_ice, D_liq)
    v_ice = p3_particle_terminal_velocity(state, D_ice, vel, ρₐ)
    v_liq = CM2.rain_particle_terminal_velocity(D_liq, vel.rain, ρₐ)
    return E * K * abs(v_ice - v_liq)
end

function dry_rate_integrand(D_ice, D_liq, ice_dist, vel, liquid_dist, ρₐ, L_liq, N_liq)  # TODO: Change from L_liq to q_liq
    ∂ₜvolume = volumetric_collision_rate_integrand(D_ice, D_liq, ice_dist, vel, ρₐ)  # [m³ / s]
    N′_liq_psd = CM2.size_distribution(liquid_dist, D_liq, L_liq / ρₐ, ρₐ, N_liq)  # [m⁻³ / m]
    m_liq = ρ_liq * CO.volume_sphere_D(D_liq)  # [kg]
    return ∂ₜvolume * m_liq * N′_liq_psd  # [kg / s / m]
end

function compute_dry_rate(D_ice, ice_dist, vel, liquid_dist, ρₐ, L_liq, N_liq)
    (dry_rate, _) = HC.hcubature(0, 1) do D_liq
        dry_rate_integrand(D_ice, D_liq, ice_dist, vel, liquid_dist, ρₐ, L_liq, N_liq)
    end
    return dry_rate  # [kg/s]
end

function compute_wet_rate(
    D_ice::FT,
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    T::FT,  # Temperature, K
    ice_dist::P3Distribution{FT},
    ρₐ::FT,  # Air density, kg/m³
    vel::CMP.Chen2022VelType,
) where {FT}
    (; T_freeze, vent) = ice_dist.state.params
    (; vent_a, vent_b) = vent
    (; ν_air, D_vapor, K_therm) = aps
    HV = TD.latent_heat_vapor(tps, T)
    Δρ = TD.q_vap_saturation(tps, T_freeze, ρₐ, TD.PhaseNonEquil{FT}) -
         TD.q_vap_saturation(tps, T,        ρₐ, TD.PhaseNonEquil{FT})
    ΔT = T - T_freeze
    HF = TD.latent_heat_fusion(tps, T)
    C = FT(4184)
    N_sc = ν_air / D_vapor
    v_terminal = ice_particle_terminal_velocity(ice_dist.state, D_ice, vel, ρₐ)
    N_Re = D_ice * v_terminal / ν_air
    a = vent_a + vent_b * N_sc^FT(1 / 3) * N_Re^FT(1 / 2)
    return 2 * FT(π) * D_ice * a * (-K_therm * ΔT + HV * D_vapor * Δρ) / (HF + C * ΔT)
end

function particle_collision_rate_with_liquid(
    ice_dist::P3Distribution{FT},
    vel::CMP.Chen2022VelType,
    cloud_dist::CMP.CloudParticlePDF_SB2006{FT},
    rain_dist::Union{CMP.RainParticlePDF_SB2006{FT}, CMP.RainParticlePDF_SB2006_limited{FT}},
    ρₐ::FT,         # Air density, kg/m³
    L_cloud::FT,    # Cloud water content, kg/m³
    N_cloud::FT,    # Cloud water number concentration, 1/m³
    L_rain::FT,     # Rain water content, kg/m³
    N_rain::FT,     # Rain water number concentration, 1/m³
    T::FT,          # Temperature, K
) where {FT}
    # This returns a function in terms of ice and liquid particle diameters.
    function ice_integrand(D_ice)
        N′_ice_psd = N′ice(ice_dist, D_ice)

        # Calculate collision-based dry growth rate [kg / s]
        ∂ₜdry_cloud = compute_dry_rate(D_ice, ice_dist, vel, cloud_dist, ρₐ, L_cloud, N_cloud)
        ∂ₜdry_rain  = compute_dry_rate(D_ice, ice_dist, vel, rain_dist,  ρₐ, L_rain,  N_rain)
        ∂ₜdry = ∂ₜdry_cloud + ∂ₜdry_rain  # [kg / s]

        # Calculate heat-transfer-limited wet growth rate [kg / s]
        ∂ₜwet = compute_wet_rate(D_ice, aps, tps, T, ice_dist, ρₐ, vel)

        # Calculate freezing and shedding components at `D_ice` [kg / s]
        freezing = min(∂ₜdry, ∂ₜwet)
        shedding = ∂ₜdry - freezing

        # Freezing:
        ∂ₜL_ice_integrand = freezing * N′_ice_psd  # [kg m⁻³ / s / m]
        ∂ₜN_ice_integrand = FT(0)  # collisions do not change the number of ice particles, they just grow
        
        # TODO: Since `mass_rain_1mm` and `ρₐ` are constants, only need to integrate `shedding * N′_ice_psd` once
        # Mass of 1-mm-sized raindrops [kg]
        mass_rain_1mm = rain_dist.ρw * CO.volume_sphere_D(FT(1e-3))  # [kg]
        # Number concentration of 1-mm-sized raindrops [1/m³]
        ∂ₜN_rain_integrand = shedding / mass_rain_1mm * N′_ice_psd  # [m⁻³ / s / m]
        # shedding / ρₐ = [kg/s / (kg/m³)] = [kg/kg m³/s]
        ∂ₜq_rain_integrand = shedding / ρₐ * N′_ice_psd  # [kg/kg / s / m]

        # Return both components as a vector
        return SVector( # Integrating over `D_ice` gives another unit of `[m]`, so `[X / s / m]` --> `[X / s]`
            # ∂ₜN_ice, ∂ₜL_ice, ∂ₜL_rim, ∂ₜB_rim, 
            ∂ₜN_ice_integrand,  # [m⁻³   / s / m]
            ∂ₜL_ice_integrand,  # [kg m⁻³ / s / m]
            # ∂ₜN_rain, ∂ₜq_rain
            ∂ₜN_rain_integrand,  # [m⁻³   / s / m]
            ∂ₜq_rain_integrand,  # [kg/kg / s / m]
        )
    end

    return integrand
end

function bulk_liquid_sink_rate_integrand(D_ice, D_liq, ice_dist, vel, liquid_dist, ρₐ, L_liq, N_liq)
    ∂ₜvolume = volumetric_collision_rate_integrand(D_ice, D_liq, ice_dist, vel, ρₐ)  # [m³ / s]
    N′_ice_psd = N′ice(ice_dist, D_ice)  # [m⁻³ / m]
    N′_liq_psd = CM2.size_distribution(liquid_dist, D_liq, L_liq / ρₐ, ρₐ, N_liq)  # [m⁻³ / m]
    volume_liq = CO.volume_sphere_D(D_liq)  # [m³]

    ∂ₜN_liq_integrand = ∂ₜvolume * N′_ice_psd * N′_liq_psd  # [m⁻³ / s / m²]
    ∂ₜq_liq_integrand = ∂ₜN_liq_integrand * volume_liq  # [kg/kg / s / m²]

    return SVector(
        ∂ₜN_liq_integrand,  # [m⁻³ / s / m²]
        ∂ₜq_liq_integrand,  # [kg/kg / s / m²]
    )
end

import CloudMicrophysics.Microphysics2M as CM2
function bulk_collision_rate_with_liquid(
    ice_dist::P3Distribution{FT},                   # Ice distribution
    vel::CMP.Chen2022VelType,                       # Velocity distribution
    cloud_dist::CMP.CloudParticlePDF_SB2006{FT},    # Cloud distribution
    rain_dist::Union{                               # Rain distribution
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
    ρₐ::FT,         # Air density, kg/m³
    L_cloud::FT,    # Cloud water content, kg/m³
    N_cloud::FT,    # Cloud water number concentration, 1/m³
    L_rain::FT,     # Rain water content, kg/m³
    N_rain::FT,     # Rain water number concentration, 1/m³
    T::FT,          # Temperature, K
) where {FT}
    # Compute the sink terms for cloud and rain
    ((∂ₜN_cloud_sink, ∂ₜq_cloud_sink), _) = HC.hcubature((0, 0), (1, 1)) do (D_ice, D_liq)
        bulk_liquid_sink_rate_integrand(D_ice, D_liq, ice_dist, vel, cloud_dist, ρₐ, L_cloud, N_cloud)
    end

    ((∂ₜN_rain_sink, ∂ₜq_rain_sink), _) = HC.hcubature((0, 0), (1, 1)) do (D_ice, D_liq)
        bulk_liquid_sink_rate_integrand(D_ice, D_liq, ice_dist, vel, rain_dist, ρₐ, L_rain, N_rain)
    end

    # Compute the source terms due to collisions with cloud and rain
    if T < T_freeze
        # Case 1: Need to diagnose wet growth
        particle_rates = particle_collision_rate_with_liquid(
            ice_dist, vel, cloud_dist, rain_dist, ρₐ, L_cloud, N_cloud, L_rain, N_rain, T
        )

        (rates, _) = HC.hcubature(particle_rates, (0, 0), (1, 1))
    else
        # Above freezing, all cloud&rain that collides with ice is shed assuming a shed drop size of 1 mm
        volume_rain_1mm = CO.volume_sphere_D(FT(1e-3))  # [kg]
        # Calculate source terms for rain

        # all collided mass becomes rain
        ∂ₜq_rain = ∂ₜq_cloud_sink + ∂ₜq_rain_sink  # [kg/kg / s]
        # the number of new rain drops calculated gives assumed size/mass
        ∂ₜN_rain = ∂ₜq_rain / volume_rain_1mm  # [m⁻³ / s]

        # ∂ₜN_ice, ∂ₜL_ice, ∂ₜN_rain, ∂ₜq_rain
        rates = SVector(FT(0), FT(0), ∂ₜN_rain, ∂ₜq_rain)
    end
    
    ∂ₜN_ice, ∂ₜL_ice, ∂ₜN_rain, ∂ₜq_rain = rates
    
    ### TODO: Implement soaking: 
    # " If wet growth conditions are diagnosed, then particles also
    #   become soaked and undergo densification with B_rim = L_rim / ρ^*, where ρ^* = 900 kg/m³.
    #   This densification is assumed to occur within one time step."
    # NOTE: This "within one time step" business doesn't sound like a good idea to me.
    ∂ₜL_rim = FT(0)
    ∂ₜB_rim = FT(0)
  
    return (;
        ∂ₜN_ice,
        ∂ₜL_ice,
        ∂ₜL_rim,
        ∂ₜB_rim,
        ∂ₜN_cloud = - ∂ₜN_cloud_sink,
        ∂ₜq_cloud = - ∂ₜq_cloud_sink,
        ∂ₜN_rain = ∂ₜN_rain - ∂ₜN_rain_sink,
        ∂ₜq_rain = ∂ₜq_rain - ∂ₜq_rain_sink,
    )
end

