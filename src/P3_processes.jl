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
    volumetric_collision_rate_integrand(state, vel, ρₐ)
    volumetric_collision_rate_integrand(state, vel, ρₐ, D_ice, D_liq)

Returns a function that computes the volumetric collision rate integrand for ice-liquid collisions [m³/s].
The returned function takes ice and liquid particle diameters as arguments.

# Arguments
- `state`: a [`P3State`](@ref) object
- `vel`: a [`Chen2022VelType`](@ref) object
- `ρₐ`: air density
- `D_ice`: ice particle diameter (only for the second method)
- `D_liq`: liquid particle diameter (only for the second method)

# Returns
A function `(D_ice, D_liq) -> E * K * |vᵢ - vₗ|` where:
- `D_ice` and `D_liq` are the (maximum) diameters of the ice and liquid particles
- `E` is the collision efficiency
- `K` is the collision cross section
- `vᵢ` and `vₗ` are the terminal velocities of ice and liquid particles

Note that `E`, `K`, `vᵢ` and `vₗ` are all, in general, functions of `D_ice` and `D_liq`.

This function is a component of integrals like

```math
∫ ∫ E * K * |vᵢ - vₗ| * N′_ice * N′_liq dDᵢ dDₗ
```

The second method signature is a convenience method to evaluate the integrand 
at specific points, `D_ice` and `D_liq`.
"""
function volumetric_collision_rate_integrand(state, vel, ρₐ)
    v_ice = ice_particle_terminal_velocity(state, vel, ρₐ)
    v_liq = CO.liquid_particle_terminal_velocity(vel, ρₐ)
    function integrand(D_ice::FT, D_liq::FT) where {FT}
        E = FT(1)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
        K = collision_cross_section_ice_liquid(state, D_ice, D_liq)    
        return E * K * abs(v_ice(D_ice) - v_liq(D_liq))
    end
    
    return integrand
end

function volumetric_collision_rate_integrand(state, vel, ρₐ, D_ice, D_liq)
    integrand = volumetric_collision_rate_integrand(state, vel, ρₐ)
    return integrand(D_ice, D_liq)
end

# function dry_rate_integrand(D_ice::FT, D_liq, state, vel, liquid_dist, ρₐ, L_liq, N_liq) where {FT}  # TODO: Change from L_liq to q_liq
#     ∂ₜvolume = volumetric_collision_rate_integrand(state, vel, ρₐ, D_ice, D_liq)  # [m³ / s]
#     m_liq = liquid_dist.ρw * CO.volume_sphere_D(D_liq)  # [kg]
#     N′_liq_psd = CM2.size_distribution(liquid_dist, L_liq / ρₐ, ρₐ, N_liq, D_liq)  # [m⁻³ / m]
#     return FT(∂ₜvolume * m_liq * N′_liq_psd)  # [kg / s / m]  # TODO: Fix type instability
# end

function get_collision_integration_bounds(::Type{FT}) where {FT}
    # TODO: Use mode/mean/std of relevant PSDs to set intermediate bounds
    mm = FT(1e-3)
    order_bnds = FT[0; 0.01mm; 0.1mm; 1mm; 10mm; 100mm; 1]
    return order_bnds
end

"""
    compute_max_freeze_rate(aps, tps, T, state, ρₐ, vel)

Computes a function that returns the maximum possible freezing rate for an 
ice particle of diameter `D_ice` [kg/s].

This represents the thermodynamic upper limit to collisional freezing, which 
occurs when the heat transfer from the ice particle to the environment is 
balanced by the latent heat of fusion.

# Arguments
- `aps`: air properties
- `tps`: thermodynamics parameters
- `T`: temperature [K]
- `state`: a [`P3State`](@ref) object
- `ρₐ`: air density [kg/m³]
- `vel`: a [`Chen2022VelType`](@ref) object

# Returns
A function that returns the maximum possible freezing rate for an ice particle 
of diameter `D_ice` [kg/s]. Evaluates to 0 if T ≥ T_freeze.
"""
function compute_max_freeze_rate(aps, tps, T, state, ρₐ, vel)  # TODO: Consistent argument order
    FT = eltype(state)
    (; T_freeze) = state.params
    if T ≥ T_freeze
        # Above freezing, no wet growth
        return Returns(FT(0))
    end
    
    (; D_vapor, K_therm) = aps
    HV = TD.latent_heat_vapor(tps, T)
    Δqᵥ_sat = TD.q_vap_saturation(tps, T_freeze, ρₐ, TD.PhaseNonEquil) -
              TD.q_vap_saturation(tps, T,        ρₐ, TD.PhaseNonEquil) # TODO: Check units w/ paper
    ΔT = T - T_freeze
    HF = TD.latent_heat_fusion(tps, T)
    C = FT(4184)
    a = ventilation_factor(state, vel, ρₐ, aps)
    max_rate(D_ice) = 2FT(π) * D_ice * a(D_ice) * (-K_therm * ΔT + HV * D_vapor * Δqᵥ_sat) / (HF + C * ΔT)
    return max_rate
end

"""
    compute_local_rime_density(T, D_ice, D_liq, v_ice, v_liq)

Provides a function that computes the local rime density [kg/m³] 

From Milbrandt and Morrison (2013), based on the laboratory measurements of 
Cober and List (1993).

# Arguments
- `state`: a [`P3State`](@ref) object
- `vel`: a [`Chen2022VelType`](@ref) object
- `ρₐ`: air density [kg/m³]
- `T`: temperature [K]

# Returns
A function that computes the local rime density [kg/m³] using the equation:

```math
ρ′_r = 0.078 + 0.184 R_i - 0.015 R_i^2
```
where
```math
R_i = ( D_{liq} ⋅ |v_{liq} - v_{ice}| ) / ( 2 T_{sfc} )
```
and ``T_{sfc}`` is the surface temperature [°C], ``D_{liq}`` is the liquid particle
diameter [μm], ``v_{liq/ice}`` is the particle terminal velocity [m/s].
So the units of ``R_i`` are [μm m s⁻¹ / °C]. The units of ``ρ′_r`` are [g cm⁻³].
This function does the appropriate unit conversions to return ``ρ′_r`` in [kg/m³].

They assume for simplicity that `T_sfc` equals `T`, the ambient air temperature.
For real graupel, `T_sfc` is slightly higher than `T` due to latent heat release 
of freezing liquid particles onto the ice particle. MM13 found little sensitivity 
to "realistic" increases in `T_sfc`.

Implementation follows MM13, "Appendix B", part c. "Riming" (p. 427), which is based on 
CL93, "4. Experimental results", part d. "Density measurement" (p. 1599)

Note that MM15 only uses this parameterization for collisions with cloud droplets.
For rain drops, they use the solid bulk ice density, ``ρ^* = 900 kg/m^3``.
"""
function compute_local_rime_density(state, vel, ρₐ, T)
    FT = eltype(state)
    # For now, use solid bulk ice density
    return Returns(FT(900))
    
    #=
    # TODO: Implement MM13 / CL93. Draft below, but units are wrong.
    (; T_freeze) = state.params
    T°C = T - T_freeze  # Convert to °C
    # TODO: Externalize these parameters
    a, b, c = FT(0.078), FT(0.184), FT(-0.015)
    ρ′ᵣ_min, ρ′ᵣ_max = FT(50), FT(900)
    paper_to_si_factor = FT(1e-9)

    v_ice = ice_particle_terminal_velocity(state, vel, ρₐ)
    v_liq = CO.liquid_particle_terminal_velocity(vel, ρₐ)
    function compute_ρ′ᵣ(D_ice, D_liq)
        v_impact = abs(v_ice(D_ice) - v_liq(D_liq))  # Eq (B4)
        Rᵢ = (D_liq * v_impact) / (2 * T°C)          # Eq (B2)
        ρ′ᵣ = a + b * Rᵢ + c * Rᵢ^2                  # Eq (B1)
        return ρ′ᵣ * paper_to_si_factor
        # return clamp(ρ′ᵣ, ρ′ᵣ_min, ρ′ᵣ_max)
    end
    return compute_ρ′ᵣ
    =#
end

function particle_collision_rate_with_liquid(
    ice_dist::P3Distribution,
    vel::CMP.Chen2022VelType,
    cloud_dist::CMP.CloudParticlePDF_SB2006,
    rain_dist::CMP.RainParticlePDF_SB2006,
    aps::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    ρₐ,         # Air density, kg/m³
    L_cloud,    # Cloud water content, kg/m³
    N_cloud,    # Cloud water number concentration, 1/m³
    L_rain,     # Rain water content, kg/m³
    N_rain,     # Rain water number concentration, 1/m³
    T,          # Temperature, K
)
    (; state) = ice_dist
    # This returns a function in terms of ice and liquid particle diameters.
    function ice_integrand(D_ice::FT) where {FT}
        N′_ice_psd = N′ice(ice_dist, D_ice)

        # Calculate collision-based dry growth rate [kg / s]
        bounds = get_collision_integration_bounds(FT)
        (∂ₜdry_cloud, _) = QGK.quadgk_segbuf(bounds...) do D_liq
            dry_rate_integrand(D_ice, D_liq, state, vel, cloud_dist, ρₐ, L_cloud, N_cloud)
        end
        (∂ₜdry_rain, _) = QGK.quadgk_segbuf(bounds...) do D_liq
            dry_rate_integrand(D_ice, D_liq, state, vel, rain_dist, ρₐ, L_rain, N_rain)
        end
        ∂ₜdry = ∂ₜdry_cloud + ∂ₜdry_rain  # [kg / s]

        # Calculate heat-transfer-limited wet growth rate [kg / s]
        ∂ₜwet = compute_max_freeze_rate(D_ice, aps, tps, T, state, ρₐ, vel)

        # Calculate freezing and shedding components at `D_ice` [kg / s]
        freezing = min(∂ₜdry, ∂ₜwet)
        shedding = ∂ₜdry - freezing

        # Freezing:
        ∂ₜL_ice_integrand = freezing * N′_ice_psd  # [kg m⁻³ / s / m]
        
        # TODO: Since `mass_rain_1mm` and `ρₐ` are constants, only need to integrate `shedding * N′_ice_psd` once
        # Mass of 1-mm-sized raindrops [kg]
        mass_rain_1mm = rain_dist.ρw * CO.volume_sphere_D(FT(1e-3))  # [kg]
        # Number concentration of 1-mm-sized raindrops [1/m³]
        ∂ₜN_rain_integrand = shedding / mass_rain_1mm * N′_ice_psd  # [m⁻³ / s / m]
        # shedding / ρₐ = [kg/s / (kg/m³)] = [kg/kg m³/s]
        ∂ₜq_rain_integrand = shedding / ρₐ * N′_ice_psd  # [kg/kg / s / m]

        # Return both components as a vector
        return SA.SVector( # Integrating over `D_ice` gives another unit of `[m]`, so `[X / s / m]` --> `[X / s]`
            # ∂ₜL_ice, ∂ₜL_rim, ∂ₜB_rim, 
            ∂ₜL_ice_integrand,  # [kg m⁻³ / s / m]
            # ∂ₜN_rain, ∂ₜq_rain
            ∂ₜN_rain_integrand,  # [m⁻³   / s / m]
            ∂ₜq_rain_integrand,  # [kg/kg / s / m]
        )
    end

    return ice_integrand
end

function bulk_liquid_sink_rate_integrand(D_ice, D_liq, ice_dist, vel, liquid_dist, ρₐ, L_liq, N_liq)
    ∂ₜvolume = volumetric_collision_rate_integrand(ice_dist.state, vel, ρₐ, D_ice, D_liq)  # [m³ / s]
    N′_ice_psd = N′ice(ice_dist, D_ice)  # [m⁻³ / m]
    N′_liq_psd = CM2.size_distribution(liquid_dist, L_liq / ρₐ, ρₐ, N_liq, D_liq)  # [m⁻³ / m]
    volume_liq = CO.volume_sphere_D(D_liq)  # [m³]

    ∂ₜN_liq_integrand = ∂ₜvolume * N′_ice_psd * N′_liq_psd  # [m⁻³ / s / m²]
    ∂ₜq_liq_integrand = ∂ₜN_liq_integrand * volume_liq  # [kg/kg / s / m²]

    FT = typeof(D_ice)
    return SA.SVector{2, FT}(
        ∂ₜN_liq_integrand,  # [m⁻³ / s / m²]
        ∂ₜq_liq_integrand,  # [kg/kg / s / m²]
    )
end

function bulk_collision_rate_with_liquid(
    ice_dist::P3Distribution{FT},                   # Ice distribution
    vel::CMP.Chen2022VelType,                       # Velocity distribution
    cloud_dist::CMP.CloudParticlePDF_SB2006,    # Cloud distribution
    rain_dist::CMP.RainParticlePDF_SB2006,      # Rain distribution
    aps::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    ρₐ,         # Air density, kg/m³
    L_cloud,    # Cloud water content, kg/m³
    N_cloud,    # Cloud water number concentration, 1/m³
    L_rain,     # Rain water content, kg/m³
    N_rain,     # Rain water number concentration, 1/m³
    T,          # Temperature, K
) where {FT}
    mm = FT(1e-3)
    (; state) = ice_dist    
    (; T_freeze) = state.params
    ρw = cloud_dist.ρw
    @assert ρw == rain_dist.ρw "Cloud and rain should have the same liquid water density"

    # Integrand components
    ∂ₜvolume = volumetric_collision_rate_integrand(state, vel, ρₐ)
    liq_mass(D_liq) = ρw * CO.volume_sphere_D(D_liq)
    D_shed = FT(1e-3)  # 1mm
    ρᵣ′ = compute_local_rime_density(state, vel, ρₐ, T)

    max_freeze_rate = compute_max_freeze_rate(aps, tps, T, state, ρₐ, vel)

    # Particle size distributions
    N′_cloud_psd = FT ∘ CM2.size_distribution(cloud_dist, L_cloud / ρₐ, ρₐ, N_cloud)
    N′_rain_psd =  FT ∘ CM2.size_distribution(rain_dist,  L_rain / ρₐ,  ρₐ, N_rain)

    # Initialize integration buffers by evaluating a representative integral
    bounds = get_collision_integration_bounds(FT)
    segbuf_liq = QGK.quadgk_segbuf(D_liq -> N′_cloud_psd(D_liq) + N′_rain_psd(D_liq), bounds...)[3]
    segbuf_ice = QGK.quadgk_segbuf(D_ice -> N′ice(ice_dist, D_ice), bounds...)[3]

    function bulk_liquid_sink_rate_at_D_ice(Nₗ′, D_ice)
        ((∂ₜq_liq_part, ∂ₜN_liq_part), _) = 
            QGK.quadgk(bounds...; eval_segbuf = segbuf_liq) do D_liq
                return (
                    # ∂ₜq_cloud = ∫ (∫ ∂ₜvol ⋅ Nₗ′ ⋅ vol ⋅ dDₗ) ⋅ Nᵢ′ ⋅ dDᵢ
                    ∂ₜvolume(D_ice, D_liq) * Nₗ′(D_liq) * CO.volume_sphere_D(D_liq),  
                    # ∂ₜN_cloud = ∫ (∫ ∂ₜvol ⋅ Nₗ′       ⋅ dDₗ) ⋅ Nᵢ′ ⋅ dDᵢ
                    ∂ₜvolume(D_ice, D_liq) * Nₗ′(D_liq),        
                )
            end
        ∂ₜliq_mass_collisions = ρw * ∂ₜN_liq_part
        return ∂ₜliq_mass_collisions, ∂ₜq_liq_part, ∂ₜN_liq_part
    end


    function ice_integrand(D_ice)
        N′_ice_psd = N′ice(ice_dist, D_ice)

        # Inner integral over liquid particle diameters
        ∂ₜcloud_mass_collisions, ∂ₜq_cloud_part, ∂ₜN_cloud_part = 
            bulk_liquid_sink_rate_at_D_ice(N′_cloud_psd, D_ice)
        ∂ₜrain_mass_collisions, ∂ₜq_rain_part, ∂ₜN_rain_part = 
            bulk_liquid_sink_rate_at_D_ice(N′_rain_psd, D_ice)
        
        # Partition the mass collisions between freezing and shedding
        ∂ₜmass_collisions = ∂ₜcloud_mass_collisions + ∂ₜrain_mass_collisions  # [kg / s]
        ∂ₜmax_mass_freezing = max_freeze_rate(D_ice)

        frac_cloud_collisions = ∂ₜcloud_mass_collisions / ∂ₜmass_collisions

        ∂ₜmass_freezing = min(∂ₜmass_collisions, ∂ₜmax_mass_freezing)
        ∂ₜmass_shedding = ∂ₜmass_collisions - ∂ₜmass_freezing
        frac_freezing = ∂ₜmass_freezing / ∂ₜmass_collisions

        # Assume that fraction of cloud/rain that is freezes/sheds is proportional to the relative proportion of cloud/rain that collides
        ∂ₜcloud_mass_freezing = ∂ₜmass_freezing * frac_cloud_collisions
        ∂ₜrain_mass_freezing  = ∂ₜmass_freezing * (1 - frac_cloud_collisions)
        
        ∂ₜcloud_mass_shedding = ∂ₜmass_shedding * frac_cloud_collisions
        ∂ₜrain_mass_shedding  = ∂ₜmass_shedding * (1 - frac_cloud_collisions)
        
        # Rate conversion factors
        shed_mass = ρw * CO.volume_sphere_D(D_shed)  # [kg]
        
        # Resulting rates
        ∂ₜL_rim_source_integrand = ∂ₜmass_freezing * N′_ice_psd               # [kg m⁻³ / s / m]
        ∂ₜN_rain_source_integrand = ∂ₜmass_shedding / shed_mass * N′_ice_psd  # [   m⁻³ / s / m]
        ∂ₜL_rain_source_integrand = ∂ₜmass_shedding * N′_ice_psd              # [kg m⁻³ / s / m]
        ∂ₜB_rim_source_integrand = ∂ₜmass_freezing / ρᵣ′
        
        # Integrating over `D_ice` gives another unit of `[m]`, so `[X / s / m]` --> `[X / s]`
        return (
            ∂ₜq_cloud_part  * N′_ice_psd,   # ∂ₜq_cloud sink (QCCOL)
            ∂ₜN_cloud_part  * N′_ice_psd,   # ∂ₜN_cloud sink (NCCOL)
            ∂ₜq_rain_part   * N′_ice_psd,   # ∂ₜq_rain sink (QRCOL)
            ∂ₜN_rain_part   * N′_ice_psd,   # ∂ₜN_rain sink (NRCOL)
            ∂ₜL_rim_source_integrand,       # ∂ₜL_rim & ∂ₜL_ice source ("QCCOL+QRCOL")
            ∂ₜN_rain_source_integrand,      # ∂ₜN_rain source (NRSHD)
            ∂ₜL_rain_source_integrand,      # ∂ₜL_rain source (QRSHD)
            ∂ₜB_rim_source_integrand,       # ∂ₜB_rim source (BIWET)
        )

    end

    (rates, _) = QGK.quadgk(ice_integrand, bounds...; eval_segbuf = segbuf_ice)
    ∂ₜq_cloud_sink, ∂ₜN_cloud_sink, ∂ₜq_rain_sink, ∂ₜN_rain_sink, ∂ₜL_rim_source, ∂ₜN_rain_source, ∂ₜL_rain_source, ∂ₜB_rim_source = rates

    # Return the rates
    return (;
        # Liquid phase
        ∂ₜq_cloud = -∂ₜq_cloud_sink,
        ∂ₜq_rain = -∂ₜq_rain_sink,
        ∂ₜN_cloud = -∂ₜN_cloud_sink,
        ∂ₜN_rain = ∂ₜN_rain_sink,
        # Ice phase
        ∂ₜL_rim = ∂ₜL_rim_source,
        ∂ₜL_ice = ∂ₜL_rim_source,
        ∂ₜN_ice = FT(0),
        ∂ₜB_rim = ∂ₜB_rim_source,
    )

    # Compute collection rates (sink terms) for cloud and rain
    ((∂ₜN_cloud_sink, ∂ₜq_cloud_sink), _) = HC.hcubature(FT[0, 0], FT[1, 1]; initdiv = 100) do (D_ice, D_liq)
        bulk_liquid_sink_rate_integrand(D_ice, D_liq, ice_dist, vel, cloud_dist, ρₐ, L_cloud, N_cloud)
    end

    ((∂ₜN_rain_sink, ∂ₜq_rain_sink), _) = HC.hcubature(FT[0, 0], FT[1, 1]; initdiv = 100) do (D_ice, D_liq)
        bulk_liquid_sink_rate_integrand(D_ice, D_liq, ice_dist, vel, rain_dist, ρₐ, L_rain, N_rain)
    end

    # Compute source terms due to collisions with cloud and rain
    if T < T_freeze
        # Case 1: Need to diagnose wet growth
        particle_rates = particle_collision_rate_with_liquid(
            ice_dist, vel, cloud_dist, rain_dist, aps, tps, ρₐ, L_cloud, N_cloud, L_rain, N_rain, T
        )
        (rates, _) = QGK.quadgk(particle_rates, bounds...; eval_segbuf = segbuf_ice)
        ∂ₜL_ice, ∂ₜN_rain, ∂ₜq_rain = rates

        # Compute rime volume mixing ratio changes
        ρ_star = FT(900)  # [kg/m³]  (solid bulk ice density)
        # TODO: Implement proper wet growth term (BIWET)
        ∂ₜB_rim = ∂ₜL_ice / ρ_star  # [m³/m³/s]
        ∂ₜL_rim = ∂ₜL_ice  # [kg/m³/s]
    else
        # Case 2: Above freezing, all collected cloud and rain mass is shed as 1mm drops
        volume_rain_1mm = CO.volume_sphere_D(FT(1e-3))  # [m³]
        ∂ₜq_rain = ∂ₜq_cloud_sink + ∂ₜq_rain_sink  # [kg/kg/s]
        ∂ₜN_rain = ∂ₜq_rain / volume_rain_1mm  # [m⁻³/s]
        ∂ₜL_ice = FT(0)  # [kg/m³/s]
        ∂ₜL_rim = FT(0)  # [kg/m³/s]
        ∂ₜB_rim = FT(0)  # [m³/m³/s]
    end

    ### TODO: Implement soaking: 
    # " If wet growth conditions are diagnosed, then particles also
    #   become soaked and undergo densification with B_rim = L_rim / ρ^*, where ρ^* = 900 kg/m³.
    #   This densification is assumed to occur within one time step."
    # NOTE: This "within one time step" business doesn't sound like a good idea to me.
    # ∂ₜL_rim = FT(0)
    # ∂ₜB_rim = FT(0)

    return (;
        ∂ₜN_ice = FT(0),  # collisions do not change the number of ice particles
        ∂ₜL_ice,
        ∂ₜL_rim,
        ∂ₜB_rim,
        ∂ₜN_cloud = -∂ₜN_cloud_sink,  # NCCOL
        ∂ₜq_cloud = -∂ₜq_cloud_sink,  # QCCOL
        ∂ₜN_rain = ∂ₜN_rain - ∂ₜN_rain_sink,  # NRSHD - NRCOL
        ∂ₜq_rain = ∂ₜq_rain - ∂ₜq_rain_sink,  # QCSHD - QRCOL
    )
end

#=
Notes:
- Source terms for the rime volume mixing ratio (B_rim), are calculated by the
    ratio of the process rate for L_rim and the appropriate density.
- Freezing of cloud water and rain and rime generated by collection of rain by ice
    are assumed to produce ice with a density near solid bulk ice ρ^* = 900 kg/m³.
- Wet growth (BIWET) represents an additional sink term for B_rim, whereby
    B_rim decreases (i.e. particles become soaked and undergo densification).
- We expect to see sink terms from collisions/collection for these terms:
    - B_rim 
        * +QRCOL: mass of rain collected by ice
        * +BIWET: wet growth of rime
    - L_rim
        * +QCCOL: mass of cloud water collected by ice (sink to L_cloud)
        * +QRCOL: mass of rain collected by ice (sink to L_rim)
    - L_ice
        * + all source terms for L_rim
    - N_ice
        * no collision/collection terms for N_ice
    - 
=#
