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

import CloudMicrophysics.Microphysics2M as CM2
function bulk_collision_rate_with_liquid(
    ice_dist::P3Distribution{FT},
    vel::CMP.Chen2022VelType,
    liquid_dist::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
        CMP.CloudParticlePDF_SB2006{FT}
    },
    ρₐ::FT,
    Lₗ::FT,
    Nₗ::FT,
) where {FT}
    (; ρ_l, ρ_i) = ice_dist.state.params

    # dN/dt = ∫ E ⋅ K ⋅ |Vᵢ - Vₗ| ⋅ Nₗ ⋅ Nᵢ dDᵢ dDₗ
    function integrand(Dᵢ, Dₗ)
        T = typeof(Dᵢ)
        E = T(1)  # TODO - Make it a function of Dᵢ and Dₗ
        K = collision_cross_section_ice_liquid(state, Dᵢ, Dₗ)
        vᵢ = p3_particle_terminal_velocity(state, Dᵢ, vel, ρₐ)
        vₗ = terminal_velocity(Dₗ, vel.rain, ρₐ)
        Nᵢ = N′ice(ice_dist, Dᵢ)
        Nₗ = CM2.size_distribution(liquid_dist, Dₗ, Lₗ / ρₐ, ρₐ, Nₗ)
        return E * K * abs(vᵢ - vₗ) * Nᵢ * Nₗ
    end
    mₗ(Dₗ) = ρ_l * CO.volume_sphere_D(Dₗ)
    
    (dNdt_liquid, _) = HC.hcubature((0, 0), (1, 1)) do (Dᵢ, Dₗ)
        - integrand(Dᵢ, Dₗ)
    end
    (dLdt_liquid, _) = HC.hcubature((0, 0), (1, 1)) do (Dᵢ, Dₗ)
        - integrand(Dᵢ, Dₗ) * mₗ(Dₗ)
    end
    
    # How this impacts other prognostic variables:
    dNdt_ice = FT(0)  # assume no splintering/breakup (i.e. no change in N_ice)
    dLdt_ice = - dLdt_liquid
    dLdt_rim = - dLdt_liquid
    dBdt_rim = dLdt_rim / ρ_i

    return (; dNdt_ice, dLdt_ice, dLdt_rim, dBdt_rim, dNdt_liquid, dLdt_liquid)

end
