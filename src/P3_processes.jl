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
