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
    ice_melt(p3, Chen2022, aps, tps, L_ice, N_ice, Tₐ, ρₐ, F_rim, ρ_rim, dt)

 - p3 - a struct containing p3 parameters
 - Chen2022 - struct containing Chen 2022 velocity parameters
 - aps - air properties
 - tps - thermodynamics parameters
 - L_ice - ice content
 - N_ice - ice number concentration
 - T - temperature (K)
 - ρ_a - air density
 - F_r - rime mass fraction (q_rim/ q_i)
 - ρ_r - rime density (q_rim/B_rim)
 - dt - model time step (for limiting the tendnecy)

Returns the melting rate of ice (QIMLT in Morrison and Mildbrandt (2015)).
"""
function ice_melt(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    L_ice::FT,
    N_ice::FT,
    Tₐ::FT,
    ρₐ::FT,
    F_rim::FT,
    ρ_rim::FT,
    dt::FT,
) where {FT}
    dLdt = FT(0)
    dNdt = FT(0)
    if L_ice > eps(FT) && N_ice > eps(FT)
        # process not dependent on F_liq
        # (we want ice core shape params)
        F_liq_ = FT(0)
        # Get constants
        (; ν_air, D_vapor, K_therm) = aps
        L_f = TD.latent_heat_fusion(tps, Tₐ)
        N_sc = ν_air / D_vapor

        # Get the P3 diameter distribution...
        th = thresholds(p3, ρ_rim, F_rim)
        (λ, N_0) = distribution_parameter_solver(
            p3,
            L_ice,
            N_ice,
            ρ_rim,
            F_rim,
            F_liq_,
        )
        N(D) = N′ice(p3, D, λ, N_0)
        # ... and D_max for the integral
        bound = get_ice_bound(p3, λ, FT(1e-6))

        # Ice particle terminal velocity
        v(D) = ice_particle_terminal_velocity(p3, D, Chen2022, ρₐ, F_rim, th)
        # Reynolds number
        N_Re(D) = D * v(D) / ν_air
        # Ventillation factor
        F_v(D) = p3.vent_a + p3.vent_b * N_sc^(1 / 3) * N_Re(D)^(1 / 2)
        dmdD(D) = p3_dmdD(p3, D, F_rim, th)

        f(D) =
            4 * K_therm / L_f * (Tₐ - p3.T_freeze) * dmdD(D) / D * F_v(D) * N(D)
        # Integrate
        (dLdt, error) = QGK.quadgk(d -> f(d), 0, bound, rtol = FT(1e-6))

        # only consider melting (not fusion)
        dLdt = max(0, dLdt)
        # compute change of N_ice proportional to change in L
        dNdt = N_ice / L_ice * dLdt

        # ... and dont exceed the available number and mass of water droplets
        dNdt = min(dNdt, N_ice / dt)
        dLdt = min(dLdt, L_ice / dt)
    end
    return (; dNdt, dLdt)
end
