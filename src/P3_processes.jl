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
    # process not dependent on F_liq
    # (we want ice core shape params)
    F_liq_ = FT(0)
    # Get constants
    (; ν_air, D_vapor, K_therm) = aps
    L_f = TD.latent_heat_fusion(tps, Tₐ)
    N_sc = ν_air / D_vapor

    # Get the P3 diameter distribution...
    th = thresholds(p3, ρ_rim, F_rim)
    (λ, N_0) =
        distribution_parameter_solver(p3, L_ice, N_ice, ρ_rim, F_rim, F_liq_)
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

    f(D) = 4 * K_therm / L_f * (Tₐ - p3.T_freeze) * dmdD(D) / D * F_v(D) * N(D)
    # Integrate
    (dLdt, error) = QGK.quadgk(d -> f(d), 0, bound, rtol = FT(1e-6))

    # only consider melting (not fusion)
    dLdt = max(0, dLdt)
    # compute change of N_ice proportional to change in L
    dNdt = N_ice / L_ice * dLdt

    # ... and dont exceed the available number and mass of water droplets
    dNdt = min(dNdt, N_ice / dt)
    dLdt = min(dLdt, L_ice / dt)
    return (; dNdt, dLdt)
end

"""
    ice_collisions(type, p3, Chen2022, ρ_a, F_r, qᵣ, qᵢ, Nᵣ, Nᵢ, ρ, ρ_r, E_ri)

 - type - defines what is colliding with the ice ("cloud" or "rain")
 - p3 - a struct with P3 scheme parameters 
 - Chen2022 - a struct with terminal velocity parameters as in Chen (2022) 
 - L_p3_tot - total p3 ice mass content [kg/m3]
 - N_ice - number mixing ratio of ice 
 - L_c - mass content of colliding species
 - N_c - number mixing ratio of colliding species
 - ρ_a - density of air
 - F_rim - rime mass fraction (L_rim / L_ice)
 - ρ_rim - rime density (L_rim / B_rim)
 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - T - temperature (in K)
 - E_ci - collision efficiency between ice and colliding species

 Returns the rate of collisions between p3 ice and liquid. 
 Equivalent to the measure of QRCOL in Morrison and Mildbrandt (2015)
"""
function ice_collisions(
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
        CMP.CloudParticlePDF_SB2006{FT},
    },
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    qᵢ::FT,
    Nᵢ::FT,
    q_c::FT,
    N_c::FT,
    ρ_a::FT,
    F_r::FT,
    ρ_r::FT,
    F_liq::FT,
    T::FT,
    E_ci = FT(1),
) where {FT}
    ρ_l = p3.ρ_l
    th = thresholds(p3, ρ_rim, F_rim)
    (λ, N_0) = distribution_parameter_solver(p3, L_p3_tot, N_ice, ρ_rim, F_liq, F_rim)
    (colliding_bound, ice_bound) =
        (CM2.get_distribution_bound(pdf, ))

    if T > p3.T_freeze
        f_warm(D_c, Dᵢ) =
            E_ci / ρ_a *
            CM2.particle_size_distribution(pdf, D_c, q_c, ρ_a, N_c) *
            N′ice(p3, Dᵢ, λ, N_0) *
            a(D_c, Dᵢ) *
            p3_mass(p3, Dᵢ, F_r, F_liq, th) *
            vel_diff(pdf, D_c, Dᵢ, p3, Chen2022, ρ_a, F_r, F_liq, th)
        (dqdt, error) = HC.hcubature(
            d -> f_warm(d[1], d[2]),
            (FT(0), FT(0)),
            (colliding_bound, ice_bound),
        )
    else
        f_cold(D_c, Dᵢ) =
            E_ci / ρ_a *
            CM2.particle_size_distribution(pdf, D_c, q_c, ρ_a, N_c) *
            N′ice(p3, Dᵢ, λ, N_0) *
            a(D_c, Dᵢ) *
            mass_s(D_c, ρ_l) *
            vel_diff(pdf, D_c, Dᵢ, p3, Chen2022, ρ_a, F_r, F_liq, th)
        (dqdt, error) = HC.hcubature(
            d -> f_cold(d[1], d[2]),
            (FT(0), FT(0)),
            (colliding_bound, ice_bound),
        )
    end

    return dqdt
end