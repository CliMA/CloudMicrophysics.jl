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