"""
A structure containing the rates of change of the specific humidities and number
densities of liquid and rain water.
"""
Base.@kwdef struct P3Rates{FT}
    "Rate of change of total ice mass content"
    dLdt_p3_tot::FT = FT(0)
    "Rate of change of rime mass content"
    dLdt_rim::FT = FT(0)
    "Rate of change of mass content of liquid on ice"
    dLdt_liq::FT = FT(0)
    "Rate of change of rime volume"
    dB_rim_dt::FT = FT(0)
    "Rate of change of ice number concentration"
    dNdt_ice::FT = FT(0)
    "Rate of change of rain mass content"
    dLdt_rai::FT = FT(0)
    "Rate of change of rain number concentration"
    dNdt_rai::FT = FT(0)
end

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
 - L_p3_tot- total ice content
 - N_ice - ice number concentration
 - T - temperature (K)
 - ρ_a - air density
 - F_r - rime mass fraction (q_rim/ q_i)
 - ρ_r - rime density (q_rim/B_rim)
 - dt - model time step (for limiting the tendnecy)

Returns the melting rate of ice (QIMLT in Morrison and Mildbrandt (2015))
modified for the liquid fraction scheme described in Cholette et al (2019).
"""
function ice_melt(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    L_p3_tot::FT,
    N_ice::FT,
    Tₐ::FT,
    ρₐ::FT,
    F_rim::FT,
    ρ_rim::FT,
    F_liq::FT,
    dt::FT,
) where {FT}
    # initialize rates and compute solid ice content
    (; dLdt_p3_tot, dLdt_rim, dLdt_liq, dLdt_rai, dNdt_ice, dNdt_rai) =
        P3Rates{FT}()
    L_ice = (1 - F_liq) * L_p3_tot
    if F_liq > 0.99
        # transfer total ice mass to rain
        dLdt_p3_tot = L_p3_tot
        dNdt_ice = N_ice
        dLdt_rim = F_rim * L_ice
        dLdt_rai = L_p3_tot
        dNdt_rai = N_ice
    elseif L_ice > eps(FT) && N_ice > eps(FT)
        # process not dependent on F_liq
        # (we want ice core shape params)
        # define temp F_liq_ = 0 for shape solver
        F_liq_ = FT(0)

        # Get constants
        (; ν_air, D_vapor, K_therm) = aps
        L_f = TD.latent_heat_fusion(tps, Tₐ)
        N_sc = ν_air / D_vapor

        # Get the P3 diameter distribution...
        D_th = D_th_helper(p3)
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

        # check if we have a non-negligible portion of 
        # the PSD above D_th
        small = ifelse(bound < D_th, true, false)

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
        if small
            (dLdt_rai, error) = QGK.quadgk(d -> f(d), 0, bound, rtol = FT(1e-6))
        else
            (dLdt_rai, error) = QGK.quadgk(d -> f(d), 0, D_th, rtol = FT(1e-6))
            (dLdt_liq, error) =
                QGK.quadgk(d -> f(d), D_th, bound, rtol = FT(1e-6))
        end

        # only consider melting (not fusion)
        dLdt_rai = max(FT(0), dLdt_rai)
        dLdt_liq = max(FT(0), dLdt_liq)

        # compute sink of solid ice
        dLdt_ice = dLdt_rai + dLdt_liq

        # normalize dLdt_rai and dLdt_liq so that their sum does not
        # exceed the available solid ice mass
        if dLdt_ice > (L_ice / dt)
            n = L_ice / (dLdt_ice)
            dLdt_rai *= n
            dLdt_liq *= n
            # set change in solid ice mass = max
            dLdt_ice = (L_ice / dt)
        end

        # compute change in L_rim is such that F_rim is unchanged
        dLdt_rim = F_rim * dLdt_ice


        # change in total ice mass is equal to mass lost to rain
        dLdt_p3_tot = dLdt_rai

        # compute change of N_ice and N_rai proportional to the change of L_p3_tot
        dNdt_ice = N_ice / L_ice * dLdt_p3_tot
        dNdt_rai = dNdt_ice

        # ... and don't exceed the available number
        dNdt_ice = max(0, min(dNdt_ice, N_ice / dt))
        dNdt_rai = max(0, min(dNdt_rai, N_ice / dt))
    end
    return P3Rates{FT}(
        dLdt_p3_tot = dLdt_p3_tot,
        dLdt_rim = dLdt_rim,
        dLdt_liq = dLdt_liq,
        dLdt_rai = dLdt_rai,
        dNdt_ice = dNdt_ice,
        dNdt_rai = dNdt_rai,
    )
end

"""
    ice_shed(p3, L_ice, N_ice, F_rim, ρ_rim, F_liq, dt)

 - p3 - a struct containing p3 parameters
 - L_ice - ice content
 - N_ice - ice number concentration
 - F_rim - rime mass fraction (L_rim / L_ice)
 - ρ_r - rime density (L_rim/B_rim)
 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - dt - model time step (for limiting the tendnecy)

Returns the sink of L_liq due to shedding—N_ice remains constant.
"""
function ice_shed(
    p3::PSP3,
    L_ice::FT,
    N_ice::FT,
    F_rim::FT,
    ρ_rim::FT,
    F_liq::FT,
    dt::FT,
) where {FT}
    dLdt = FT(0)
    if L_ice > eps(FT) && N_ice > eps(FT)
        # process dependent on F_liq
        # (we want whole particle shape params)
        # Get constants

        # Get the P3 diameter distribution...
        th = thresholds(p3, ρ_rim, F_rim)
        (λ, N_0) =
            distribution_parameter_solver(p3, L_ice, N_ice, ρ_rim, F_rim, F_liq)
        N(D) = N′ice(p3, D, λ, N_0)

        # ... and D_max for the integral
        bound = get_ice_bound(p3, λ, FT(1e-6))

        # critical size for shedding
        shed_bound = FT(9e-3)

        # liquid mass
        m_liq(D) = F_liq * mass_s(D, p3.ρ_l)
        # integrand (mass shed is a function of F_rim)
        f(D) = F_rim * m_liq(D) * N(D)


        # if we have no particles that are big
        # enough to undergo shedding, we return 0
        # TODO - maybe there is a better way to
        # toggle between shedding and no shedding
        # Integrate
        (dLdt, error) = ifelse(
            shed_bound >= bound,
            (FT(0), FT(0)),
            QGK.quadgk(d -> f(d), shed_bound, bound, rtol = FT(1e-6)),
        )

        # ... don't exceed the available liquid mass content
        dLdt = min(dLdt, F_liq * L_ice / dt)
    end
    # return rates struct, with dNdt
    # assuming raindrops with D = 1 mm
    return P3Rates{FT}(
        dLdt_p3_tot = dLdt,
        dLdt_liq = dLdt,
        dLdt_rai = dLdt,
        dNdt_rai = dLdt / mass_s(FT(1e-3), p3.ρ_l),
    )
end
