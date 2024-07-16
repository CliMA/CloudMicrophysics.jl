"""
    a(D1, D2) 

 - D1 - maximum dimension of first particle 
 - D2 - maximum dimension of second particle 

 Returns the collision kernel (assumed to be of the form π(r1 + r2)^2) for the 
 two colliding particles
"""
function a(D1::FT, D2::FT) where {FT}
    # TODO make this more accurate for non-spherical particles 
    return π * (D1 / 2 + D2 / 2)^2
end

"""
    ice_collisions(type, p3, Chen2022, ρ_a, F_r, qᵣ, qᵢ, Nᵣ, Nᵢ, ρ, ρ_r, E_ri)

 - type - defines what is colliding with the ice ("cloud" or "rain")
 - p3 - a struct with P3 scheme parameters 
 - Chen2022 - a struct with terminal velocity parameters as in Chen (2022) 
 - qᵢ - mass mixing ratio of ice 
 - Nᵢ - number mixing ratio of ice 
 - q_c - mass mixing ratio of colliding species 
 - N_c - number mixing ratio of colliding species
 - ρ_a - density of air 
 - F_r - rime mass fraction (q_rim/ q_i) 
 - ρ_r - rime density (q_rim/B_rim) 
 - T - temperature (in K)
 - E_ci - collision efficiency between ice and colliding species 

 Returns the rate of collisions between cloud ice and rain 
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
    T::FT,
    E_ci = FT(1),
) where {FT}
    ρ_l = p3.ρ_l
    th = thresholds(p3, ρ_r, F_r)
    (λ, N_0) = distribution_parameter_solver(p3, qᵢ, Nᵢ, ρ_r, F_r)
    (colliding_bound, ice_bound) =
        integration_bounds(pdf, p3, eps(FT), λ, q_c * ρ_a, N_c, ρ_a)

    if T > p3.T_freeze
        f_warm(D_c, Dᵢ) =
            E_ci / ρ_a *
            CM2.particle_size_distribution(pdf, D_c, q_c, ρ_a, N_c) *
            N′ice(p3, Dᵢ, λ, N_0) *
            a(D_c, Dᵢ) *
            p3_mass(p3, Dᵢ, F_r, th) *
            vel_diff(pdf, D_c, Dᵢ, p3, Chen2022, ρ_a, F_r, th)
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
            vel_diff(pdf, D_c, Dᵢ, p3, Chen2022, ρ_a, F_r, th)
        (dqdt, error) = HC.hcubature(
            d -> f_cold(d[1], d[2]),
            (FT(0), FT(0)),
            (colliding_bound, ice_bound),
        )
    end

    return dqdt
end

"""
    dmdt_mass(p3, D, F_r, th) 

 - p3 - a struct containing p3 parameters 
 - D - maximum dimension of the particle 
 - F_r - rime mass fraction (q_rim/ q_i) 
 - th - thresholds as calculated by thresholds()

Returns the value equivalent to dm(D)/dt * 4 / D for each P3 regime
    4 / D comes from dD/dt
"""
function dmdt_mass(p3, D, F_r, th)
    D_th = D_th_helper(p3)
    if D_th > D
        return 2 * π * p3.ρ_i * D
    elseif F_r == 0
        return 4 * α_va_si(p3) * p3.β_va * D^(p3.β_va - 2)
    elseif th.D_gr > D >= D_th
        return 4 * α_va_si(p3) * p3.β_va * D^(p3.β_va - 2)
    elseif th.D_cr > D >= th.D_gr
        return 2 * π * th.ρ_g * D
    elseif D >= th.D_cr
        return 4 * α_va_si(p3) / (1 - F_r) * p3.β_va * D^(p3.β_va - 2)
    end
end

"""
    p3_melt(p3, Chen2022, aps, tps, q, N, ρ, T, ρ_a, F_r, ρ_r)

 - p3 - a struct containing p3 parameters 
 - Chen2022 - struct containing Chen 2022 velocity parameters 
 - aps - air properties
 - tps - thermodynamics parameters
 - q - mass mixing ratio of ice 
 - N - number mixing ratio of ice 
 - T - temperature (K) 
 - ρ_a - air density
 - F_r - rime mass fraction (q_rim/ q_i)
 - ρ_r - rime density (q_rim/B_rim) 

 Returns the calculated melting rate of ice
 Equivalent to the measure of QIMLT in Morrison and Mildbrandt (2015) 
"""
function p3_melt(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::FT,
    N::FT,
    T::FT,
    ρ_a::FT,
    F_r::FT,
    ρ_r::FT,
) where {FT}
    # Get constants
    (; ν_air, D_vapor, K_therm) = aps
    a = p3.vent_a
    b = p3.vent_b
    L_f = TD.latent_heat_fusion(tps, T)
    T_freeze = p3.T_freeze
    N_sc = ν_air / D_vapor
    ρ_l = p3.ρ_l

    # Get distribution values  
    th = thresholds(p3, ρ_r, F_r)
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, F_r)

    # Get bound 
    ice_bound = get_ice_bound(p3, λ, eps(FT))

    # Define function pieces
    N_re(D) =           # TODO: What is this for non-spherical particles?
        D * TV.velocity_chen(
            D,
            Chen2022.snow_ice,
            ρ_a,
            p3_mass(p3, D, F_r, th),
            p3_area(p3, D, F_r, th),
            p3.ρ_i,
        ) / ν_air
    F(D) = a + b * N_sc^(1 / 3) * N_re(D)^(1 / 2)
    dmdt(D) =
        dmdt_mass(p3, D, F_r, th) / ρ_l * K_therm / L_f * (T - T_freeze) * F(D)
    f(D) = 1 / (2 * ρ_a) * dmdt(D) * N′ice(p3, D, λ, N_0)

    # Integrate 
    (dqdt, error) = QGK.quadgk(d -> f(d), FT(0), 2 * ice_bound)

    return dqdt
end

"""
    p3_het_freezing(mass, tps, q, N, T, ρ_a, qᵥ, aero_type)

 - mass - true if calculating change in mass, false for change in number
 - tps - thermodynamics parameters
 - q - mass mixing ratio of rain
 - N - number mixing ratio of rain
 - T - temperature in K 
 - ρ_a - density of air 
 - qᵥ - mixing ratio of water vapor 
 - aero_type - type of aerosols present 

 Returns the rate of hetergoeneous freezing within rain or cloud water 
    If mass false corresponds to NRHET in Morrison and Mildbrandt (2015)
    If mass true corresponds to QRHET in Morrison and Mildbrandt (2015)
"""
function p3_rain_het_freezing(
    mass::Bool,
    pdf_r::CMP.RainParticlePDF_SB2006{FT},
    p3::PSP3,
    tps::TDP.ThermodynamicsParameters{FT},
    q::FT,
    N::FT,
    T::FT,
    ρ_a::FT,
    qᵥ::FT,
    aero_type,
) where {FT}
    ρ_w = p3.ρ_l

    Rₐ = TD.gas_constant_air(tps, TD.PhasePartition(qᵥ))
    R_v = TD.Parameters.R_v(tps)
    e = qᵥ * ρ_a * R_v / Rₐ

    a_w = CO.a_w_eT(tps, e, T)
    a_w_ice = CO.a_w_ice(tps, T)
    Δa_w = a_w - a_w_ice
    J_immersion = CM.HetIceNucleation.ABIFM_J(aero_type, Δa_w)

    bound = CM2.get_distribution_bound(pdf_r, q, N, ρ_a, eps(FT))
    if mass
        f_rain_mass(D) =
            J_immersion *
            CM2.particle_size_distribution(pdf_r, D, q * ρ_a, ρ_a, N) *
            mass_s(D, ρ_w) *
            π *
            D^2
        dqdt, = QGK.quadgk(d -> f_rain_mass(d), FT(0), 2 * bound)
        return dqdt
    else
        f_rain(D) =
            J_immersion *
            CM2.particle_size_distribution(pdf_r, D, q * ρ_a, ρ_a, N) *
            π *
            D^2
        dNdt, = QGK.quadgk(d -> f_rain(d), FT(0), 2 * bound)
        return dNdt
    end
end