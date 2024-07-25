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
 - qᵢ - mass mixing ratio of ice (includes)
 - Nᵢ - number mixing ratio of ice 
 - q_c - mass mixing ratio of colliding species 
 - N_c - number mixing ratio of colliding species
 - ρ_a - density of air 
 - F_r - rime mass fraction (q_rim/ q_i,ice) 
 - ρ_r - rime density (q_rim/B_rim)
 - F_liq - liquid fraction (q_liq/q_i,tot)
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
    F_liq::FT,
    T::FT,
    E_ci = FT(1),
) where {FT}
    ρ_l = p3.ρ_l
    th = thresholds(p3, ρ_r, F_r)
    (λ, N_0) = distribution_parameter_solver(p3, qᵢ, Nᵢ, ρ_r, F_liq, F_r)
    (colliding_bound, ice_bound) =
        integration_bounds(pdf, p3, eps(FT), λ, q_c * ρ_a, N_c, ρ_a)

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

"""
    dmdD_(p3, D, ρ, F_r)
 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension [m]
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m3]
 - F_r - rime mass fraction (q_rim/q_i)

Returns dm/dD * (4 / D)
for each particle regime
"""
# for spherical ice:
# small (D < D_th) or completely rimed (D_gr <= D < D_th)
# or for liquid on the particle: F_liq != 0
dmdD_s(D::FT, ρ::FT) where {FT} = 2 * π * ρ * D
# for unrimed nonspherical and dense nonspherical ice
dmdD_nl(p3::PSP3, D::FT) where {FT} =
    4 * α_va_si(p3) * p3.β_va * D^(p3.β_va - 2)
# for partially rimed ice
dmdD_r(p3::PSP3, D::FT, F_r::FT) where {FT} =
    4 * α_va_si(p3) / (1 - F_r) * p3.β_va * D^(p3.β_va - 2)

"""
    dmdD_mass(p3, D, F_r, F_liq, th) 
 - p3 - a struct containing p3 parameters 
 - D - maximum dimension of the particle 
 - F_r - rime mass fraction (q_rim/ q_i)
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - th - thresholds as calculated by thresholds()
Returns the value equivalent to dm(D)/dt * 4 / D for each P3 regime
    4 / D is a term from dD/dt which is useful to
    include before the other terms
"""
function dmdD_mass(p3, D, F_r, F_liq, th)
    D_th = D_th_helper(p3)
    if D_th > D
        return (1 - F_liq) * dmdD_s(D, p3.ρ_i) + F_liq * dmdD_s(D, p3.ρ_l)
    elseif F_r == 0
        return (1 - F_liq) * dmdD_nl(p3, D) + F_liq * dmdD_s(D, p3.ρ_l)
    elseif th.D_gr > D >= D_th
        return (1 - F_liq) * dmdD_nl(p3, D) + F_liq * dmdD_s(D, p3.ρ_l)
    elseif th.D_cr > D >= th.D_gr
        return (1 - F_liq) * dmdD_s(D, th.ρ_g) + F_liq * dmdD_s(D, p3.ρ_l)
    elseif D >= th.D_cr
        return (1 - F_liq) * dmdD_r(p3, D, F_r) + F_liq * dmdD_s(D, p3.ρ_l)
    end
end

"""
    p3_melt(p3, Chen2022, aps, tps, q, N, ρ, T, ρ_a, F_r, ρ_r, F_liq)
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
 - F_liq - liquid fraction (q_liq/q_i,tot)
 Returns the calculated melting rate of ice
 Equivalent to the measure of QIMLT in Morrison and Mildbrandt (2015)

Returns the total melted ice mass as dqdt; returns the melted ice mass
transferred directly to rain (D < D_th) as dqdt_r; and the melted ice mass which
accumulates as liquid on ice particles (D > D_th) as dqdt_i.
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
    F_liq::FT,
) where {FT}
    # Get constants
    (; ν_air, D_vapor, K_therm) = aps
    a = p3.vent_a
    b = p3.vent_b
    L_f = TD.latent_heat_fusion(tps, T)
    T_freeze = p3.T_freeze
    N_sc = ν_air / D_vapor
    ρ_l = p3.ρ_l

    # use non liquid mass
    # and feed F_liq = 0 to shape solver
    # to get ice core size distribution
    q = (1 - F_liq) * q

    # Get distribution values
    D_th = D_th_helper(p3)
    th = thresholds(p3, ρ_r, F_r)
    (λ, N_0) = distribution_parameter_solver(p3, q, N, ρ_r, FT(0), F_r) # F_liq = 0

    # Get bound 
    ice_bound = get_ice_bound(p3, λ, eps(FT))

    # check whether all mass should be transferred to rain
    # (i.e. whether "all" particles have D < D_th)
    small_d = false
    if (2 * ice_bound) < D_th
        small_d = true
    end

    # Define function pieces
    N_re(D) =           # TODO: What is this for non-spherical particles?
        D * mixed_phase_particle_velocity(
            p3,
            D,
            Chen2022,
            ρ_a,
            F_r,
            F_liq, # F_liq = nonzero because we want the real velocity
            th,
        ) / ν_air
    F(D) = a + b * N_sc^(1 / 3) * N_re(D)^(1 / 2)
    # F_liq = 0 for dmdt since we only care about the change in ice core mass
    dmdt(D) =
        dmdD_mass(p3, D, F_r, FT(0), th) / ρ_l * K_therm / L_f *
        (T - T_freeze) *
        F(D)
    f(D) = 1 / (2 * ρ_a) * dmdt(D) * N′ice(p3, D, λ, N_0)

    # Integrate
    if !small_d
        (dqdt_r, error) = QGK.quadgk(d -> f(d), FT(0), D_th)
        (dqdt_i, error) = QGK.quadgk(d -> f(d), D_th, 2 * ice_bound)
        dqdt = dqdt_i + dqdt_r
        return dqdt, dqdt_r, dqdt_i
    else
        (dqdt, error) = QGK.quadgk(d -> f(d), FT(0), 2 * ice_bound)
        return dqdt, dqdt, FT(0)
    end
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
