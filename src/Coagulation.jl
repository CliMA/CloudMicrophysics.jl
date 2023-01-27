module Coagulation

using Cubature
using CLIMAParameters

include("CoagCorrectionFactors.jl")
import .CoagCorrectionFactors as CCF



# Returns a log-normal distribution for the given aerosol mode.
function lognormal_dist(am)
    # Transform geometric mean and stdev for lognormal distribution
    mu = log(2 * am.r_dry)
    sigma = log(am.stdev)
    return x ->
        1 / (x * sigma * sqrt(2 * pi)) *
        exp(-(log(x) - mu)^2 / (2 * sigma^2))
end

"""
    quadrature(ad, air_pressure, temp, params)

 - `ad`: AerosolDistribution, a tuple of aerosol modes. Currently only MAM3 is supported.
 - `particle_density_ait`: Particle density for the aitken mode (kg/m^3)
 - `particle_density_acc`: Particle density for the accumulation mode (kg/m^3)
 - `air_pressure`: Ambient air pressure (pa)
 - `temp`: Ambient air temperature (K)
 - `params`: NamedTuple of relevant parameters
Calculates the rates of change to the 0th, 2nd, and 3rd moments of the Aitken
and Accumulation modes due to inter- and intramodal coagulation. This is done
Gaussian quadrature over the log-normal distributions using Cubature.jl.
"""
function quadrature(
    ad,
    particle_density_ait,
    particle_density_acc,
    air_pressure,
    temp,
    params,
)
    # Set up lognormal distributions
    aitken = ad.Modes[1]
    accum = ad.Modes[2]
    aitken_distribution = lognormal_dist(aitken)
    accumulation_distribution = lognormal_dist(accum)
    k_fm = sqrt(
        6 * params.k_Boltzmann * temp /
        (particle_density_ait + particle_density_acc)
    )
    # u.s. standard atmosphere 1962 page 14 expression for dynamic viscosity is:
    #   1.458e-6 * t * sqrt(t) / ( t + 110.4)
    gas_viscosity = 1.458e-6 * temp^1.5 / ( temp + 110.4 )
    k_nc = 2 * params.k_Boltzmann * temp / (3 * gas_viscosity)
    # Coag coefficient
    p0 = params.MSLP        #  standard surface pressure (pa)
    t0 = params.T_surf_ref  #  standard surface temperature (K)
    # From CAM: 
    # 6.6328e-8 is the sea level value given in table i.2.8
    # on page 10 of u.s. standard atmosphere 1962
    lambda = 6.6328e-8 * p0 * temp / (t0 * air_pressure)
    # Binkowski and Shankar, 95 - a5
    beta_fm(dp1, dp2) = k_fm * sqrt(1 / dp1^3 + 1 / dp2^3) * (dp1 + dp2)^2
    # Binkowski and Shankar, 95 - a6
    beta_nc(dp1, dp2) =
        k_nc *
        (dp1 + dp2) *
        (1 / dp1 + 1 / dp2 + 2.492 * lambda * (1 / dp1^2 + 1 / dp2^2))

    # Calculate quadrature
    # Intramodal
    intramodal_rates = intracoag_quadrature(beta_fm, beta_nc, aitken_distribution, accumulation_distribution)
    # Intermodal
    intermodal_rates = intercoag_quadrature(beta_fm, beta_nc, aitken_distribution, accumulation_distribution)
    return (intramodal_rates, intermodal_rates)
end

# Wrapper function for computing Gaussian quadrature. Takes an integrand formatted for Cubature.jl
function cubature(integrand)
    start = eps()
    stop = 1e2
    (result, err) = Cubature.hcubature(integrand, [start, start], [stop, stop])
    return result
end

"""
    intracoag_quadrature(beta_fm, beta_nc, distribution)
 - `beta_fm`: Coagulation coefficient for free-molecule Knudsen regime
 - `beta_nc`: Coagulation coefficient for near-continuum Knudsen regime
 - `distribution`: A function for the given lognormal aerosol distribution
Helper function for `coagulation_quadrature`, calculates 
intramodal coagulation rates for the given modal distribution.
"""
function intracoag_quadrature(beta_fm, beta_nc, aitken_dist, accum_dist)
    integrand1(moment, beta, distribution) =
        dp ->
            dp[1]^moment *
            dp[2]^moment *
            beta(dp[1], dp[2]) *
            distribution(dp[1]) *
            distribution(dp[2])
    integrand2(moment, beta, distribution) =
        dp ->
            (dp[1]^3 + dp[2]^3)^(moment/3) *
            beta(dp[1], dp[2]) *
            distribution(dp[1]) *
            distribution(dp[2])
    aitken_rates = Vector{Float64}(undef, 3)
    # Aitken
    for moment in 0:2
        rate_fm = 
            0.5 * cubature(integrand1(moment, beta_fm, aitken_dist))
            - cubature(integrand2(moment, beta_fm, aitken_dist))
        rate_nc = 
            0.5 * cubature(integrand1(moment, beta_nc, aitken_dist))
            - cubature(integrand2(moment, beta_nc, aitken_dist))
        aitken_rates[moment+1] = rate_fm * rate_nc / (rate_fm + rate_nc + eps())
    end
    # Accumulation
    accum_rates = Vector{Float64}(undef, 3)
    for moment in 0:2
        rate_fm = 
            0.5 * cubature(integrand1(moment, beta_fm, accum_dist))
            - cubature(integrand2(moment, beta_fm, accum_dist))
        rate_nc = 
            0.5 * cubature(integrand1(moment, beta_nc, accum_dist))
            - cubature(integrand2(moment, beta_nc, accum_dist))
        accum_rates[moment+1] = rate_fm * rate_nc / (rate_fm + rate_nc + eps())
    end
    return (aitken_rates, accum_rates)
end

"""
    intracoag_quadrature(beta_fm, beta_nc, distribution)
 - `beta_fm`: Coagulation coefficient for free-molecule Knudsen regime
 - `beta_nc`: Coagulation coefficient for near-continuum Knudsen regime
 - `aitken_dist`: A function for the aitken mode log-normal distribution
 - `accum_dist`: A function for the accumulation mode log-normal distribution
Helper function for `coagulation_quadrature`, calculates 
intermodal coagulation rates between the aitken and accumulation modes.
"""
function intercoag_quadrature(beta_fm, beta_nc, aitken_dist, accum_dist)
    integrand1(moment, beta, dist1, dist2) =
        dp ->
            dp[1]^moment *
            beta(dp[1], dp[2]) *
            dist1(dp[1]) *
            dist2(dp[2])
    
    integrand2(moment, beta, dist1, dist2) =
        dp ->
            (dp[1]^3 + dp[2]^3)^(moment/3) *
            beta(dp[1], dp[2]) *
            dist1(dp[1]) *
            dist2(dp[2])

    # Aitken
    aitken_rates = Vector{Float64}(undef, 4)
    for moment in 0:3
        rate_fm = -cubature(integrand1(moment, beta_fm, aitken_dist, accum_dist))
        rate_nc = -cubature(integrand1(moment, beta_nc, aitken_dist, accum_dist))
        aitken_rates[moment+1] = rate_fm * rate_nc / (rate_fm + rate_nc + eps())
    end
    # Accumulation:
    accum_rates = Vector{Float64}(undef, 4)
    for moment in 0:3
        rate_fm = 
            cubature(integrand2(moment, beta_fm, accum_dist, aitken_dist))
            -cubature(integrand1(moment, beta_fm, accum_dist, aitken_dist))
        rate_nc = 
        cubature(integrand2(moment, beta_nc, accum_dist, aitken_dist))
            -cubature(integrand1(moment, beta_nc, accum_dist, aitken_dist))
        accum_rates[moment+1] = rate_fm * rate_nc / (rate_fm + rate_nc + eps())
    end
    return (aitken_rates, accum_rates)
end

"""
    whitby_coagulation(ad, air_pressure, temp, parameter_file)

    - `ad`: AerosolDistribution, a tuple of aerosol modes. Currently only MAM3 is supported.
    - `particle_density_ait`: Particle density for the aitken mode (kg/m^3)
    - `particle_density_acc`: Particle density for the accumulation mode (kg/m^3)
    - `air_pressure`: Ambient air pressure (pa)
    - `temp`: Ambient air temperature (K)
    - `parameter_file`: TOML file to read in parameters. This may be changed once the coagulation rate is actually applied to the aerosol data struct.
Calculates and applies changes from coagulation to the aerosol modes.
This method follows CAM5 (10.5194/gmd-5-709-2012), which is largely based off of 
Whitby et al., 1991. 
"""
function whitby_coagulation(
    ad,
    air_pressure,
    temp,
    particle_density_ait,
    particle_density_acc,
    params
)
    surface_pressure = params.MSLP        #  standard surface pressure (pa)
    surface_temp = params.T_surf_ref  #  standard surface temperature (K)
    k_Boltzmann = params.k_Boltzmann

    lambda =
        6.6328e-8 * surface_pressure * temp / (surface_temp * air_pressure)
    gas_viscosity = 1.458e-6 * temp^1.5 / ( temp + 110.4 )
    aitken = ad.Modes[1]
    accum = ad.Modes[2]
    Kn_ait = lambda / aitken.r_dry
    Kn_acc = lambda / accum.r_dry
    ESG_ait = exp(1 / 8 * log(aitken.stdev)^2)
    ESG_acc = exp(1 / 8 * log(accum.stdev)^2)
    # Issue: CAM5 and Whitby don't match for K_nc and K_fm
    # Currently implementing Whitby, should be cam
    K_fm_ait = sqrt(3 * k_Boltzmann * temp / particle_density_ait)
    K_fm_acc = sqrt(3 * k_Boltzmann * temp / particle_density_acc)
    K_fm_both =
        sqrt(6 * k_Boltzmann * temp / (particle_density_ait + particle_density_acc))
    K_nc = sqrt(2 * k_Boltzmann * temp / (3 * gas_viscosity))

    # Call coags for aitken and accumulation modes:
    (aitken_m0_intracoag, aitken_m2_intracoag) =
        whitby_intramodal_coag(aitken, K_fm_ait, K_nc, Kn_ait, ESG_ait)
    (accum_m0_intracoag, accum_m2_intracoag) =
        whitby_intramodal_coag(accum, K_fm_acc, K_nc, Kn_acc, ESG_acc)
    (aitken_m0_intercoag, aitken_m2_intercoag, accum_m2_intercoag, aitken_m3_intercoag) = whitby_intermodal_coag(
        ad,
        K_fm_both,
        K_nc,
        Kn_ait,
        Kn_acc,
        ESG_ait,
        ESG_acc,
    )

    # aitken.N += intra_m_0_ait + inter_m_0_ait
    # accum.N += intra_m_0_acc
    return (
        aitken_m0_intracoag,
        accum_m0_intracoag,
        aitken_m2_intracoag,
        accum_m2_intracoag,
        aitken_m0_intercoag,
        aitken_m2_intercoag,
        accum_m2_intercoag,
        aitken_m3_intercoag,
        -aitken_m3_intercoag
    )
end


"""
    whitby_intramodal_coag(ad, K_fm, K_nc, Kn_g, ESG)

 - `am` - AerosolMode
 - `K_fm` - Non-size dependent term for free-molecule coagulation coefficient (m^2.5/s)
 - `K_nc` - Non-size dependent term for continuum/near-continuum coagulation coefficient (m^3/s)
 - `Kn_g` - Knudsen number for the aitken mode
 - `ESG` - exp(1/8 * ln^2(aitken.stdev))
Calculates the intramodal coagulation rates for the aitken and accumulation modes, as specified in Whitby et al., 1991, Appendix H.
Larger modes are not treated, as specified in Liu et al., 2012.
The equations below are analytical expressions for coagulation integrals for the 0-th and 6-th moments,
as found in Whitby et al., 1991.
For more information on how these expressions are produced, see the documentation page.
"""
function whitby_intramodal_coag(am, K_fm, K_nc, Kn_g, ESG)
    A = 1.43 * Kn_g^0.0814
    # Fixed in CAM5, fix it here
    A = 1.246
    sqrt_diam = sqrt(2 * am.r_dry)

    # Whitby H.11a 
    m_0_fm =
        -am.N^2 *
        K_fm *
        CCF.intramodal_m0_correction(am.stdev) *
        sqrt_diam *
        (ESG + ESG^25 + 2 * ESG^5)
    # Whitby H.12a
    m_0_nc = -am.N^2 * K_nc * (1 + ESG^8 + A * Kn_g * (ESG^20 + ESG^4))

    m_2_fm =
        am.N^2 *
        K_fm *
        sqrt_diam^5CCF.intramodal_m2_correction(am.stdev) *
        (ESG^25 + 2 * ESG^13 + ESG^17 + ESG^73 + 2 * ESG^37 + ESG^17)

    m_2_nc =
        am.N^2 *
        K_nc *
        (2 * am.r_dry)^2 *
        (2 * ESG^16 + ESG^8 + ESG^40 + A * Kn_g * (2 * ESG^4 + ESG^20 + ESG^52))

    # Harmonic mean
    m_0 = m_0_fm * m_0_nc / (m_0_fm + m_0_nc)
    m_2 = m_2_fm * m_2_nc / (m_2_fm + m_2_nc)
    return m_0, m_2
end

"""
    whitby_intermodal_coag(ad, K_fm, K_nc, Kn_ait, Kn_acc, ESG_ait, ESG_acc)

 - `ad` - AerosolDistribution, a tuple of aerosol modes
 - `K_fm` - Non-size dependent term for free-molecule coagulation coefficient (m^2.5/s)
 - `K_nc` - Non-size dependent term for continuum/near-continuum coagulation coefficient (m^3/s)
 - `Kn_ait` - Knudsen number for the aitken mode
 - `Kn_acc` - Knudsen number for the accumulation mode
 - `ESG_ait` - exp(1/8 * ln^2(aitken.stdev))
 - `ESG_acc` - exp(1/8 * ln^2(accum.stdev))
Calculates the intermodal coagulation rates for the aitken and accumulation modes, as specified in Whitby et al., 1991, Appendix H
The equations below are analytical expressions for coagulation integrals for the 0-th and 6-th moments,
as found in Whitby et al., 1991.
"""
function whitby_intermodal_coag(
    ad,
    K_fm,
    K_nc,
    Kn_ait,
    Kn_acc,
    ESG_ait,
    ESG_acc,
)
    aitken = ad.Modes[1]
    accum = ad.Modes[2]
    sqrt_diam_ait = sqrt(2 * aitken.r_dry)
    sqrt_diam_acc = sqrt(2 * accum.r_dry)
    R = sqrt_diam_acc / sqrt_diam_ait
    R_2 = R^2
    R_3 = R^3
    R_4 = R^4
    R_6 = R^6
    # Issue: a is hardcoded to one value in CAM5? 
    A = 1.246
    A_ait = A
    A_acc = A
    # Free molecular form:
    # Whitby H.7a
    m_0_ait_fm =
        -aitken.N *
        accum.N *
        K_fm *
        CCF.intermodal_m0_correction(R_2, aitken.stdev, accum.stdev) *
        sqrt_diam_ait *
        (
            ESG_ait +
            R * ESG_acc +
            2 * R_2 * ESG_ait * ESG_acc^4 +
            R^4 * ESG_ait^9 * ESG_acc^16 +
            (1 / R_3) * ESG_ait^16 * ESG_acc^9 +
            2 * (1 / R) * ESG_ait^4 * ESG_acc
        )
    # From CAM5 code, Whitby doesn't specify 2nd moment integral approximation
    # Not sure about units - should this be multiplied by number concentration?
    m_2_ait_fm =
        aitken.N *
        accum.N *
        K_fm *
        sqrt_diam_ait^5 *
        CCF.intermodal_m2_ait_correction(R_2, aitken.stdev, accum.stdev) *
        (
            ESG_ait^25 +
            2 * R_2 * ESG_ait^9 * ESG_acc^04 +
            R_4 * ESG_ait * ESG_acc^16 +
            (1 / R_3) * ESG_ait^64 * ESG_acc^9 +
            2 * (1 / R) * ESG_ait^36 * ESG_acc +
            R * ESG_ait^16 * ESG_acc
        )
    # Whitby H.7b
    m_3_fm =
        -aitken.N *
        accum.N *
        K_fm *
        CCF.intermodal_m3_correction(R_2, aitken.stdev, accum.stdev) *
        sqrt_diam_ait^7 *
        (
            ESG_ait^49 +
            R * ESG_ait^36 * ESG_acc +
            2 * R_2 * ESG_ait^25 * ESG_acc^4 +
            R_4 * ESG_ait^9 * ESG_acc^16 +
            (1 / R_3) * ESG_ait^100 * ESG_acc^9 +
            2 * (1 / R) * ESG_ait^4 * ESG_acc
        )
    # Near-continuum:
    # Whitby H.10a
    m_0_ait_nc =
        -aitken.N * accum.N * K_nc +
        A_ait * Kn_ait * (ESG_ait^4 + R_2 * ESG_ait^16 * ESG_acc^4) +
        A_acc +
        Kn_acc * (ESG_acc^4 + (1 / R_2) * ESG_acc^16 * ESG_ait^4) +
        (R_2 + (1 / R_2) * ESG_ait^4 * ESG_acc^4)
    # Unsure about A_ait or A_acc for this case
    m_2_ait_nc =
        K_nc *
        (2 * aitken.r_dry)^2 *
        (
            2 * ESG_ait^16 +
            R_2 * ESG_ait^4 * ESG_acc^4 +
            A_ait *
            Kn_ait *
            (
                ESG_ait^4 +
                (1 / R_2) * ESG_ait^16 * ESG_acc^4 +
                (1 / R_4) * ESG_ait^36 * ESG_acc^16 +
                R_2 * ESG_acc^4
            )
        )
    # Whitby H.10b
    m_3_nc =
        -aitken.N *
        accum.N *
        K_nc *
        aitken.r_dry^3 *
        (
            2 * ESG_ait^36 +
            A_ait * Kn_ait * (ESG_ait^16 + R_2 * ESG_ait^4 * ESG_acc^4) +
            A_acc *
            Kn_acc *
            (ESG_ait^36 * ESG_acc^4 + (1 / R_2) * ESG_ait^64 * ESG_acc^16) +
            R_2 * ESG_ait^16 * ESG_acc^4 +
            (1 / R_2) * ESG_ait^64 * ESG_acc^4
        )
    # Harmonic mean
    m_0_ait = m_0_ait_fm * m_0_ait_nc / (m_0_ait_fm + m_0_ait_nc)
    m_2_ait = m_2_ait_fm * m_2_ait_nc / (m_2_ait_fm + m_2_ait_nc)
    m_3 = m_3_fm * m_3_nc / (m_3_fm + m_3_nc)
    # From CAM5 :  coagacat2 = ( ( one + r6 ) ** two3rds - rx4 ) * i1
    m_2_acc =
        (((1 + R_6)^(2 / 3) - R_4) * m_2_ait) *
        CCF.intermodal_m2_acc_correction(R_2, aitken.stdev, accum.stdev)
    return (m_0_ait, m_2_ait, m_2_acc, m_3)
end

end

