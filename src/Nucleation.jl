module Nucleation

using CLIMAParameters

"""
    cloud_h2so4_nucleation_rate(h2so4_conc, nh3_conc, negative_ion_conc, temp)

 - `h2so4_conc` - Concentration of h2so4 (1/cm3) - This will be changed to (1/m3)
 - `nh3_conc` - Concentration of nh3 (1/cm3)
 - `negative_ion_conc` - Concentration of negative ions (1/cm3)
 - `temp` - Temperature (K)
Calculates  the rate of binary H2SO4-H2O and ternary H2SO4-H2O-NH3 nucleation for a single timestep.
The particle formation rate is parameterized using data from the CLOUD experiment, through neutral and ion-induced channels.
For more background, see Nucleation documentation page or doi:10.1126/science.aaf2649 Appendix 8-10
"""
function cloud_h2so4_nucleation_rate(
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
    params,
)
    k(T, u, v, w) = exp(u - exp(v * (T / 1000 - w)))
    f_y(h2so4, nh3, p_t_y, p_A_y, a_y) = nh3 / (a_y + h2so4^p_t_y / nh3^p_A_y)
    k_b_n = k(temp, params.u_b_n, params.v_b_n, params.w_b_n)
    k_b_i = k(temp, params.u_b_i, params.v_b_i, params.w_b_i)
    k_t_n = k(temp, params.u_t_n, params.v_t_n, params.w_t_n)
    k_t_i = k(temp, params.u_t_i, params.v_t_i, params.w_t_i)
    f_n = f_y(h2so4_conc, nh3_conc, params.p_t_n, params.p_A_n, params.a_n)
    f_i = f_y(h2so4_conc, nh3_conc, params.p_t_i, params.p_A_i, params.a_i)
    # Calculate rates
    binary_rate =
        k_b_n * h2so4_conc^params.p_b_n +
        k_b_i * h2so4_conc^params.p_b_i * negative_ion_conc
    ternary_rate =
        k_t_n * f_n * h2so4_conc^params.p_t_n +
        k_t_i * f_i * h2so4_conc^params.p_t_i * negative_ion_conc
    return binary_rate, ternary_rate
end

"""
    cloud_organic_nucleation_rate(negative_ion_conc, monoterpene_conc, O3_conc, OH_conc, temp, condensation_sink, params)

 - `negative_ion_conc` - Concentration of negative ions (1/cm3)
 - `monoterpene_conc` - Concentration of monoterpenes (1/cm3)
 - `O3_conc` - Concentration of O3 (1/cm3)
 - `OH_conc` - Concentration of OH (1/cm3)
 - `temp` - Temperature (K)
 - `condensation_sink` - Condensation sink (1/s?)
 - `params` - Parameters from the TOML file
"""
function cloud_organic_nucleation_rate(
    negative_ion_conc,
    monoterpene_conc,
    O3_conc,
    OH_conc,
    temp,
    condensation_sink,
    params,
)
    # Y_ params from Dunne et al. 2016
    a_1 = params.a_1
    a_2 = params.a_2
    a_3 = params.a_3
    a_4 = params.a_4
    a_5 = params.a_5
    Y_MTO3 = params.Y_MTO3
    Y_MTOH = params.Y_MTOH
    k_MTO3 = params.k_MTO3 * exp(444 / temp)
    k_MTOH = params.k_MTOH * exp(444 / temp)
    HOM_conc =
        (
            Y_MTO3 * k_MTO3 * monoterpene_conc * O3_conc +
            Y_MTOH * k_MTOH * monoterpene_conc * OH_conc
        ) / condensation_sink
    rate =
        a_1 +
        HOM_conc^(a_2 + a_5 / HOM_conc) +
        a_3 * HOM_conc^(a_4 + a_5 / HOM_conc) * negative_ion_conc
    return rate
end

"""
    cloud_organic_and_h2so4_nucleation_rate(h2so4_conc monoterpene_conc, OH_conc, temp, condensation_sink, params)

- `h2so4_conc` - Concentration of sulfuric acid (1/cm3)
- `monoterpene_conc` - Concentration of monoterpenes (1/cm3)
- `OH_conc` - Concentration of OH (1/cm3)
- `temp` - Temperature (K)
- `condensation_sink` - Condensation sink (1/s?)
- `params` - Parameters from the TOML file
"""
function cloud_organic_and_h2so4_nucleation_rate(
    h2so4_conc,
    monoterpene_conc,
    OH_conc,
    temp,
    condensation_sink,
    params,
)
    k_H2SO4org = params.k_H2SO4org
    k_MTOH = params.k_MTOH * exp(444 / temp)
    bioOxOrg = k_MTOH * monoterpene_conc * OH_conc / condensation_sink
    rate = 0.5 * k_H2SO4org * h2so4_conc^2 * bioOxOrg
    return rate
end

"""
    vehkamaki_nucleation_timestep(rh, temp, so4)

 - `rh` - Relative humidity, expressed as a percentage between 0 and 1
 - `temp` - Temperature (K)
 - `so4` - Total gas phase concentration of sulfuric acid (m⁻³)

Calculates the rate of binary H2SO4-H2O nucleation for a single timestep under the given conditions.
The parameterization is valid for the temperature range 230.15–300.15 K, 
relative humidities 0.01–100% and total sulfuric acid concentrations 1e10-1e17 m⁻³.

All newly formed particles will be placed into the Aitken mode.
The initial nucleation rate comes from Vehkamaki et al, 2002 doi:10.1029/2002JD002184
Adjustment factor comes from Kerminen and Kulmala, 2002 doi:10.1029/2002JD002184
"""
function vehkamaki_nucleation_timestep(rh, temp, so4)
    so4 = so4 / 1e6 # convert to cm⁻³
    ln_so4 = log(so4)
    log_rh = log(rh)
    x_crit = (
        0.740997 - 0.00266379 * temp - 0.00349998 * ln_so4 +
        0.0000504022 * temp * ln_so4 +
        0.00201048 * log_rh - 0.000183289 * temp * log_rh +
        0.00157407 * log_rh^2 - 0.0000179059 * temp * log_rh^2 +
        0.000184403 * log_rh^3 - 1.50345e-6 * temp * log_rh^3
    )

    a_T = (
        0.14309 + 2.21956 * temp - 0.0273911 * temp^2 +
        0.0000722811 * temp^3 +
        5.91822 / x_crit
    )
    b_T = (
        0.117489 + 0.462532 * temp - 0.0118059 * temp^2 +
        0.0000404196 * temp^3 +
        15.7963 / x_crit
    )
    c_T = (
        -0.215554 - 0.0810269 * temp + 0.00143581 * temp^2 -
        4.7758e-6 * temp^3 - 2.91297 / x_crit
    )
    d_T = (
        -3.58856 + 0.049508 * temp - 0.00021382 * temp^2 + 3.10801e-7 * temp^3 - 0.0293333 / x_crit
    )
    e_T = (
        1.14598 - 0.600796 * temp + 0.00864245 * temp^2 -
        0.0000228947 * temp^3 - 8.44985 / x_crit
    )
    f_T = (
        2.15855 + 0.0808121 * temp - 0.000407382 * temp^2 -
        4.01957e-7 * temp^3 + 0.721326 / x_crit
    )
    g_T = (
        1.6241 - 0.0160106 * temp +
        0.0000377124 * temp^2 +
        3.21794e-8 * temp^3 - 0.0113255 / x_crit
    )
    h_T = (
        9.71682 - 0.115048 * temp +
        0.000157098 * temp^2 +
        4.00914e-7 * temp^3 +
        0.71186 / x_crit
    )
    i_T = (
        -1.05611 + 0.00903378 * temp - 0.0000198417 * temp^2 +
        2.46048e-8 * temp^3 - 0.0579087 / x_crit
    )
    j_T = (
        -0.148712 + 0.00283508 * temp - 9.24619e-6 * temp^2 +
        5.00427e-9 * temp^3 - 0.0127081 / x_crit
    )

    # Nucleation rate 1/(m^3 s)
    J =
        exp(
            a_T +
            b_T * log_rh +
            c_T * log_rh^2 +
            d_T * log_rh^3 +
            e_T * ln_so4 +
            f_T * log_rh * ln_so4 +
            g_T * log_rh^2 * ln_so4 +
            h_T * ln_so4^2 +
            i_T * log_rh * ln_so4^2 +
            j_T * ln_so4^3,
        ) * 1e6
    return J
end

end
