module Nucleation

"""
    apparent_nucleation_rate(
        output_diam,
        nucleation_rate,
        condensation_growth_rate,
        coag_sink,
        coag_sink_input_diam,
        input_diam,
    )

 - output_diam: The smallest diameter in the aerosol model (nm).
 - nucleation_rate: The formation rate of particles with diameter input_diam.
 - condensation_growth_rate: The growth rate by condensation.
 - coag_sink: The coagulation rate of particles with diameter less than output_diam.
 - coag_sink_input_diam: The coagulation rate of particles with diameter less than input_diam.
 - input_diam: The diameter of freshly nucleated particles (nm). Currently defaults to 1.7 nm.
Computes the apparent nucleation rate from Lehtinen et al., 2007 (doi: 10.1016/j.jaerosci.2007.06.009)
This is a rough approximation accounting for the effect of coagulation as newly formed particles
grow to the output diameter.
"""
function apparent_nucleation_rate(
    output_diam,
    nucleation_rate,
    condensation_growth_rate,
    coag_sink,
    coag_sink_input_diam,
    input_diam = 1.7,
)
    m = log(coag_sink / coag_sink_input_diam) / (log(output_diam / input_diam))
    γ = 1 / (m + 1) * ((output_diam / input_diam)^(m + 1) - 1)
    J_x =
        nucleation_rate *
        exp(-γ * input_diam * coag_sink_input_diam / condensation_growth_rate)
    return J_x
end

"""
    h2so4_nucleation_rate(h2so4_conc, nh3_conc, negative_ion_conc, temp)

 - `h2so4_conc` - Concentration of h2so4 (1/m³)
 - `nh3_conc` - Concentration of nh3 (1/m³)
 - `negative_ion_conc` - Concentration of negative ions (1/m³)
 - `temp` - Temperature (K)
 - `params` - NamedTuple parameter set obtained from CLIMAParameters. 
Calculates the rate of binary H2SO4-H2O and ternary H2SO4-H2O-NH3 nucleation for a single timestep (1/m³/s).
The particle formation rate is parameterized using data from the CLOUD experiment, through neutral and ion-induced channels.
This is an implementation of Dunne et al 1016 doi:10.1126/science.aaf2649 Appendix 8-10
"""
function h2so4_nucleation_rate(
    h2so4_conc,
    nh3_conc,
    negative_ion_conc,
    temp,
    params,
)
    (;
        u_b_n,
        v_b_n,
        w_b_n,
        u_b_i,
        v_b_i,
        w_b_i,
        u_t_n,
        v_t_n,
        w_t_n,
        u_t_i,
        v_t_i,
        w_t_i,
        p_t_n,
        p_A_n,
        a_n,
        p_t_i,
        p_A_i,
        a_i,
        p_b_n,
        p_b_i,
    ) = params
    # Change units from 1/m³ to 1/cm³
    h2so4_conc *= 1e-6
    nh3_conc *= 1e-6

    # Reference concentration for h2so4 and nh3 (Units: 1e6/cm³)
    ref_conc = 1e6

    # Calculate factors
    k(T, u, v, w) = exp(u - exp(v * (T / 1000 - w)))
    f_y(h2so4, nh3, p_t_y, p_A_y, a_y) =
        (nh3 / ref_conc) /
        (a_y + (h2so4 / ref_conc)^p_t_y / (nh3 / ref_conc)^p_A_y)
    k_b_n = k(temp, u_b_n, v_b_n, w_b_n)
    k_b_i = k(temp, u_b_i, v_b_i, w_b_i)
    k_t_n = k(temp, u_t_n, v_t_n, w_t_n)
    k_t_i = k(temp, u_t_i, v_t_i, w_t_i)
    f_n = f_y(h2so4_conc, nh3_conc, p_t_n, p_A_n, a_n)
    f_i = f_y(h2so4_conc, nh3_conc, p_t_i, p_A_i, a_i)

    # Calculate rates 1/cm³/s
    binary_rate =
        k_b_n * (h2so4_conc / ref_conc)^p_b_n +
        k_b_i * (h2so4_conc / ref_conc)^p_b_i * negative_ion_conc
    ternary_rate =
        k_t_n * f_n * (h2so4_conc / ref_conc)^p_t_n +
        k_t_i * f_i * (h2so4_conc / ref_conc)^p_t_i * negative_ion_conc
    # Convert to 1/m³/s
    return (;
        binary_rate = binary_rate * 1e6,
        ternary_rate = ternary_rate * 1e6,
    )
end

"""
    organic_nucleation_rate(negative_ion_conc, monoterpene_conc, O3_conc, OH_conc, temp, condensation_sink, params)

 - `negative_ion_conc` - Concentration of negative ions (1/m³)
 - `monoterpene_conc` - Concentration of monoterpenes (1/m³)
 - `O3_conc` - Concentration of O3 (1/m³)
 - `OH_conc` - Concentration of OH (1/m³)
 - `temp` - Temperature (K)
 - `condensation_sink` - Condensation sink (1/s)
 - `params`  - NamedTuple parameter set obtained from CLIMAParameters. 

Returns nucleation rate of pure biogenic particles (1/m³/s)
The parameterization is an implementation of Kirkby at al 2016 doi.org/10.1038/nature17953
"""
function organic_nucleation_rate(
    negative_ion_conc,
    monoterpene_conc,
    O3_conc,
    OH_conc,
    temp,
    condensation_sink,
    params,
)
    # Convert units from 1/m³ to 1/cm³ - assume units are meant to be in 1/cm³
    negative_ion_conc *= 1e-6
    monoterpene_conc *= 1e-6
    O3_conc *= 1e-6
    OH_conc *= 1e-6

    # Y_* params from Dunne et al. 2016
    (; Y_MTO3, Y_MTOH, k_MTO3, k_MTOH, exp_MTO3, exp_MTOH) = params
    k_MTO3 = k_MTO3 * exp(exp_MTO3 / temp)
    k_MTOH = k_MTOH * exp(exp_MTOH / temp)
    HOM_conc =
        (
            Y_MTO3 * k_MTO3 * monoterpene_conc * O3_conc +
            Y_MTOH * k_MTOH * monoterpene_conc * OH_conc
        ) / condensation_sink
    return organic_nucleation_rate(negative_ion_conc, HOM_conc, params)
end

function organic_nucleation_rate(negative_ion_conc, HOM_conc, params)
    (; a_1, a_2, a_3, a_4, a_5) = params
    # HOM reference concentration: 1e7/cm³
    ref_conc = 1e7
    rate =
        a_1 * (HOM_conc / ref_conc)^(a_2 + a_5 / (HOM_conc / ref_conc)) +
        a_3 *
        (HOM_conc / ref_conc)^(a_4 + a_5 / (HOM_conc / ref_conc)) *
        negative_ion_conc
    # Convert from (1/cm³/s) to (1/m³/s)
    return rate * 1e6
end

"""
    organic_and_h2so4_nucleation_rate(h2so4_conc monoterpene_conc, OH_conc, temp, condensation_sink, params)

- `h2so4_conc` - Concentration of sulfuric acid (1/m³)
- `monoterpene_conc` - Concentration of monoterpenes (1/m³)
- `OH_conc` - Concentration of OH (1/m³)
- `temp` - Temperature (K)
- `condensation_sink` - Condensation sink (1/s)
- `params` - NamedTuple parameter set obtained from CLIMAParameters. 

Returns nucleation rate of sulfuric acid and oxidized organic molecules (1/m³/s).
This is an implementation of Riccobono et al 2014 (https://doi.org/10.1126/science.1243527)
"""
function organic_and_h2so4_nucleation_rate(
    h2so4_conc,
    monoterpene_conc,
    OH_conc,
    temp,
    condensation_sink,
    params,
)
    (; k_MTOH, exp_MTOH) = params
    k_MTOH = k_MTOH * exp(exp_MTOH / temp) # Units: 1/cm³/s
    bioOxOrg = k_MTOH * monoterpene_conc * OH_conc / condensation_sink
    # Convert from 1/cm³ to 1/m³
    bioOxOrg *= 1e6
    return organic_and_h2so4_nucleation_rate(h2so4_conc, bioOxOrg, params)
end

function organic_and_h2so4_nucleation_rate(h2so4_conc, bioOxOrg, params)
    (; k_H2SO4org) = params
    # Convert bioOxOrg and k_H2SO4org from 1/m³ to 1/cm³
    k_H2SO4org *= 1e-6
    bioOxOrg *= 1e-6
    # Convert from 1/m³ to 10⁶/cm³
    # h2so4_conc *= 1e-6
    rate = 0.5 * k_H2SO4org * h2so4_conc^2 * bioOxOrg
    # Convert from 1/cm³/s to 1/m³/s
    return rate * 1e6
end

end
