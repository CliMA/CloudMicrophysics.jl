module Coagulation

import ..AerosolModel as AM
include("CoagCorrectionFactors.jl")
import .CoagCorrectionFactors as CCF

"""
    whitby_coagulation(ad, temp, particle_density, gas_viscosity, K_b)

 - `ad` - AerosolDistribution, a tuple of aerosol modes. Currently only MAM3 is supported.
 - `temp` - Air temperature (K)
 - `particle_density` - Particle density (kg/m^3)
 - `gas_viscosity` - Gas viscosity (kg/(m s))
 - `K_b` - Boltzmann constant. TODO: Obtain this from climaparameters
Calculates and applies changes from coagulation to the aerosol modes.
This method follows CAM5 (10.5194/gmd-5-709-2012), which is largely based off of 
Whitby et al., 1991. 
Current notable differences from CAM5:
 - K_nc and K_fm are calculated from Whitby, not CAM
 - Need to add 2nd moment calculations, but methods aren't backed up by Whitby or Binkowski
"""
function whitby_coagulation(
        ad,
        temp,
        particle_density,
        gas_viscosity,
        K_b
    )
    # TODO: add these values to climaparameters - this comes from CAM5
    surface_temp = 288.15
    surface_pressure = 101325.0
    mean_free_path = 6.6328e-8 * surface_pressure * temp / (surface_temp * pressure)
    aitken = ad.Modes[1]
    accum = ad.Modes[2]
    Kn_ait = mean_free_path /  aitken.r_dry
    Kn_acc = mean_free_path /  accum.r_dry
    sqrt_diam_ait = sqrt(atiken.r_dry)
    sqrt_diam_acc = sqrt(accum.r_dry)
    ESG_ait = exp(1 / 8 * log(aitken.stdev)^2)
    ESG_acc = exp(1 / 8 * log(accum.stdev)^2)
    # Issue: CAM5 and Whitby don't match for K_nc and K_fm
    # Currently implementing Whitby, should be cam
    K_fm = sqrt(3 * K_b * T / particle_density)
    K_nc = sqrt(2 * K_b * T / (3 * gas_viscosity))
    
    # Call coags for aitken and accumulation modes:
    intra_m_0_ait = whitby_intramodal_coag(aitken, K_fm, K_nc, Kn_ait, ESG_ait)
    intra_m_0_acc = whitby_intramodal_coag(accum, K_fm, K_nc, Kn_acc, ESG_acc)
    (inter_m_0_ait, m_3) = whitby_intermodal_coag(ad, K_fm, K_nc, Kn_ait, Kn_acc, ESG_ait, ESG_acc)

    aitken.N += intra_m_0_ait + inter_m_0_ait
    accum.N += intra_m_0_acc


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
function whitby_intramodal_coag(
    am,
    K_fm,
    K_nc,
    Kn_g,
    ESG
    )
    A_nc = 1.43 * Kn_g^0.0814

    # Whitby H.11a 
    m_0_fm = -am.N^2 * K_fm * CCF.intramodal_correction(am.stdev) * sqrt_diameter * (ESG + ESG^25 + 2 * ESG^5)
    # Whitby H.11b
    m_6_fm = -2 * am.N^2 * K_fm * CCF.intramodal_correction(am.stdev) * sqrt_diameter^6 * (ESG^85 + 2*ESG^89 + ESG^109)
    # Whitby H.12a
    m_0_nc = -am.N^2 * K_nc * (1 + ESG^8 + A_nc * Kn_g * (ESG^20 + ESG^4))
    #Whitby H.12b
    m_6_nc = -2 * am.N^2 * K_nc * (2*am.r_dry)^6 * (ESG^72 + ESG^80 + A_nc * Kn_g * (ESG^52 + ESG^68))

    return m_0_fm * m_0_nc / (m_0_fm + m_0_nc)
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
    ESG_acc
    )
    aitken = ad.Modes[1]
    accum = ad.Modes[2]

    R = sqrt(2*accum.r_dry) / sqrt(2*aitken.r_dry)
    R_2 = R^2
    R_3 = R^3
    R_4 = R^4
    A_ait = 1.43 * Kn_ait^0.0814
    A_acc = 1.43 * Kn_acc^0.0814
    # Free molecular form:
    # Whitby H.7a
    m_0_ait_fm = 
        - aitken.N * accum.N * K_fm * CCF.intermodal_m0_correction(R_2, aitken.stdev, accum.stdev) * sqrt_diam_ait * 
        (ESG_ait + R * ESG_acc + 2 * R_2 * ESG_ait * ESG_acc^4 +
         R^4 * ESG_ait^9 * ESG_acc^16 + 
         (1/R_3) * ESG_ait^16 * ESG_acc^9 +
         2 * (1/R) * ESG_ait^4 * ESG_acc)
    # Whitby H.7b
    m_3_fm = 
        - aitken.N * accum.N * K_fm * CCF.intermodal_m3_correction(R_2, aitken.stdev, accum.stdev) * sqrt_diam_ait^7 *
        (ESG_ait^49 + R * ESG^36 * ESG_acc + 2 * R_2 * ESG_ait^25 * ESG_acc^4 +
         R_4 * ESG_ait^9 * ESG_acc^16 + (1/R_3) * ESG_ait^100 * ESG_acc^9 +
         2 * (1/R) * ESG_ait^4 * ESG_acc)
    # Whitby H.7c
    # m_6_ait_fm = 
    #     - aitken.N * accum.N * K_fm * b_6_ait * sqrt_diam_ait^6 * sqrt_diam_acc^7 *
    #     (ESG_ait^169 + R * ESG_ait^144 * ESG_acc + 2 * R_2 * ESG_ait^121 * ESG_acc^4 +
    #      R_4 * ESG_ait^81 * ESG_acc^16 + (1/R_3) * ESG_ait^256 * ESG_acc^9 +
    #      2 * (1/R) * ESG_ait^196 * ESG_acc)
    # # Whitby H.7d
    # m_6_acc_fm = 
    #     aitken.N * accum.N * K_fm * b_6_acc * sqrt_diam_ait^7 * sqrt_diam_acc^6 *
    #     (ESG_ait^49 * ESG_acc^36 + R * ESG^36 * ESG_acc^49 + 2 * R_2 * ESG_ait^25 * ESG_acc^64 +
    #     R_4 * ESG_ait^9 * ESG_acc^100 + (1/R_3) * ESG_ait^100 * ESG_acc^9 +
    #     2 * (1/R) * ESG_ait^64 * ESG_acc^25)
    # Near-continuum:
    # Whitby H.10a
    m_0_ait_nc =
        - aitken.N * accum.N * K_nc +
        A_ait * Kn_ait * (ESG_ait^4 + R_2 * ESG_ait^16 * ESG_acc^4) +
        A_acc + Kn_acc * (ESG_acc^4 + (1/R_2) * ESG_acc^16 * ESG_ait^4) +
        (R_2 + (1/R_2) * ESG_ait^4 * ESG_acc^4)
    # Whitby H.10b
    m_3_nc =
        - aitken.N * accum.N * K_nc * aitken.r_dry^3 *
        (2 * ESG_ait^36 + A_ait * Kn_ait * (ESG_ait^16 + R_2 * ESG_ait^4 * ESG_acc^4) +
        A_acc * Kn_acc * (ESG_ait^36 * ESG_acc^4 + (1/R_2) * ESG_ait^64 * ESG_acc^16) +
        R_2 * ESG_ait^16 * ESG_acc^4 + (1/R_2) * ESG_ait^64 * ESG_acc^4)
    # Whitby H.10c
    m_6_ait_nc =
        - aitken.N * accum.N * K_nc * aitken.r_dry^6 *
        (2 * ESG_ait^144 + A_ait * Kn_ait * (ESG_ait^100 + R_2 * ESG_ait^64 * ESG_acc^4) +
        A_acc * Kn_acc * (ESG_ait^144 * ESG_acc^4 + (1/R_2) * ESG_ait^196 * ESG_acc^16) +
        R_2 * ESG_ait^100 * ESG_acc^4 + (1/R_2) * ESG_ait^196 * ESG_acc^4)        
    # Whitby H.10d
    m_6_acc_nc =
        2 * aitken.N * accum.N * K_nc * aitken.r_dry^3 * accum.r_dry^3 *
        (2 * ESG_ait^36 * ESG_acc^36 +
        A_ait * Kn_ait * (ESG_ait^16 * ESG_acc^36 + R_2 *ESG_ait^4 * ESG_acc^64) +
        A_acc * Kn_acc * (ESG_ait^36 * ESG_acc^16 + (1/R_2) * ESG_ait^64 * ESG_acc^4) +
        R_2 * ESG_ait^16 * ESG_acc^64 + (1/R_2) * ESG_ait^64 * ESG_acc^16)
    # Harmonic mean
    m_3 = m_3_fm * m_3_nc / (m_3_fm + m_3_nc)
    m_0_ait = m_0_ait_fm * m_0_ait_nc / (m_0_ait_fm + m_0_ait_nc)
    return (m_0_ait, m_3)
end

end

