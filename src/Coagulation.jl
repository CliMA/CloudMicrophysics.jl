module Coagulation

"""
    whitby_intramodal_coag(ad, temp, rho_p)

 - `ad` - AerosolDistribution, a tuple of aerosol modes
 - `temp` - Temperature (K)
 - `rho_p` - Particle density (kg/m^3)
Calculates the change of the 0th and 6th moments for the Aitken and Accumulation modes
due to intramodal coagulation. Larger modes are not treated, as specified in Liu et al., 2012.

The equations below are analytical expressions for coagulation integrals for the 0-th and 6-th moments,
as found in Whitby et al., 1991.
For more information on how these expressions are produced, see the documentation page.
"""
function whitby_intramodal_coag(
    am,
    K_fm,
    K_nc,
    ESG, 
)
    SQD_gn = sqrt(D_gn)  # D_gn is the Geometric mean size of the number weighted lognormal distribution - aka 2*r_dry?

    # Whitby H.11a
    m_0_fm(am) = -am.N^2 * K_fm * b_0 * SQD_gn * (ESG + ESG^25 + 2 * ESG^5)
    # Whitby H.11b
    m_6_fm(am) = -2 * am.N^2 * K_fm * b_6 * SQD_gn^6 * (ESG^85 + 2*ESG^89 + ESG^109)
    # Whitby H.12a
    m_0_nc(am) = -am.N^2 * K_nc * (1 + ESG^8 + A * Kn_g * (ESG^20 + ESG^4))
    #Whitby H.12b
    m_6_nc(am) = -2 * am.N^2 * K_nc * D_gn^6 * (ESG^72 + ESG^80 + A * Kn_g * (ESG^52 + ESG^68))

    return (
        m_0_fm(am),
        m_6_fm(am),
        m_0_nc(am),
        m_6_nc(am)
    )
end

end
