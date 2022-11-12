module Coagulation
include("Parameters.jl")
include("Common.jl")
include("CommonTypes.jl")
include("AerosolModel.jl")
include("CoagCorrectionFactors.jl")
import .CoagCorrectionFactors as CCF

import .AerosolModel as AM

function whitby_coagulation(ad)
    # TODO: add these values to climaparameters - this comes from CAM5
    surface_temp = 288.15
    surface_pressure = 101325.0
    mean_free_path = 6.6328e-8 * surface_pressure * temp / (surface_temp * pressure)

    # Call coags for aitken and accumulation modes:
    
end


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
    mean_free_path
)
    # Knudsen number (m)
    Kn_g = mean_free_path /  am.r_dry
    # Near-continuum factor
    A_nc = 1.43 * Kn_g^0.0814
    sqrt_diameter = sqrt(2 * am.r_dry)  

    # Whitby H.11a 
    m_0_fm(am) = -am.N^2 * K_fm * CCF.intramodal_correction(am.stdev) * sqrt_diameter * (ESG + ESG^25 + 2 * ESG^5)
    # Whitby H.11b
    m_6_fm(am) = -2 * am.N^2 * K_fm * CCF.intramodal_correction(am.stdev) * sqrt_diameter^6 * (ESG^85 + 2*ESG^89 + ESG^109)
    # Whitby H.12a
    m_0_nc(am) = -am.N^2 * K_nc * (1 + ESG^8 + A_nc * Kn_g * (ESG^20 + ESG^4))
    #Whitby H.12b
    m_6_nc(am) = -2 * am.N^2 * K_nc * (2*am.r_dry)^6 * (ESG^72 + ESG^80 + A_nc * Kn_g * (ESG^52 + ESG^68))

    return (
        m_0_fm(am),
        m_6_fm(am),
        m_0_nc(am),
        m_6_nc(am)
    )
end


function whitby_intermodal_aitken_coag(
    am,
    K_fm,
    K_nc,
    ESG,
    mean_free_path
    )

end

function whitby_intermodal_accumulation_coag(
    am,
    K_fm,
    K_nc,
    ESG,
    mean_free_path
    )

end

end
