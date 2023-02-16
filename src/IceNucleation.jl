"""
Parameterization for heterogenous ice nucleation. 
"""

"""
Working on diffusion coefficent, ventilation factor.
"""
function deposition_activated_fraction(
    ice_sat_ratio::NTuple{N, T},
    T::FT,              # Temp in Kelvins
    dustType::Int32     # i.e. 1 would correspond to Arizona Test Dust, 2-->desert dust
) where {FT <: Real}
    if dustType == 1    # Arizona Test Dust
        S0 = 1.05;      # threshold ice saturatio ratio
        if T > 220
            a = 4.7;    # average at 223K
        else
            a = 9.2;    # average at 210K
        end
    end
    if dustType == 2    # desert dust
        if T > 220
            a = 0.5;    # "less than 1"
            S0 = 1.14;  # threshold ice saturatio ratio
        else
            a = 2;
            S0 = 1.05;
        end
    end
    for (S) in ice_sat_ratio
        AF = exp(a * (S - S0)) - 1;
    end

    return AF
end

function immersion_activated_fraction

end

function karcher_coeffs(
    T::FT,
    P::FT,
    Si::FT
    )

    # constants
    R = 8.314;                  # universal ideal gas constant, J*K^-1*mol^-1
    k_B = 1.38*10^(-23);        # boltzmann constant, m^2*kg*s^-2*K^-1
    avogNum = 6.02*10^(-23);    # avogadro's number
    grav = 9.81;                # acceleration of gravity, m*s^-2
    depCoeff = 0.5;             # deposition coefficient of water molecules impinging on ice surface
    capacitance = 1;            # capacitance factor (accounts for geometry of ice crystal)

    # water
    L_subl = 2939*10^3          # latent heat of sublimation of water, J/kg
    m_water = 3*10^(-23)/1000;  # mass of one water molecule, kg/molecule
    MW_water = 18/1000;         # molar mass of water, kg/mol
    cp_water = 4.186*1000;      # specific heat capacity of water (constant P), J*kg^-1*K^-1
    specVol_water = 0.001/1000*MW_water/avogNum; # specific volume of water, m^3/molecule
    v_th_water = sqrt(8/pi*k_B/m_water*T);       # thermal speed of water, temp dependent, m/s

    # air
    MW_air = 28.97/1000;        # molar mass of air, kg/mol 

    # other
    n_sat = P_sat/T/k_B;   # water vapor number density @ ice saturation, P_sat and temp dependent, from IG law, molecules/m^3
    # define P_sat w/ Antoine's eqn? !!!!!!!!!!!!!!!!!!!!!!!
    diffCoeff = 9999999;   # diffusion coefficient of water in air !!FIX ME FIX ME FIX ME!!
    
    # Coefficients from Karcher et al (2006)
    a1 = ((L_subl*MW_water*grav)/(cp_water*R*T^2)) - ((MW_air*grav)/(R*T)); # temperature dependence, m^-1
    a2 = (n_sat)^-1;      # m^3/molecules
    a3 = ((L_subl^2*MW_water*m_water) / (cp_water*MW_air)) / (T*P); # temperature & pressure dependence, m^3/molecule

    b1 = specVol_water*depCoeff*v_th_water*n_sat*((Si-1)/4);    # S_i dependence
    b2 = depCoeff*v_th_water/4/capacitance/diffCoeff;           

    return a1, a2, a3, b1, b2
end