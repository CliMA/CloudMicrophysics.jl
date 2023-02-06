"""
Parameterization for heterogenous ice nucleation. 
"""

"""
Working on deposition_activated_fraction function. Check Smax and have AF be constant after reaching it?
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

