"""
Parameterization for heterogenous ice nucleation. 
"""

"""
Working on mohler_activated_fraction function. Need to define a, S0, and Smax
"""
function mohler_activated_fraction(
    S0::FT, # threshold ice saturatio ratio
    Smax::FT, # max ice saturation ratio. maybe make a separate function for this?
    ice_sat_ratio::NTuple{N,T}
) where {FT <: Real}

    # at S >= S_max, there is a constant f_max
    for (S) in ice_sat_ratio
        AF = exp(a*(S - S0)) - 1;
    end

    return AF
end



