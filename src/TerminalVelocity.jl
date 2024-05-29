"""
Terminal velocity calculations including: 
- Chen 2022 parametrization for ice and rain 
"""
module TerminalVelocity

import CloudMicrophysics.Common as CO
import CloudMicrophysics.Parameters as CMP

export velocity_chen

"""
    ϕ_ice(mass, area, ρ_i) 

 - mass - mass of particle 
 - area - area of particle 
 - ρ_i - density of ice

 Returns the aspect ratio (ϕ) for a particle with given mass, area, and ice density
"""
function ϕ_ice(mass::FT, area::FT, ρ_i::FT) where {FT}
    return 16 * ρ_i^2 * area^3 / (9 * π * mass^2)
end

"""
    velocity_chen(D, Chen2022, ρ_a, mass, area, ρ_i)

 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
 - ρ_a - density of air  
 - mass - mass of particle 
 - area - area of particle 
 - ρ_i - density of ice

 Returns the terminal velocity of ice at given particle dimension using Chen 2022 parametrizations
"""
function velocity_chen(
    D::FT,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    ρ_a::FT,
    mass::FT,
    area::FT,
    ρ_i::FT,
) where {FT}
    if D <= Chen2022.cutoff
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)
    else
        (ai, bi, ci) = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
    end

    κ = FT(-1 / 6)

    v = sum(@. sum(ai * D^bi * exp(-ci * D)))

    return ϕ_ice(mass, area, ρ_i)^κ * v
end

"""
    velocity_chen(p3, Chen2022, ρ_a, D)
  
 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
 - ρ_a - density of air
 
 Returns the terminal velocity of rain given Chen 2022 velocity parametrizations
"""
function velocity_chen(
    D::FT,
    Chen2022::CMP.Chen2022VelTypeRain,
    ρ_a::FT,
) where {FT}
    (ai, bi, ci) = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)

    v = sum(@. sum(ai * D^bi * exp(-ci * D)))

    return v
end

end
