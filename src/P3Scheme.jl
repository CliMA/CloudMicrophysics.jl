"""
Predicted particle properties scheme (P3) for cloud ice, which includes:
 - m(D) regime
 - A(D) regime
 - particle fall speed regime
 - add more as we go! (TODO)
"""
module P3Scheme

import SpecialFunctions as SF

import Thermodynamics as TD

import ..CommonTypes as CT
import ..Common as CO
import ..Parameters as CMP

const APS = CMP.AbstractCloudMicrophysicsParameters
const ρ_i::Float64 = 917
const α_va::Float64 = 7.38e-11
const β_va::Float64 = 1.9
const D_th::Float64 = ((Float64(π) * 6)/(6 * α_va))^(1/(β_va - 3))
const γ::Float64 = 0.2285
const σ::Float64 = 1.88

# TODO: should I define the constants like ρ_i and α_va here?
#       polish code
#       figure out what B_rim is
#       figure out how to compute ρ_d, ρ_g, D_cr, D_gr:
#           they "form a closed set of equations that can be solved by iteration"
#           so do we use IterativeSolvers.jl here? Or do we have an in-house solver?
#       finish adding A(D), V(D) regimes
"""
ρ_r(q_rim, B_rim)

 - q_rim: rime mass mixing ratio
 - B_rim: TODO???

Predicted rime density
"""
function ρ_r(q_rim::FT, B_rim::FT) where {FT <: Real}
    return q_rim/B_rim
end

"""
m_s(D, ρ)

 - D: maximum particle dimension
 - ρ: bulk ice density (note to self: use ρ_i for small ice and ρ_g for graupel)

m(D) relation for spherical ice (small ice or completely rimed ice)
"""
function m_s(D::FT, ρ::FT) where {FT <: Real}
    return (FT(π)/6) * ρ * D^3
end

"""
m_nl(D)

 - D: maximum particle dimension


m(D) relation for large, nonspherical ice (used for unrimed and dense types)
"""
function m_nl(D::FT) where {FT <: Real}
    return α_va * D^β_va
end

"""
m_r(D, q_rim, q_i)

 - D: maximum particle dimension
 - q_rim: rime mass mixing ratio
 - q_i: total ice mass mixing ratio
 - NOTE: if calculating rime mass fraction F_r
  iteratively, we can just replace the 
  mixing ratios with F_r = q_rim/q_i

m(D) relation for partially rimed ice
"""
function m_r(D::FT, q_rim::FT, q_i::FT) where {FT <: Real}
    F_r = q_rim/q_i
    return (α_va/(1-F_r)) * D^β_va
end

"""
m(D, q_rim, q_i)

 - D: maximum particle dimension
 - q_rim: rime mass mixing ratio
 - q_i: total ice mass mixing ratio

 m(D) regime,
 which computes breakpoints, classifies particles, and returns mass
"""
function m(D::FT, q_rim::FT, q_i::FT, B_rim::FT) where {FT <: Real}
    if D <= D_th:
        return m_s(D, ρ_i) # small spherical ice
    elseif q_rim == 0:
        return m_nl(D) # large, nonspherical, unrimed ice
    else:
        # find variable densities:
        ρ_r = q_rim/B_rim
        ρ_d = 1 # unsure how to compute this, depends on D_cr/gr which depend on it
        ρ_g = 1 # ^
        D_cr = 1 # ^
        D_gr = 2 # ^
        if D <= D_cr:
            return m_r(D, q_rim, q_i) # partially rimed ice
        elseif {D > D_cr} & {D >= D_gr}:
            return m_s(D, ρ_g) # graupel
        elseif D < D_gr:
            return m_nl(D) # dense nonspherical ice
        end
    end
end

"""
A_s(D)


 - D: maximum particle dimension
 
Particle projected area relation for assumed spherical particles
"""
function A_s(D::FT) where {FT <: Real}
    return (FT(π)/4) * D^2
end

"""
A_ns(D)

 - D: maximum particle dimension

Particle projected area relation for assumed nonspherical particles,
    from Mitchell 1996
"""
function A_ns(D::FT) where {FT <: Real}
    return γ * D^σ
end

"""
A(D, q_rim, q_i)

 - D: maximum particle dimension
 - q_rim: rime mass mixing ratio
 - q_i: total ice mass mixing ratio
 - NOTE: if calculating rime mass fraction F_r
  iteratively, we can just replace the 
  mixing ratios with F_r = q_rim/q_i

Assumed particle projected area regime,
which computes breakpoints, classifies particles, and returns area
"""
function A(D::FT, q_rim::FT, q_i::FT) where {FT <: Real}
    # placeholder values until iterative solver is implemented:
    D_cr = 1
    ρ_d = 1
    ρ_g = 1
    D_gr = 1
    F_r = q_rim/q_i
    if D < D_th:
        return A_s(D) # small, spherical ice
    elseif q_rim == 0:
        return A_ns(D) # large, unrimed ice
    else:
        if D <= D_cr:
            return (F_r * A_s(D)) + (1 - F_r) * A_ns(D) # partially rimed ice
        elseif {D > D_cr} & {D >= D_gr}:
            return A_s(D) # graupel
        elseif D < D_gr:
            return A_ns(D) # dense nonspherical ice
        end
    end
end