"""
Predicted particle properties scheme (P3) for cloud ice, which includes:
 - m(D) regime
 - A(D) regime
 - fall speed regime
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

# TODO: add more constants... should I define the constants like ρ_i and α_va here?
#       code helper functions for m(D) for various types
#       add m(D) that utilizes helper functions

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
m_va(D)

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
  iteratively, we can just replace the qs with F_r = q_rim/q_i

m(D) relation for partially rimed ice
"""
function m_r(D::FT, q_rim::FT, q_i::FT) where {FT <: Real}
    F_r = q_rim/q_i
    return (α_va/(1-F_r)) * D^β_va
end