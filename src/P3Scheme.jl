"""
Predicted particle properties scheme (P3) for ice, which includes:
 - m(D) regime
 - A(D) regime
 - particle fall speed regime
 - add more (source/sink) as we go! (TODO)
"""
module P3Scheme

import CloudMicrophysics as CM
import Thermodynamics as TD
import NonlinearSolve as NLS

const CT = CM.CommonTypes
const CO = CM.Common
const CMP = CM.Parameters
const APS = CMP.AbstractCloudMicrophysicsParameters

const FT = Float64

const ρ_i::FT = 917.0
const β_va::FT = 1.9
const α_va::FT = (7.38e-11) * 10^((6 * β_va) - 3)
const D_th::FT = ((FT(π) * ρ_i) / (6 * α_va))^(1 / (β_va - 3))
const γ::FT = 0.2285
const σ::FT = 1.88
const δ_o::FT = 5.83 # UNSURE about these four constants below
const C_o::FT = 0.6
const a_o::FT = 1.7e-3
const b_o::FT = 0.8

"""
breakpoints(ρ_r, F_r)

- ρ_r: predicted rime density (q_rim/B_rim)
- F_r: rime mass fraction (q_rim/q_i)

Solves the nonlinear system consisting of D_cr, D_gr, ρ_g, ρ_d
 for a given predicted rime density and rime mass fraction.      
"""
function breakpoints(ρ_r::FT, F_r::FT) where {FT <: Real}
    if ρ_r == 0.0 || F_r == 0.0
        return [NaN64, NaN64, NaN64, NaN64]
    
    else
        # Let u[1] = D_cr, u[2] = D_gr, u[3] = ρ_g, u[3] = ρ_d,
        # and let each corresponding component function of F        
        # be defined F[i] = x[i] - a[i] where x[i] = a[i]
        # such that F[x] = 0:
        function f(u, p)
            # implementation of function with domain shift exp(u)
            return [
                (exp(u[1])) -
                (
                    (1 / (1 - F_r)) * ((6 * α_va) / (FT(π) * exp(u[3])))
                )^(1 / (3 - β_va)),
                (exp(u[2])) -
                (((6 * α_va) / (FT(π) * (exp(u[3]))))^(1 / (3 - β_va))),
                (exp(u[3])) - (ρ_r * F_r) - ((1 - F_r) * (exp(u[4]))),
                (exp(u[4])) - (
                    (
                        (6 * α_va) *
                        ((exp(u[1])^(β_va - 2)) - ((exp(u[2]))^(β_va - 2)))
                    ) / (
                        FT(π) *
                        (β_va - 2) *
                        (max((exp(u[1])) - (exp(u[2])), 1e-16))
                    )
                ),
            ]
        end
        
        u0 = [-5.0, -6.0, 6.0, 8.0] # guess for solver
        p = [0.0] # (no parameters)
        prob_obj = NLS.NonlinearProblem(f, u0, p)
        sol = NLS.solve(prob_obj, NLS.NewtonRaphson(), reltol = 1e-9)
        D_cr, D_gr, ρ_g, ρ_d = exp.(sol) # shift back into desired domain space
        
        return [D_cr, D_gr, ρ_g, ρ_d]
    end
end

"""
m_s(D, ρ)

 - D: maximum particle dimension
 - ρ: bulk ice density (note to self: use ρ_i for small ice and ρ_g for graupel)

m(D) relation for spherical ice (small ice or completely rimed ice)
"""
function m_s(D::FT, ρ::FT) where {FT <: Real}
    return (FT(π) / 6) * ρ * D^3
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
m_r(D, F_r)

 - D: maximum particle dimension
 - F_r: rime mass fraction (q_rim/q_i)

m(D) relation for partially rimed ice
"""
function m_r(D::FT, F_r::FT) where {FT <: Real}
    return (α_va / (1 - F_r)) * D^β_va
end

"""
m(D, ρ_r, F_r)

 - D: maximum particle dimension
 - breakpoints: result of breakpoints() function, i.e.
  a vector containing D_cr, D_gr, ρ_g, ρ_d, in that order
 - F_r: rime mass fraction (q_rim/q_i)

 m(D) regime,
 which computes breakpoints, classifies particles, and returns mass;
 used to create figures for the docs page.
"""
function m(D::FT, breakpoints::Vector{FT}, F_r::FT) where {FT <: Real}
    if D <= D_th
        return m_s(D, ρ_i) # small spherical ice
    elseif F_r == 0
        return m_nl(D) # large, nonspherical, unrimed ice
    else
        if D >= breakpoints[1]
            return m_r(D, F_r) # partially rimed ice
        elseif D < breakpoints[1]
            if D >= breakpoints[2]
                return m_s(D, breakpoints[3]) # graupel
            elseif D < breakpoints[2]
                return m_nl(D) # dense nonspherical ice
            end
        end
    end
end

end
