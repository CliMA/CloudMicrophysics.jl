"""
Predicted particle properties scheme (P3) for ice, which includes:
 - m(D) regime
 - A(D) regime
 - particle fall speed regime
 - add more (source/sink) as we go! (TODO)
"""
module P3Scheme

import SpecialFunctions as SF

import Thermodynamics as TD

import ..CommonTypes as CT
import ..Common as CO
import ..Parameters as CMP
import NLsolve as NLS

const FT = Float64
const APS = CMP.AbstractCloudMicrophysicsParameters
const ρ_i::FT = 917
const β_va::FT = 1.9
const α_va::FT = 1e-3 * (7.38e-11) * 10^(6 * β_va)
const D_th::FT = ((FT(π) * 6) / (6 * α_va))^(1 / (β_va - 3))
const γ::FT = 0.2285
const σ::FT = 1.88
const δ_o::FT = 5.83 # UNSURE about these four constants below
const C_o::FT = 0.6
const a_o::FT = 1e-5
const b_o::FT = 1.0

# more constants for now until we add constants to CliMA parameters:
include(joinpath(pkgdir(CloudMicrophysics), "test", "create_parameters.jl"))
toml_dict = CLIMAParameters.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)
thermo_params = CloudMicrophysics.Parameters.thermodynamics_params(param_set)
const ρ_air::FT = param_set.ρ0_SB2006
const ν_air::FT = param_set.ν_air

# TODO: should I define the constants like ρ_i and α_va here?
#       polish code
#       figure out what B_rim is
#       figure out how to compute ρ_d, ρ_g, D_cr, D_gr:
#           they "form a closed set of equations that can be solved by iteration"
#           so do we use IterativeSolvers.jl here? Or do we have an in-house solver?
#       finish adding A(D), V(D) regimes
#       instead of using if/elif/else statements, would it be better to define
#           new types (i.e. small spherical ice, graupel, etc) ?

"""
breakpoints(ρ_r, F_r)

- ρ_r: predicted rime density (q_rim/B_rim)
- F_r: rime mass fraction (q_rim/q_i)

Computes breakpoints (D_cr, D_th, D_gr) of the system
 which govern m(D) and A(D) regimes.
"""
function breakpoints(ρ_r::FT, F_r::FT) where {FT <: Real}
    D_th = ((FT(π) * 6) / (6 * α_va))^(1 / (β_va - 3))
    # Let x[1] = D_cr, x[2] = D_gr, x[3] = ρ_g, x[3] = ρ_d,
    # and let each corresponding component function of F
    # be defined F[i] = x[i] - a[i] where x[i] = a[i]
    # such that F[x] = 0:
    function f!(F, x)
        F[1] =
            x[1] - ((1 / (1 - F_r)) *
            ((6 * α_va) /
            (FT(π) * x[3]))) ^ (1 / (3 - β_va))
        F[2] =
            x[2] - (((6 * α_va) /
            (FT(π) * x[3])) ^ (1 / (3 - β_va)))
        F[3] =
            x[3] - (ρ_r * F_r) - ((1 - F_r) * x[4])
        F[4] =
            x[4] -
            (((6 * α_va) *
            ((x[1]^(β_va - 2)) -
            (x[2]^(β_va - 2)))) /
            (FT(π) *
            (β_va - 2) *
            (x[1] - x[2])))
    end
    sol = NLS.nlsolve(f!,
        [0.001; 0.0003; 1000; 1000],
        autodiff = :forward)
    D_cr, D_gr, ρ_g, ρ_d = sol.zero
    return D_th, D_gr, D_cr

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
 - ρ_r: predicted rime density (q_rim/B_rim)
 - F_r: rime mass fraction (q_rim/q_i)

 m(D) regime,
 which computes breakpoints, classifies particles, and returns mass
"""
function m(D::FT, ρ_r::FT, F_r::FT) where {FT <: Real}
    D_th, D_gr, D_cr = breakpoints(ρ_r, F_r)
    if D <= D_th
        return m_s(D, ρ_i) # small spherical ice
    elseif q_rim == 0
        return m_nl(D) # large, nonspherical, unrimed ice
    else
        if D <= D_cr
            return m_r(D, F_r) # partially rimed ice
        elseif D > D_cr & D >= D_gr
            return m_s(D, ρ_g) # graupel
        elseif D < D_gr
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
    return (FT(π) / 4) * D^2
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
A(D, ρ_r, F_r)

 - D: maximum particle dimension
 - ρ_r: predicted rime density (q_rim/B_rim)
 - F_r: rime mass fraction (q_rim/q_i)

Assumed particle projected area regime,
which computes breakpoints, classifies particles, and returns area
"""
function A(D::FT, ρ_r::FT, F_r::FT) where {FT <: Real}
    D_th, D_gr, D_cr = breakpoints(ρ_r, F_r)
    if D < D_th
        return A_s(D) # small, spherical ice
    elseif q_rim == 0
        return A_ns(D) # large, unrimed ice
    else
        if D <= D_cr
            return (F_r * A_s(D)) + (1 - F_r) * A_ns(D) # partially rimed ice
        elseif D > D_cr & D >= D_gr
            return A_s(D) # graupel
        elseif D < D_gr
            return A_ns(D) # dense nonspherical ice
        end
    end
end

"""
V(param_set, D, ρ, ρ_mod="n")

 - param_set: abstract set with Earth parameters
 - D: maximum particle dimension
 - ρ: particle density
 - ρ_mod="n": optional keyword argument specifying whether
    the density modification is to be applied
    (default is no modification)

Assumed particle fall speed regime
"""
function V(param_set::APS, D::FT, ρ::FT, ρ_mod="n") where {FT <: Real}
    ν_air = param_set.ν_air
    g = param_set.grav
    ρ_air = param_set.ρ0_SB2006
    μ_air = 1.67e-5
    X = (2 * α_va * g * ρ_air * D^(β_va - 2 + σ)) /
     (γ * μ_air^2)
    C_1 = 4 / (δ_o^2 * C_o^0.5)
    C_2 = (δ_o^2) / 4
    b_1 = ((C_1 * X^0.5) /
     (2 * (((1 + (C_1 * X^0.5))^0.5) - 1) *
     (1 + (C_1 * X^0.5))^0.5)) -
     ((a_o * b_o * X^b_o) /
     (C_2 *
     (((1 + (C_1 * X^0.5))^0.5) - 1)^2))
    a_1 = ((C_2 *
     (((1 + (C_1 * X^0.5))^0.5) - 1)^2) -
     (a_o * X^b_o)) \
     X^b_1
    X, C_1, C_2, a_1. b_1
    a = a_1 * 
     ν_air^(1 - 2 * b_1) *
     ((2 * α_va * g) /
     (ρ_air * γ))^b_1
    b = (b_1 * (2 - σ) - 1) \
     (1 - b_1)
     V = a * D^b
     if ρ_mod == "n"
        V = a * D^b
        return V
    elseif ρ_mod == "y"
        return (ρ / ρ_air)^0.54 * V
    end
end

end
