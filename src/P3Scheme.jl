"""
Predicted particle properties scheme (P3) for ice, which includes:
 - threshold solver
 - m(D) regime
 - TODO

Implementation of Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

Note: Particle size is defined as its maximum length (i.e. max dimesion).
"""
module P3Scheme

import NonlinearSolve as NLS
import CLIMAParameters as CP
import ..Parameters as CMP

const PSP3 = CMP.CloudMicrophysicsParametersP3

export thresholds

"""
    α_va_si(p3)

 - p3 - a struct with P3 scheme parameters

Returns `α_va` coefficient for the assumed particle mass(size) relation for
large unrimed ice and dense nonspherical ice, in base SI units: kg m^(-β_va).
`β_va` is another coefficient of the mass(size) relation.
From measurements of mass grown by vapor diffusion and aggregation
in midlatitude cirrus by Brown and Francis (1995)
doi: 10.1175/1520-0426(1995)012<0410:IMOTIW>2.0.CO;2
"""
α_va_si(p3::PSP3) = p3.α_va * 10^(6 * p3.β_va - 3)

"""
    D_th_helper(p3)

 - p3 - a struct with P3 scheme parameters

Returns the critical size separating spherical and nonspherical ice, in meters.
Eq. 8 in Morrison and Milbrandt (2015).
"""
D_th_helper(p3::PSP3) = (π * p3.ρ_i / 6 / α_va_si(p3))^(1 / (p3.β_va - 3))

"""
    D_cr_helper(p3, F_r, ρ_g)

 - p3 - a struct with P3 scheme parameters
 - F_r - rime mass fraction (q_rim/q_i) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the size of equal mass for graupel and partially rimed ice, in meters.
Eq. 14 in Morrison and Milbrandt (2015).
"""
function D_cr_helper(p3::PSP3, F_r, ρ_g)
    α_va = α_va_si(p3)
    return (1 / (1 - F_r) * 6 * α_va / π / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    D_gr_helper(p3, ρ_g)

 - p3 - a struct with P3 scheme parameters
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the size of equal mass for graupel and unrimed ice, in meters.
Eq. 15 in Morrison and Milbrandt (2015).
"""
function D_gr_helper(p3::PSP3, ρ_g)
    α_va = α_va_si(p3)
    return (6 * α_va / π / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    ρ_g_helper(ρ_r, F_r, ρ_d)

 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the density of total (deposition + rime) ice mass for graupel, in kg/m3
Eq. 16 in Morrison and Milbrandt (2015).
"""
ρ_g_helper(ρ_r, F_r, ρ_d) = F_r * ρ_r + (1 - F_r) * ρ_d

"""
    ρ_d_helper(p3, D_cr, D_gr)

 - p3 - a struct with P3 scheme parameters
 - D_cr - is the size of equal mass for graupel and partially rimed ice, in meters
 - D_gr - the size of equal mass for graupel and unrimed ice, in meters

Returns the density of unrimed ice mass, in kg/m3
Eq. 17 in Morrison and Milbrandt (2015).
"""
function ρ_d_helper(p3::PSP3, D_cr, D_gr)
    α_va = α_va_si(p3)
    β_m2 = p3.β_va - 2
    FT = eltype(α_va)
    return 6 * α_va * (D_cr^β_m2 - D_gr^β_m2) / π / β_m2 /
           max(D_cr - D_gr, eps(FT))
end

"""
    thresholds(p3, ρ_r, F_r, u0)

 - p3 - a struct with P3 scheme parameters
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i) [-]
 - u0 - initial guess for the solver.
   Best option is to pass the natural logarithm of the solution from the previous time step.
   The default value is set to around log.(thresholds(400.0, 0.5, log.([0.00049, 0.00026, 306.668, 213.336])))

Solves the nonlinear system consisting of D_cr, D_gr, ρ_g, ρ_d
for a given rime density and rime mass fraction, where:
 - D_cr - is the threshold size separating partially rimed ice and graupel [m],
 - D_gr - is the threshold size separating graupel and dense nonspherical ice [m],
 - ρ_g - is the effective density of a spherical graupel particle [kg/m3],
 - ρ_d - is the density of the unrimed portion of the particle [kg/m3],
"""
function thresholds(
    p3::PSP3,
    ρ_r::FT,
    F_r::FT,
    u0::Vector{FT} = [FT(-7.6), FT(-8.2), FT(5.7), FT(5.4)], # TODO - Vectors won't work on GPU
) where {FT <: Real}

    @assert ρ_r > FT(0)   # rime density must be positive ...
    @assert ρ_r <= p3.ρ_l # ... and as a bulk ice density can't exceed the density of water
    @assert F_r > FT(0)   # rime mass fraction must be positive ...
    @assert F_r < FT(1)   # ... and there must always be some unrimed part
    p = (; p3, ρ_r, F_r)

    # Domain shift exp(u) to constrain the solutions to be positive.
    # We are solving a set of four equations: u[i] - _helper(u[i], p) = 0
    function f(u, p)
        _D_cr = exp(u[1])
        _D_gr = exp(u[2])
        _ρ_g = exp(u[3])
        _ρ_d = exp(u[4])
        return [
            _D_cr - D_cr_helper(p.p3, p.F_r, _ρ_g),
            _D_gr - D_gr_helper(p.p3, _ρ_g),
            _ρ_g - ρ_g_helper(p.ρ_r, p.F_r, _ρ_d),
            _ρ_d - ρ_d_helper(p.p3, _D_cr, _D_gr),
        ]
    end

    P3_prob = NLS.NonlinearProblem(f, u0, p)
    sol = NLS.solve(P3_prob, NLS.NewtonRaphson(), abstol = eps(FT))
    D_cr, D_gr, ρ_g, ρ_d = exp.(sol) # shift back into desired domain space
    return [D_cr, D_gr, ρ_g, ρ_d]
end

end
