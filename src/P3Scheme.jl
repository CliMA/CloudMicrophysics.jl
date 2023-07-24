"""
Predicted particle properties scheme (P3) for ice, which includes:
 - m(D) regime
 - A(D) regime
 - particle fall speed regime
 - add more (source/sink) as we go! (TODO)
"""
module P3Scheme

import NonlinearSolve as NLS
import ..Parameters as CMP

const APS = CMP.AbstractCloudMicrophysicsParameters
FT = Float64

# exponent in power law from Brown and Francis 1995 for mass grown by
# vapor diffusion and aggregation in midlatitude cirrus: (unitless I think?)
const β_va::FT = 1.9
# coefficient in power law modified from Brown and Francis 1995 for mass grown by
# vapor diffusion and aggregation in midlatitude cirrus: (units of kg m^(-β_va) I think?)
const α_va::FT = (7.38e-11) * 10^((6 * β_va) - 3)

"""
D_th(param_set::APS)

 - prs: abstract set with Earth parameters

Returns the threshold particle dimension from p. 292 of Morrison and Milbrandt 2015 
between small spherical and large, nonspherical unrimed ice, D_th (m),
which is a constant function of α_va, β_va, and ρ_i (bulk density of ice (kg m^{-3})), with no D-dependency.
"""
function D_th(param_set::APS) where {FT <: Real}
    ρ_i::FT = CMP.ρ_cloud_ice(param_set)
    # add alpha and beta once they get into CloudMicrophysicsParameters
    return ((FT(π) * ρ_i) / (6 * α_va))^(1 / (β_va - 3))
end

"""
thresholds(ρ_r, F_r, u0)

- ρ_r: predicted rime density (q_rim/B_rim)
    - [ρ_r] = ``kg m^{-3}``
- F_r: rime mass fraction (q_rim/q_i)
    - [F_r] = none, unitless
- u0: initial guess for solver:
    - Ideally, one passes the natural logarithm of the
    values returned from the previous time step for u0
    - However, if this is not possible, the default value is set to around
    log.(thresholds(400.0, 0.5, log.([0.00049, 0.00026, 306.668, 213.336])))

Solves the nonlinear system consisting of D_cr, D_gr, ρ_g, ρ_d
for a given predicted rime density and rime mass fraction, where:
    - D_cr, defined on p. 292 of Morrison and Milbrandt 2015,
    is the threshold particle dimension separating partially rimed ice and graupel,
    given here in meters
    - D_gr, defined on p. 293 of Morrison and Milbrandt 2015,
    is the threshold particle dimension separating graupel and dense nonspherical ice,
    given here in meters
    - ρ_g, defined on pp. 292-293 of Morrison and Milbrandt 2015,
    is the effective density of an assumed spherical graupel particle,
    given here in ``kg m^{-3}``
    - ρ_d, defined on p. 293 of Morrison and Milbrandt 2015,
    is the density of the unrimed portion of the particle,
    given here in ``kg m^{-3}``
"""
function thresholds(
    ρ_r::FT,
    F_r::FT,
    u0::Vector{FT} = [-7.6, -8.2, 5.7, 5.4],
) where {FT <: Real}
    if ρ_r == 0.0
        throw(
            DomainError(
                ρ_r,
                "D_cr, D_gr, ρ_g, ρ_d are not physically relevant when no rime is present.",
            ),
        )
    elseif F_r == 0.0
        throw(
            DomainError(
                ρ_r,
                "D_cr, D_gr, ρ_g, ρ_d are not physically relevant when no rime is present.",
            ),
        )
    elseif ρ_r > 997
        throw(
            DomainError(
                ρ_r,
                "Predicted rime density ρ_r, being a density of bulk ice, cannot exceed the density of water.",
            ),
        )
    elseif F_r == 1.0 || F_r > 1.0
        throw(
            DomainError(
                F_r,
                "The rime mass fraction F_r is not physically defined for values greater than or equal to 1 because some fraction of the total mass must always consist of the mass of the unrimed portion of the particle.",
            ),
        )
    elseif F_r < 0
        throw(DomainError(F_r, "Rime mass fraction F_r cannot be negative."))
    elseif ρ_r < 0
        throw(
            DomainError(ρ_r, "Predicted rime density ρ_r cannot be negative."),
        )
    else
        # Let u[1] = D_cr, u[2] = D_gr, u[3] = ρ_g, u[4] = ρ_d,
        # and let each corresponding component function of F        
        # be defined F[i] = x[i] - a[i] where x[i] = a[i]
        # such that F[x] = 0:
        function f(u, p)
            # Implementation of function with domain shift exp(u):
            # This domain shift is necessary because it constrains
            # the solver to search only for positive values, preventing
            # a DomainError with complex exponentiation.
            # The domain restriction is reasonable in any case 
            # because the quantities of interest, [D_cr, D_gr, ρ_g, ρ_d],
            # should all be positive regardless of ρ_r and F_r.
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

        p = Nothing # (no parameters)
        prob_obj = NLS.NonlinearProblem(f, u0, p)
        sol = NLS.solve(prob_obj, NLS.NewtonRaphson(), reltol = 1e-9)
        D_cr, D_gr, ρ_g, ρ_d = exp.(sol) # shift back into desired domain space

        return [D_cr, D_gr, ρ_g, ρ_d]
    end
end

end