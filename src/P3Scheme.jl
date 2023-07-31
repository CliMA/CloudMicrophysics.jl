"""
Predicted particle properties scheme (P3) for ice, which includes:
 - threshold solver
 - m(D) regime
"""
module P3Scheme

import NonlinearSolve as NLS
import NCDatasets as NC


# THINGS TO ADD TO PARAMETERS
# exponent in power law from Brown and Francis 1995 for mass grown by
# vapor diffusion and aggregation in midlatitude cirrus: (unitless I think?)
# const β_va::FT = 1.9
# coefficient in power law modified from Brown and Francis 1995 for mass grown by
# vapor diffusion and aggregation in midlatitude cirrus: (units of kg m^(-β_va) I think?)
# const α_va::FT = (7.38e-11) * 10^((6 * β_va) - 3)
# const ρ_i::FT = 916.7
# the threshold particle dimension from p. 292 of Morrison and Milbrandt 2015 
# between small spherical and large, nonspherical unrimed ice, D_th (m),
# which is a constant function of α_va, β_va, and ρ_i with no D-dependency
# const D_th::FT = ((FT(π) * ρ_i) / (6 * α_va))^(1 / (β_va - 3))
# ρ_i::FT = CMP.ρ_cloud_ice(param_set) # ρ_i (bulk density of ice (kg m^{-3}))

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
    u0::Vector{FT} = [FT(-7.6), FT(-8.2), FT(5.7), FT(5.4)],
) where {FT <: Real}

    β_va::FT = 1.9
    α_va::FT = (7.38e-11) * 10^((6 * β_va) - 3)

    if ρ_r == FT(0.0)
        throw(
            DomainError(
                ρ_r,
                "D_cr, D_gr, ρ_g, ρ_d are not physically relevant when no rime is present.",
            ),
        )
    elseif F_r == FT(0.0)
        throw(
            DomainError(
                F_r,
                "D_cr, D_gr, ρ_g, ρ_d are not physically relevant when no rime is present.",
            ),
        )
    elseif ρ_r > FT(997.0)
        throw(
            DomainError(
                ρ_r,
                "Predicted rime density ρ_r, being a density of bulk ice, cannot exceed the density of water (997 kg m^-3).",
            ),
        )
    elseif F_r == FT(1.0) || F_r > FT(1.0)
        throw(
            DomainError(
                F_r,
                "The rime mass fraction F_r is not physically defined for values greater than or equal to 1 because some fraction of the total mass must always consist of the mass of the unrimed portion of the particle.",
            ),
        )
    elseif F_r < FT(0)
        throw(DomainError(F_r, "Rime mass fraction F_r cannot be negative."))
    elseif ρ_r < FT(0)
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

"""
generate_threshold_table(ρ_r_dim, F_r_dim)

- ρ_r_dim: dimension of ρ_r vector in look-up table
- F_r_dim: dimension of F_r vector in look-up table

Employs thresholds(ρ_r, F_r, u0) to solve the nonlinear system
consisting of D_cr, D_gr, ρ_g, ρ_d for a range of predicted
rime density and rime mass fraction values.
generate_threshold_table() then returns a NetCDF file
look-up table which can be read as arrays containing
D_cr, D_gr, ρ_g, ρ_d without the longer run time
and memory which are byproducts of thresholds()'s dependency
on NonlinearSolve.jl.
"""
function generate_threshold_table(ρ_r_dim::T, F_r_dim::T, FT::Type) where {T <: Integer}
    # generate ranges of values based on the dimensions given in the function arguments:
    F_r_range = collect(range(start = 1e-10, stop = 0.995, length = F_r_dim))
    ρ_r_range = collect(range(start = 50, stop = 996, length = ρ_r_dim))

    # initialize arrays and indexes which will allow us to iterate over the arrays
    # and populate them with look-up values provided by thresholds()
    D_cr_vals = Array{FT}(undef, ρ_r_dim, F_r_dim)
    D_gr_vals = Array{FT}(undef, ρ_r_dim, F_r_dim)
    ρ_g_vals = Array{FT}(undef, ρ_r_dim, F_r_dim)
    ρ_d_vals = Array{FT}(undef, ρ_r_dim, F_r_dim)
    ρi = 0
    Fi = 0

    # populate the arrays with look-up values using a nested for loop
    for ρ_r in ρ_r_range
        ρi += 1
        println(string(ρi) * " " * string(ρ_r))
        for F_r in F_r_range
            Fi += 1
            println(string(Fi) * " " * string(F_r))
            D_cr_vals[ρi, Fi] = thresholds(ρ_r, F_r)[1]
            D_gr_vals[ρi, Fi] = thresholds(ρ_r, F_r)[2]
            ρ_g_vals[ρi, Fi] = thresholds(ρ_r, F_r)[3]
            ρ_d_vals[ρi, Fi] = thresholds(ρ_r, F_r)[4]
        end
        Fi = 0
    end

    # save as NetCDF: create NetCDF file
    table = NC.NCDataset("./test.nc", "c")

    # define dimensions ρ_r, F_r
    NC.defDim(table, "ρ_r", ρ_r_dim)
    NC.defDim(table, "F_r", F_r_dim)

    # define and write variables stored in the look-up table
    D_cr = NC.defVar(table, "D_cr", FT, ("ρ_r", "F_r"))
    D_cr[:,:] = D_cr_vals
    D_gr = NC.defVar(table, "D_gr", FT, ("ρ_r", "F_r"))
    D_gr[:,:] = D_gr_vals
    ρ_g = NC.defVar(table, "ρ_g", FT, ("ρ_r", "F_r"))
    ρ_g[:,:] = ρ_g_vals
    ρ_d = NC.defVar(table, "ρ_d", FT, ("ρ_r", "F_r"))
    ρ_d[:,:] = ρ_d_vals

    # add file attributes

    close(table)
end

"""
read_threshold_table(path, heatmap = "n")

- path: path to threshold table
- heatmap = "n" : option that toggles creation of a heatmap graph


Reads in the NetCDF file created by generate_threshold_table(),
then returns arrays which can be indexed at each time step to generate
D_cr, D_gr, ρ_g, ρ_d without:
(i) the longer run time and memory which are byproducts of thresholds()'s dependency
on NonlinearSolve.jl and
(ii) opening and closing the NetCDF file containing the lookup table at each timestep 
"""
function read_threshold_table(path::String = "/Users/rowan/Desktop/p3_scheme_work/CloudMicrophysics.jl/test.nc", heatmap::String = "n")
    table = NC.NCDataset(path, "r")
    D_cr = table["D_cr"]
    D_gr = table["D_gr"]
    ρ_g = table["ρ_g"]
    ρ_d = table["ρ_d"]
    D_cr_vals = D_cr[:,:]
    D_gr_vals = D_gr[:,:]
    ρ_g_vals = ρ_g[:,:]
    ρ_d_vals = ρ_d[:,:]

    return D_cr_vals, D_gr_vals, ρ_g_vals, ρ_d_vals
end

"""
lookup_threshold(ρ_r, F_r, vals; opt)

 - ρ_r: predicted rime density (q_rim/B_rim, kg/m3)
 - F_r: rime mass fraction (q_rim/q_i, --)
 - vals: matrix containing lookup values
 - opt: kwarg specifying which of the four look-up quantities is desired

Uses a matrix generated by read_threshold_table()
to look-up one or more quantities. This function can be used at each time step
because it does not open and close NetCDF objects; rather, it uses
arrays stored in local memory to generate values. If values
lie between lookup table points, a linear weighting between table points
is used to generate the returned values.
"""
function lookup_threshold(ρ_r::FT, F_r::FT, vals::Matrix{FT}; opt::String) where {FT <: Real}
end

end
