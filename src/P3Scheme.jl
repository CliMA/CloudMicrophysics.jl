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
        sol = NLS.solve(prob_obj, NLS.NewtonRaphson(), abstol = eps(FT))
        D_cr, D_gr, ρ_g, ρ_d = exp.(sol) # shift back into desired domain space
        return [D_cr, D_gr, ρ_g, ρ_d]
    end
    sol = NLS.nlsolve(f!,
        [0.001; 0.0003; 1000; 1000],
        autodiff = :forward)
    D_cr, D_gr, ρ_g, ρ_d = sol.zero
    return D_th, D_gr, D_cr
end

"""
generate_threshold_table(ρ_r_axis, F_r_axis)

- ρ_r_axis: range of ρ_r values used in the table
- F_r_axis: range of F_r values used in the table

Employs thresholds(ρ_r, F_r, u0) to solve the nonlinear system
consisting of D_cr, D_gr, ρ_g, ρ_d for a range of predicted
rime density and rime mass fraction values.
generate_threshold_table() then returns a NetCDF file
look-up table.
"""
function generate_threshold_table(;
    ρ_r_axis::Vector{FT} = exp10.(range(start = 2, stop = 2.965, length = 100)),
    F_r_axis::Vector{FT} = exp10.(
        range(start = -2, stop = -0.02, length = 100)
    ),
) where {FT <: Real}
    # generate ranges of values based on the dimensions given in the function arguments:
    F_r_vals = collect(F_r_axis)
    ρ_r_vals = collect(ρ_r_axis)

    # initialize arrays and indexes which will allow us to iterate over the arrays
    # and populate them with look-up values provided by thresholds()
    D_cr_vals = Array{FT}(undef, length(ρ_r_axis), length(F_r_axis))
    D_gr_vals = Array{FT}(undef, length(ρ_r_axis), length(F_r_axis))
    ρ_g_vals = Array{FT}(undef, length(ρ_r_axis), length(F_r_axis))
    ρ_d_vals = Array{FT}(undef, length(ρ_r_axis), length(F_r_axis))
    ρi = 0
    Fi = 0

    # populate the arrays with look-up values using a nested for loop
    for ρ_r in ρ_r_vals
        ρi += 1
        for F_r in F_r_vals
            Fi += 1
            D_cr_vals[ρi, Fi] = thresholds(ρ_r, F_r)[1]
            D_gr_vals[ρi, Fi] = thresholds(ρ_r, F_r)[2]
            ρ_g_vals[ρi, Fi] = thresholds(ρ_r, F_r)[3]
            ρ_d_vals[ρi, Fi] = thresholds(ρ_r, F_r)[4]
        end
        Fi = 0
    end

    # save as NetCDF: create NetCDF file
    table = NC.NCDataset("./p3_lookup.nc", "c")

    # define dimensions ρ_r, F_r
    NC.defDim(table, "ρ_r", length(ρ_r_axis))
    NC.defDim(table, "F_r", length(F_r_axis))

    # define and write variables stored in the look-up table
    D_cr = NC.defVar(table, "D_cr", FT, ("ρ_r", "F_r"))
    D_cr[:, :] = D_cr_vals
    D_gr = NC.defVar(table, "D_gr", FT, ("ρ_r", "F_r"))
    D_gr[:, :] = D_gr_vals
    ρ_g = NC.defVar(table, "ρ_g", FT, ("ρ_r", "F_r"))
    ρ_g[:, :] = ρ_g_vals
    ρ_d = NC.defVar(table, "ρ_d", FT, ("ρ_r", "F_r"))
    ρ_d[:, :] = ρ_d_vals

    close(table)
end

"""
read_threshold_table(path)

 - path: path to threshold table file

Reads in the NetCDF file created by generate_threshold_table(),
then returns a tuple of 4 matrices (containing look-up values for
D_cr, D_gr, ρ_g, ρ_d, respectively)
which can be indexed to generate
D_cr, D_gr, ρ_g, ρ_d.
"""
function read_threshold_table(path::String = "../p3_lookup.nc")
    # open file object
    table = NC.NCDataset(path, "r")

    # define quantities from variable names in the NetCDF file
    D_cr = table["D_cr"]
    D_gr = table["D_gr"]
    ρ_g = table["ρ_g"]
    ρ_d = table["ρ_d"]

    # grab values from the NetCDF file and write them to matrices
    D_cr_vals = D_cr[:, :]
    D_gr_vals = D_gr[:, :]
    ρ_g_vals = ρ_g[:, :]
    ρ_d_vals = ρ_d[:, :]

    # close table
    close(table)

    return D_cr_vals, D_gr_vals, ρ_g_vals, ρ_d_vals
end

"""
lookup_threshold(ρ_r, F_r, vals, ρ_r_axis, F_r_axis)

 - ρ_r: predicted rime density (q_rim/B_rim, kg/m3)
 - F_r: rime mass fraction (q_rim/q_i, --)
 - vals: tuple of matrices containing lookup values (generated by read_threshold_table())
 - ρ_r_axis, F_r_axis: ranges of ρ_r, F_r values which should match those passed to generate_threshold_table() -- see above

Uses a tuple of 4 matrices generated by read_threshold_table()
to look-up one or more quantities.
If values lie between lookup table points,
a linear interpolation between table points
is used to generate the returned values.
"""
function lookup_threshold(
    ρ_r::FT,
    F_r::FT,
    vals::NTuple{4, Matrix{FT}};
    ρ_r_axis::Vector{FT} = exp10.(range(start = 2, stop = 2.965, length = 100)),
    F_r_axis::Vector{FT} = exp10.(
        range(start = -2, stop = -0.02, length = 100)
    ),
) where {FT <: Real}
    ρ_r_vals = collect(ρ_r_axis)
    F_r_vals = collect(F_r_axis)
    # if out of bounds of the axes, throw an error
    if ρ_r < ρ_r_vals[1] ||
       F_r < F_r_vals[1] ||
       ρ_r > ρ_r_vals[length(ρ_r_vals)] ||
       F_r > F_r_vals[length(F_r_vals)]
        throw(
            ArgumentError(
                "ρ_r and/or F_r is out of bounds. Please select value(s) in the ranges of the ρ_r_axis and F_r_axis range object.",
            ),
        )
    else # proceed with the look-up
        # let ρ_i, F_i be the indices of the inputs (ρ_r and F_r)
        # in the vectors ρ_r_vals and F_r_vals;
        # that is, if ρ_r ∈ ρ_r_vals and F_r ∈ F_r_vals
        # if the inputs ρ_r, F_r are not in the vector of
        # representative values used in the look-up table
        # axis, indexin() will return nothing

        ρ_i = indexin(ρ_r, ρ_r_axis)[1]
        F_i = indexin(F_r, F_r_axis)[1]

        # define a function that can be used to index into the look-up
        # arrays, given coordinates corresponding to certain ρ_r ∈ ρ_r_vals
        # and F_r ∈ F_r_vals which returns a 4-element vector (result)
        # containing D_cr, D_gr, ρ_g, ρ_d
        function lookup(ρ_index, F_index)
            result = []
            for quantity in vals
                push!(result, quantity[ρ_index, F_index])
            end
            return result
        end

        # we can simply use lookup() if ρ_r and F_r align with values
        # in the dimensions of the look-up table 
        if ρ_r ∈ ρ_r_vals && F_r ∈ F_r_vals
            return lookup(ρ_i, F_i)

            # for the case in which ρ_r, F_r, or both do not align directly with values
            # in the dimensions of the look-up table:
        else
            # define a function which uses binary search to generate the
            # indices corresponding to values immediately less than (left)
            # and above (right) the input values for ρ_r and F_r, called lookup_val
            # in the vector representing a dimension of the look-up table,
            # called table_axis_vals
            function get_left_right(table_axis_vals::Vector, lookup_val::FT)
                # grab value at center of range and initialize indices
                i = 0
                j = length(table_axis_vals)

                # execute binary search while loop
                while i < j
                    n = floor(Int64, (i + j) / 2)
                    if table_axis_vals[n] < lookup_val
                        i = n + 1
                    else
                        j = n
                    end
                end

                return i - 1, i # "left" and "right" indices
            end

            # define a function which, given an input x (ρ_r or F_r) which is not
            # an element of the vector representing the corresponding dimension
            # of the look-up table, the elements immediately less than (x_left)
            # and greater than (x_right) this input x, and their corresponding
            # look-up 4-vectors containing D_cr, D_gr, ρ_g, ρ_d,
            # returns the linearly interpolated look-up 4-vector (y)
            # between y_left and y_right:
            function linear_interpolate(x, x_left, y_left, x_right, y_right)
                return y_left .+ (
                    (x_left - x) .* ((y_right .- y_left) ./ (x_right - x_left))
                )
            end

            # if ρ_r is not an element of the vector corresponding to the look-up table axis ρ_r_vals
            # but F_r is, we generate a linear interpolation between the values of ρ_r immediately above
            # and below the input ρ_r
            if ρ_r ∉ ρ_r_vals && F_r ∈ F_r_vals
                # get left and right indices
                ρ_i_left, ρ_i_right = get_left_right(ρ_r_vals, ρ_r)

                # get "y_left" and "y_right" for the linear_interpolate() function
                lookup_left = lookup(ρ_i_left, F_i)
                lookup_right = lookup(ρ_i_right, F_i)

                # return linear interpolation
                return linear_interpolate(
                    ρ_r,
                    ρ_r_vals[ρ_i_left],
                    lookup_left,
                    ρ_r_vals[ρ_i_right],
                    lookup_right,
                )

                # similarly for F_r:
            elseif F_r ∉ F_r_vals && ρ_r ∈ ρ_r_vals
                F_i_left, F_i_right = get_left_right(F_r_vals, F_r)
                lookup_left = lookup(ρ_i, F_i_left)
                lookup_right = lookup(ρ_i, F_i_right)
                return linear_interpolate(
                    F_r,
                    F_r_vals[F_i_left],
                    lookup_left,
                    F_r_vals[F_i_right],
                    lookup_right,
                )

                # if both input values are not elements of their corresponding dimension vectors,
                # we must perform a two-dimensional linear interpolation:
            elseif F_r ∉ F_r_vals && ρ_r ∉ ρ_r_vals
                # get left and right indices
                ρ_i_left, ρ_i_right = get_left_right(ρ_r_vals, ρ_r)
                F_i_left, F_i_right = get_left_right(F_r_vals, F_r)

                # get two sets of "y_left" and "y_right" for two linear interpolations:
                # keeping ρ constant, generate linear interpolations of F_r
                lookup_ρleft_Fleft = lookup(ρ_i_left, F_i_left)
                lookup_ρleft_Fright = lookup(ρ_i_left, F_i_right)
                lookup_ρleft = linear_interpolate(
                    F_r,
                    F_r_vals[F_i_left],
                    lookup_ρleft_Fleft,
                    F_r_vals[F_i_right],
                    lookup_ρleft_Fright,
                )

                lookup_ρright_Fleft = lookup(ρ_i_right, F_i_left)
                lookup_ρright_Fright = lookup(ρ_i_right, F_i_right)
                lookup_ρright = linear_interpolate(
                    F_r,
                    F_r_vals[F_i_left],
                    lookup_ρright_Fleft,
                    F_r_vals[F_i_right],
                    lookup_ρright_Fright,
                )

                # return the linear interpolation in ρ_r of the two linear interpolations in F_r:
                return linear_interpolate(
                    ρ_r,
                    ρ_r_vals[ρ_i_left],
                    lookup_ρleft,
                    ρ_r_vals[ρ_i_right],
                    lookup_ρright,
                )
            end
        end
    end
end

end

end # module P3Scheme.jl
