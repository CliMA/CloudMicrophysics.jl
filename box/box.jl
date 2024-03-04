import OrdinaryDiffEq as ODE
import Random as RD
import Distributions as DS

import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het

"""
    ODE problem definitions
"""
function box_model(dY, Y, p, t)

    # Get simulation parameters
    (; tps, A_aero, aerosol, cooling_rate) = p
    # Numerical precision used in the simulation
    FT = eltype(Y)

    # Our state vector
    T = Y[1]          # temperature
    N_liq = Y[3]      # number concentration of existing water droplets
    N_ice = Y[4]      # number concentration of activated ice crystals

    Δa = FT(1) - CMO.a_w_ice(tps, T)
    J_immer = CMI_het.ABIFM_J(aerosol, Δa) # m^-2 s^-1

    # Update the tendecies
    dT_dt = -cooling_rate
    dN_ice_dt = FT(0)
    dN_liq_dt = FT(0)
    if N_liq > 0
        dN_ice_dt = J_immer * N_liq * A_aero
        dN_liq_dt = -dN_ice_dt
    end

    # Set tendencies
    dY[1] = dT_dt          # temperature
    dY[2] = FT(0)          # not used here
    dY[3] = dN_liq_dt      # number concentration of droplets
    dY[4] = dN_ice_dt      # number concentration of ice crystals
end

"""
    ODE problem definitions
"""
function box_model_with_probability(dY, Y, p, t)

    # Get simulation parameters
    (; tps, const_dt, A_aero, aerosol, cooling_rate, Aj_sorted, N₀) = p
    # Numerical precision used in the simulation
    FT = eltype(Y)

    # Our state vector
    T = Y[1]                      # temperature
    A = Y[2]                      # available surface area for freezing
    N_liq = max(FT(0), Y[3])      # number concentration of existing water droplets
    N_ice = max(FT(0), Y[4])      # number concentration of activated ice crystals

    # immersion freezing rate
    Δa = FT(1) - CMO.a_w_ice(tps, T)
    J_immer = CMI_het.ABIFM_J(aerosol, Δa)

    # Update the tendecies
    dT_dt = -cooling_rate
    dAj_dt = FT(0)
    dN_ice_dt = FT(0)
    dN_liq_dt = FT(0)
    if N_liq > 0
        n_frz = 0
        for j in range(start = 1, stop = N₀, step = 1)
            # Sample from the surface area distribution
            Aj = Aj_sorted[j]
            # Compute the freezing probability (eq 10)
            Pj_frz = FT(1) - exp(-Aj * J_immer * const_dt)
            # Sample from bimodal distribution for if a droplet is frozen or not (eq 5)
            freeze_event = RD.rand(DS.Binomial(1, Pj_frz), 1)[1]
            # Sum up all the droplets that froze
            n_frz += freeze_event
            if Bool(freeze_event)
                Aj_sorted[j] = FT(0)
            end
        end
        Aj_sum = sum(Aj_sorted)
        dAj_dt = (Aj_sum - A) / const_dt
        dN_ice_dt = max(FT(0), n_frz / const_dt)
        dN_liq_dt = -dN_ice_dt
    end

    # Set tendencies
    dY[1] = dT_dt          # temperature
    dY[2] = dAj_dt         # total Aj
    dY[3] = dN_liq_dt      # mumber concentration of droplets
    dY[4] = dN_ice_dt      # number concentration of ice activated particles

end

"""
    run_box(IC, t_0, t_end, p)

Returns the solution of an ODE probelm defined by the box model.

Inputs:
 - IC - A vector with the initial conditions for
   [T, A_sum, N_liq, N_ice]
 - t_0 - simulation start time
 - t_end - simulation end time
 - p - a named tuple with simulation parameters.

Initial condition contains (all in base SI units):
 - T - temperature
 - A_sum - available surface area for freezing
 - N_liq - cloud droplet number concnentration
 - N_ice - ice crystal number concentration

The named tuple p should contain:
 - tps - a struct with free parameters for Thermodynamics package,
 - const_dt - simulation timestep,
 - aerosol - insoluble aerosol parameters
 - cooling_rate - cooling rate
 - A_aero - assumed surface area for freezing (when assuming constant A)
 - Aj_sorted - a vector with available surface area (when variable)
 - N₀ - initial number of liquid droplets
"""
function run_box(IC, t_0, t_end, p)

    FT = eltype(IC)
    (; const_dt, flag) = p

    problem = if flag
        ODE.ODEProblem(box_model_with_probability, IC, (FT(t_0), FT(t_end)), p)
    else
        ODE.ODEProblem(box_model, IC, (FT(t_0), FT(t_end)), p)
    end
    sol = ODE.solve(
        problem,
        ODE.Euler(),
        dt = const_dt,
        reltol = 10 * eps(FT),
        abstol = 10 * eps(FT),
    )
    return sol
end
