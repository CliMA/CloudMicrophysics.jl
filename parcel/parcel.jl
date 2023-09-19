import Thermodynamics as TD
import CloudMicrophysics as CM

import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Common as CMO
import CloudMicrophysics.HetIceNucleation as CMI_het
import CloudMicrophysics.Parameters as CMP

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

#! format: off
"""
    ODE problem definitions
"""
function parcel_model(dY, Y, p, t)

    # Get simulation parameters
    (; prs, const_dt, r_nuc, w, α_m) = p
    (; ice_nucleation_modes, growth_modes) = p
    # Numerical precision used in the simulation
    FT = eltype(Y)

    # Our state vector
    # TODO - We are solving for both q_ and supersaturation.
    # We should decide if we want to solve for S (as in the cirrus box model)
    # or for q_ which is closer to what ClimaAtmos is doing.
    S_liq = Y[1]      # saturation ratio over liquid water
    p_a = Y[2]        # pressure
    T = Y[3]          # temperature
    q_vap = Y[4]      # vapor specific humidity
    q_liq = Y[5]      # liquid water specific humidity
    q_ice = Y[6]      # ice specific humidity
    N_aer = Y[7]      # number concentration of interstitial aerosol
    N_liq = Y[8]      # number concentration of existing water droplets
    N_ice = Y[9]      # number concentration of activated ice crystals
    x_sulph = Y[10]   # percent mass sulphuric acid

    # Constants
    R_v = CMP.R_v(prs)
    grav = CMP.grav(prs)
    ρ_ice = CMP.ρ_cloud_ice(prs)
    ρ_liq = CMP.ρ_cloud_liq(prs)

    # Get thermodynamic parameters, phase partition and create thermo state.
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_liq + q_ice, q_liq, q_ice)
    ts = TD.PhaseNonEquil_pTq(thermo_params, p_a, T, q)

    # Saturation ratio over ice
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    e_sl = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    ξ = e_sl / e_si
    S_i = ξ * S_liq

    # Constants and variables that depend on the moisture content
    R_a = TD.gas_constant_air(thermo_params, q)
    cp_a = TD.cp_m(thermo_params, q)
    L_subl = TD.latent_heat_sublim(thermo_params, T)
    L_vap = TD.latent_heat_vapor(thermo_params, T)
    ρ_air = TD.air_density(thermo_params, ts)

    # Adiabatic parcel coefficients
    a1 = L_vap * grav / cp_a / T^2 / R_v - grav / R_a / T
    a2 = 1 / q_vap
    a3 = L_vap^2 / R_v / T^2 / cp_a
    a4 = L_vap * L_subl / R_v / T^2 / cp_a

    # TODO - we should zero out all tendencies and augemnt them
    # TODO - add immersion, homogeneous, ...
    #dN_act_dt_immersion = FT(0)
    #dN_act_dt_homogeneous = FT(0)

    dN_act_dt_depo = FT(0)
    if "DustDeposition" in ice_nucleation_modes
        AF = CMI_het.dust_activated_number_fraction(
            prs,
            S_i,
            T,
            CMT.DesertDustType(),
        )
        dN_act_dt_depo = max(FT(0), AF * N_aer - N_ice) / const_dt
    end
    dN_ice_dt = dN_act_dt_depo
    dN_aer_dt = -dN_act_dt_depo
    dqi_dt_new_depo = dN_act_dt_depo * 4 / 3 * π * r_nuc^3 * ρ_ice / ρ_air

    # Growth
    dqi_dt_depo = FT(0)
    if "Deposition" in growth_modes && N_ice > 0
        # Deposition on existing crystals (assuming all are the same...)
        G_i = CMO.G_func(prs, T, TD.Ice())
        r_i = cbrt(q_ice / N_ice / (4 / 3 * π) / ρ_ice * ρ_air)
        C_i = r_i
        dqi_dt_depo = 4 * π / ρ_air * (S_i - 1) * G_i * r_i * N_ice
        # TODO - check if r_i is non-zero if initial conditions say q_ice = 0
    end

    dql_dt_cond = FT(0)
    if "Condensation" in growth_modes && N_liq > 0
        # Condensation on existing droplets (assuming all are the same)
        G_l = CMO.G_func(prs, T, TD.Liquid())
        r_l = cbrt(q_liq / N_liq / (4 / 3 * π) / ρ_liq * ρ_air)
        C_l = r_l
        dql_dt_cond = 4 * π / ρ_air * (S_liq - 1) * G_l * r_l * N_liq
    end

    # Update the tendecies
    dq_ice_dt = dqi_dt_new_depo + dqi_dt_depo
    dq_liq_dt = dql_dt_cond
    dS_l_dt = a1 * w * S_liq - (a2 + a3) * S_liq * dq_liq_dt - (a2 + a4) * S_liq * dq_ice_dt
    dp_a_dt = -p_a * grav / R_a / T * w
    dT_dt = -grav / cp_a * w + L_vap / cp_a * dq_liq_dt + L_subl / cp_a * dq_ice_dt
    dq_vap_dt = -dq_ice_dt - dq_liq_dt

    # Set tendencies
    dY[1] = dS_l_dt        # saturation ratio over liquid water
    dY[2] = dp_a_dt        # pressure
    dY[3] = dT_dt          # temperature
    dY[4] = dq_vap_dt      # vapor specific humidity
    dY[5] = dq_liq_dt      # liquid water specific humidity
    dY[6] = dq_ice_dt      # ice specific humidity
    dY[7] = dN_aer_dt      # number concentration of interstitial aerosol
    dY[8] = FT(0)          # mumber concentration of droplets
    dY[9] = dN_ice_dt      # number concentration of activated particles
    dY[10] = FT(0)         # sulphuric acid concentration

    # TODO - add diagnostics output (radius, S, etc)

end
#! format: on

"""
    run_parcel(IC, t_0, t_end, p)

Returns the solution of an ODE probelm defined by the parcel model.

Inputs:
 - IC - A vector with the initial conditions for
   [S_l, p_a, T, q_vap, q_liq, q_ice, N_aer, N_liq, N_ice, x_sulph]
 - t_0 - simulation start time
 - t_end - simulation end time
 - p - a named tuple with simulation parameters.

Initial condition contains (all in base SI units):
 - S_l - saturation ratio over liquid water,
 - p_a - atmospheric pressure
 - T - temperature
 - q_vap - water vapor specific humidity
 - q_liq - cloud liquid water specific humidity
 - q_ice - cloud ice specific humidity
 - N_aer - aerosol number concnetration
 - N_liq - cloud droplet number concnentration
 - N_ice - ice crystal number concentration
 - x_sulph - sulphuric acid concentration

The named tuple p should contain:
 - prs - a struct with free parameters for CloudMicrophysics package,
 - const_dt - simulation timestep,
 - r_nuc - assumed radius of newly nucleated ice crystals,
 - w - vertical velocity,
 - α_m - accomodation coefficient
 - ice_nucleation_modes - a vector with enabled ice nucleation paths. Possible options: ("DustDeposition",)
 - growth_modes - a vector with enabled growth modes. Possible options: ("Condensation", "Deposition")
"""
function run_parcel(IC, t_0, t_end, p)

    FT = eltype(IC)
    (; const_dt, ice_nucleation_modes, growth_modes) = p

    println(" ")
    println("Running parcel model with: ")
    print("Ice nucleation modes: ")
    if "DustDeposition" in ice_nucleation_modes
        print("Deposition on dust particles ")
    end
    print("\n")
    print("Growth modes: ")
    if "Condensation" in growth_modes
        print("Condensation on liquid droplets ")
    end
    if "Deposition" in growth_modes
        print("Deposition on ice crystals ")
    end
    print("\n")
    println(" ")

    problem = ODE.ODEProblem(parcel_model, IC, (FT(t_0), FT(t_end)), p)
    sol = ODE.solve(
        problem,
        ODE.Euler(),
        dt = const_dt,
        reltol = 10 * eps(FT),
        abstol = 10 * eps(FT),
    )
    return sol
end
