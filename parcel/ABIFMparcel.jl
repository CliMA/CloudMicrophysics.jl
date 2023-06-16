import OrdinaryDiffEq as ODE
import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const CMO = CM.Common
const CMI = CM.HetIceNucleation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

"""
    ODE problem definition
"""
function cirrus_box(dY, Y, p, t)

    # Get simulation parameters
    (; prs, const_dt, r_nuc, w, α_m) = p
    # Numerical precision used in the simulation
    FT = eltype(Y)

    # Our state vector
    S_i = Y[1]        # supersaturation over ice
    N_act = Y[2]      # number concentration of activated particles
    p_a = Y[3]        # pressure
    T = Y[4]          # temperature
    q_vap = Y[5]      # vapor specific humidity
    q_ice = Y[6]      # ice specific humidity
    N_aerosol = Y[7]  # number concentration of interstitial aerosol
    J_immer_t = Y[8]  # J per unit area only
    P_ice_t = Y[9]    # ice produced
    x_sulph = Y[10]   # percent mass sulphuric acid

    # Constants
    R_v = CMP.R_v(prs)
    grav = CMP.grav(prs)
    ρ_ice = CMP.ρ_cloud_liq(prs)

    # Get thermodynamic parameters, phase partition and create thermo state.
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_ice, FT(0), q_ice)
    ts = TD.PhaseNonEquil_pTq(thermo_params, p_a, T, q)

    # Constants and variables that depend on the moisture content
    R_a = TD.gas_constant_air(thermo_params, q)
    cp_a = TD.cp_m(thermo_params, q)
    L_subl = TD.latent_heat_sublim(thermo_params, T)
    L_fus = TD.latent_heat_fusion(thermo_params, T)
    ρ = TD.air_density(thermo_params, ts)

    # We are solving for both q_ and supersaturation. This is unneccessary.
    # I'm checking here if we are consistent. We should decide
    # if we want to solve for S (as in the cirrucs box model) or for q_
    # which is closer to what ClimaAtmos is doing.
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    e = S_i * e_si
    e_check = q_vap * p_a * R_v / R_a
    @assert isapprox(e, e_check; rtol = 0.01)

    # Adiabatic parcel coefficients
    a1 = L_subl * grav / cp_a / T^2 / R_v - grav / R_a / T
    a2 = p_a / e_si * R_v / R_a
    a3 = L_subl^2 / R_v / T^2 / cp_a
    a4 = L_subl * L_fus / R_v / T^2 / cp_a

    # Activating new crystals
    # AF = CMI.dust_activated_number_fraction(prs, S_i, T, CMT.DesertDustType())
    J_immer = CMI.ABIFM_J(prs, x_sulph, T) * 10000 # converting cm^-2 s^-1 to m^-2 s^-1
    P_ice = J_immer * 4*π*r_nuc^2 * N_aerosol # per sec
    τ_relax = const_dt
    dN_act_dt = max(FT(0), P_ice * τ_relax)
    dN_aerosol_dt = -dN_act_dt
    dqi_dt_new_particles = dN_act_dt * 4 / 3 * π * r_nuc^3 * ρ_ice / ρ

    # Growing existing crystals (assuming all are the same...)
    G = CMO.G_func(prs, T, TD.Ice())
    r = N_act > 0 ? cbrt(q_ice / N_act / (4 / 3 * π) / ρ_ice * ρ) : 0
    C = r
    dqi_dt_deposition = 1 / ρ * N_act * α_m * 4 * π * C * (S_i - 1) * G
    dqi_dt_deposition = FT(0.0) # Use this if you don't want to grow

    # Sum of all phase changes
    dqi_dt = dqi_dt_new_particles + dqi_dt_deposition
    dqw_dt = FT(0.0)
    # TODO - update dqw_dt when implementing homo. and immersion freezing

    # Update the tendecies
    dS_i_dt = a1 * w * S_i - (a2 + a3 * S_i) * dqi_dt - (a2 + a4 * S_i) * dqw_dt
    dp_a_dt = -p_a * grav / R_a / T * w
    dT_dt = -grav / cp_a * w + L_subl / cp_a * dqi_dt
    dq_vap_dt = -dqi_dt
    dq_ice_dt = dqi_dt
    x_sulph_dt = FT(0.0)
    # dq_liq_dt = dqw_dt # Use this when introducing liquid water

    # Set tendencies
    dY[1] = dS_i_dt        # supersaturation over ice
    dY[2] = dN_act_dt      # number concentration of activated particles
    dY[3] = dp_a_dt        # pressure
    dY[4] = dT_dt          # temperature
    dY[5] = dq_vap_dt      # vapor specific humidity
    dY[6] = dq_ice_dt      # ice specific humidity
    dY[7] = dN_aerosol_dt  # number concentration of interstitial aerosol
    dY[8] = J_immer        # nucleation rate coefficient per unit area per unit time
    dY[9] = P_ice          # ice production rate
    dY[10] = x_sulph_dt    # %wt. sulphuric acid
    # add dY state for dq_liq_dt when introducing liquid

    # TODO - add diagnostics output (radius, S, etc)
end

"""
    Wrapper for initial condition
"""
function get_initial_condition(
    prs,
    N_act,
    p_a,
    T,
    q_vap,
    q_liq,
    q_ice,
    N_aerosol,
    x_sulph
)
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_liq + q_ice, q_liq, q_ice)
    R_a = TD.gas_constant_air(thermo_params, q)
    R_v = CMP.R_v(prs)
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    e = q_vap * p_a * R_v / R_a
    S_i = e / e_si
    J_immer_t = 0.0
    P_ice_t = 0.0
    x_sulph = x_sulph

    return [S_i, N_act, p_a, T, q_vap, q_ice, N_aerosol, J_immer_t, P_ice_t, x_sulph]
end

"""
    Wrapper for running the simulation. Rising parcel under cirrus initial conditions
        undergoing immersion freezing and growth by deposition only. Conservation of water 
        has not been fixed to account for immersion freezing. 
"""
function run_parcel(FT)

    # Boiler plate code to have access to model parameters and constants
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)
    thermo_params = CMP.thermodynamics_params(prs)

    # Initial conditions for 1st period
    N_aerosol = FT(500*1000) # m^-3
    N_0 = FT(0)
    p_0 = FT(20000)
    T_0 = FT(230)
    q_vap_0 = FT(0.0003345)
    q_liq_0 = FT(0)
    q_ice_0 = FT(0)
    x_sulph = FT(0.1)


    # Simulation time
    t_max = 30 * 60

    # Simulation parameters passed into ODE solver
    r_nuc = FT(1e-6) # assumed size of nucleated particles
    w = FT(0.1) # updraft speed, m/s
    α_m = FT(0.5) # accomodation coefficient
    const_dt = 0.1 # model timestep
    p = (; prs, const_dt, r_nuc, w, α_m)

    # Simulation 1
    IC1 = get_initial_condition(
        prs,
        N_0,
        p_0,
        T_0,
        q_vap_0,
        q_liq_0,
        q_ice_0,
        N_aerosol,
        x_sulph
    )
    prob1 = ODE.ODEProblem(cirrus_box, IC1, (FT(0), t_max), p)
    sol1 = ODE.solve(
        prob1,
        ODE.Euler(),
        dt = const_dt,
        reltol = 10 * eps(FT),
        abstol = 10 * eps(FT),
    )

    # Plot results
    fig = MK.Figure(resolution = (1000, 800))
    ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [-]")
    ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
    ax3 = MK.Axis(fig[2, 1], ylabel = "N act [1/dm3]", yscale = log10)
    ax4 = MK.Axis(fig[2, 2], ylabel = "N areo [1/dm3]")
    ax5 = MK.Axis(fig[3, 1], ylabel = "q_vap [g/kg]")
    ax6 = MK.Axis(fig[3, 2], ylabel = "q_ice [g/kg]")
    ax7 = MK.Axis(fig[4, 1], ylabel = "J [cm^-2 s^-1]", xlabel = "Height [m]", yscale = log10)
    ax8 = MK.Axis(fig[4, 2], ylabel = "P_ice [min^-1]", xlabel = "Height [m]", yscale = log10)

    MK.ylims!(ax1, 1.0, 1.5)
    MK.ylims!(ax3, 3, 2e3)
    MK.ylims!(ax5, 0.3, 0.4)
    MK.ylims!(ax6, 0, 0.013)
    MK.ylims!(ax7, 10, 10e8)
    MK.ylims!(ax8, 10^(-4), 10e3)

    MK.lines!(ax1, sol1.t * w, sol1[1, :])
    MK.lines!(ax2, sol1.t * w, sol1[4, :])
    MK.lines!(ax3, sol1.t * w, sol1[2, :] * 1e-3)
    MK.lines!(ax4, sol1.t * w, sol1[7, :] * 1e-3)
    MK.lines!(ax5, sol1.t * w, sol1[5, :] * 1e3)
    MK.lines!(ax6, sol1.t * w, sol1[6, :] * 1e3)
    MK.lines!(ax7, sol1.t * w, sol1[8, :] * 1e-4)
    MK.lines!(ax8, sol1.t * w, sol1[9, :] * 60)

    MK.save("cirrus_box.svg", fig)

    println(sol1[2,1:15])

end

run_parcel(Float64)
