import OrdinaryDiffEq as ODE
import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

const CMT = CM.CommonTypes
const CMO = CM.Common
const CMI_het = CM.HetIceNucleation
const CMI_hom = CM.HomIceNucleation
const CMP = CM.Parameters

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

"""
    ODE problem definition
"""
function cirrus_box(dY, Y, p, t)

    # Get simulation parameters
    (; prs, const_dt, r_nuc, w, α_m, freeze_mode, deposition_growth) = p
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
    x_sulph = Y[8]   # percent mass sulphuric acid

    # Constants
    R_v = CMP.R_v(prs)
    grav = CMP.grav(prs)
    ρ_ice = CMP.ρ_cloud_liq(prs)

    # Get thermodynamic parameters, phase partition and create thermo state.
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_ice, FT(0), q_ice) # TODO add liquid
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
    τ_relax = const_dt
    if freeze_mode == "deposition"

        AF = CMI_het.dust_activated_number_fraction(
            prs,
            S_i,
            T,
            CMT.DesertDustType(),
        )

        dN_act_dt = max(FT(0), AF * N_aerosol - N_act) / τ_relax
    elseif freeze_mode == "ABIFM"

        Delta_a_w =
            T > FT(185) && T < FT(235) ? CMO.Delta_a_w(prs, x_sulph, T) : FT(0)
        J_immer =
            T > FT(185) && T < FT(235) ?
            CMI_het.ABIFM_J(CMT.DesertDustType(), Delta_a_w) : FT(0)
        P_ice = J_immer * 4 * π * r_nuc^2 * N_aerosol

        dN_act_dt = max(FT(0), P_ice)
    elseif freeze_mode == "homogeneous"

        Delta_a_w =
            T > FT(185) && T < FT(235) ? CMO.Delta_a_w(prs, x_sulph, T) : FT(0)
        J_homogeneous =
            Delta_a_w > 0.26 && Delta_a_w < 0.34 ?
            CMI_hom.homogeneous_J(Delta_a_w) : FT(0)
        P_ice = J_homogeneous * 4 / 3 * π * r_nuc^3 * N_aerosol

        dN_act_dt = max(FT(0), P_ice)
    else
        @warn "Invalid freezing mode in run_parcel argument. Running without freezing.\nPlease choose between \"deposition\", \"ABIFM\", or \"homogeneous\" (include quotation marks)."
        dN_act_dt = FT(0)
    end
    dN_aerosol_dt = -dN_act_dt
    dqi_dt_new_particles = dN_act_dt * 4 / 3 * π * r_nuc^3 * ρ_ice / ρ

    # Growing existing crystals (assuming all are the same...)
    G = CMO.G_func(prs, T, TD.Ice())
    r = N_act > 0 ? cbrt(q_ice / N_act / (4 / 3 * π) / ρ_ice * ρ) : 0
    C = r
    dqi_dt_deposition =
        deposition_growth == true ?
        1 / ρ * N_act * α_m * 4 * π * C * (S_i - 1) * G : FT(0)

    # Sum of all phase changes
    dqi_dt = dqi_dt_new_particles + dqi_dt_deposition
    dqw_dt = freeze_mode == "deposition" ? FT(0) : FT(0)
    # TODO - update dqw_dt when implementing homo. and immersion freezing

    # Update the tendecies
    dS_i_dt = a1 * w * S_i - (a2 + a3 * S_i) * dqi_dt - (a2 + a4 * S_i) * dqw_dt
    dp_a_dt = -p_a * grav / R_a / T * w
    dT_dt = -grav / cp_a * w + L_subl / cp_a * dqi_dt
    dq_vap_dt = -dqi_dt
    dq_ice_dt = dqi_dt
    # dq_liq_dt = dqw_dt # Use this when introducing liquid water
    x_sulph_dt = FT(0)

    # Set tendencies
    dY[1] = dS_i_dt        # supersaturation over ice
    dY[2] = dN_act_dt      # number concentration of activated particles
    dY[3] = dp_a_dt        # pressure
    dY[4] = dT_dt          # temperature
    dY[5] = dq_vap_dt      # vapor specific humidity
    dY[6] = dq_ice_dt      # ice specific humidity
    dY[7] = dN_aerosol_dt  # number concentration of interstitial aerosol
    dY[8] = x_sulph_dt        # nucleation rate coefficient per unit area per unit time
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
    x_sulph,
)
    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_liq + q_ice, q_liq, q_ice)
    R_a = TD.gas_constant_air(thermo_params, q)
    R_v = CMP.R_v(prs)
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    e = q_vap * p_a * R_v / R_a
    S_i = e / e_si
    x_sulph = x_sulph

    return [S_i, N_act, p_a, T, q_vap, q_ice, N_aerosol, x_sulph]
end

"""
    Wrapper for running the simulation following the same framework as in
    Tully et al 2022 (10.5194/egusphere-2022-1057)
     - the simulation consists of 3 periods mimicking 3 large scale model steps
     - each period is 30 minutes long
     - each period is run with user specified constant timestep
     - the large scale initial conditions between each period are different
    Possible freeze_mode inputs are the following strings
     - "deposition"
"""
function run_parcel(FT, freeze_mode, deposition_growth = true)

    # Boiler plate code to have access to model parameters and constants
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)
    thermo_params = CMP.thermodynamics_params(prs)

    # Initial conditions for 1st period
    N_aerosol = FT(2000 * 1e3)
    N_0 = FT(0)
    p_0 = FT(20000)
    T_0 = FT(230)
    q_vap_0 = FT(0.0003345)
    q_liq_0 = FT(0)
    q_ice_0 = FT(0)
    x_sulph = (0.1)
    # Initial conditions for the 2nd period
    T2 = FT(229.25)
    q_vap2 = FT(0.00034)
    # Initial conditions for the 3rd period
    T3 = FT(228.55)
    q_vap3 = FT(0.000345)

    # Simulation time
    t_max = 30 * 60

    # Simulation parameters passed into ODE solver
    r_nuc = FT(0.5 * 1.e-6) # assumed size of nucleated particles, meters
    w = FT(3.5 * 1e-2) # updraft speed
    α_m = FT(0.5) # accomodation coefficient
    const_dt = 0.1 # model timestep
    p = (; prs, const_dt, r_nuc, w, α_m, freeze_mode, deposition_growth)

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
        x_sulph,
    )
    prob1 = ODE.ODEProblem(cirrus_box, IC1, (FT(0), t_max), p)
    sol1 = ODE.solve(
        prob1,
        ODE.Euler(),
        dt = const_dt,
        reltol = 10 * eps(FT),
        abstol = 10 * eps(FT),
    )

    # Simulation 2
    # (alternatively set T and take q_vap from the previous simulation)
    #IC2 = get_initial_condition(sol1[2, end], sol1[3, end], T2, sol1[5, end], 0.0, sol1[6, end], sol1[7, end])
    IC2 = get_initial_condition(
        prs,
        sol1[2, end],
        sol1[3, end],
        sol1[4, end],
        q_vap2,
        q_liq_0,
        sol1[6, end],
        sol1[7, end],
        x_sulph,
    )
    prob2 = ODE.ODEProblem(
        cirrus_box,
        IC2,
        (FT(sol1.t[end]), sol1.t[end] + t_max),
        p,
    )
    sol2 = ODE.solve(
        prob2,
        ODE.Euler(),
        dt = const_dt,
        reltol = 10 * eps(FT),
        abstol = 10 * eps(FT),
    )

    # Simulation 3
    # (alternatively set T and take q_vap from the previous simulation)
    #IC3 = get_initial_condition(sol2[2, end], sol2[3, end], T3, sol2[5, end], 0.0, sol2[6, end], sol2[7, end])
    IC3 = get_initial_condition(
        prs,
        sol2[2, end],
        sol2[3, end],
        sol2[4, end],
        q_vap3,
        q_liq_0,
        sol2[6, end],
        sol2[7, end],
        x_sulph,
    )
    prob3 = ODE.ODEProblem(
        cirrus_box,
        IC3,
        (FT(sol2.t[end]), sol2.t[end] + t_max),
        p,
    )
    sol3 = ODE.solve(
        prob3,
        ODE.Euler(),
        dt = const_dt,
        reltol = 10 * eps(FT),
        abstol = 10 * eps(FT),
    )

    # Plot results
    fig = MK.Figure(resolution = (800, 600))
    ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [-]")
    ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
    ax3 = MK.Axis(fig[2, 1], ylabel = "N act [1/dm3]", yscale = log10)
    ax4 = MK.Axis(fig[2, 2], ylabel = "N areo [1/dm3]")
    ax5 = MK.Axis(fig[3, 1], ylabel = "q_vap [g/kg]", xlabel = "Height [m]")
    ax6 = MK.Axis(fig[3, 2], ylabel = "q_ice [g/kg]", xlabel = "Height [m]")

    MK.ylims!(ax1, 1.0, 1.5)
    MK.ylims!(ax3, 3, 2e3)

    MK.lines!(ax1, sol1.t * w, sol1[1, :])
    MK.lines!(ax1, sol2.t * w, sol2[1, :])
    MK.lines!(ax1, sol3.t * w, sol3[1, :])

    MK.lines!(ax2, sol1.t * w, sol1[4, :])
    MK.lines!(ax2, sol2.t * w, sol2[4, :])
    MK.lines!(ax2, sol3.t * w, sol3[4, :])

    MK.lines!(ax3, sol1.t * w, sol1[2, :] * 1e-3)
    MK.lines!(ax3, sol2.t * w, sol2[2, :] * 1e-3)
    MK.lines!(ax3, sol3.t * w, sol3[2, :] * 1e-3)

    MK.lines!(ax4, sol1.t * w, sol1[7, :] * 1e-3)
    MK.lines!(ax4, sol2.t * w, sol2[7, :] * 1e-3)
    MK.lines!(ax4, sol3.t * w, sol3[7, :] * 1e-3)

    MK.lines!(ax5, sol1.t * w, sol1[5, :] * 1e3)
    MK.lines!(ax5, sol2.t * w, sol2[5, :] * 1e3)
    MK.lines!(ax5, sol3.t * w, sol3[5, :] * 1e3)

    MK.lines!(ax6, sol1.t * w, sol1[6, :] * 1e3)
    MK.lines!(ax6, sol2.t * w, sol2[6, :] * 1e3)
    MK.lines!(ax6, sol3.t * w, sol3[6, :] * 1e3)

    MK.save("cirrus_box.svg", fig)
end

run_parcel(Float64, "deposition")
