import OrdinaryDiffEq as ODE
import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# boilerplate code to get free parameter values
include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))

"""
    Wrapper for initial condition
"""
function get_initial_condition(
    prs,
    p_a,
    T,
    q_vap,
    q_liq,
    q_ice,
    N_aer,
    N_liq,
    N_ice,
    x_sulph,
)

    thermo_params = CMP.thermodynamics_params(prs)
    q = TD.PhasePartition(q_vap + q_liq + q_ice, q_liq, q_ice)
    R_a = TD.gas_constant_air(thermo_params, q)
    R_v = CMP.R_v(prs)
    e_si = TD.saturation_vapor_pressure(thermo_params, T, TD.Ice())
    e = q_vap * p_a * R_v / R_a
    S_i = e / e_si

    return [S_i, p_a, T, q_vap, q_liq, q_ice, N_aer, N_liq, N_ice, x_sulph]
end

"""
    Wrapper for running the simulation following the same framework as in
    Tully et al 2023 (https://doi.org/10.5194/gmd-16-2957-2023)
     - the simulation consists of 3 periods mimicking 3 large scale model steps
     - each period is 30 minutes long
     - each period is run with user specified constant timestep
     - the large scale initial conditions between each period are different
"""
function Tully_et_al_2023(FT)

    # Boiler plate code to have access to model parameters and constants
    toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
    prs = cloud_microphysics_parameters(toml_dict)
    thermo_params = CMP.thermodynamics_params(prs)

    # Initial conditions for 1st period
    N_aerosol = FT(2000 * 1e3)
    N_droplets = FT(0)
    N_0 = FT(0)
    p_0 = FT(20000)
    T_0 = FT(230)
    q_vap_0 = FT(0.0003345)
    q_liq_0 = FT(0)
    q_ice_0 = FT(0)
    x_sulph = FT(0)
    # Initial conditions for the 2nd period
    T2 = FT(229.25)
    q_vap2 = FT(0.00034)
    # Initial conditions for the 3rd period
    T3 = FT(228.55)
    q_vap3 = FT(0.000345)

    # Simulation time
    t_max = 30 * 60

    # Simulation parameters passed into ODE solver
    r_nuc = FT(0.5 * 1.e-4 * 1e-6)             # assumed size of nucleated particles
    w = FT(3.5 * 1e-2)                         # updraft speed
    α_m = FT(0.5)                              # accomodation coefficient
    const_dt = 0.1                             # model timestep
    ice_nucleation_modes = ["DustDeposition",] # switch on deposition on dust
    growth_modes = ["Deposition",]             # switch on deposition growth
    p = (; prs, const_dt, r_nuc, w, α_m, ice_nucleation_modes, growth_modes)

    # Simulation 1
    IC1 = get_initial_condition(
        prs,
        p_0,
        T_0,
        q_vap_0,
        q_liq_0,
        q_ice_0,
        N_aerosol,
        N_droplets,
        N_0,
        x_sulph,
    )
    sol1 = run_parcel(IC1, 0, t_max, p)

    # Simulation 2
    # (alternatively set T and take q_vap from the previous simulation)
    #IC2 = get_initial_condition(sol1[2, end], sol1[3, end], T2, sol1[5, end], 0.0, sol1[6, end], sol1[7, end])
    IC2 = get_initial_condition(
        prs,
        sol1[2, end],
        sol1[3, end],
        q_vap2,
        q_liq_0,
        sol1[6, end],
        sol1[7, end],
        sol1[8, end],
        sol1[9, end],
        x_sulph,
    )
    sol2 = run_parcel(IC2, sol1.t[end], sol1.t[end] + t_max, p)

    # Simulation 3
    # (alternatively set T and take q_vap from the previous simulation)
    #IC3 = get_initial_condition(sol2[2, end], sol2[3, end], T3, sol2[5, end], 0.0, sol2[6, end], sol2[7, end])
    IC3 = get_initial_condition(
        prs,
        sol2[2, end],
        sol2[3, end],
        q_vap3,
        q_liq_0,
        sol2[6, end],
        sol2[7, end],
        sol2[8, end],
        sol2[9, end],
        x_sulph,
    )
    sol3 = run_parcel(IC3, sol2.t[end], sol2.t[end] + t_max, p)

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

    MK.lines!(ax2, sol1.t * w, sol1[3, :])
    MK.lines!(ax2, sol2.t * w, sol2[3, :])
    MK.lines!(ax2, sol3.t * w, sol3[3, :])

    MK.lines!(ax3, sol1.t * w, sol1[9, :] * 1e-3)
    MK.lines!(ax3, sol2.t * w, sol2[9, :] * 1e-3)
    MK.lines!(ax3, sol3.t * w, sol3[9, :] * 1e-3)

    MK.lines!(ax4, sol1.t * w, sol1[7, :] * 1e-3)
    MK.lines!(ax4, sol2.t * w, sol2[7, :] * 1e-3)
    MK.lines!(ax4, sol3.t * w, sol3[7, :] * 1e-3)

    MK.lines!(ax5, sol1.t * w, sol1[4, :] * 1e3)
    MK.lines!(ax5, sol2.t * w, sol2[4, :] * 1e3)
    MK.lines!(ax5, sol3.t * w, sol3[4, :] * 1e3)

    MK.lines!(ax6, sol1.t * w, sol1[6, :] * 1e3)
    MK.lines!(ax6, sol2.t * w, sol2[6, :] * 1e3)
    MK.lines!(ax6, sol3.t * w, sol3[6, :] * 1e3)

    MK.save("cirrus_box.svg", fig)
end

Tully_et_al_2023(Float64)
