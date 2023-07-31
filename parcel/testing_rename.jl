import OrdinaryDiffEq as ODE
import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import CLIMAParameters as CP

# boilerplate code to get free parameter values
include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "parcel.jl"))

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

    return [S_i, N_act, p_a, T, q_vap, q_ice, N_aerosol, x_sulph]
end

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
    x_sulph = FT(0.1)

    # Simulation time
    t_max = 30 * 60 * 30

    # Simulation parameters passed into ODE solver
    r_nuc = FT(0.5 * 1.e-4 * 1e-6) # assumed size of nucleated particles
    w = FT(1.0 * 1e-2) # updraft speed
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
    prob1 = ODE.ODEProblem(parcel_model, IC1, (FT(0), t_max), p)
    sol1 = ODE.solve(
        prob1,
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

    MK.lines!(ax2, sol1.t * w, sol1[4, :])

    MK.lines!(ax3, sol1.t * w, sol1[2, :] * 1e-3)

    MK.lines!(ax4, sol1.t * w, sol1[7, :] * 1e-3)

    MK.lines!(ax5, sol1.t * w, sol1[5, :] * 1e3)

    MK.lines!(ax6, sol1.t * w, sol1[6, :] * 1e3)

    MK.save("cirrus_box.svg", fig)
end

run_parcel(Float64, "ABIFM")
