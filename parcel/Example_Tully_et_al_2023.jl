import OrdinaryDiffEq as ODE
import CairoMakie as MK

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

# definition of the ODE problem for parcel model
include(joinpath(pkgdir(CM), "parcel", "Parcel.jl"))

"""
    Wrapper for initial condition
"""
function get_initial_condition(tps, p_air, T, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC)
    q = TD.PhasePartition(qᵥ + qₗ + qᵢ, qₗ, qᵢ)
    R_a = TD.gas_constant_air(tps, q)
    R_v = TD.Parameters.R_v(tps)
    e_sl = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    e = eᵥ(qᵥ, p_air, R_a, R_v)
    Sₗ = e / e_sl

    return [Sₗ, p_air, T, qᵥ, qₗ, qᵢ, Nₐ, Nₗ, Nᵢ, ln_INPC]
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

    # get free parameters
    tps = TD.Parameters.ThermodynamicsParameters(FT)
    # Initial conditions for 1st period
    N_aerosol = FT(2000 * 1e3)
    N_droplets = FT(0)
    N_0 = FT(0)
    p_0 = FT(20000)
    T_0 = FT(230)
    q_vap_0 = FT(3.3e-4)
    q_liq_0 = FT(0)
    q_ice_0 = FT(0)
    ln_INPC = FT(0)
    # Initial conditions for the 2nd period
    T2 = FT(229.25)
    q_vap2 = FT(3.3e-4)
    # Initial conditions for the 3rd period
    T3 = FT(228.55)
    q_vap3 = FT(3.3e-4)

    # Simulation time
    t_max = 30 * 60

    # Simulation parameters passed into ODE solver
    r_nuc = FT(0.5 * 1.e-4 * 1e-6)                    # assumed size of nucleated particles
    w = FT(3.5 * 1e-2)                                # updraft speed
    const_dt = 0.1                                    # model timestep
    aerosol = CMP.DesertDust(FT)                      # aerosol type
    ice_nucleation_modes = ["MohlerAF", "MohlerRate"] # ice nucleation modes
    deposition_growth = "Deposition"                  # switch on deposition growth
    size_distribution = "Monodisperse"

    # Plots
    fig = MK.Figure(size = (800, 600))
    ax1 = MK.Axis(fig[1, 1], ylabel = "Supersaturation [-]")
    ax2 = MK.Axis(fig[1, 2], ylabel = "Temperature [K]")
    ax3 = MK.Axis(fig[2, 1], ylabel = "N act [1/dm3]", yscale = log10)
    ax4 = MK.Axis(fig[2, 2], ylabel = "N areo [1/dm3]")
    ax5 = MK.Axis(fig[3, 1], ylabel = "q_vap [g/kg]", xlabel = "Height [m]")
    ax6 = MK.Axis(fig[3, 2], ylabel = "q_ice [g/kg]", xlabel = "Height [m]")
    MK.ylims!(ax3, 3, 2e3)

    for mode in ice_nucleation_modes
        params = parcel_params{FT}(
            const_dt = const_dt,
            r_nuc = r_nuc,
            w = w,
            aerosol = aerosol,
            deposition = mode,
            deposition_growth = deposition_growth,
            liq_size_distribution = size_distribution,
        )
        # Simulation 1
        IC1 = get_initial_condition(
            tps,
            p_0,
            T_0,
            q_vap_0,
            q_liq_0,
            q_ice_0,
            N_aerosol,
            N_droplets,
            N_0,
            ln_INPC,
        )
        sol1 = run_parcel(IC1, 0, t_max, params)
        if mode == "MohlerAF"
            # Simulation 2
            # (alternatively set T and take q_vap from the previous simulation)
            #IC2 = get_initial_condition(sol1[2, end], sol1[3, end], T2, sol1[5, end], 0.0, sol1[6, end], sol1[7, end])
            IC2 = get_initial_condition(
                tps,
                sol1[2, end],
                #sol1[3, end],
                T2,
                q_vap2,
                q_liq_0,
                sol1[6, end],
                sol1[7, end],
                sol1[8, end],
                sol1[9, end],
                ln_INPC,
            )
            sol2 = run_parcel(IC2, sol1.t[end], sol1.t[end] + t_max, params)

            # Simulation 3
            # (alternatively set T and take q_vap from the previous simulation)
            #IC3 = get_initial_condition(sol2[2, end], sol2[3, end], T3, sol2[5, end], 0.0, sol2[6, end], sol2[7, end])
            IC3 = get_initial_condition(
                tps,
                sol2[2, end],
                #sol2[3, end],
                T3,
                q_vap3,
                q_liq_0,
                sol2[6, end],
                sol2[7, end],
                sol2[8, end],
                sol2[9, end],
                ln_INPC,
            )
            sol3 = run_parcel(IC3, sol2.t[end], sol2.t[end] + t_max, params)

            # Plot results
            sol = [sol1, sol2, sol3]
            clr = ["blue", "orange", "green"]
            #! format: off
            for it in [1, 2, 3]
                MK.lines!(ax1, sol[it].t * w, sol[it][1, :] .- 1, color = clr[it])
                MK.lines!(ax2, sol[it].t * w, sol[it][3, :], color = clr[it])
                MK.lines!(ax3, sol[it].t * w, sol[it][9, :] * 1e-3, color = clr[it])
                MK.lines!(ax4, sol[it].t * w, sol[it][7, :] * 1e-3, color = clr[it])
                MK.lines!(ax5, sol[it].t * w, sol[it][4, :] * 1e3, color = clr[it])
                MK.lines!(ax6, sol[it].t * w, sol[it][6, :] * 1e3, color = clr[it])
                MK.lines!(
                    ax1,
                    sol[it].t * w,
                    S_i.(tps, sol[it][3, :], sol[it][1, :]) .- 1,
                    linestyle = :dash,
                    color = clr[it],
                )
            end
            #! format: on

        elseif mode == "MohlerRate"
            MK.lines!(ax1, sol1.t * w, sol1[1, :] .- 1, color = :lightblue)
            MK.lines!(ax2, sol1.t * w, sol1[3, :], color = :lightblue)
            MK.lines!(ax3, sol1.t * w, sol1[9, :] * 1e-3, color = :lightblue)
            MK.lines!(ax4, sol1.t * w, sol1[7, :] * 1e-3, color = :lightblue)
            MK.lines!(ax5, sol1.t * w, sol1[4, :] * 1e3, color = :lightblue)
            MK.lines!(ax6, sol1.t * w, sol1[6, :] * 1e3, color = :lightblue)

            MK.lines!(
                ax1,
                sol1.t * w,
                S_i.(tps, sol1[3, :], sol1[1, :]) .- 1,
                linestyle = :dash,
                color = :lightblue,
            )
        end
    end
    MK.save("cirrus_box.svg", fig)
end
Tully_et_al_2023(Float32)
