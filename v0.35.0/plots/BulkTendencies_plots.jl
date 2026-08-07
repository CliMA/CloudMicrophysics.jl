import CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.BulkMicrophysicsTendencies as BMT
using CairoMakie

"""
Integrate the original nonlinear microphysics system with explicit Euler substepping.

This is used as a reference trajectory for the nonlinear system.
"""
function integrate_bulk_microphysics_reference(
    cm::BMT.Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T0,
    q_tot,
    q_lcl0,
    q_icl0,
    q_rai0,
    q_sno0,
    t_end;
    dt = 0.1,
    N_lcl = zero(ρ),
)
    FT = typeof(q_tot)
    nsteps = Int(round(t_end / dt))
    times = collect(FT, range(zero(FT), step = FT(dt), length = nsteps + 1))

    q_lcl = similar(times)
    q_icl = similar(times)
    q_rai = similar(times)
    q_sno = similar(times)
    T = similar(times)

    q_lcl[1] = q_lcl0
    q_icl[1] = q_icl0
    q_rai[1] = q_rai0
    q_sno[1] = q_sno0
    T[1] = T0

    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    for i in 1:nsteps
        rates = BMT.bulk_microphysics_tendencies(
            cm,
            mp,
            tps,
            ρ,
            T[i],
            q_tot,
            q_lcl[i],
            q_icl[i],
            q_rai[i],
            q_sno[i],
            N_lcl,
        )

        q_lcl[i + 1] = q_lcl[i] + FT(dt) * rates.dq_lcl_dt
        q_icl[i + 1] = q_icl[i] + FT(dt) * rates.dq_icl_dt
        q_rai[i + 1] = q_rai[i] + FT(dt) * rates.dq_rai_dt
        q_sno[i + 1] = q_sno[i] + FT(dt) * rates.dq_sno_dt

        dT_dt =
            Lv_over_cp * (rates.dq_lcl_dt + rates.dq_rai_dt) +
            Ls_over_cp * (rates.dq_icl_dt + rates.dq_sno_dt)

        T[i + 1] = T[i] + FT(dt) * dT_dt
    end

    return (; times, q_lcl, q_icl, q_rai, q_sno, T)
end

"""
Integrate the linearized implicit microphysics model over one full step `t_end`,
using `nsub` internal linearized substeps.

Returns the intermediate state after each internal substep, so for `nsub = 1`
the returned trajectory contains exactly two points.
"""
function integrate_bulk_microphysics_linearized_one_step(
    cm::BMT.Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T0,
    q_tot,
    q_lcl0,
    q_icl0,
    q_rai0,
    q_sno0,
    t_end;
    nsub = 1,
    N_lcl = zero(ρ),
)
    FT = typeof(q_tot)

    times = collect(FT, range(zero(FT), stop = FT(t_end), length = nsub + 1))

    q_lcl = similar(times)
    q_icl = similar(times)
    q_rai = similar(times)
    q_sno = similar(times)
    T = similar(times)

    q_lcl[1] = q_lcl0
    q_icl[1] = q_icl0
    q_rai[1] = q_rai0
    q_sno[1] = q_sno0
    T[1] = T0

    Δt_sub = FT(t_end) / FT(nsub)

    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    for i in 1:nsub
        rates = BMT._average_bulk_microphysics_tendencies(
            cm,
            mp,
            tps,
            ρ,
            T[i],
            q_tot,
            q_lcl[i],
            q_icl[i],
            q_rai[i],
            q_sno[i],
            Δt_sub,
            N_lcl,
        )

        q_lcl[i + 1] = q_lcl[i] + Δt_sub * rates.dq_lcl_dt
        q_icl[i + 1] = q_icl[i] + Δt_sub * rates.dq_icl_dt
        q_rai[i + 1] = q_rai[i] + Δt_sub * rates.dq_rai_dt
        q_sno[i + 1] = q_sno[i] + Δt_sub * rates.dq_sno_dt

        dT_dt =
            Lv_over_cp * (rates.dq_lcl_dt + rates.dq_rai_dt) +
            Ls_over_cp * (rates.dq_icl_dt + rates.dq_sno_dt)

        T[i + 1] = T[i] + Δt_sub * dT_dt
    end

    return (; times, q_lcl, q_icl, q_rai, q_sno, T)
end

"""
Advance the system over one full step `t_end` using the instantaneous nonlinear
tendency evaluated only at the initial state.

This is just one explicit Euler step, included for comparison.
"""
function integrate_bulk_microphysics_instantaneous_one_step(
    cm::BMT.Microphysics1Moment,
    mp::CMP.Microphysics1MParams,
    tps,
    ρ,
    T0,
    q_tot,
    q_lcl0,
    q_icl0,
    q_rai0,
    q_sno0,
    t_end;
    N_lcl = zero(ρ),
)
    FT = typeof(q_tot)

    times = FT[zero(FT), FT(t_end)]

    q_lcl = similar(times)
    q_icl = similar(times)
    q_rai = similar(times)
    q_sno = similar(times)
    T = similar(times)

    q_lcl[1] = q_lcl0
    q_icl[1] = q_icl0
    q_rai[1] = q_rai0
    q_sno[1] = q_sno0
    T[1] = T0

    rates = BMT.bulk_microphysics_tendencies(
        cm,
        mp,
        tps,
        ρ,
        T0,
        q_tot,
        q_lcl0,
        q_icl0,
        q_rai0,
        q_sno0,
        N_lcl,
    )

    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)

    q_lcl[2] = q_lcl0 + FT(t_end) * rates.dq_lcl_dt
    q_icl[2] = q_icl0 + FT(t_end) * rates.dq_icl_dt
    q_rai[2] = q_rai0 + FT(t_end) * rates.dq_rai_dt
    q_sno[2] = q_sno0 + FT(t_end) * rates.dq_sno_dt

    dT_dt =
        Lv_over_cp * (rates.dq_lcl_dt + rates.dq_rai_dt) +
        Ls_over_cp * (rates.dq_icl_dt + rates.dq_sno_dt)

    T[2] = T0 + FT(t_end) * dT_dt

    return (; times, q_lcl, q_icl, q_rai, q_sno, T)
end

"""
Plot the nonlinear reference solution together with linearized one-step solutions
using different numbers of internal substeps.

The nonlinear reference is integrated with many small explicit steps.
Each linearized curve uses the full `t_end` as one step, split internally into
`nsub` linearized implicit substeps.
Also shows the single explicit update from the instantaneous initial tendency as
a dashed line.
"""
function plot_bulk_microphysics_linearized_convergence(;
    FT = Float64,
    t_end = 60.0,
    dt_ref = 0.1,
    nsubs = [1, 2, 5, 10],
)
    tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
    mp = CMP.Microphysics1MParams(FT)
    cm = BMT.Microphysics1Moment()

    # Example initial condition
    ρ = FT(1.0)
    T0 = TDI.T_freeze(tps) + FT(5.0)
    q_lcl0 = FT(1e-3)
    q_icl0 = FT(5e-4)
    q_rai0 = FT(1e-3)
    q_sno0 = FT(5e-4)

    q_vap0 = FT(0.01)
    q_tot = q_vap0 + q_lcl0 + q_icl0 + q_rai0 + q_sno0

    ref = integrate_bulk_microphysics_reference(
        cm, mp, tps, ρ, T0, q_tot, q_lcl0, q_icl0, q_rai0, q_sno0, FT(t_end);
        dt = FT(dt_ref),
    )

    nolin = integrate_bulk_microphysics_reference(
        cm, mp, tps, ρ, T0, q_tot, q_lcl0, q_icl0, q_rai0, q_sno0, FT(t_end);
        dt = FT(t_end / nsubs[end]),
    )

    inst = integrate_bulk_microphysics_instantaneous_one_step(
        cm, mp, tps, ρ, T0, q_tot, q_lcl0, q_icl0, q_rai0, q_sno0, FT(t_end),
    )

    runs = Dict{Int, Any}()
    for nsub in nsubs
        runs[nsub] = integrate_bulk_microphysics_linearized_one_step(
            cm, mp, tps, ρ, T0, q_tot, q_lcl0, q_icl0, q_rai0, q_sno0, FT(t_end);
            nsub = nsub,
        )
    end

    fig = Figure(size = (1100, 800))

    ax1 = Axis(fig[1, 1], title = "Cloud liquid", xlabel = "time [s]", ylabel = "q_lcl [g/kg]")
    ax2 = Axis(fig[1, 2], title = "Cloud ice", xlabel = "time [s]", ylabel = "q_icl [g/kg]")
    ax3 = Axis(fig[2, 1], title = "Rain", xlabel = "time [s]", ylabel = "q_rai [g/kg]")
    ax4 = Axis(fig[2, 2], title = "Snow", xlabel = "time [s]", ylabel = "q_sno [g/kg]")

    # Reference
    lines!(ax1, ref.times, ref.q_lcl * 1e3, linewidth = 3, label = "nonlinear reference (dt = $dt_ref s)")
    lines!(ax2, ref.times, ref.q_icl * 1e3, linewidth = 3)
    lines!(ax3, ref.times, ref.q_rai * 1e3, linewidth = 3)
    lines!(ax4, ref.times, ref.q_sno * 1e3, linewidth = 3)

    # One-step explicit line from instantaneous initial tendency
    lines!(
        ax1,
        inst.times,
        inst.q_lcl * 1e3,
        linestyle = :dash,
        linewidth = 2,
        label = "instantaneous tendencies, nsubsteps = 1",
    )
    lines!(ax2, inst.times, inst.q_icl * 1e3, linestyle = :dash, linewidth = 2)
    lines!(ax3, inst.times, inst.q_rai * 1e3, linestyle = :dash, linewidth = 2)
    lines!(ax4, inst.times, inst.q_sno * 1e3, linestyle = :dash, linewidth = 2)

    # No linearization, nsubsteps = nsubs[end]
    nsub_end = nsubs[end]
    lines!(
        ax1,
        nolin.times,
        nolin.q_lcl * 1e3,
        linestyle = :dashdot,
        linewidth = 2,
        label = "instantaneous tendencies, nsubsteps = $nsub_end",
    )
    lines!(ax2, nolin.times, nolin.q_icl * 1e3, linestyle = :dashdot, linewidth = 2)
    lines!(ax3, nolin.times, nolin.q_rai * 1e3, linestyle = :dashdot, linewidth = 2)
    lines!(ax4, nolin.times, nolin.q_sno * 1e3, linestyle = :dashdot, linewidth = 2)

    q_lcl_max = maximum(ref.q_lcl)
    q_icl_max = maximum(ref.q_icl)
    q_rai_max = maximum(ref.q_rai)
    q_sno_max = maximum(ref.q_sno)
    # Linearized runs
    for nsub in nsubs
        sol = runs[nsub]
        label = "nsubsteps = $nsub"
        l1 = lines!(ax1, sol.times, sol.q_lcl * 1e3, label = label)
        l2 = lines!(ax2, sol.times, sol.q_icl * 1e3, label = label)
        l3 = lines!(ax3, sol.times, sol.q_rai * 1e3, label = label)
        l4 = lines!(ax4, sol.times, sol.q_sno * 1e3, label = label)
        scatter!(ax1, sol.times, sol.q_lcl * 1e3, color = l1.color, markersize = 10)
        scatter!(ax2, sol.times, sol.q_icl * 1e3, color = l2.color, markersize = 10)
        scatter!(ax3, sol.times, sol.q_rai * 1e3, color = l3.color, markersize = 10)
        scatter!(ax4, sol.times, sol.q_sno * 1e3, color = l4.color, markersize = 10)
        q_lcl_max = max(q_lcl_max, maximum(sol.q_lcl))
        q_icl_max = max(q_icl_max, maximum(sol.q_icl))
        q_rai_max = max(q_rai_max, maximum(sol.q_rai))
        q_sno_max = max(q_sno_max, maximum(sol.q_sno))
    end

    # Set limits based on reference (before plotting dashed line)
    ylims!(ax1, FT(-1e-2), 1.2 * q_lcl_max * 1e3)
    ylims!(ax2, FT(-1e-2), 1.2 * q_icl_max * 1e3)
    ylims!(ax3, FT(-1e-2), 1.2 * q_rai_max * 1e3)
    ylims!(ax4, FT(-1e-2), 1.2 * q_sno_max * 1e3)
    for ax in (ax1, ax2, ax3, ax4)
        ax.xgridvisible = true
        ax.ygridvisible = true

        ax.leftspinevisible = true
        ax.rightspinevisible = true
        ax.topspinevisible = true
        ax.bottomspinevisible = true
    end
    axislegend(ax1, position = :lt, framevisible = false)
    rowgap!(fig.layout, 10)
    colgap!(fig.layout, 10)

    return fig
end

colors = Makie.wong_colors()[[1, 2, 7, 3, 4, 5, 6]]
set_theme!(palette = (color = colors,))
fig = plot_bulk_microphysics_linearized_convergence(
    FT = Float64,
    t_end = 60.0,
    dt_ref = 0.1,
    nsubs = [1, 2, 5, 10],
)
save("bulk_microphysics_linearized_convergence.svg", fig)
