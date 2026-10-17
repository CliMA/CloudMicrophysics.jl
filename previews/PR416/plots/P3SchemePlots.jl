import CairoMakie as CMK
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)
F_liq = FT(0)



function define_axis(fig, row_range, col_range, title, ylabel, yticks, aspect)
    return CMK.Axis(
        fig[row_range, col_range],
        title = title,
        xlabel = CMK.L"$D$ (mm)",
        ylabel = ylabel,
        xscale = CMK.log10,
        yscale = CMK.log10,
        xminorticksvisible = true,
        xminorticks = CMK.IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = CMK.IntervalsBetween(3),
        xticks = [0.01, 0.1, 1, 10],
        yticks = yticks,
        aspect = aspect,
        limits = ((0.02, 10.0), nothing),
    )
end

#! format: off
function p3_relations_plot()

    D_range = range(3e-5, stop = 1e-2, length = Int(1e4))
    logocolors = CMK.Colors.JULIA_LOGO_COLORS
    cl = [logocolors.blue, logocolors.green, logocolors.red, logocolors.purple]
    lw = 3

    fig = CMK.Figure(size=(1400, 900))

    # define plot axis
    #[row, column]
    ax1 = define_axis(fig, 1:7,  1:9,   CMK.L"m(D) regime for $ρ_r = 400 kg m^{-3}$", CMK.L"$m$ (kg)", [1e-10, 1e-8, 1e-6], 1.9)
    ax3 = define_axis(fig, 1:7,  10:18, CMK.L"m(D) regime for $F_r = 0.95$",          CMK.L"$m$ (kg)", [1e-10, 1e-8, 1e-6], 1.8)
    ax2 = define_axis(fig, 8:15, 1:9,   CMK.L"A(D) regime for $ρ_r = 400 kg m^{-3}$", CMK.L"$A$ ($m^2$)", [1e-8, 1e-6, 1e-4], 1.7)
    ax4 = define_axis(fig, 8:15, 10:18, CMK.L"A(D) regime for $F_r = 0.95$", CMK.L"$A$ ($m^2$)", [1e-8, 1e-6, 1e-4], 1.6)

    # Get thresholds
    sol4_5 = P3.thresholds(p3, 400.0, 0.5)
    sol4_8 = P3.thresholds(p3, 400.0, 0.8)
    # m(D)
    fig1_0 = CMK.lines!(ax1, D_range * 1e3, [P3.p3_mass(p3, D, 0.0, F_liq        ) for D in D_range], color = cl[1], linewidth = lw)
    fig1_5 = CMK.lines!(ax1, D_range * 1e3, [P3.p3_mass(p3, D, 0.5, F_liq, sol4_5) for D in D_range], color = cl[2], linewidth = lw)
    fig1_8 = CMK.lines!(ax1, D_range * 1e3, [P3.p3_mass(p3, D, 0.8, F_liq, sol4_8) for D in D_range], color = cl[3], linewidth = lw)
    # a(D)
    fig2_0 = CMK.lines!(ax2, D_range * 1e3, [P3.p3_area(p3, D, 0.0, F_liq        ) for D in D_range], color = cl[1], linewidth = lw)
    fig2_5 = CMK.lines!(ax2, D_range * 1e3, [P3.p3_area(p3, D, 0.5, F_liq, sol4_5) for D in D_range], color = cl[2], linewidth = lw)
    fig2_8 = CMK.lines!(ax2, D_range * 1e3, [P3.p3_area(p3, D, 0.8, F_liq, sol4_8) for D in D_range], color = cl[3], linewidth = lw)
    # plot verical lines
    for ax in [ax1, ax2]
        global d_tha  = CMK.vlines!(ax, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
        global d_cr_5 = CMK.vlines!(ax, sol4_5[1]          * 1e3, linestyle = :dot,  color = cl[2], linewidth = lw)
        global d_cr_8 = CMK.vlines!(ax, sol4_8[1]          * 1e3, linestyle = :dot,  color = cl[3], linewidth = lw)
        global d_gr_5 = CMK.vlines!(ax, sol4_5[2]          * 1e3, linestyle = :dash, color = cl[2], linewidth = lw)
        global d_gr_8 = CMK.vlines!(ax, sol4_8[2]          * 1e3, linestyle = :dash, color = cl[3], linewidth = lw)
    end

    # get thresholds
    sol_2 = P3.thresholds(p3, 200.0, 0.95)
    sol_4 = P3.thresholds(p3, 400.0, 0.95)
    sol_8 = P3.thresholds(p3, 800.0, 0.95)
    # m(D)
    fig3_200 = CMK.lines!(ax3, D_range * 1e3, [P3.p3_mass(p3, D, 0.95, F_liq, sol_2) for D in D_range], color = cl[1], linewidth = lw)
    fig3_400 = CMK.lines!(ax3, D_range * 1e3, [P3.p3_mass(p3, D, 0.95, F_liq, sol_4) for D in D_range], color = cl[2], linewidth = lw)
    fig3_800 = CMK.lines!(ax3, D_range * 1e3, [P3.p3_mass(p3, D, 0.95, F_liq, sol_8) for D in D_range], color = cl[3], linewidth = lw)
    # a(D)
    fig3_200 = CMK.lines!(ax4, D_range * 1e3, [P3.p3_area(p3, D, 0.5, F_liq, sol_2) for D in D_range], color = cl[1], linewidth = lw)
    fig3_400 = CMK.lines!(ax4, D_range * 1e3, [P3.p3_area(p3, D, 0.5, F_liq, sol_4) for D in D_range], color = cl[2], linewidth = lw)
    fig3_800 = CMK.lines!(ax4, D_range * 1e3, [P3.p3_area(p3, D, 0.5, F_liq, sol_8) for D in D_range], color = cl[3], linewidth = lw)
    # plot verical lines
    for ax in [ax3, ax4]
        global d_thb    = CMK.vlines!(ax, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
        global d_cr_200 = CMK.vlines!(ax, sol_2[1] * 1e3,           linestyle = :dot,  color = cl[1], linewidth = lw)
        global d_cr_400 = CMK.vlines!(ax, sol_4[1] * 1e3,           linestyle = :dot,  color = cl[2], linewidth = lw)
        global d_cr_800 = CMK.vlines!(ax, sol_8[1] * 1e3,           linestyle = :dot,  color = cl[3], linewidth = lw)
        global d_gr_200 = CMK.vlines!(ax, sol_2[2] * 1e3,           linestyle = :dash, color = cl[1], linewidth = lw)
        global d_gr_400 = CMK.vlines!(ax, sol_4[2] * 1e3,           linestyle = :dash, color = cl[2], linewidth = lw)
        global d_gr_800 = CMK.vlines!(ax, sol_8[2] * 1e3,           linestyle = :dash, color = cl[3], linewidth = lw)
    end
    # add legend
    CMK.Legend(fig[16:17, 1], [fig1_0, fig1_5, fig1_8], [CMK.L"$F_{r} = 0.0$", CMK.L"$F_{r} = 0.5$", CMK.L"$F_{r} = 0.8$"], framevisible = false)
    CMK.Legend(fig[16:17, 3], [d_tha], [CMK.L"$D_{th}$"], framevisible = false)
    CMK.Legend(fig[16:17, 7], [d_cr_5, d_cr_8], [CMK.L"$D_{cr}$ for $F_{r} = 0.5$", CMK.L"$D_{cr}$ for $F_{r} = 0.8$"], framevisible = false)
    CMK.Legend(fig[16:17, 5], [d_gr_5, d_gr_8], [CMK.L"$D_{gr}$ for $F_{r} = 0.5$", CMK.L"$D_{gr}$ for $F_{r} = 0.8$"], framevisible = false)

    CMK.Legend(fig[16:17, 13], [fig3_200, fig3_400, fig3_800],    [CMK.L"$\rho_{r} = 200.0 kg m^{-3}$", CMK.L"$\rho_{r} = 400.0 kg m^{-3}$", CMK.L"$\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    CMK.Legend(fig[16:17, 14], [d_thb],                           [CMK.L"$D_{th}$"], framevisible = false)
    CMK.Legend(fig[16:17, 17], [d_cr_200, d_cr_400, d_cr_800],    [CMK.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", CMK.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", CMK.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    CMK.Legend(fig[16:17, 16], [d_gr_200, d_gr_400, d_gr_800],    [CMK.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", CMK.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", CMK.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)

    CMK.resize_to_layout!(fig)
    CMK.save("P3Scheme_relations.svg", fig)
end
#! format: on

p3_relations_plot()
