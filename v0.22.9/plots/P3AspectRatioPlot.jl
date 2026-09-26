import CairoMakie as CMK
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

function define_axis(
    fig,
    row_range,
    col_range,
    title,
    ylabel,
    yticks,
    aspect;
    logscale = true,
)
    return CMK.Axis(
        fig[row_range, col_range],
        title = title,
        xlabel = CMK.L"$D$ (mm)",
        ylabel = ylabel,
        xscale = CMK.log10, #ifelse(logscale, CMK.log10, CMK.identity),
        yscale = ifelse(logscale, CMK.log10, CMK.identity),
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
function p3_aspect()

    D_range = range(3e-5, stop = 1e-2, length = Int(1e4))
    logocolors = CMK.Colors.JULIA_LOGO_COLORS
    cl = [logocolors.blue, logocolors.green, logocolors.red, logocolors.purple]
    lw = 3

    fig = CMK.Figure()

    # define plot axis
    #[row, column]
    ax1 = define_axis(fig, 1:7,  1:9,   CMK.L"ϕᵢ(D) regime for $ρ_r = 400 kg m^{-3}$", CMK.L"$ϕᵢ$ (-)", [0.0, 0.5, 1.0], 1.9, logscale = false)
    ax3 = define_axis(fig, 1:7,  10:18, CMK.L"ϕᵢ(D) regime for $F_rim = 0.4$",          CMK.L"$ϕᵢ$ (-)", [0.0, 0.5, 1.0], 1.8, logscale = false)

    # Get thresholds
    sol4_0 = P3.thresholds(p3, 400.0, 0.0)
    sol4_5 = P3.thresholds(p3, 400.0, 0.5)
    sol4_8 = P3.thresholds(p3, 400.0, 0.8)
    # ϕᵢ(D)
    fig1_0 = CMK.lines!(ax1, D_range * 1e3, [P3.ϕᵢ(p3, D, 0.0, sol4_0) for D in D_range], color = cl[1], linewidth = lw)
    fig1_5 = CMK.lines!(ax1, D_range * 1e3, [P3.ϕᵢ(p3, D, 0.5, sol4_5) for D in D_range], color = cl[2], linewidth = lw)
    fig1_8 = CMK.lines!(ax1, D_range * 1e3, [P3.ϕᵢ(p3, D, 0.8, sol4_8) for D in D_range], color = cl[3], linewidth = lw)
    for ax in [ax1]
        global d_tha  = CMK.vlines!(ax, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
        global d_cr_5 = CMK.vlines!(ax, sol4_5[1]          * 1e3, linestyle = :dot,  color = cl[2], linewidth = lw)
        global d_cr_8 = CMK.vlines!(ax, sol4_8[1]          * 1e3, linestyle = :dot,  color = cl[3], linewidth = lw)
        global d_gr_5 = CMK.vlines!(ax, sol4_5[2]          * 1e3, linestyle = :dash, color = cl[2], linewidth = lw)
        global d_gr_8 = CMK.vlines!(ax, sol4_8[2]          * 1e3, linestyle = :dash, color = cl[3], linewidth = lw)
    end

    # get thresholds
    sol_2 = P3.thresholds(p3, 200.0, 0.4)
    sol_4 = P3.thresholds(p3, 400.0, 0.4)
    sol_8 = P3.thresholds(p3, 800.0, 0.4)
    # ϕᵢ(D)
    fig2_200 = CMK.lines!(ax3, D_range * 1e3, [P3.ϕᵢ(p3, D, 0.4, sol_2) for D in D_range], color = cl[1], linewidth = lw)
    fig2_400 = CMK.lines!(ax3, D_range * 1e3, [P3.ϕᵢ(p3, D, 0.4, sol_4) for D in D_range], color = cl[2], linewidth = lw)
    fig2_800 = CMK.lines!(ax3, D_range * 1e3, [P3.ϕᵢ(p3, D, 0.4, sol_8) for D in D_range], color = cl[3], linewidth = lw)
    # plot verical lines
    for ax in [ax3]
        global d_thb    = CMK.vlines!(ax, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
        global d_cr_200 = CMK.vlines!(ax, sol_2[1] * 1e3,           linestyle = :dot,  color = cl[1], linewidth = lw)
        global d_cr_400 = CMK.vlines!(ax, sol_4[1] * 1e3,           linestyle = :dot,  color = cl[2], linewidth = lw)
        global d_cr_800 = CMK.vlines!(ax, sol_8[1] * 1e3,           linestyle = :dot,  color = cl[3], linewidth = lw)
        global d_gr_200 = CMK.vlines!(ax, sol_2[2] * 1e3,           linestyle = :dash, color = cl[1], linewidth = lw)
        global d_gr_400 = CMK.vlines!(ax, sol_4[2] * 1e3,           linestyle = :dash, color = cl[2], linewidth = lw)
        global d_gr_800 = CMK.vlines!(ax, sol_8[2] * 1e3,           linestyle = :dash, color = cl[3], linewidth = lw)
    end
    # add legend
    CMK.Legend(fig[8:9, 1], [fig1_0, fig1_5, fig1_8], [CMK.L"$F_{rim} = 0.0$", CMK.L"$F_{rim} = 0.5$", CMK.L"$F_{rim} = 0.8$"], framevisible = false)
    CMK.Legend(fig[8:9, 3], [d_tha], [CMK.L"$D_{th}$"], framevisible = false)
    CMK.Legend(fig[8:9, 7], [d_cr_5, d_cr_8], [CMK.L"$D_{cr}$ for $F_{rim} = 0.5$", CMK.L"$D_{cr}$ for $F_{rim} = 0.8$"], framevisible = false)
    CMK.Legend(fig[8:9, 5], [d_gr_5, d_gr_8], [CMK.L"$D_{gr}$ for $F_{rim} = 0.5$", CMK.L"$D_{gr}$ for $F_{rim} = 0.8$"], framevisible = false)

    CMK.Legend(fig[8:9, 13], [fig2_200, fig2_400, fig2_800],    [CMK.L"$\rho_{r} = 200.0 kg m^{-3}$", CMK.L"$\rho_{r} = 400.0 kg m^{-3}$", CMK.L"$\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    CMK.Legend(fig[8:9, 14], [d_thb],                           [CMK.L"$D_{th}$"], framevisible = false)
    CMK.Legend(fig[8:9, 17], [d_cr_200, d_cr_400, d_cr_800],    [CMK.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", CMK.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", CMK.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    CMK.Legend(fig[8:9, 16], [d_gr_200, d_gr_400, d_gr_800],    [CMK.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", CMK.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", CMK.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)

    CMK.resize_to_layout!(fig)
    CMK.save("P3Scheme_aspect_ratio.svg", fig)
end
#! format: on

p3_aspect()
