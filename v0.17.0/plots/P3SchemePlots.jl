import CairoMakie as Plt
import CloudMicrophysics as CM
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

"""
    A_(p3, D)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension

Returns particle projected area as a function of size for different particle regimes
"""
# for spherical particles
A_s(D::FT) = FT(π) / 4 * D^2
# for nonspherical particles
A_ns(p3::PSP3, D::FT) = p3.γ * D^p3.σ
# partially rimed ice
A_r(p3::PSP3, F_r::FT, D::FT) = F_r * A_s(D) + (1 - F_r) * A_ns(p3, D)

"""
    area(p3, D, F_r, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_r - rime mass fraction (q_rim/q_i)
 - th - P3Scheme nonlinear solve output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns area(D), used to create figures for the documentation.
"""
function area(
    p3::PSP3,
    D::FT,
    F_r::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT <: Real}
    # Area regime:
    if P3.D_th_helper(p3) > D
        return A_s(D)                      # small spherical ice
    end
    if F_r == 0
        return A_ns(p3, D)                 # large nonspherical unrimed ice
    end
    if th.D_gr > D >= P3.D_th_helper(p3)
        return A_ns(p3, D)                 # dense nonspherical ice
    end
    if th.D_cr > D >= th.D_gr
        return A_s(D)                      # graupel
    end
    if D >= th.D_cr
        return A_r(p3, F_r, D)             # partially rimed ice
    end

end

# Plot Fig. 1a and 1b from Morrison and Milbrandt 2015.
function p3_mass_plot()

    D_range = range(3e-5, stop = 1e-2, length = Int(1e4))
    logocolors = Plt.Colors.JULIA_LOGO_COLORS
    cl = [logocolors.blue, logocolors.green, logocolors.red, logocolors.purple]
    lw = 3

    fig1_a = Plt.Figure()
    ax1_a = Plt.Axis(
        fig1_a[1:7, 1:8],
        title = Plt.L"m(D) regime for $ρ_r = 400 kg m^{-3}$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$m$ (kg)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 1.75,
        limits = ((0.02, 10.0), nothing),
    )

    sol_5 = P3.thresholds(p3, 400.0, 0.5)
    sol_8 = P3.thresholds(p3, 400.0, 0.8)

    #! format: off
    fig1_a_0 = Plt.lines!(ax1_a, D_range * 1e3, [P3.p3_mass(p3, D, 0.0       ) for D in D_range], color = cl[1], linewidth = lw)
    fig1_a_5 = Plt.lines!(ax1_a, D_range * 1e3, [P3.p3_mass(p3, D, 0.5, sol_5) for D in D_range], color = cl[2], linewidth = lw)
    fig1_a_8 = Plt.lines!(ax1_a, D_range * 1e3, [P3.p3_mass(p3, D, 0.8, sol_8) for D in D_range], color = cl[3], linewidth = lw)

    d_tha  = Plt.vlines!(ax1_a, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
    d_cr_5 = Plt.vlines!(ax1_a, sol_5[1]           * 1e3, linestyle = :dot,  color = cl[2], linewidth = lw)
    d_cr_8 = Plt.vlines!(ax1_a, sol_8[1]           * 1e3, linestyle = :dot,  color = cl[3], linewidth = lw)
    d_gr_5 = Plt.vlines!(ax1_a, sol_5[2]           * 1e3, linestyle = :dash, color = cl[2], linewidth = lw)
    d_gr_8 = Plt.vlines!(ax1_a, sol_8[2]           * 1e3, linestyle = :dash, color = cl[3], linewidth = lw)

    leg1_a = Plt.Legend(fig1_a[8:9, 1], [fig1_a_0, fig1_a_5, fig1_a_8], [Plt.L"$F_{r} = 0.0$", Plt.L"$F_{r} = 0.5$", Plt.L"$F_{r} = 0.8$"], framevisible = false)
    leg1_a_dth = Plt.Legend(fig1_a[8:9, 3], [d_tha], [Plt.L"$D_{th}$"], framevisible = false)
    leg1_a_dcr = Plt.Legend(fig1_a[8:9, 7], [d_cr_5, d_cr_8], [Plt.L"$D_{cr}$ for $F_{r} = 0.5$", Plt.L"$D_{cr}$ for $F_{r} = 0.8$"], framevisible = false)
    leg1_a_dgr = Plt.Legend(fig1_a[8:9, 5], [d_gr_5, d_gr_8], [Plt.L"$D_{gr}$ for $F_{r} = 0.5$", Plt.L"$D_{gr}$ for $F_{r} = 0.8$"], framevisible = false)

    #! format: on
    Plt.save("MorrisonandMilbrandtFig1a.svg", fig1_a)

    fig1_b = Plt.Figure()
    ax1_b = Plt.Axis(
        fig1_b[1:10, 1:11],
        title = Plt.L"m(D) regime for $F_r = 0.95$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$m$ (kg)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 1.67,
        limits = ((0.02, 10.0), nothing),
    )

    sol_2 = P3.thresholds(p3, 200.0, 0.95)
    sol_4 = P3.thresholds(p3, 400.0, 0.95)
    sol_8 = P3.thresholds(p3, 800.0, 0.95)

    #! format: off
    fig1_b200 = Plt.lines!(ax1_b, D_range * 1e3, [P3.p3_mass(p3, D, 0.95, sol_2) for D in D_range], color = cl[1], linewidth = lw)
    fig1_b400 = Plt.lines!(ax1_b, D_range * 1e3, [P3.p3_mass(p3, D, 0.95, sol_4) for D in D_range], color = cl[2], linewidth = lw)
    fig1_b800 = Plt.lines!(ax1_b, D_range * 1e3, [P3.p3_mass(p3, D, 0.95, sol_8) for D in D_range], color = cl[3], linewidth = lw)

    d_thb    = Plt.vlines!(ax1_b, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
    d_cr_200 = Plt.vlines!(ax1_b, sol_2[1] * 1e3,           linestyle = :dot,  color = cl[1], linewidth = lw)
    d_cr_400 = Plt.vlines!(ax1_b, sol_4[1] * 1e3,           linestyle = :dot,  color = cl[2], linewidth = lw)
    d_cr_800 = Plt.vlines!(ax1_b, sol_8[1] * 1e3,           linestyle = :dot,  color = cl[3], linewidth = lw)
    d_gr_200 = Plt.vlines!(ax1_b, sol_2[2] * 1e3,           linestyle = :dash, color = cl[1], linewidth = lw)
    d_gr_400 = Plt.vlines!(ax1_b, sol_4[2] * 1e3,           linestyle = :dash, color = cl[2], linewidth = lw)
    d_gr_800 = Plt.vlines!(ax1_b, sol_8[2] * 1e3,           linestyle = :dash, color = cl[3], linewidth = lw)

    leg1_b     = Plt.Legend(fig1_b[11:12, 4], [fig1_b200, fig1_b400, fig1_b800], [Plt.L"$\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    leg1_b_dth = Plt.Legend(fig1_b[11:12, 5], [d_thb],                           [Plt.L"$D_{th}$"], framevisible = false)
    leg1_b_dcr = Plt.Legend(fig1_b[11:12, 8], [d_cr_200, d_cr_400, d_cr_800],    [Plt.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    leg1_b_dgr = Plt.Legend(fig1_b[11:12, 6], [d_gr_200, d_gr_400, d_gr_800],    [Plt.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)

    #! format: on
    Plt.save("MorrisonandMilbrandtFig1b.svg", fig1_b)
end

# Plot area(size) for different regimes for the documentation.
function p3_area_plot()

    D_range = range(3e-5, stop = 1e-2, length = Int(1e4))
    logocolors = Plt.Colors.JULIA_LOGO_COLORS
    cl = [logocolors.blue, logocolors.green, logocolors.red, logocolors.purple]
    lw = 3

    fig = Plt.Figure()
    ax = Plt.Axis(
        fig[1:7, 1:8],
        title = Plt.L"A(D) regime for $ρ_r = 400 kg m^{-3}$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$A$ ($m^2$)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 1.75,
        limits = ((0.02, 10.0), nothing),
    )
    sol_5 = P3.thresholds(p3, 400.0, 0.5)
    sol_8 = P3.thresholds(p3, 400.0, 0.8)

    #! format: off
    fig_0 = Plt.lines!(ax, D_range * 1e3, [area(p3, D, 0.0       ) for D in D_range], color = cl[1], linewidth = lw)
    fig_5 = Plt.lines!(ax, D_range * 1e3, [area(p3, D, 0.5, sol_5) for D in D_range], color = cl[2], linewidth = lw)
    fig_8 = Plt.lines!(ax, D_range * 1e3, [area(p3, D, 0.8, sol_8) for D in D_range], color = cl[3], linewidth = lw)

    d_th =   Plt.vlines!(ax, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
    d_cr_5 = Plt.vlines!(ax, sol_5[1] * 1e3,           linestyle = :dot,  color = cl[2], linewidth = lw)
    d_cr_8 = Plt.vlines!(ax, sol_8[1] * 1e3,           linestyle = :dot,  color = cl[3], linewidth = lw)
    d_gr_5 = Plt.vlines!(ax, sol_5[2] * 1e3,           linestyle = :dash, color = cl[2], linewidth = lw)
    d_gr_8 = Plt.vlines!(ax, sol_8[2] * 1e3,           linestyle = :dash, color = cl[3], linewidth = lw)

    leg     = Plt.Legend(fig[8:9, 1], [fig_0, fig_5, fig_8], [Plt.L"$F_{r} = 0.0$", Plt.L"$F_{r} = 0.5$", Plt.L"$F_{r} = 0.8$"], framevisible = false)
    leg_dth = Plt.Legend(fig[8:9, 3], [d_th],                [Plt.L"$D_{th}$"], framevisible = false)
    leg_dcr = Plt.Legend(fig[8:9, 7], [d_cr_5, d_cr_8],      [Plt.L"$D_{cr}$ for $F_{r} = 0.5$", Plt.L"$D_{cr}$ for $F_{r} = 0.8$"], framevisible = false)
    leg_dgr = Plt.Legend(fig[8:9, 5], [d_gr_5, d_gr_8],      [Plt.L"$D_{gr}$ for $F_{r} = 0.5$", Plt.L"$D_{gr}$ for $F_{r} = 0.8$"], framevisible = false)
    #! format:on
    Plt.save("P3Scheme_Area_1.svg", fig)

    fig = Plt.Figure()
    ax = Plt.Axis(
        fig[1:10, 1:11],
        title = Plt.L"A(D) regime for $F_r = 0.95$",
        xlabel = Plt.L"$D$ (mm)",
        ylabel = Plt.L"$A$ ($m^2$)",
        xscale = Plt.log10,
        yscale = Plt.log10,
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(3),
        xminorticksvisible = true,
        xminorticks = Plt.IntervalsBetween(5),
        xticks = [0.01, 0.1, 1, 10],
        aspect = 1.67,
        limits = ((0.02, 10.0), nothing),
    )
    sol_2 = P3.thresholds(p3, 200.0, 0.95)
    sol_4 = P3.thresholds(p3, 400.0, 0.95)
    sol_8 = P3.thresholds(p3, 800.0, 0.95)

    #! format: off
    fig200 = Plt.lines!(ax, D_range * 1e3, [area(p3, D, 0.5, sol_2) for D in D_range], color = cl[1], linewidth = lw)
    fig400 = Plt.lines!(ax, D_range * 1e3, [area(p3, D, 0.5, sol_4) for D in D_range], color = cl[2], linewidth = lw)
    fig800 = Plt.lines!(ax, D_range * 1e3, [area(p3, D, 0.5, sol_8) for D in D_range], color = cl[3], linewidth = lw)

    d_thb    = Plt.vlines!(ax, P3.D_th_helper(p3) * 1e3, linestyle = :dash, color = cl[4], linewidth = lw)
    d_cr_200 = Plt.vlines!(ax, sol_2[1] * 1e3,           linestyle = :dot,  color = cl[1], linewidth = lw)
    d_cr_400 = Plt.vlines!(ax, sol_4[1] * 1e3,           linestyle = :dot,  color = cl[2], linewidth = lw)
    d_cr_800 = Plt.vlines!(ax, sol_8[1] * 1e3,           linestyle = :dot,  color = cl[3], linewidth = lw)
    d_gr_200 = Plt.vlines!(ax, sol_2[2] * 1e3,           linestyle = :dash, color = cl[1], linewidth = lw)
    d_gr_400 = Plt.vlines!(ax, sol_4[2] * 1e3,           linestyle = :dash, color = cl[2], linewidth = lw)
    d_gr_800 = Plt.vlines!(ax, sol_8[2] * 1e3,           linestyle = :dash, color = cl[3], linewidth = lw)

    leg1_b     = Plt.Legend(fig[11:12, 3], [fig200, fig400, fig800],       [Plt.L"$\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    leg1_b_dth = Plt.Legend(fig[11:12, 4], [d_thb],                        [Plt.L"$D_{th}$"], framevisible = false)
    leg1_b_dcr = Plt.Legend(fig[11:12, 7], [d_cr_200, d_cr_400, d_cr_800], [Plt.L"$D_{cr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$D_{cr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    leg1_b_dgr = Plt.Legend(fig[11:12, 5], [d_gr_200, d_gr_400, d_gr_800], [Plt.L"$D_{gr}$ for $\rho_{r} = 200.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 400.0 kg m^{-3}$", Plt.L"$D_{gr}$ for $\rho_{r} = 800.0 kg m^{-3}$",], framevisible = false)
    #! format: on
    Plt.save("P3Scheme_Area_2.svg", fig)
end

p3_mass_plot()
p3_area_plot()
