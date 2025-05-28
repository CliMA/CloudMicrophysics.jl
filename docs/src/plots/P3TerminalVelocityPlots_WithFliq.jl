import CairoMakie as MK

import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP

FT = Float64

p3 = CMP.ParametersP3(FT)

# Testing terminal velocity with liquid fraction

function get_values(
    p3,
    Chen2022::CMP.Chen2022VelType,
    L::FT,
    N::FT,
    F_liq::FT,
    ρ_a::FT,
    x_resolution::Int,
    y_resolution::Int,
) where {FT <: Real}
    F_rims = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    ρ_rs = range(FT(25), stop = FT(975), length = y_resolution)

    V_m = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)
    use_aspect_ratio = true

    for i in 1:x_resolution
        for j in 1:y_resolution
            F_rim = F_rims[i]
            ρ_r = ρ_rs[j]

            V_m[i, j] = P3.ice_terminal_velocity(
                p3,
                Chen2022,
                L,
                N,
                ρ_r,
                F_rim,
                F_liq,
                ρ_a,
                use_aspect_ratio,
            )[2]
            # get D_m in mm for plots
            D_m[i, j] = 1e3 * P3.D_m(p3, L, N, ρ_r, F_rim, F_liq)
        end
    end
    return (; F_rims, ρ_rs, V_m, D_m)
end

function make_axis(fig, row, col, title)
    return MK.Axis(
        fig[row, col],
        xlabel = "F_rim",
        ylabel = "ρ_r",
        title = title,
        width = 350,
        height = 350,
        limits = (0, 1.0, 25, 975),
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    )
end
#! format: off
function figure_2()

    Chen2022 = CMP.Chen2022VelType(FT)
    # density of air in kg/m^3
    ρ_a = FT(1.2) #FT(1.293)
    # small D_m
    L_s_0 = FT(0.0008)
    L_s_5 = FT(0.00091)
    L_s_9 = FT(0.00096)
    N_s = FT(1e6)
    # medium D_m
    L_m_0 = FT(0.22)
    L_m_5 = FT(0.555)
    L_m_9 = FT(0.635)
    N_m = FT(1e6)
    # large D_m
    L_l_0 = FT(0.7)
    L_l_5 = FT(2.6)
    L_l_9 = FT(3)
    N_l = FT(1e6)
    # get V_m and D_m
    xres = 100
    yres = 100

    (F_rims_0, ρ_rs_0, V_ms_0, D_ms_0) = get_values(p3, Chen2022, L_s_0, N_s, FT(0), ρ_a, xres, yres)
    (F_rimm_0, ρ_rm_0, V_mm_0, D_mm_0) = get_values(p3, Chen2022, L_m_0, N_m, FT(0), ρ_a, xres, yres)
    (F_riml_0, ρ_rl_0, V_ml_0, D_ml_0) = get_values(p3, Chen2022, L_l_0, N_l, FT(0), ρ_a, xres, yres)

    (F_rims_5, ρ_rs_5, V_ms_5, D_ms_5) = get_values(p3, Chen2022, L_s_5, N_s, FT(0.5), ρ_a, xres, yres)
    (F_rimm_5, ρ_rm_5, V_mm_5, D_mm_5) = get_values(p3, Chen2022, L_m_5, N_m, FT(0.5), ρ_a, xres, yres)
    (F_riml_5, ρ_rl_5, V_ml_5, D_ml_5) = get_values(p3, Chen2022, L_l_5, N_l, FT(0.5), ρ_a, xres, yres)

    (F_rims_9, ρ_rs_9, V_ms_9, D_ms_9) = get_values(p3, Chen2022, L_s_9, N_s, FT(0.9), ρ_a, xres, yres)
    (F_rimm_9, ρ_rm_9, V_mm_9, D_mm_9) = get_values(p3, Chen2022, L_m_9, N_m, FT(0.9), ρ_a, xres, yres)
    (F_riml_9, ρ_rl_9, V_ml_9, D_ml_9) = get_values(p3, Chen2022, L_l_9, N_l, FT(0.9), ρ_a, xres, yres)

    fig = MK.Figure()

    # Plot velocities as in Fig 2 in Morrison and Milbrandt 2015

    args = (color = :black, labels = true, levels = 4, linewidth = 1.5, labelsize = 18)

    # set the colorscale to be the same for all figures
    crange_small(len) = range(0.39, 0.62, length = len)
    crange_med(len) = range(3, 7.5, length = len)
    crange_large(len) = range(5, 10, length = len)

    ax1 = make_axis(fig, 1, 1, "Vₘ with small Dₘ, F_liq = 0")
    hm = MK.contourf!(ax1, F_rims_0, ρ_rs_0, V_ms_0, levels = crange_small(10), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax1, F_rims_0, ρ_rs_0, D_ms_0; args...)
    # MK.Colorbar(fig[2, 1], hm, vertical = false)

    ax2 = make_axis(fig, 1, 2, "Vₘ with medium Dₘ, F_liq = 0")
    hm = MK.contourf!(ax2, F_rimm_0, ρ_rm_0, V_mm_0, levels = crange_med(10), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax2, F_rimm_0, ρ_rm_0, D_mm_0; args...)
    # MK.Colorbar(fig[2, 2], hm, vertical = false)

    ax3 = make_axis(fig, 1, 3, "Vₘ with large Dₘ, F_liq = 0")
    hm = MK.contourf!(ax3, F_riml_0, ρ_rl_0, V_ml_0, levels = crange_large(10), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax3, F_riml_0, ρ_rl_0, D_ml_0; args...)
    # MK.Colorbar(fig[2, 3], hm, vertical = false)

    MK.linkxaxes!(ax1, ax2)
    MK.linkxaxes!(ax2, ax3)
    MK.linkyaxes!(ax1, ax2)
    MK.linkyaxes!(ax2, ax3)

    ## Plot F_liq = 0.5 as second row of comparisons

    ax4 = make_axis(fig, 3, 1, "Vₘ with small Dₘ, F_liq = 0.5")
    hm = MK.contourf!(ax4, F_rims_5, ρ_rs_5, V_ms_5, levels = crange_small(25), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax4, F_rims_5, ρ_rs_5, D_ms_5; args...)
    # MK.Colorbar(fig[4, 1], hm, vertical = false)

    ax5 = make_axis(fig, 3, 2, "Vₘ with medium Dₘ, F_liq = 0.5")
    hm = MK.contourf!(ax5, F_rimm_5, ρ_rm_5, V_mm_5, levels = crange_med(25), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax5, F_rimm_5, ρ_rm_5, D_mm_5; args...)
    # MK.Colorbar(fig[4, 2], hm, vertical = false)

    ax6 = make_axis(fig, 3, 3, "Vₘ with large Dₘ, F_liq = 0.5")
    hm = MK.contourf!(ax6, F_riml_5, ρ_rl_5, V_ml_5, levels = crange_large(30), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax6, F_riml_5, ρ_rl_5, D_ml_5; args...)
    # MK.Colorbar(fig[4, 3], hm, vertical = false)

    MK.linkxaxes!(ax1, ax4)
    MK.linkxaxes!(ax4, ax5)
    MK.linkxaxes!(ax5, ax6)
    MK.linkyaxes!(ax1, ax4)
    MK.linkyaxes!(ax4, ax5)
    MK.linkyaxes!(ax5, ax6)

    ## Plot F_liq = 0.9 as second row of comparisons

    ax7 = make_axis(fig, 5, 1, "Vₘ with small Dₘ, F_liq = 0.9")
    hm = MK.contourf!(ax7, F_rims_9, ρ_rs_9, V_ms_9, levels = crange_small(45), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax7, F_rims_9, ρ_rs_9, D_ms_9; args...)
    MK.Colorbar(fig[6, 1], hm, vertical = false)

    ax8 = make_axis(fig, 5, 2, "Vₘ with medium Dₘ, F_liq = 0.9")
    hm = MK.contourf!(ax8, F_rimm_9, ρ_rm_9, V_mm_9, levels = crange_med(45), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax8, F_rimm_9, ρ_rm_9, D_mm_9; args...)
    MK.Colorbar(fig[6, 2], hm, vertical = false)

    ax9 = make_axis(fig, 5, 3, "Vₘ with large Dₘ, F_liq = 0.9")
    hm = MK.contourf!(ax9, F_riml_9, ρ_rl_9, V_ml_9, levels = crange_large(45), extendlow = :auto, extendhigh = :auto)
    MK.contour!(ax9, F_riml_9, ρ_rl_9, D_ml_9; args...)
    MK.Colorbar(fig[6, 3], hm, vertical = false)

    MK.linkxaxes!(ax1, ax7)
    MK.linkxaxes!(ax7, ax8)
    MK.linkxaxes!(ax8, ax9)
    MK.linkyaxes!(ax1, ax7)
    MK.linkyaxes!(ax7, ax8)
    MK.linkyaxes!(ax8, ax9)


    MK.resize_to_layout!(fig)
    # MK.display(fig)
    MK.save("MorrisonandMilbrandtFig2_with_F_liq.svg", fig)
end

figure_2()
