import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)

# Testing terminal velocity with liquid fraction

function get_values(
    p3::PSP3,
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
    aspect_ratio = false

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
                aspect_ratio,
            )[2]
            # get D_m in mm for plots
            D_m[i, j] = 1e3 * P3.D_m(p3, L, N, ρ_r, F_rim, F_liq)
        end
    end
    return (; F_rims, ρ_rs, V_m, D_m)
end

function make_axis_top(fig, col, title)
    return Plt.Axis(
        fig[1, col],
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
function make_axis_middle(fig, col, title)
    return Plt.Axis(
        fig[3, col],
        height = 350,
        width = 350,
        xlabel = "F_rim",
        ylabel = "ρ_r",
        title = title,
        limits = (0, 1.0, 25, 975),
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    )
end
function make_axis_bottom(fig, col, title)
    return Plt.Axis(
        fig[5, col],
        height = 350,
        width = 350,
        xlabel = "F_rim",
        ylabel = "ρ_r",
        title = title,
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
    L_s_33 = FT(0.00091)
    L_s_67 = FT(0.00096)
    N_s = FT(1e6)
    # medium D_m
    L_m_0 = FT(0.22)
    L_m_33 = FT(0.555)
    L_m_67 = FT(0.635)
    N_m = FT(1e6)
    # large D_m
    L_l_0 = FT(0.7)
    L_l_33 = FT(2.6)
    L_l_67 = FT(3)
    N_l = FT(1e6)
    # get V_m and D_m
    xres = 100
    yres = 100

    (F_rims_0, ρ_rs_0, V_ms_0, D_ms_0) = get_values(p3, Chen2022, L_s_0, N_s, FT(0), ρ_a, xres, yres)
    (F_rimm_0, ρ_rm_0, V_mm_0, D_mm_0) = get_values(p3, Chen2022, L_m_0, N_m, FT(0), ρ_a, xres, yres)
    (F_riml_0, ρ_rl_0, V_ml_0, D_ml_0) = get_values(p3, Chen2022, L_l_0, N_l, FT(0), ρ_a, xres, yres)

    (F_rims_33, ρ_rs_33, V_ms_33, D_ms_33) = get_values(p3, Chen2022, L_s_33, N_s, FT(0.33), ρ_a, xres, yres)
    (F_rimm_33, ρ_rm_33, V_mm_33, D_mm_33) = get_values(p3, Chen2022, L_m_33, N_m, FT(0.33), ρ_a, xres, yres)
    (F_riml_33, ρ_rl_33, V_ml_33, D_ml_33) = get_values(p3, Chen2022, L_l_33, N_l, FT(0.33), ρ_a, xres, yres)

    (F_rims_67, ρ_rs_67, V_ms_67, D_ms_67) = get_values(p3, Chen2022, L_s_67, N_s, FT(0.67), ρ_a, xres, yres)
    (F_rimm_67, ρ_rm_67, V_mm_67, D_mm_67) = get_values(p3, Chen2022, L_m_67, N_m, FT(0.67), ρ_a, xres, yres)
    (F_riml_67, ρ_rl_67, V_ml_67, D_ml_67) = get_values(p3, Chen2022, L_l_67, N_l, FT(0.67), ρ_a, xres, yres)

    fig = Plt.Figure()

    # Plot velocities as in Fig 2 in Morrison and Milbrandt 2015

    args = (color = :black, labels = true, levels = 4, linewidth = 1.5, labelsize = 18)

    # set the colorscale to be the same for all figures
    crange_small = range(0.39, 0.6, length = 20)
    crange_med(len) = range(4, 8, length = len)
    crange_large(len) = range(5, 10, length = len)

    ax1 = make_axis_top(fig, 1, "Small Dₘ, F_liq = 0")
    hm = Plt.contourf!(ax1, F_rims_0, ρ_rs_0, V_ms_0, levels = crange_small, extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax1, F_rims_0, ρ_rs_0, D_ms_0; args...)
    Plt.Colorbar(fig[6, 1], hm, vertical = false)

    ax2 = make_axis_top(fig, 2, "Medium Dₘ, F_liq = 0")
    hm = Plt.contourf!(ax2, F_rimm_0, ρ_rm_0, V_mm_0, levels = crange_med(10), extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax2, F_rimm_0, ρ_rm_0, D_mm_0; args...)
    Plt.Colorbar(fig[6, 2], hm, vertical = false)

    ax3 = make_axis_top(fig, 3, "Large Dₘ, F_liq = 0")
    hm = Plt.contourf!(ax3, F_riml_0, ρ_rl_0, V_ml_0, levels = crange_large(10), extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax3, F_riml_0, ρ_rl_0, D_ml_0; args...)
    Plt.Colorbar(fig[6, 3], hm, vertical = false)

    Plt.linkxaxes!(ax1, ax2)
    Plt.linkxaxes!(ax2, ax3)
    Plt.linkyaxes!(ax1, ax2)
    Plt.linkyaxes!(ax2, ax3)

    ## Plot F_liq = 0.33 as second row of comparisons

    ax4 = make_axis_middle(fig, 1, "Small Dₘ, F_liq = 0.33")
    hm = Plt.contourf!(ax4, F_rims_33, ρ_rs_33, V_ms_33, levels = crange_small, extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax4, F_rims_33, ρ_rs_33, D_ms_33; args...)

    ax5 = make_axis_middle(fig, 2, "Medium Dₘ, F_liq = 0.33")
    hm = Plt.contourf!(ax5, F_rimm_33, ρ_rm_33, V_mm_33, levels = crange_med(20), extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax5, F_rimm_33, ρ_rm_33, D_mm_33; args...)

    ax6 = make_axis_middle(fig, 3, "Large Dₘ, F_liq = 0.33")
    hm = Plt.contourf!(ax6, F_riml_33, ρ_rl_33, V_ml_33, levels = crange_large(30), extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax6, F_riml_33, ρ_rl_33, D_ml_33; args...)

    Plt.linkxaxes!(ax1, ax4)
    Plt.linkxaxes!(ax4, ax5)
    Plt.linkxaxes!(ax5, ax6)
    Plt.linkyaxes!(ax1, ax4)
    Plt.linkyaxes!(ax4, ax5)
    Plt.linkyaxes!(ax5, ax6)

    ## Plot F_liq = 0.67 as second row of comparisons

    ax7 = make_axis_bottom(fig, 1, "Small Dₘ, F_liq = 0.67")
    hm = Plt.contourf!(ax7, F_rims_67, ρ_rs_67, V_ms_67, levels = crange_small, extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax7, F_rims_67, ρ_rs_67, D_ms_67; args...)

    ax8 = make_axis_bottom(fig, 2, "Medium Dₘ, F_liq = 0.67")
    hm = Plt.contourf!(ax8, F_rimm_67, ρ_rm_67, V_mm_67, levels = crange_med(45), extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax8, F_rimm_67, ρ_rm_67, D_mm_67; args...)

    ax9 = make_axis_bottom(fig, 3, "Large Dₘ, F_liq = 0.67")
    hm = Plt.contourf!(ax9, F_riml_67, ρ_rl_67, V_ml_67, levels = crange_large(45), extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax9, F_riml_67, ρ_rl_67, D_ml_67; args...)

    Plt.linkxaxes!(ax1, ax7)
    Plt.linkxaxes!(ax7, ax8)
    Plt.linkxaxes!(ax8, ax9)
    Plt.linkyaxes!(ax1, ax7)
    Plt.linkyaxes!(ax7, ax8)
    Plt.linkyaxes!(ax8, ax9)


    Plt.resize_to_layout!(fig)
    Plt.save("MorrisonandMilbrandtFig2_with_F_liq.svg", fig)
end

figure_2()
