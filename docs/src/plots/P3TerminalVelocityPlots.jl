import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)
# keep this fig the same for the velocity docs section
# (i.e. F_liq = 0)
F_liq = FT(0)

function get_values(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    L::FT,
    N::FT,
    ρ_a::FT,
    x_resolution::Int,
    y_resolution::Int,
) where {FT}
    F_rims = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    ρ_rs = range(FT(25), stop = FT(975), length = y_resolution)

    V_m = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)
    aspect_ratio = false
    F_liq = FT(0)

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
function make_axis_bottom(fig, col, title)
    return Plt.Axis(
        fig[3, col],
        height = 350,
        width = 350,
        xlabel = "F_rim",
        ylabel = "ρ_r",
        title = title,
    )
end

#! format: off
function figure_2()
    Chen2022 = CMP.Chen2022VelType(FT)
    # density of air in kg/m^3
    ρ_a = FT(1.2) #FT(1.293)
    # small D_m
    L_s = FT(0.0008)
    N_s = FT(1e6)
    # medium D_m
    L_m = FT(0.22)
    N_m = FT(1e6)
    # large D_m
    L_l = FT(0.7)
    N_l = FT(1e6)
    # get V_m and D_m
    xres = 100
    yres = 100
    (F_rims, ρ_rs, V_ms, D_ms) = get_values(p3, Chen2022, L_s, N_s, ρ_a, xres, yres)
    (F_rimm, ρ_rm, V_mm, D_mm) = get_values(p3, Chen2022, L_m, N_m, ρ_a, xres, yres)
    (F_riml, ρ_rl, V_ml, D_ml) = get_values(p3, Chen2022, L_l, N_l, ρ_a, xres, yres)

    fig = Plt.Figure()

    # Plot velocities as in Fig 2 in Morrison and Milbrandt 2015

    args = (color = :black, labels = true, levels = 3, linewidth = 1.5, labelsize = 18)

    ax1 = make_axis_top(fig, 1, "Small Dₘ")
    hm = Plt.contourf!(ax1, F_rims, ρ_rs, V_ms)
    Plt.contour!(ax1, F_rims, ρ_rs, D_ms; args...)
    Plt.Colorbar(fig[2, 1], hm, vertical = false)

    ax2 = make_axis_top(fig, 2, "Medium Dₘ")
    hm = Plt.contourf!(ax2, F_rimm, ρ_rm, V_mm)
    Plt.contour!(ax2, F_rimm, ρ_rm, D_mm; args...)
    Plt.Colorbar(fig[2, 2], hm, vertical = false)

    ax3 = make_axis_top(fig, 3, "Large Dₘ")
    hm = Plt.contourf!(ax3, F_riml, ρ_rl, V_ml)
    Plt.contour!(ax3, F_riml, ρ_rl, D_ml; args...)
    Plt.Colorbar(fig[2, 3], hm, vertical = false)

    Plt.linkxaxes!(ax1, ax2)
    Plt.linkxaxes!(ax2, ax3)
    Plt.linkyaxes!(ax1, ax2)
    Plt.linkyaxes!(ax2, ax3)

    ## Plot D_m as second row of comparisons

    ax4 = make_axis_bottom(fig, 1, "Small Dₘ vs F_rim and ρ_r")
    hm = Plt.contourf!(ax4, F_rims, ρ_rs, D_ms)
    Plt.Colorbar(fig[4, 1], hm, vertical = false)

    ax5 = make_axis_bottom(fig, 2, "Medium Dₘ vs F_rim and ρ_r")
    hm = Plt.contourf!(ax5, F_rimm, ρ_rm, D_mm)
    Plt.Colorbar(fig[4, 2], hm, vertical = false)

    ax6 = make_axis_bottom(fig, 3, "Large Dₘ vs F_rim and ρ_r")
    hm = Plt.contourf!(ax6, F_riml, ρ_rl, D_ml)
    Plt.Colorbar(fig[4, 3], hm, vertical = false)

    Plt.linkxaxes!(ax1, ax4)
    Plt.linkxaxes!(ax4, ax5)
    Plt.linkxaxes!(ax5, ax6)
    Plt.linkyaxes!(ax1, ax4)
    Plt.linkyaxes!(ax4, ax5)
    Plt.linkyaxes!(ax5, ax6)

    Plt.resize_to_layout!(fig)
    Plt.save("MorrisonandMilbrandtFig2.svg", fig)
end
#! format: on

# Terminal Velocity figure
figure_2()
