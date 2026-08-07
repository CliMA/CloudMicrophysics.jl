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
    Chen2022_ice::CMP.Chen2022VelTypeSnowIce,
    Chen2022_rain::CMP.Chen2022VelTypeRain,
    q::FT,
    N::FT,
    F_liq::FT,
    ρ_a::FT,
    x_resolution::Int,
    y_resolution::Int,
) where {FT <: Real}
    F_rs = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    ρ_rs = range(FT(25), stop = FT(975), length = y_resolution)

    V_m = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)

    for i in 1:x_resolution
        for j in 1:y_resolution
            F_r = F_rs[i]
            ρ_r = ρ_rs[j]

            V_m[i, j] = P3.terminal_velocity_tot(
                p3,
                Chen2022_ice,
                Chen2022_rain,
                q,
                N,
                ρ_r,
                F_liq,
                F_r,
                ρ_a,
            )[2]
            # # get D_m in mm for plots
            D_m[i, j] = 1e3 * P3.D_m(p3, q, N, ρ_r, F_r, F_liq)
        end
    end
    return (; F_rs, ρ_rs, V_m, D_m)
end

function make_axis_top(fig, col, title)
    return Plt.Axis(
        fig[1, col],
        xlabel = "F_r",
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
        xlabel = "F_r",
        ylabel = "ρ_r",
        title = title,
    )
end

#! format: off
function figure_2(F_liq)
    Chen2022 = CMP.Chen2022VelType(FT)
    # density of air in kg/m^3
    ρ_a = FT(1.2) #FT(1.293)
    # small D_m
    q_s = FT(0.0008)
    N_s = FT(1e6)
    # medium D_m
    q_m = FT(0.22)
    N_m = FT(1e6)
    # large D_m
    q_l = FT(0.7)
    N_l = FT(1e6)
    # get V_m and D_m
    xres = 100
    yres = 100

    (F_rs, ρ_rs, V_ms, D_ms) = get_values(p3, Chen2022.snow_ice, Chen2022.rain, q_s, N_s, F_liq, ρ_a, xres, yres)
    (F_rm, ρ_rm, V_mm, D_mm) = get_values(p3, Chen2022.snow_ice, Chen2022.rain, q_m, N_m, F_liq, ρ_a, xres, yres)
    (F_rl, ρ_rl, V_ml, D_ml) = get_values(p3, Chen2022.snow_ice, Chen2022.rain, q_l, N_l, F_liq, ρ_a, xres, yres)

    fig = Plt.Figure()

    # Plot velocities as in Fig 2 in Morrison and Milbrandt 2015

    args = (color = :black, labels = true, levels = 3, linewidth = 1.5, labelsize = 18)

    # set the colorscale to be the same for all figures
    crange_small = range(0, 0.6, length = 20)
    crange_med_large = range(0, 9, length = 20)

    ax1 = make_axis_top(fig, 1, "Small Dₘ, F_liq = $F_liq")
    hm = Plt.contourf!(ax1, F_rs, ρ_rs, V_ms, levels = crange_small, extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax1, F_rs, ρ_rs, D_ms; args...)
    Plt.Colorbar(fig[2, 1], hm, vertical = false)

    ax2 = make_axis_top(fig, 2, "Medium Dₘ, F_liq = $F_liq")
    hm = Plt.contourf!(ax2, F_rm, ρ_rm, V_mm, levels = crange_med_large, extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax2, F_rm, ρ_rm, D_mm; args...)
    Plt.Colorbar(fig[2, 2], hm, vertical = false)

    ax3 = make_axis_top(fig, 3, "Large Dₘ, F_liq = $F_liq")
    hm = Plt.contourf!(ax3, F_rl, ρ_rl, V_ml, levels = crange_med_large, extendlow = :auto, extendhigh = :auto)
    Plt.contour!(ax3, F_rl, ρ_rl, D_ml; args...)
    Plt.Colorbar(fig[2, 3], hm, vertical = false)

    Plt.linkxaxes!(ax1, ax2)
    Plt.linkxaxes!(ax2, ax3)
    Plt.linkyaxes!(ax1, ax2)
    Plt.linkyaxes!(ax2, ax3)

    ## Plot D_m as second row of comparisons
    crange_small = range(0, 0.2, length = 20)
    crange_med_large = range(0, 9, length = 20)

    ax4 = make_axis_bottom(fig, 1, "Small Dₘ vs F_r and ρ_r, F_liq = $F_liq")
    hm = Plt.contourf!(ax4, F_rs, ρ_rs, D_ms, levels = crange_small, extendlow = :auto, extendhigh = :auto)
    Plt.Colorbar(fig[4, 1], hm, vertical = false)

    ax5 = make_axis_bottom(fig, 2, "Medium Dₘ vs F_r and ρ_r, F_liq = $F_liq")
    hm = Plt.contourf!(ax5, F_rm, ρ_rm, D_mm, levels = crange_med_large, extendlow = :auto, extendhigh = :auto)
    Plt.Colorbar(fig[4, 2], hm, vertical = false)

    ax6 = make_axis_bottom(fig, 3, "Large Dₘ vs F_r and ρ_r, F_liq = $F_liq")
    hm = Plt.contourf!(ax6, F_rl, ρ_rl, D_ml, levels = crange_med_large, extendlow = :auto, extendhigh = :auto)
    Plt.Colorbar(fig[4, 3], hm, vertical = false)

    Plt.linkxaxes!(ax1, ax4)
    Plt.linkxaxes!(ax4, ax5)
    Plt.linkxaxes!(ax5, ax6)
    Plt.linkyaxes!(ax1, ax4)
    Plt.linkyaxes!(ax4, ax5)
    Plt.linkyaxes!(ax5, ax6)


    Plt.resize_to_layout!(fig)
    Plt.save("MorrisonandMilbrandtFig2_$(F_liq).svg", fig)
end
#! format: on

# Terminal Velocity figures
F_liq = [FT(0), FT(0.33), FT(0.67), FT(0.99)]
for i in F_liq
    figure_2(i)
end
