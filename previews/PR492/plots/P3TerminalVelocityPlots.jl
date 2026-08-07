import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt
# import ColorSchemes as CLR

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)
F_liq = FT(0) # set to 0 for now

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
    V_m_ϕ = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)
    D_m_regimes = zeros(x_resolution, y_resolution)
    ϕᵢ = zeros(x_resolution, y_resolution)

    for i in 1:x_resolution
        for j in 1:y_resolution
            F_rim = F_rims[i]
            ρ_r = ρ_rs[j]
            use_aspect_ratio = true
            do_not_use_aspect = false
            th = P3.thresholds(p3, ρ_r, F_rim)

            V_m[i, j] = P3.ice_terminal_velocity(
                p3,
                Chen2022,
                L,
                N,
                ρ_r,
                F_rim,
                F_liq,
                ρ_a,
                do_not_use_aspect,
            )[2]

            V_m_ϕ[i, j] = P3.ice_terminal_velocity(
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

            D_m[i, j] = P3.D_m(p3, L, N, ρ_r, F_rim, F_liq)
            D_m_regimes[i, j] = D_m[i, j]
            ϕᵢ[i, j] = P3.ϕᵢ(p3, D_m[i, j], F_rim, th)

            # plot the regimes
            D_th = P3.D_th_helper(p3)
            if D_th > D_m[i, j]
                # small spherical ice
                D_m_regimes[i, j] = 0
            elseif F_rim == 0
                # large nonspherical unrimed ice
                D_m_regimes[i, j] = 0
            elseif th.D_gr > D_m[i, j] >= D_th
                # dense nonspherical ice
                D_m_regimes[i, j] = 1
            elseif th.D_cr > D_m[i, j] >= th.D_gr
                # graupel
                D_m_regimes[i, j] = 2
            else #elseif D >= th.D_cr
                # partially rimed ice
                D_m_regimes[i, j] = 3
            end
        end
    end
    D_m *= 1e3
    return (; F_rims, ρ_rs, D_m_regimes, D_m, ϕᵢ, V_m, V_m_ϕ)
end

function make_axis(fig, row, col, title)
    return Plt.Axis(
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

    # D_th
    D_th = P3.D_th_helper(p3)

    # get V_m and D_m
    xres = 100
    yres = 100

    (F_rims, ρ_rs, D_m_regimes_s, D_m_s, ϕᵢ_s, V_m_s, V_m_ϕ_s) =
        get_values(p3, Chen2022, L_s, N_s, ρ_a, xres, yres)
    (F_rimm, ρ_rm, D_m_regimes_m, D_m_m, ϕᵢ_m, V_m_m, V_m_ϕ_m) =
        get_values(p3, Chen2022, L_m, N_m, ρ_a, xres, yres)
    (F_riml, ρ_rl, D_m_regimes_l, D_m_l, ϕᵢ_l, V_m_l, V_m_ϕ_l) =
        get_values(p3, Chen2022, L_l, N_l, ρ_a, xres, yres)

    fig = Plt.Figure()

    # Plot velocities as in Fig 2 in Morrison and Milbrandt 2015

    contour_args = (
        color = :black,
        labels = true,
        levels = 3,
        linewidth = 1.5,
        labelsize = 18,
    )

    ax1 = make_axis(fig, 1, 1, "Particle regimes with small Dₘ")
    hm = Plt.contourf!(
        ax1,
        F_rims,
        ρ_rs,
        D_m_regimes_s,
        levels = 3,
        colormap = Plt.cgrad(:PuBuGn_3, 3, categorical = true),
    )
    Plt.Colorbar(
        fig[2, 1],
        hm,
        vertical = false,
        label = "dense nonspherical ice (1), graupel (2), partially rimed ice (3)",
        ticks = [1, 2, 3],
        labelsize = 14,
    )

    ax2 = make_axis(fig, 1, 2, "Particle regimes with medium Dₘ")
    hm = Plt.contourf!(
        ax2,
        F_rimm,
        ρ_rm,
        D_m_regimes_m,
        levels = 3,
        colormap = Plt.cgrad(:PuBuGn_3, 3, categorical = true),
    )
    Plt.Colorbar(
        fig[2, 2],
        hm,
        vertical = false,
        label = "graupel (2), partially rimed ice (3)",
        ticks = [1, 2, 3],
        labelsize = 14,
    )

    ax3 = make_axis(fig, 1, 3, "Particle regimes with large Dₘ")
    hm = Plt.contourf!(
        ax3,
        F_riml,
        ρ_rl,
        D_m_regimes_l,
        levels = 3,
        colormap = Plt.cgrad(:PuBuGn_3, 3, categorical = true),
    )
    Plt.Colorbar(
        fig[2, 3],
        hm,
        vertical = false,
        label = "graupel (2), partially rimed ice (3)",
        ticks = [1, 2, 3],
        labelsize = 14,
    )

    ax4 = make_axis(fig, 3, 1, "Dₘ (mm) (L = 8e-4 kgm⁻³, N = 1e6 m⁻³)")
    hm = Plt.contourf!(ax4, F_rims, ρ_rs, D_m_s)
    Plt.Colorbar(fig[4, 1], hm, vertical = false)

    ax5 = make_axis(fig, 3, 2, "Dₘ (mm) (L = 0.22 kgm⁻³, N = 1e6 m⁻³)")
    hm = Plt.contourf!(ax5, F_rimm, ρ_rm, D_m_m)
    Plt.Colorbar(fig[4, 2], hm, vertical = false)

    ax6 = make_axis(fig, 3, 3, "Dₘ (mm) (L = 0.7 kgm⁻³, N = 1e6 m⁻³)")
    hm = Plt.contourf!(ax6, F_riml, ρ_rl, D_m_l)
    Plt.Colorbar(fig[4, 3], hm, vertical = false)

    ax7 = make_axis(fig, 5, 1, "ϕᵢ with small Dₘ")
    hm = Plt.contourf!(ax7, F_rims, ρ_rs, ϕᵢ_s, levels = 20)
    Plt.Colorbar(fig[6, 1], hm, vertical = false)
    Plt.contour!(ax7, F_rims, ρ_rs, D_m_s, label = "D_m"; contour_args...)

    ax8 = make_axis(fig, 5, 2, "ϕᵢ with medium Dₘ")
    hm = Plt.contourf!(ax8, F_rimm, ρ_rm, ϕᵢ_m)
    Plt.Colorbar(fig[6, 2], hm, vertical = false)
    Plt.contour!(ax8, F_rimm, ρ_rm, D_m_m, label = "D_m"; contour_args...)

    ax9 = make_axis(fig, 5, 3, "ϕᵢ with large Dₘ")
    hm = Plt.contourf!(ax9, F_riml, ρ_rl, ϕᵢ_l)
    Plt.Colorbar(fig[6, 3], hm, vertical = false)
    Plt.contour!(ax9, F_riml, ρ_rl, D_m_l, label = "D_m"; contour_args...)

    ax10 = make_axis(fig, 7, 1, "Vₘ (ϕᵢ = 1) with small Dₘ")
    hm = Plt.contourf!(ax10, F_rims, ρ_rs, V_m_s)
    Plt.contour!(ax10, F_rims, ρ_rs, D_m_s; contour_args...)
    Plt.Colorbar(fig[8, 1], hm, vertical = false)

    ax11 = make_axis(fig, 7, 2, "Vₘ (ϕᵢ = 1) with medium Dₘ")
    hm = Plt.contourf!(ax11, F_rimm, ρ_rm, V_m_m)
    Plt.contour!(ax11, F_rimm, ρ_rm, D_m_m; contour_args...)
    Plt.Colorbar(fig[8, 2], hm, vertical = false)

    ax12 = make_axis(fig, 7, 3, "Vₘ  (ϕᵢ = 1) with large Dₘ")
    hm = Plt.contourf!(ax12, F_riml, ρ_rl, V_m_l)
    Plt.contour!(ax12, F_riml, ρ_rl, D_m_l; contour_args...)
    Plt.Colorbar(fig[8, 3], hm, vertical = false)

    ax13 = make_axis(fig, 9, 1, "Vₘ (using ϕᵢ) with small Dₘ")
    hm = Plt.contourf!(ax13, F_rims, ρ_rs, V_m_ϕ_s)
    Plt.contour!(ax13, F_rims, ρ_rs, D_m_s; contour_args...)
    Plt.Colorbar(fig[10, 1], hm, vertical = false)

    ax14 = make_axis(fig, 9, 2, "Vₘ (using ϕᵢ) with medium Dₘ")
    hm = Plt.contourf!(ax14, F_rimm, ρ_rm, V_m_ϕ_m)
    Plt.contour!(ax14, F_rimm, ρ_rm, D_m_m; contour_args...)
    Plt.Colorbar(fig[10, 2], hm, vertical = false)

    ax15 = make_axis(fig, 9, 3, "Vₘ  (using ϕᵢ) with large Dₘ")
    hm = Plt.contourf!(ax15, F_riml, ρ_rl, V_m_ϕ_l)
    Plt.contour!(ax15, F_riml, ρ_rl, D_m_l; contour_args...)
    Plt.Colorbar(fig[10, 3], hm, vertical = false)

    Plt.linkxaxes!(ax1, ax2)
    Plt.linkxaxes!(ax2, ax3)
    Plt.linkyaxes!(ax1, ax2)
    Plt.linkyaxes!(ax2, ax3)
    Plt.linkxaxes!(ax1, ax4)
    Plt.linkxaxes!(ax4, ax5)
    Plt.linkxaxes!(ax5, ax6)
    Plt.linkyaxes!(ax1, ax4)
    Plt.linkyaxes!(ax4, ax5)
    Plt.linkyaxes!(ax5, ax6)
    Plt.linkxaxes!(ax1, ax7)
    Plt.linkxaxes!(ax7, ax8)
    Plt.linkxaxes!(ax5, ax9)
    Plt.linkyaxes!(ax1, ax7)
    Plt.linkyaxes!(ax7, ax8)
    Plt.linkyaxes!(ax8, ax9)
    Plt.linkxaxes!(ax1, ax10)
    Plt.linkxaxes!(ax10, ax11)
    Plt.linkxaxes!(ax11, ax12)
    Plt.linkyaxes!(ax1, ax10)
    Plt.linkyaxes!(ax10, ax11)
    Plt.linkyaxes!(ax11, ax12)
    Plt.linkxaxes!(ax1, ax13)
    Plt.linkxaxes!(ax13, ax14)
    Plt.linkxaxes!(ax14, ax15)
    Plt.linkyaxes!(ax1, ax13)
    Plt.linkyaxes!(ax13, ax14)
    Plt.linkyaxes!(ax14, ax15)

    Plt.resize_to_layout!(fig)
    Plt.save("MorrisonandMilbrandtFig2.svg", fig)
end
#! format: on

# Terminal Velocity figure
figure_2()
