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
    ρ_a::FT,
    x_resolution::Int,
    y_resolution::Int,
) where {FT <: Real}
    F_rims = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    F_liqs = range(FT(0), stop = FT(1), length = y_resolution)
    ρ_r = FT(900)
    use_aspect_ratio = true

    V_m = zeros(x_resolution, y_resolution)
    D_m_regimes = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)

    for i in 1:x_resolution
        for j in 1:y_resolution
            F_rim = F_rims[i]
            F_liq = F_liqs[j]
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
                use_aspect_ratio,
            )[2]
            D_m[i, j] = P3.D_m(p3, L, N, ρ_r, F_rim, F_liq)
            D_m_regimes[i, j] = D_m[i, j]

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
    return (; F_rims, F_liqs, V_m, D_m, D_m_regimes)
end

function make_axis(fig, row, col, title)
    return Plt.Axis(
        fig[row, col],
        xlabel = "F_rim",
        ylabel = "F_liq",
        title = title,
        width = 350,
        height = 350,
        limits = (0, 1.0, 0, 1.0),
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
    )
end

function fig1()
    Chen2022 = CMP.Chen2022VelType(FT)
    ρ_a = FT(1.2)

    # small D_m
    L_s = FT(0.0008)
    N_s = FT(1e6)
    # medium D_m
    L_m = FT(0.22)
    N_m = FT(1e6)
    # large D_m
    L_l = FT(0.7)
    N_l = FT(1e6)

    crange_small = range(0.39, 0.6, 20)
    crange_med = range(3, 6, 20)
    crange_large = range(4, 7, 20)

    # get V_m and D_m
    xres = 100
    yres = 100

    (F_rims, F_liqs, V_ms, D_ms, D_ms_regimes) =
        get_values(p3, Chen2022, L_s, N_s, ρ_a, xres, yres)
    (F_rimm, F_liqm, V_mm, D_mm, D_mm_regimes) =
        get_values(p3, Chen2022, L_m, N_m, ρ_a, xres, yres)
    (F_riml, F_liql, V_ml, D_ml, D_ml_regimes) =
        get_values(p3, Chen2022, L_l, N_l, ρ_a, xres, yres)

    args = (
        color = :black,
        labels = true,
        levels = 3,
        linewidth = 1.5,
        labelsize = 18,
    )


    fig = Plt.Figure()

    ax1_ = make_axis(fig, 1, 1, "Particle regimes with small Dₘ")
    hm = Plt.contourf!(
        ax1_,
        F_rims,
        F_liqs,
        D_ms_regimes,
        levels = 3,
        colormap = Plt.cgrad(:PuBuGn_3, 3, categorical = true),
    )
    Plt.Colorbar(
        fig[2, 1],
        hm,
        vertical = false,
        label = "graupel (2), partially rimed ice (3)",
        ticks = [1, 2, 3],
        labelsize = 14,
    )

    ax2_ = make_axis(fig, 1, 2, "Particle regimes with medium Dₘ")
    hm = Plt.contourf!(
        ax2_,
        F_rimm,
        F_liqm,
        D_mm_regimes,
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

    ax3_ = make_axis(fig, 1, 3, "Particle regimes with large Dₘ")
    hm = Plt.contourf!(
        ax3_,
        F_riml,
        F_liql,
        D_ml_regimes,
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

    ax1 = make_axis(fig, 3, 1, "Vₘ with small Dₘ")
    hm = Plt.contourf!(
        ax1,
        F_rims,
        F_liqs,
        V_ms,
        levels = crange_small,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(fig[4, 1], hm, vertical = false)
    Plt.contour!(ax1, F_rims, F_liqs, D_ms; args...)


    ax2 = make_axis(fig, 3, 2, "Vₘ with medium Dₘ")
    hm = Plt.contourf!(
        ax2,
        F_rimm,
        F_liqm,
        V_mm,
        levels = crange_med,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(fig[4, 2], hm, vertical = false)
    Plt.contour!(ax2, F_rimm, F_liqm, D_mm; args...)


    ax3 = make_axis(fig, 3, 3, "Vₘ with large Dₘ")
    hm = Plt.contourf!(
        ax3,
        F_riml,
        F_liql,
        V_ml,
        levels = crange_large,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(fig[4, 3], hm, vertical = false)
    Plt.contour!(ax3, F_riml, F_liql, D_ml; args...)


    # Plt.linkaxes!(ax1_, ax2_, ax3_)
    Plt.linkaxes!(ax1, ax2, ax3)

    Plt.resize_to_layout!(fig)
    Plt.save("P3TerminalVelocity_F_liq_rim.svg", fig)
end

fig1()
