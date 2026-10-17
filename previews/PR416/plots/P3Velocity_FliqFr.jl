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
    ρ_a::FT,
    x_resolution::Int,
    y_resolution::Int,
) where {FT <: Real}
    F_rs = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    F_liqs = range(FT(0), stop = FT(1), length = y_resolution)
    ρ_r = FT(900)

    V_m = zeros(x_resolution, y_resolution)

    for i in 1:x_resolution
        for j in 1:y_resolution
            F_r = F_rs[i]
            F_liq = F_liqs[j]

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
        end
    end
    return (; F_rs, F_liqs, V_m)
end

function make_axis(fig, col, title)
    return Plt.Axis(
        fig[1, col],
        xlabel = "F_r",
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

    (F_rs, F_liqs, V_ms) = get_values(
        p3,
        Chen2022.snow_ice,
        Chen2022.rain,
        q_s,
        N_s,
        ρ_a,
        xres,
        yres,
    )
    (F_rm, F_liqm, V_mm) = get_values(
        p3,
        Chen2022.snow_ice,
        Chen2022.rain,
        q_m,
        N_m,
        ρ_a,
        xres,
        yres,
    )
    (F_rl, F_liql, V_ml) = get_values(
        p3,
        Chen2022.snow_ice,
        Chen2022.rain,
        q_l,
        N_l,
        ρ_a,
        xres,
        yres,
    )

    fig = Plt.Figure()

    #F_r = 0
    ax1 = make_axis(fig, 1, "Vₘ with small Dₘ")
    hm = Plt.contourf!(
        ax1,
        F_rs,
        F_liqs,
        V_ms,
        levels = 0.39:0.0175:0.61,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(
        fig[2, 1],
        hm,
        vertical = false,
        ticks = [0.4, 0.45, 0.5, 0.55, 0.6],
    )

    ax2 = make_axis(fig, 2, "Vₘ with medium Dₘ")
    hm = Plt.contourf!(
        ax2,
        F_rm,
        F_liqm,
        V_mm,
        levels = 3:0.5:8,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(fig[2, 2], hm, vertical = false, ticks = [3, 4, 5, 6, 7, 8])

    ax3 = make_axis(fig, 3, "Vₘ with large Dₘ")
    hm = Plt.contourf!(
        ax3,
        F_rl,
        F_liql,
        V_ml,
        levels = 2:0.8:10,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(fig[2, 3], hm, vertical = false, ticks = [2, 4, 6, 8, 10])


    Plt.linkaxes!(ax1, ax2, ax3)

    Plt.resize_to_layout!(fig)
    Plt.save("P3Velocity_FliqFr.svg", fig)
end

fig1()
