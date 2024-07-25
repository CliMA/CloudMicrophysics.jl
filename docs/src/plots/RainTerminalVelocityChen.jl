import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.TerminalVelocity as TV
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)

function get_values(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    ρ_a::FT,
    F_r::FT,
    th,
    start::FT,
    stop::FT,
    res::Int,
) where {FT <: Real}

    # particle dimension from 0 to 0.5 cm
    D_ps = range(start, stop = stop, length = res)

    V_r = zeros(res)
    V_i = zeros(res)

    for i in 1:res
        V_r[i] = TV.velocity_chen(D_ps[i], Chen2022.rain, ρ_a)
        V_i[i] = TV.velocity_chen(D_ps[i], Chen2022.snow_ice, ρ_a, P3.p3_mass(p3, D_ps[i], F_r, FT(0), th),
            P3.p3_area(p3, D_ps[i], F_r, FT(0), th), p3.ρ_i)
    end
    return (D_ps, V_r, V_i)
end

function fig1()

    # small, med, large D_m
    start = FT(0)
    stops = (FT(500e-6), FT(5e-3), FT(5e-2))
    # try F_r = 0.5
    F_r = FT(0.5)
    ρ_r = FT(900)
    ρ_a = FT(1.2)
    Chen2022 = CMP.Chen2022VelType(FT)
    th = P3.thresholds(p3, ρ_r, F_r)

    # get values
    res = 100

    fig = Plt.Figure()

    labels = ["rain", "snow"]
    colors = [:deepskyblue2, :lightskyblue]

    fig = Plt.Figure(size = (1200, 400))

    #small
    ax1 = Plt.Axis(
        fig[1:7, 1:9],
        title = "Small particles",
        xlabel = "D (m)",
        ylabel = "V (m)",
    )

    Ds, V_r, V_i = get_values(
        p3,
        Chen2022,
        ρ_a,
        F_r,
        th,
        start,
        stops[1],
        res,
    )
    Plt.lines!(ax1, Ds, V_r, color = colors[1], label = labels[1], linewidth = 2)
    Plt.lines!(ax1, Ds, V_i, color = colors[2], label = labels[2], linewidth = 2)

    #med
    ax2 = Plt.Axis(
        fig[1:7, 10:18],
        title = "Medium particles",
        xlabel = "D (m)",
        ylabel = "V (m)",
    )

    Ds, V_r2, V_i2 = get_values(
        p3,
        Chen2022,
        ρ_a,
        F_r,
        th,
        start,
        stops[2],
        res,
    )
    Plt.lines!(ax2, Ds, V_r2, color = colors[1], label = labels[1], linewidth = 2)
    Plt.lines!(ax2, Ds, V_i2, color = colors[2], label = labels[2], linewidth = 2)

    #large
    ax3 = Plt.Axis(
        fig[1:7, 19:27],
        title = "Large particles",
        xlabel = "D (m)",
        ylabel = "V (m)",
    )

    Ds, V_r3, V_i3 = get_values(
        p3,
        Chen2022,
        ρ_a,
        F_r,
        th,
        start,
        stops[1],
        res,
    )
    Plt.lines!(ax3, Ds, V_r3, color = colors[1], label = labels[1], linewidth = 2)
    Plt.lines!(ax3, Ds, V_i3, color = colors[2], label = labels[2], linewidth = 2)

    Plt.Legend(
        fig[4, 28:31],
        [
            Plt.LineElement(color = :deepskyblue2),
            Plt.LineElement(color = :lightskyblue),
        ],
        labels,
        framevisible = false,
    )

    # Plt.linkaxes!(ax1, ax2, ax3)
    Plt.resize_to_layout!(fig)
    Plt.display(fig)
    # Plt.save("P3VvsD.svg", fig)
end

fig1()