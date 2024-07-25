import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)

function get_values(
    p3::PSP3,
    q::FT,
    N::FT,
    F_r::FT,
    res::Int,
) where {FT <: Real}

    # particle dimension from 0 to 0.5 cm
    F_liqs = range(FT(0), stop = FT(1), length = res)

    # ρ_r = 900
    ρ_r = FT(900)

    D = zeros(res)

    for i in 1:res
        D[i] = P3.D_m(p3, q, N, ρ_r, F_r, F_liqs[i])
    end
    return (F_liqs, D)
end

function fig1()

    # small, med, large D_m
    qs = (FT(0.0008), FT(0.22), FT(0.7))
    Ns = (FT(1e6), FT(1e6), FT(1e6))
    F_rs = (FT(0), FT(0.5), FT(0.8))

    # get values
    res = 100

    fig = Plt.Figure()

    labels = ["Fᵣᵢₘ = 0", "Fᵣᵢₘ = 0.5", "Fᵣᵢₘ = 0.8"]
    colors = [:lightskyblue, :deepskyblue2, :deepskyblue4]

    fig = Plt.Figure(size = (1200, 400))

    #small
    ax1 = Plt.Axis(
        fig[1:7, 1:9],
        title = "Small particles",
        xlabel = "F_liq",
        ylabel = "Dₘ (m)",
    )

    for i in 1:3
        F_liqs, D_1 = get_values(
            p3,
            qs[1],
            Ns[1],
            F_rs[i],
            res,
        )
        Plt.lines!(ax1, F_liqs, D_1, color = colors[i], label = labels[i], linewidth = 2)
    end

    #med
    ax2 = Plt.Axis(
        fig[1:7, 10:18],
        title = "Medium particles",
        xlabel = "F_liq",
        ylabel = "Dₘ (m)",
    )

    for i in 1:3
        F_liqs, D_2 = get_values(
            p3,
            qs[2],
            Ns[2],
            F_rs[i],
            res,
        )
        Plt.lines!(ax2, F_liqs, D_2, color = colors[i], label = labels[i], linewidth = 2)
    end

    #large
    ax3 = Plt.Axis(
        fig[1:7, 19:27],
        title = "Large particles",
        xlabel = "F_liq",
        ylabel = "Dₘ (m)",
    )

    for i in 1:3
        F_liqs, D_3 = get_values(
            p3,
            qs[3],
            Ns[3],
            F_rs[i],
            res,
        )
        Plt.lines!(ax3, F_liqs, D_3, color = colors[i], label = labels[i], linewidth = 2)
    end

    Plt.Legend(
        fig[4, 28:31],
        [
            Plt.LineElement(color = :lightskyblue),
            Plt.LineElement(color = :deepskyblue2),
            Plt.LineElement(color = :deepskyblue4),
        ],
        labels,
        framevisible = false,
    )

    # Plt.linkaxes!(ax1, ax2, ax3)
    Plt.resize_to_layout!(fig)
    Plt.display(fig)
    # Plt.save("P3DvsFliq.svg", fig)
end

fig1()
