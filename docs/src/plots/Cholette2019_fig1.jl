import CairoMakie as MK

import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Parameters.ParametersP3 as PSP3

FT = Float64

p3 = CMP.ParametersP3(FT)

# Testing terminal velocity with liquid fraction

function get_values(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelType,
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
    th,
    res::Int,
) where {FT <: Real}

    # particle dimension from 0 to 1 cm
    D_ps = range(FT(0), stop = FT(0.01), length = res)

    # ρ_r = 900
    ρ_r = FT(900)

    V = zeros(res)
    Vϕ = zeros(res)
    do_not_use_aspect = false
    use_aspect_ratio = true

    for i in 1:res
        Vϕ[i] = P3.p3_particle_terminal_velocity(
            p3,
            D_ps[i],
            Chen2022,
            ρ_a,
            F_r,
            F_liq,
            th,
            use_aspect_ratio,
        )
        V[i] = P3.p3_particle_terminal_velocity(
            p3,
            D_ps[i],
            Chen2022,
            ρ_a,
            F_r,
            F_liq,
            th,
            do_not_use_aspect,
        )
    end
    return (D_ps, V, Vϕ)
end

function fig1()
    Chen2022 = CMP.Chen2022VelType(FT)
    ρ_a = FT(1.2)
    ρ_r = FT(900)

    F_liqs = (FT(0), FT(0.33), FT(0.67), FT(1))
    F_rs = (FT(0), FT(1 - eps(FT)))

    labels =
        ["F_liq = 0", "F_liq = 0.33", "F_liq = 0.67", "F_liq = 1", "ϕᵢ = 1"]
    colors = [:black, :blue, :red, :green]
    res = 50

    fig = MK.Figure(size = (1200, 400))

    #F_r = 0
    ax1 = MK.Axis(
        fig[1:7, 1:9],
        title = "Fᵣ = 0",
        xlabel = "Dₚ (m)",
        ylabel = "V (m s⁻¹)",
    )

    for i in 1:4
        D_ps, V_0, V_0_ϕ = get_values(
            p3,
            Chen2022,
            F_liqs[i],
            F_rs[1],
            ρ_a,
            P3.thresholds(p3, ρ_r, F_rs[1]),
            res,
        )
        MK.lines!(ax1, D_ps, V_0_ϕ, color = colors[i], label = labels[i])
        MK.lines!(
            ax1,
            D_ps,
            V_0,
            color = colors[i],
            label = labels[i],
            linestyle = :dash,
        )
    end

    # F_r = 1
    ax2 = MK.Axis(
        fig[1:7, 10:18],
        title = "Fᵣ = 1",
        xlabel = "Dₚ (m)",
        ylabel = "V (m s⁻¹)",
    )
    for i in 1:4
        D_ps, V_1, V_1_ϕ = get_values(
            p3,
            Chen2022,
            F_liqs[i],
            F_rs[2],
            ρ_a,
            P3.thresholds(p3, ρ_r, F_rs[2]),
            res,
        )
        MK.lines!(ax2, D_ps, V_1_ϕ, color = colors[i], label = labels[i])
        MK.lines!(
            ax2,
            D_ps,
            V_1,
            color = colors[i],
            label = labels[i],
            linestyle = :dash,
        )
    end

    MK.Legend(
        fig[4, 19:22],
        [
            MK.LineElement(color = :black),
            MK.LineElement(color = :blue),
            MK.LineElement(color = :red),
            MK.LineElement(color = :green),
            MK.LineElement(color = :black, linestyle = :dash),
        ],
        labels,
        framevisible = false,
    )

    MK.linkaxes!(ax1, ax2)
    ax1.xticks = (
        [0, 0.002, 0.004, 0.006, 0.008, 0.01],
        string.([0, 0.002, 0.004, 0.006, 0.008, 0.01]),
    )
    ax2.xticks = (
        [0, 0.002, 0.004, 0.006, 0.008, 0.01],
        string.([0, 0.002, 0.004, 0.006, 0.008, 0.01]),
    )


    MK.resize_to_layout!(fig)
    MK.save("Choletteetal2019_fig1.svg", fig)
end

fig1()
