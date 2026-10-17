import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)

# Testing terminal velocity with liquid fraction

# We need terminal velocity functions that take single
# particle dimension values rather than a whole PSD:

"""
    terminal_velocity_Dhelper_snowice(p3, Chen2022, D, ρ_r, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
 - D - particle dimension
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

 Returns fall speed for a single particle of dimension D, F_liq = 0
"""
function terminal_velocity_Dhelper_snowice(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeSnowIce,
    D::FT,
    ρ_r::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    # Get the free parameters for terminal velocities of small
    # and large particles
    small = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)
    large = CO.Chen2022_vel_coeffs_large(Chen2022, ρ_a)
    get_p(prs, it) = (prs[1][it], prs[2][it], prs[3][it])

    # Get the thresholds for different particles regimes
    (; D_cr, D_gr, ρ_g, ρ_d) = P3.thresholds(p3, ρ_r, F_r)
    D_th = P3.D_th_helper(p3)
    D_ct = Chen2022.cutoff

    # TODO: Change when each value used depending on type of particle
    # TODO: or keep fixed and add to ClimaParams...?
    κ = FT(-1 / 6) #FT(1/3)
    # Redefine α_va to be in si units
    α_va = P3.α_va_si(p3)

    spheres_n(a, b, c) = (a, b, c)

    aₙₛ(a) = a * (16 * p3.ρ_i^2 * p3.γ^3 / (9 * FT(π) * α_va^2))^κ
    bₙₛ(b) = b + κ * (3 * p3.σ - 2 * p3.β_va)

    non_spheres_n(a, b, c) = (aₙₛ(a), bₙₛ(b), c)

    aᵣₛ(a) = a * (p3.ρ_i / ρ_g)^(2 * κ)

    rimed_n(a, b, c) = (aᵣₛ(a), b, c)

    v_n_D_cr(D, a, b, c) =
        a *
        D^(b) *
        exp((-c) * D) *
        (
            16 * p3.ρ_i^2 * (F_r * π / 4 * D^2 + (1 - F_r) * p3.γ * D^p3.σ)^3 /
            (9 * π * (α_va / (1 - F_r) * D^p3.β_va)^2)
        )^κ

    v_snow_ice = 0

    v(a, b, c) = a * D^b * exp(-c * D)
    if D < D_ct
        size = small
    else
        size = large
    end
    for i in 1:2
        if F_r == 0
            if D < D_th
                v_snow_ice += v(spheres_n(get_p(size, i)...)...)
            else
                v_snow_ice += v(non_spheres_n(get_p(size, i)...)...)
            end
        else
            if D <= D_th
                v_snow_ice += v(spheres_n(get_p(size, i)...)...)
            elseif D <= D_gr
                v_snow_ice += v(non_spheres_n(get_p(size, i)...)...)
            elseif D <= D_cr
                v_snow_ice += v(rimed_n(get_p(size, i)...)...)
            else
                v_snow_ice += v_n_D_cr(D, get_p(size, i)...)
            end
        end
    end
    return v_snow_ice
end

"""
    terminal_velocity_Dhelper_liq(p3, Chen2022, D, ρ_r, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
 - D - particle dimension
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

 Returns fall speed for a single particle of dimension D, F_liq != 0
"""
function terminal_velocity_Dhelper_liq(
    p3::PSP3,
    Chen2022::CMP.Chen2022VelTypeRain,
    D::FT,
    ρ_r::FT,
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    # Get the free parameters for terminal velocities of small particles
    # (No large particles for rain?)
    small = CO.Chen2022_vel_coeffs_small(Chen2022, ρ_a)
    get_p(prs, it) = (prs[1][it], prs[2][it], prs[3][it])
    κ = 0
    v(D, a, b, c) = a * D^b * exp(-c * D)
    v_liq = 0
    for i in 1:3
        v_liq += v(D, get_p(small, i)...)
    end
    return v_liq
end

"""
    terminal_velocity_Dhelper(p3, Chen2022, D, ρ_r, F_liq, F_r, ρ_a)

 - p3 - a struct with P3 scheme parameters
 - Chen2022 - a struct with terminal velocity parameters as in Chen(2022)
            - pass (Chen2022VelTypeSnowIce, Chen2022VelTypeRain)
 - D - particle dimension
 - ρ_r - rime density (q_rim/B_rim) [kg/m^3]
 - F_liq - liquid fraction (q_liq/q_i,tot)
 - F_r - rime mass fraction (q_rim/q_i)
 - ρ_a - density of air

 Returns fall speed for a single particle of dimension D
"""
function terminal_velocity_Dhelper(
    p3::PSP3,
    Chen2022_ice::CMP.Chen2022VelTypeSnowIce,
    Chen2022_rain::CMP.Chen2022VelTypeRain,
    D::FT,
    ρ_r::FT,
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
) where {FT}
    return (1 - F_liq) * terminal_velocity_Dhelper_snowice(
        p3,
        Chen2022_ice,
        D,
        ρ_r,
        F_r,
        ρ_a,
    ) +
           F_liq * terminal_velocity_Dhelper_liq(
        p3,
        Chen2022_rain,
        D,
        ρ_r,
        F_liq,
        F_r,
        ρ_a,
    )
end

function get_values(
    p3::PSP3,
    Chen2022_ice::CMP.Chen2022VelTypeSnowIce,
    Chen2022_rain::CMP.Chen2022VelTypeRain,
    F_liq::FT,
    F_r::FT,
    ρ_a::FT,
    res::Int,
) where {FT <: Real}

    # particle dimension from 0 to 1 cm
    D_ps = range(FT(0), stop = FT(0.01), length = res)

    # ρ_r = 900
    ρ_r = FT(900)

    V_m = zeros(res)

    for i in 1:res
        V_m[i] = terminal_velocity_Dhelper(
            p3,
            Chen2022_ice,
            Chen2022_rain,
            D_ps[i],
            ρ_r,
            F_liq,
            F_r,
            ρ_a,
        )
    end
    return (D_ps, V_m)
end

function fig1()
    Chen2022 = CMP.Chen2022VelType(FT)
    ρ_a = FT(1.2)

    F_liqs = (FT(0), FT(0.33), FT(0.67), FT(1))
    F_rs = (FT(0), FT(1 - eps(FT)))

    labels = ["F_liq = 0", "F_liq = 0.33", "F_liq = 0.67", "F_liq = 1"]
    colors = [:black, :blue, :red, :green]
    res = 50

    fig = Plt.Figure(size = (1200, 400))

    #F_r = 0
    ax1 = Plt.Axis(
        fig[1:7, 1:9],
        title = "Fᵣ = 0",
        xlabel = "Dₚ (m)",
        ylabel = "V (m s⁻¹)",
    )

    for i in 1:4
        D_ps, V_0 = get_values(
            p3,
            Chen2022.snow_ice,
            Chen2022.rain,
            F_liqs[i],
            F_rs[1],
            ρ_a,
            res,
        )
        Plt.lines!(ax1, D_ps, V_0, color = colors[i], label = labels[i])
    end

    # F_r = 1
    ax2 = Plt.Axis(
        fig[1:7, 10:18],
        title = "Fᵣ = 1",
        xlabel = "Dₚ (m)",
        ylabel = "V (m s⁻¹)",
    )
    for i in 1:4
        D_ps, V_1 = get_values(
            p3,
            Chen2022.snow_ice,
            Chen2022.rain,
            F_liqs[i],
            F_rs[2],
            ρ_a,
            res,
        )
        Plt.lines!(ax2, D_ps, V_1, color = colors[i], label = labels[i])
    end

    Plt.Legend(
        fig[4, 19:22],
        [
            Plt.LineElement(color = :black),
            Plt.LineElement(color = :blue),
            Plt.LineElement(color = :red),
            Plt.LineElement(color = :green),
        ],
        labels,
        framevisible = false,
    )

    Plt.linkaxes!(ax1, ax2)
    ax1.xticks = (
        [0, 0.002, 0.004, 0.006, 0.008, 0.01],
        string.([0, 0.002, 0.004, 0.006, 0.008, 0.01]),
    )
    ax2.xticks = (
        [0, 0.002, 0.004, 0.006, 0.008, 0.01],
        string.([0, 0.002, 0.004, 0.006, 0.008, 0.01]),
    )


    Plt.resize_to_layout!(fig)
    Plt.save("Choletteetal2019_fig1.svg", fig)
end

fig1()
