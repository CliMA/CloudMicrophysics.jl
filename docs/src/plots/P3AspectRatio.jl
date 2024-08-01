import ClimaParams
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64

p3 = CMP.ParametersP3(FT)

function get_values(
    p3::PSP3,
    F_liq::FT,
    F_rim::FT,
    th,
    x_resolution::Int,
    y_resolution::Int,
) where {FT <: Real}
    F_rims = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    F_liqs = range(FT(0), stop = FT(1), length = y_resolution)
    ρ_r = FT(900)

    V_m = zeros(x_resolution, y_resolution)
    ϕᵢ = zeros(res)
    F_rim = zeros(xres)
    F_liq = 
    for i in 1:xres
        for j in 1:yres
        ϕᵢ[i, j] = P3.ϕᵢ(
            P3.p3_mass(p3,
            D_ps[i],
            F_rim,
            F_liq,
            th),
            P3.p3_area(p3,
            D_ps[i],
            F_rim,
            F_liq,
            th),
            p3.ρ_i)
    end
    return (D_ps, ϕᵢ)
end

function fig1()
    ρ_a = FT(1.2)
    ρ_r = FT(900)

    res = 200

    fig = Plt.Figure(size = (1200, 400))

    ax1 = make_axis(fig, 1, "ϕᵢ")
    hm = Plt.contourf!(
        ax1,
        F_rims,
        F_liqs,
        ϕᵢ,
        levels = crange_small,
        extendlow = :auto,
        extendhigh = :auto,
    )
    Plt.Colorbar(
        fig[2, 1],
        hm,
        vertical = false,
    )

    Plt.resize_to_layout!(fig)
    Plt.display(fig)
    #Plt.save("P3AspectRatio.svg", fig)
end

fig1()