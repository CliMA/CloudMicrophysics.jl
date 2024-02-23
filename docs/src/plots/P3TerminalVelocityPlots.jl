import CLIMAParameters
import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CairoMakie as Plt

const PSP3 = CMP.ParametersP3

FT = Float64 

function get_values(p3::PSP3, Chen2022::CMP.Chen2022VelTypeSnowIce, q::FT, N::FT, ρ_a::FT, x_resolution::Int, y_resolution::Int) where {FT} 
    F_rs = range(FT(0), stop = FT(1 - eps(FT)), length = x_resolution)
    ρ_rs = range(FT(25), stop = FT(975), length = y_resolution)

    V_m = zeros(x_resolution, y_resolution)
    D_m = zeros(x_resolution, y_resolution)

    for i in 1:x_resolution 
        for j in 1:y_resolution 
            F_r = F_rs[i]
            ρ_r = ρ_rs[j]

            V_m[i, j] = P3.terminal_velocity_mass(p3, Chen2022, q, N, ρ_r, F_r, ρ_a) 
            D_m[i, j] = P3.D_m(p3, q, N, ρ_r, F_r)
        end 
    end

    minV = minimum(v for v in V_m if !isnan(v))
    maxV = maximum(v for v in V_m if !isnan(v))

    return (; F_rs, ρ_rs, V_m, D_m, minV, maxV)
end

function figure_2()
    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    # density of air in kg/m^3
    ρ_a = FT(1.293)

    f = Plt.Figure()

    # small D_m
    q_s = FT(0.001)
    N_s = FT(1e6)     
    
    Plt.Axis(
        f[1, 1], 
        xlabel = "F_r", 
        ylabel = "ρ_r", 
        title = "Small Dₘ", 
        width = 300, 
        height = 300,
        limits = (0, 1.0, 25, 975), 
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    )

    (F_rs, ρ_rs, V_m, D_m, minV, maxV) = get_values(
        p3, 
        Chen2022.snow_ice, 
        q_s, 
        N_s, 
        ρ_a, 
        50, 
        50)

    println(D_m[1:10])

    Plt.heatmap!(
        F_rs, 
        ρ_rs, 
        V_m, 
        colormap = Plt.reverse(Plt.cgrad(:Spectral))
        )
    Plt.contour!(
        F_rs, 
        ρ_rs, 
        D_m, 
        color = :black, 
        labels = true, 
        levels = 3
        )
    Plt.Colorbar(
        f[1, 2], 
        limits = (minV, maxV), 
        colormap = Plt.reverse(Plt.cgrad(:Spectral)),
        flipaxis = true,
    )
    

    Plt.resize_to_layout!(f) 
    Plt.save("MorrisonandMilbrandtFig2.svg", f)
end

println("start") 
figure_2()
println("done")