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
            # get D_m in mm for plots 
            D_m[i, j] = 1e3 * P3.D_m(p3, q, N, ρ_r, F_r)
        end 
    end
    return (; F_rs, ρ_rs, V_m, D_m)
end

function figure_2()
    p3 = CMP.ParametersP3(FT)
    Chen2022 = CMP.Chen2022VelType(FT)
    # density of air in kg/m^3
    ρ_a = FT(1.293)

    f = Plt.Figure()
    xres = 10
    yres = 10
    min = FT(0)
    max = FT(10)

    # small D_m
    q_s = FT(0.00074047)
    N_s = FT(1e6)    
    
    # medium D_m 
    q_m = FT(0.0028151) 
    N_m = FT(1e6)

    # large D_m
    q_l = FT(0.02494021) 
    N_l = FT(1e6)
    
    # get V_m and D_m
    (F_rs, ρ_rs, V_ms, D_ms) = get_values(p3, Chen2022.snow_ice, q_s, N_s, ρ_a, xres, yres)
    (F_rm, ρ_rm, V_mm, D_mm) = get_values(p3, Chen2022.snow_ice, q_m, N_m, ρ_a, xres, yres)
    (F_rl, ρ_rl, V_ml, D_ml) = get_values(p3, Chen2022.snow_ice, q_l, N_l, ρ_a, xres, yres)
    
    println("small: min: ", minimum(d for d in D_ms if !isnan(d))," max: ", maximum(d for d in D_ms if !isnan(d)))
    println("medium: min: ", minimum(d for d in D_mm if !isnan(d))," max: ", maximum(d for d in D_mm if !isnan(d)))
    println("large: min: ", minimum(d for d in D_ml if !isnan(d))," max: ", maximum(d for d in D_ml if !isnan(d)))

    #println(D_mm)

    println("") 
    #println(D_ml)

    Plt.Axis(
        f[1, 1], 
        xlabel = "F_r", 
        ylabel = "ρ_r", 
        title = "Small Dₘ", 
        width = 350, 
        height = 350,
        limits = (0, 1.0, 25, 975), 
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    )

    Plt.heatmap!(F_rs, ρ_rs, V_ms, colormap = Plt.reverse(Plt.cgrad(:Spectral)), colorrange = (min, max),)
    Plt.contour!(F_rs, ρ_rs, D_ms, color = :black, labels = true, levels = 3,)

    Plt.Axis(
        f[1, 2], 
        xlabel = "F_r", 
        ylabel = "ρ_r", 
        title = "Medium Dₘ", 
        width = 350, 
        height = 350,
        limits = (0, 1.0, 25, 975), 
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    )

    Plt.heatmap!(F_rm, ρ_rm, V_mm, colormap = Plt.reverse(Plt.cgrad(:Spectral)), colorrange = (min, max),)
    Plt.contour!(F_rm, ρ_rm, D_mm, color = :black, labels = true, levels = 3)

    Plt.Axis(
        f[1, 3], 
        xlabel = "F_r", 
        ylabel = "ρ_r", 
        title = "Large Dₘ", 
        width = 350, 
        height = 350,
        limits = (0, 1.0, 25, 975), 
        xticks = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0],
        yticks = [200, 400, 600, 800],
    )

    Plt.heatmap!(F_rl, ρ_rl, V_ml, colormap = Plt.reverse(Plt.cgrad(:Spectral)), colorrange = (min, max),)
    Plt.contour!(F_rl, ρ_rl, D_ml, color = :black, labels = true, levels = 3)

    Plt.Colorbar(
        f[2, 2], 
        limits = (min, max), 
        colormap = Plt.reverse(Plt.cgrad(:Spectral)),
        flipaxis = true,
        vertical = false, 
    )

    # Attempt to plot D_m for clearer information than the contour plots and debugging 

    Plt.Axis(f[3, 1], height = 350, width = 350)
    Plt.heatmap!(
        F_rs, 
        ρ_rs, 
        D_ms,
        colormap = Plt.reverse(Plt.cgrad(:Spectral))
    )
    Plt.Colorbar(
        f[4, 1], 
        limits = (minimum(d for d in D_ms if !isnan(d)), maximum(d for d in D_ms if !isnan(d))), 
        colormap = Plt.reverse(Plt.cgrad(:Spectral)),
        flipaxis = true,
        vertical = false, 
    )
    Plt.Axis(f[3, 2], height = 350, width = 350)
    Plt.heatmap!(
        F_rm, 
        ρ_rm, 
        D_mm,
        colormap = Plt.reverse(Plt.cgrad(:Spectral))
    )
    Plt.Colorbar(
        f[4, 2], 
        limits = (minimum(d for d in D_mm if !isnan(d)), maximum(d for d in D_mm if !isnan(d))), 
        colormap = Plt.reverse(Plt.cgrad(:Spectral)),
        flipaxis = true,
        vertical = false, 
    )
    Plt.Axis(f[3, 3], height = 350, width = 350)
    Plt.heatmap!(
        F_rl, 
        ρ_rl, 
        D_ml,
        colormap = Plt.reverse(Plt.cgrad(:Spectral))
    )
    Plt.Colorbar(
        f[4, 3], 
        limits = (minimum(d for d in D_ml if !isnan(d)), maximum(d for d in D_ml if !isnan(d))), 
        colormap = Plt.reverse(Plt.cgrad(:Spectral)),
        flipaxis = true,
        vertical = false, 
    )
    
    Plt.resize_to_layout!(f) 
    Plt.save("MorrisonandMilbrandtFig2.svg", f)
end

println("start") 
#figure_2()
#println("done")
#println(P3.D_th_helper(p3))
#println(P3.thresholds(p3, FT(500), FT(0.5)))
#println("")
#println("small = ", P3.q_gamma(p3, FT(0.5), FT(1e7), FT(log(4.9 * 10^2)), P3.thresholds(p3, FT(500), FT(0.5))))


import RootSolvers as RS
ρ_r = FT(500) 
F_r = FT(0.8)
N = FT(1e6) 
Dₘ = 0.006
println("started")
x = q <= 0 ? FT(0) : log(q)
shape_problem(x) = Dₘ - P3.D_m(p3, x, N, ρ_r, F_r)
x =
        RS.find_zero(
            shape_problem,
            RS.SecantMethod(log(0.01), log(0.015)),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            5,
        ).root

println("q_solved = ", exp(x)) 

function plotDvsq() 
    n = 100
    log10qs = range(FT(-4), stop = FT(0), length = n)
    qs = [FT(10^q) for q in log10qs]

    F_r = range(FT(0), stop = 0.9, length = n)

    N = FT(2 * 1e7)
    ρ_r = FT(700)

    Dms = zeros(n, n) 
    for i in 1:n 
        for j in 1:n 
            Dms[i, j] = P3.D_m(p3, qs[i], N, ρ_r, F_r[j])
        end
    end

    println(log10qs)
    println(qs) 
    println(Dms)

    min = 0 
    max = 0.006

    f = Plt.Figure()
    Plt.Axis(f[1, 1], xlabel = "F_r", ylabel = "log10(q)", title = "D_m") 
    Plt.heatmap!(F_r, log10qs, Dms, colormap = Plt.reverse(Plt.cgrad(:Spectral)), colorrange = (min, max))
    Plt.Colorbar(
        f[2, 1], 
        limits = (min, max), 
        colormap = Plt.reverse(Plt.cgrad(:Spectral)),
        flipaxis = true,
        vertical = false, 
    )

    Plt.save("testFigure.svg", f)

end

#plotDvsq()