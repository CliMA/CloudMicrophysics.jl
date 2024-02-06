import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import CLIMAParameters as CP
import Integrals as IN
import SpecialFunctions as SF
import RootSolvers as RS
import CairoMakie as Plt

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

function λ_diff(F_r::FT, ρ_r::FT, N::FT, λ_ex::FT, p3::PSP3) 
    # Find the P3 scheme  thresholds
    th = P3.thresholds(p3, ρ_r, F_r)
    # Convert λ to ensure it remains positive
    x = log(λ_ex)
    # Compute mass density based on input shape parameters
    q_calc = P3.q_gamma(p3, F_r, N, x, th)  
    (λ_calculated, N_0, converged) = P3.distribution_parameter_solver(p3, q_calc, N, ρ_r, F_r,)
    return (diff = abs(λ_ex - λ_calculated), converged = converged)
end

function find_gaps_in_solver(F_r::FT, ρ_r::FT, λ_min::FT, λ_max::FT, p3::PSP3, numSteps)
    N = FT(1e8)
    
    prev = FT(0)

    λ_range = range(λ_min, stop = λ_max, length = numSteps)

    xs = Vector{FT}()
    ysFr = Vector{FT}()
    ysρr = Vector{FT}()

    for λ in λ_range 
        # println("λ = ", λ, " F_r = ", F_r, " ρ_r = ", ρ_r)

        next = λ_diff(F_r, ρ_r, N, λ, p3)
        p = isequal(prev, NaN)
        n = isequal(next, NaN)
        #println("λ = ", λ, " diff = ", next)
        if n > p
            # println("min = ", λ)
            append!(xs, λ)
            append!(ysFr, F_r)
            append!(ysρr, ρ_r)
        elseif p > n 
            # println("     max = ", λ)
            append!(xs, λ)
            append!(ysFr, F_r)
            append!(ysρr, ρ_r)
        end
        prev = next

        #println("done, λ_difference = ", next)
        #println(" ")
    end
    if isequal(prev, NaN)
        append!(xs, λ_max)
        append!(ysFr, F_r)
        append!(ysρr, ρ_r)
    end
    return (xs, ysFr, ysρr) 
end

function plot_gaps(F_r_input::FT, ρ_r::FT, λ_min::FT, λ_max::FT, p3::PSP3, numSteps)
    fig1 = Plt.Figure()
    ax1 = Plt.Axis(
        fig1[1, 1],
        title = Plt.L" Working λ ",
        xlabel = Plt.L" λ ",
        ylabel = Plt.L" F_r ",
        #xscale = Plt.log10,
        # yscale = Plt.log10,
        xminorticksvisible = true,
        #xminorticks = Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        #yminorticks = Plt.IntervalsBetween(10),
        #yticks = [1e-11, 1e-10, 1e-9, 1e-8],
        #aspect = 1.75,
        #limits = ((0.02, 10.0), nothing),
    )

    Plt.xlims!(low = 0)
    Plt.xlims!(high = λ_max)
    Plt.ylims!(low = 0) 
    Plt.ylims!(high = 1)

    F_rs = range(FT(0), stop = 1-eps(FT), length = 1001)   # 1-eps(FT)
    for F_r in F_rs

        (xs, ysFr, ysρr) = find_gaps_in_solver(F_r, ρ_r, λ_min, λ_max, p3, numSteps)
        println("xs = ", xs)

        Plt.linesegments!(ax1, xs, ysFr, color = :black)

    end

    # Plt.hlines!(ax1, F_r, xmin = min, xmax = max)

    Plt.save("LambdaTesting4.svg", fig1) 


    fig2 = Plt.Figure()
    ax2 = Plt.Axis(
        fig2[1, 1],
        title = Plt.L" Working λ for F_r = 0.4",
        xlabel = Plt.L" λ ",
        ylabel = Plt.L" ρ_r ",
        #xscale = Plt.log10,
        # yscale = Plt.log10,
        xminorticksvisible = true,
        #xminorticks = Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        #yminorticks = Plt.IntervalsBetween(10),
        #yticks = [1e-11, 1e-10, 1e-9, 1e-8],
        #aspect = 1.75,
        #limits = ((0.02, 10.0), nothing),
    )
    Plt.ylims!(low = 0)  
    Plt.ylims!(high = 1000)

    ρ_rs = range(FT(100), stop = FT(900), length = 500)
    for ρ_r in ρ_rs
        (xs, ysFr, ysρr) = find_gaps_in_solver(F_r_input, ρ_r, λ_min, λ_max, p3, numSteps)

        Plt.linesegments!(ax2, xs, ysρr, color = :black)
    end 

    Plt.save("LambdaTesting3.svg", fig2) 
end

function plot_all_gaps(λ_min::FT, λ_max::FT, F_r_min::FT, F_r_max::FT, ρ_r_min::FT, ρ_r_max::FT, p3::PSP3, λSteps, F_rSteps, numPlots) where {FT}
    f = Plt.Figure()
    
    λs = range(λ_min, stop = λ_max, length = λSteps)
    F_rs = range(F_r_min, stop = F_r_max, length = F_rSteps)
    #ρ_rs = (FT(100), FT(200), FT(300), FT(400), FT(500), FT(600), FT(700), FT(800))
    #l = length(ρ_rs)
    ρ_rs = range(ρ_r_min, stop = ρ_r_max, length = numPlots)

    for i in 1:numPlots
        ρ_r = ρ_rs[i]

        j = Int(ceil(i/2))
        k = Int(mod(i - 1, 2) + 1)
        Plt.Axis(f[i, 1], xlabel = "λ", ylabel = "F_r", title = string("ρ_r = " , string(ρ_r)), width = 400, height = 300)
        
        for F_r in F_rs
            (xs, ysFr,) = find_gaps_in_solver(F_r, ρ_r, λ_min, λ_max, p3, λSteps)

            Plt.linesegments!(f[i, 1], xs, ysFr, color = :black)
        end
    end
    Plt.resize_to_layout!(f)
    Plt.save("allFrsbig.svg", f)
end

function rel_error(λ, F_r)
    ρ_r = FT(400) 
    N = FT(1e8) 
    p3 = CMP.ParametersP3(FT)

    er = log(λ_diff(F_r, ρ_r, N, λ, p3)/λ)
    if isequal(er, NaN)
        er = Inf
    end
    return er
end

function plot_relerrors(λ_min::FT, λ_max::FT, p3::PSP3, λSteps, F_rSteps, numPlots) where{FT}
    N = FT(1e8) 

    #λs = range(λ_min, stop = λ_max, length = λSteps)
    #F_rs = range(FT(0), stop = 1-eps(FT), length = F_rSteps)
    # ρ_rs = range(FT(10), stop = 990, length = numPlots)
    ρ_r = FT(400)

    #= xs = Vector{FT}()
    ys = Vector{FT}()
    zs = Vector{FT}()

    for λ in λs 
        for F_r in F_rs
            append!(xs, λ)
            append!(ys, F_r)
            er = log(λ_diff(F_r, ρ_r, N, λ, p3)/λ)
            if isequal(er, NaN)
                er = Inf
            end
            append!(zs, er)
        end
    end =#

    f = Plt.Figure()
    Plt.Axis(f[1, 1])

    (λs, F_rs, E, min, max) = get_errors(λ_min, λ_max, p3, ρ_r, N, λSteps, F_rSteps)

    println(min)
    println(max)

    Plt.heatmap!(λs, F_rs, E)
    Plt.Colorbar(f[1, 2], limits = (min, max), colormap = :viridis, flipaxis = false)

    # Plt.Colorbar(f[1, 2], label = "error", limits = FT.(range_val))

    Plt.save("ContourAttempt1.svg", f)

end

function get_errors(λ_min::FT, λ_max::FT, p3::PSP3, ρ_r::FT, N::FT, λSteps, F_rSteps) where{FT}
    λs = range(λ_min, stop = λ_max, length = λSteps)
    F_rs = range(FT(0), stop = 1-eps(FT), length = F_rSteps)
    E = zeros(λSteps, F_rSteps)
    min = Inf
    max = -Inf

    for i in 1:λSteps 
        for j in 1:F_rSteps 
            λ = λs[i] 
            F_r = F_rs[j] 

            er = log(λ_diff(F_r, ρ_r, N, λ, p3)/λ)
            #= if isequal(er, NaN)
                er = Inf
            end =#

            E[i, j] = er

            if er > max && er < Inf
                max = er
            end
            if er < min && er > -Inf
                min = er
            end

        end
    end
    return (λs = λs, F_rs = F_rs, E = E, min = min, max = max)
end

function buckets(λ_min::FT, λ_max::FT, F_r_min::FT, F_r_max::FT, p3::PSP3, ρ_r::FT, N::FT, λSteps, F_rSteps) where{FT}
    λs = range(λ_min, stop = λ_max, length = λSteps)
    F_rs = range(F_r_min, stop = F_r_max, length = F_rSteps)
    E = zeros(λSteps, F_rSteps)
    min_error = eps(FT) 
    max_error = sqrt(eps(FT))

    for i in 1:λSteps 
        for j in 1:F_rSteps 
            λ = λs[i] 
            F_r = F_rs[j] 

            (error, converged) = λ_diff(F_r, ρ_r, N, λ, p3)

            er = log(error/λ)

            if converged
                if er <= min_error
                    E[i, j] = FT(0)
                elseif er <= max_error 
                    E[i, j] = FT(1)
                elseif er <= FT(1) 
                    E[i, j] = FT(2)
                else 
                    E[i, j] = FT(3)
                end
            else 
                if er <= min_error
                    E[i, j] = FT(0.5)
                elseif er <= max_error 
                    E[i, j] = FT(1.5)
                elseif er <= FT(1) 
                    E[i, j] = FT(2.5)
                else 
                    E[i, j] = FT(3.5)
                end
            end
        end
    end
    return (λs = λs, F_rs = F_rs, E = E)
end

function plot_buckets(λ_min::FT, λ_max::FT, F_r_min::FT, F_r_max::FT, ρ_r_min::FT, ρ_r_max::FT,p3::PSP3, λSteps, F_rSteps, numPlots) where{FT}
    N = FT(1e8) 
    ρ_rs = range(ρ_r_min, stop = ρ_r_max, length = numPlots)

    f = Plt.Figure()
    for i in 1:numPlots
        ρ_r = ρ_rs[i]

        Plt.Axis(f[i, 1], xlabel = "λ", ylabel = "F_r", title = string("ρ_r = " , string(ρ_r)), width = 400, height = 300)

        (λs, F_rs, E) = buckets(λ_min, λ_max, F_r_min, F_r_max, p3, ρ_r, N, λSteps, F_rSteps)

        Plt.heatmap!(λs, F_rs, E)
        Plt.Colorbar(f[i, 2], limits = (0, 3.5), colormap = :viridis, flipaxis = false)

        # Plt.Colorbar(f[1, 2], label = "error", limits = FT.(range_val))
    end

    Plt.resize_to_layout!(f)
    Plt.save("ContourAttempt2.svg", f)

end

F_r = FT(0.4)
ρ_r = FT(400)
λ_min = FT(1)       # anything above 400 giving errors that λ <= FT(0)
λ_max = FT(1e5)
F_r_min = FT(0.0)
F_r_max = FT(1-eps(FT))
ρ_r_min = FT(100)
ρ_r_max = FT(900)

# plot_gaps(F_r, ρ_r, λ_min, λ_max, p3, 250)

#constant_Fr_ρr(FT(0.8), FT(400), FT(0), FT(2000))

#plot_all_gaps(λ_min, λ_max, F_r_min, F_r_max, ρ_r_min, ρ_r_max, p3, 200, 200, 10)

#plot_relerrors(λ_min, λ_max, p3, 200, 200, 1)

F = FT(0.34) 
ρ = FT(650)
λ = FT(35000)
N = FT(1e8) 

#diff = λ_diff(F, ρ, N, λ, p3)

#println(λ_diff(FT(0.5), FT(400), FT(1e8), FT(15000), p3))
plot_buckets(λ_min, λ_max, F_r_min, F_r_max, ρ_r_min, ρ_r_max, p3, 250, 250, 9)