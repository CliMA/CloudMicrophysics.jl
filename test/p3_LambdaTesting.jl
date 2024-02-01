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

function λ_diff(F_r::FT, ρ_r::FT, N::FT, λ::FT, p3::PSP3) 
    μ = P3.DSD_μ(p3, λ) 
    N_0 = P3.DSD_N₀(p3, N, λ)
    th = P3.thresholds(p3, ρ_r, F_r)
    q_calc = FT(P3.q_gamma(p3, F_r, N_0, log(λ), th))
    (λ_calculated, ) = P3.distribution_parameter_solver(p3, q_calc, N, ρ_r, F_r)
    return abs(λ - λ_calculated)
end

function constant_Fr_ρr(F_r::FT, ρ_r::FT, λ_min::FT, λ_max::FT)
    # lambda from 10 to 1e5 
    λ_range = range(λ_min, stop = λ_max, length = Int(abs((λ_min - λ_max)/100)))
    N = FT(1e8)
    lw = 3

    fig1 = Plt.Figure()
    ax1 = Plt.Axis(
        fig1[1, 1],
        title = Plt.L"λ errors in shape solver",
        xlabel = Plt.L"λ",
        ylabel = Plt.L"error",
        #xscale = Plt.log10,
        # yscale = Plt.log10,
        xminorticksvisible = true,
        #xminorticks = Plt.IntervalsBetween(5),
        yminorticksvisible = true,
        yminorticks = Plt.IntervalsBetween(10),
        #yticks = [1e-11, 1e-10, 1e-9, 1e-8],
        #aspect = 1.75,
        #limits = ((0.02, 10.0), nothing),
    )

    # for each lambda calculate the error given F_r and ρ_r 
    # fig1_0 = Plt.lines!(ax1, λ_range, [(first(P3.distribution_parameter_solver(p3, FT(P3.q_gamma(p3, F_r, P3.N_0_helper(N, λ, P3.μ_calc(λ)), log(λ), P3.μ_calc(λ), th)), N, ρ_r, F_r)) - λ) for λ in λ_range])
    fig1_0 = Plt.lines!(λ_range, [λ_diff(F_r, ρ_r, N, λ, p3) for λ in λ_range], linewidth = lw)

    println(λ_diff(F_r, ρ_r, N, FT(λ_min), p3))

    Plt.save("LambdaTesting1.svg", fig1)

end

function find_gaps_in_solver(F_r::FT, ρ_r::FT, λ_min::FT, λ_max::FT, p3::PSP3)
    N = FT(1e8)
    
    prev = FT(0)
    
    λ_range = range(λ_min, stop = λ_max, length = Int(abs((λ_min - λ_max)/100)))

    xs = Vector{FT}()
    ysFr = Vector{FT}()
    ysρr = Vector{FT}()

    for λ in λ_range 
        next = λ_diff(F_r, ρ_r, N, λ, p3)
        p = isequal(prev, NaN)
        n = isequal(next, NaN)
        # println("λ = ", λ, " diff = ", next)
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
    end
    if isequal(prev, NaN)
        append!(xs, λ_max)
        append!(ysFr, F_r)
        append!(ysρr, ρ_r)
    end
    return (xs, ysFr, ysρr) 
end

function plot_gaps(F_r_input::FT, ρ_r::FT, λ_min::FT, λ_max::FT, p3::PSP3)
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

    F_rs = range(FT(0), stop = FT(1-eps(FT)), length = 1000)
    for F_r in F_rs

        (xs, ysFr, ysρr) = find_gaps_in_solver(F_r, ρ_r, λ_min, λ_max, p3)

        Plt.linesegments!(ax1, xs, ysFr, color = :black)

    end

    # Plt.hlines!(ax1, F_r, xmin = min, xmax = max)

    Plt.save("LambdaTesting1.svg", fig1)


    fig2 = Plt.Figure()
    ax2 = Plt.Axis(
        fig2[1, 1],
        title = Plt.L" Working λ ",
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

    ρ_rs = range(FT(100), stop = FT(900), length = 1000)
    for ρ_r in ρ_rs
        (xs, ysFr, ysρr) = find_gaps_in_solver(F_r_input, ρ_r, λ_min, λ_max, p3)

        Plt.linesegments!(ax2, xs, ysρr, color = :black)
    end 

    Plt.save("LambdaTesting2.svg", fig2)
end

F_r = FT(0.8)
ρ_r = FT(400)
λ_min = eps(FT)
λ_max = FT(100000)

plot_gaps(F_r, ρ_r, λ_min, λ_max, p3)

#constant_Fr_ρr(FT(0.8), FT(400), FT(0), FT(2000))