import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP
import SpecialFunctions as SF
import RootSolvers as RS
import CairoMakie as CMK

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)
F_liq = FT(0) # preserving original P3
# TODO: investigate for F_liq != 0

function λ_diff(F_rim::FT, ρ_r::FT, N::FT, λ_ex::FT, p3::PSP3) where {FT}

    # Find the P3 scheme  thresholds
    th = P3.thresholds(p3, ρ_r, F_rim)
    # Get μ corresponding to λ
    μ = P3.DSD_μ(p3, λ_ex)
    # Convert λ to ensure it remains positive
    x = log(λ_ex)
    # Compute mass density based on input shape parameters
    L_calc = N * P3.L_over_N_gamma(p3, F_rim, F_liq, x, μ, th)

    (λ_calculated,) =
        P3.distribution_parameter_solver(p3, L_calc, N, ρ_r, F_rim, F_liq)
    return abs(λ_ex - λ_calculated)
end

function get_errors(
    p3::PSP3,
    log10_λ_min::FT,
    log10_λ_max::FT,
    F_rim_min::FT,
    F_rim_max::FT,
    ρ_r::FT,
    N::FT,
    λSteps::Int,
    F_rimSteps::Int,
) where {FT}

    logλs = range(FT(log10_λ_min), stop = log10_λ_max, length = λSteps)
    λs = [10^logλ for logλ in logλs]
    F_rims = range(F_rim_min, stop = F_rim_max, length = F_rimSteps)
    E = zeros(λSteps, F_rimSteps)

    for i in 1:λSteps
        for j in 1:F_rimSteps
            λ = λs[i]
            F_rim = F_rims[j]

            er = log(λ_diff(F_rim, ρ_r, N, λ, p3) / λ)
            er = er == Inf ? 9999.99 : er
            er = er == -Inf ? -9999.99 : er
            E[i, j] = er
        end
    end
    return (; λs, F_rims, E)
end

#! format: off
function plot_relerrors(
    N::FT,
    log10_λ_min::FT,
    log10_λ_max::FT,
    F_rim_min::FT,
    F_rim_max::FT,
    ρ_r_min::FT,
    ρ_r_max::FT,
    λSteps::Int,
    F_rimSteps::Int,
    numPlots::Int,
    p3::PSP3,
) where {FT}

    ρ_rs = range(ρ_r_min, stop = ρ_r_max, length = numPlots)

    f = CMK.Figure()
    x = 1
    y = 1
    for i in 1:numPlots
        i = 1

        ρ = ρ_rs[i]
        (λs, F_rims, E) = get_errors(p3, log10_λ_min, log10_λ_max, F_rim_min, F_rim_max, ρ, N, λSteps, F_rimSteps)

        CMK.Axis(f[x, y], xlabel = "λ", ylabel = "F_rim", title = string("log(relative error calculated λ) for ρ_r = ", string(ρ)), width = 400, height = 300, xscale = log10)
        hm = CMK.heatmap!(λs, F_rims, E, colormap = CMK.cgrad(:viridis, 20, categorical=true),  colorrange = (-10, 0), highclip = :red, lowclip = :indigo)
        CMK.Colorbar(f[x, y + 1], hm)

        y = y + 2
        if (y > 6)
            x = x + 1
            y = 1
        end
    end
    CMK.resize_to_layout!(f)
    CMK.save("P3LambdaHeatmap.png", f)
end
#! format: on

function μ_approximation_effects(F_rim::FT, ρ_r::FT) where {FT}

    f = CMK.Figure(size = (800, 500))

    ax1 = CMK.Axis(
        f[1, 1],
        xlabel = "L/N",
        ylabel = "μ",
        title = string("μ vs L/N for F_rim = ", F_rim, " ρ_r = ", ρ_r),
        width = 400,
        height = 300,
        xscale = log10,
    )

    ax2 = CMK.Axis(
        f[1, 2],
        xlabel = "λ",
        ylabel = "μ",
        title = string("μ vs λ for F_rim = ", F_rim, " ρ_r = ", ρ_r),
        width = 400,
        height = 300,
        xscale = log10,
    )

    ax3 = CMK.Axis(
        f[1, 3],
        xlabel = "λ",
        ylabel = "L/N",
        title = string("L/N vs λ for F_rim = ", F_rim, " ρ_r = ", ρ_r),
        width = 400,
        height = 300,
        xscale = log10,
        yscale = log10,
    )

    CMK.linkxaxes!(ax2, ax3)
    CMK.linkyaxes!(ax1, ax2)

    numpts = 100

    # Set up vectors
    th = P3.thresholds(p3, ρ_r, F_rim)
    log_λs = range(FT(3.6), stop = FT(4.6), length = numpts)
    λs = [10^log_λ for log_λ in log_λs]
    μs = [P3.DSD_μ(p3, λ) for λ in λs]

    μs_approx = [FT(0) for λ in λs]
    Ls = [FT(0) for λ in λs]
    λ_solved = [FT(0) for λ in λs]

    for i in 1:numpts
        L_over_N = P3.L_over_N_gamma(p3, F_rim, F_liq, log(λs[i]), μs[i], th)
        Ls[i] = L_over_N
        N = FT(1e6)
        (L, N) = P3.distribution_parameter_solver(
            p3,
            L_over_N * N,
            N,
            ρ_r,
            F_rim,
            F_liq,
        )
        λ_solved[i] = L
        μs_approx[i] = P3.DSD_μ_approx(p3, N * L_over_N, N, ρ_r, F_rim, F_liq)
    end

    # Plot
    CMK.lines!(ax3, λs, Ls, label = "true distribution")
    CMK.lines!(ax3, λ_solved, Ls, label = "approximated")

    CMK.lines!(ax1, Ls, μs, label = "true distribution")
    CMK.lines!(ax1, Ls, μs_approx, label = "approximated")
    CMK.lines!(ax2, λs, μs, label = "true distribution")
    CMK.lines!(ax2, λ_solved, μs_approx, label = "approximated")

    CMK.axislegend(ax1, position = :lb)
    CMK.axislegend(ax2, position = :lt)
    CMK.axislegend(ax3, position = :lb)

    CMK.resize_to_layout!(f)
    CMK.save("MuApprox.svg", f)

end

# Define variables for heatmap relative error plots:

log10_λ_min = FT(2)
log10_λ_max = FT(6)
F_rim_min = FT(0)
F_rim_max = FT(0.9)
ρ_r_min = FT(100)
ρ_r_max = FT(900)
N = FT(1e8)

λ_Steps = 40
F_rim_Steps = 40
NumPlots = 9

plot_relerrors(
    N,
    log10_λ_min,
    log10_λ_max,
    F_rim_min,
    F_rim_max,
    ρ_r_min,
    ρ_r_max,
    λ_Steps,
    F_rim_Steps,
    NumPlots,
    p3,
)

F_rim = FT(0.5)
ρ_r = FT(500)

μ_approximation_effects(F_rim, ρ_r)
