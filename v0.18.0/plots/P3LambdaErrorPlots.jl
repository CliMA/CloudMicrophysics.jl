import CloudMicrophysics as CM
import CloudMicrophysics.P3Scheme as P3
import CloudMicrophysics.Parameters as CMP
import ClimaParams as CP
import SpecialFunctions as SF
import RootSolvers as RS
import CairoMakie as Plt

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

function λ_diff(F_r::FT, ρ_r::FT, N::FT, λ_ex::FT, p3::PSP3) where {FT}

    # Find the P3 scheme  thresholds
    th = P3.thresholds(p3, ρ_r, F_r)
    # Convert λ to ensure it remains positive
    x = log(λ_ex)
    # Compute mass density based on input shape parameters
    q_calc = P3.q_gamma(p3, F_r, N, x, th)

    (λ_calculated,) = P3.distribution_parameter_solver(p3, q_calc, N, ρ_r, F_r)
    return abs(λ_ex - λ_calculated)
end

function get_errors(
    p3::PSP3,
    λ_min::FT,
    λ_max::FT,
    F_r_min::FT,
    F_r_max::FT,
    ρ_r::FT,
    N::FT,
    λSteps::Int,
    F_rSteps::Int,
) where {FT}
    λs = range(FT(λ_min), stop = λ_max, length = λSteps)
    F_rs = range(F_r_min, stop = F_r_max, length = F_rSteps)
    E = zeros(λSteps, F_rSteps)
    min = Inf
    max = -Inf

    for i in 1:λSteps
        for j in 1:F_rSteps
            λ = λs[i]
            F_r = F_rs[j]

            diff = λ_diff(F_r, ρ_r, N, λ, p3)
            er = log(diff / λ)

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

function plot_relerrors(
    N::FT,
    λ_min::FT,
    λ_max::FT,
    F_r_min::FT,
    F_r_max::FT,
    ρ_r_min::FT,
    ρ_r_max::FT,
    λSteps::Int,
    F_rSteps::Int,
    numPlots::Int,
    p3::PSP3,
) where {FT}

    ρ_rs = range(ρ_r_min, stop = ρ_r_max, length = numPlots)

    f = Plt.Figure()

    x = 1
    y = 1
    for i in 1:numPlots

        ρ = ρ_rs[i]

        Plt.Axis(
            f[x, y],
            xlabel = "λ",
            ylabel = "F_r",
            title = string(
                "log(relative error calculated λ) for ρ_r = ",
                string(ρ),
            ),
            width = 400,
            height = 300,
        )

        (λs, F_rs, E, min, max) = get_errors(
            p3,
            λ_min,
            λ_max,
            F_r_min,
            F_r_max,
            ρ,
            N,
            λSteps,
            F_rSteps,
        )

        Plt.heatmap!(λs, F_rs, E)
        Plt.Colorbar(
            f[x, y + 1],
            limits = (min, max),
            colormap = :viridis,
            flipaxis = false,
        )

        y = y + 2
        if (y > 6)
            x = x + 1
            y = 1
        end
    end

    Plt.resize_to_layout!(f)
    Plt.save("P3LambdaHeatmap.svg", f)
end

# Define variables for heatmap relative error plots:

λ_min = FT(1e2)
λ_max = FT(1e6)
F_r_min = FT(0)
F_r_max = FT(1 - eps(FT))
ρ_r_min = FT(100)
ρ_r_max = FT(900)
N = FT(1e8)

λ_Steps = 100
F_r_Steps = 100
NumPlots = 9

plot_relerrors(
    N,
    λ_min,
    λ_max,
    F_r_min,
    F_r_max,
    ρ_r_min,
    ρ_r_max,
    λ_Steps,
    F_r_Steps,
    NumPlots,
    p3,
)
