import CairoMakie as Plt
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)

function guess_value(λ::FT, p1::FT, p2::FT, q1::FT, q2::FT)
    return q1 * (λ / p1)^((log(q1) - log(q2)) / (log(p1) - log(p2)))
end

function lambda_guess_plot(F_r::FT, ρ_r::FT) where {FT}
    N = FT(1e8)

    λs = FT(1e2):FT(1e2):FT(1e6 + 1)
    th = P3.thresholds(p3, ρ_r, F_r)
    qs = [P3.q_gamma(p3, F_r, N, log(λ), th) for λ in λs]

    guesses = [guess_value(λ, λs[1], last(λs), qs[1], last(qs)) for λ in λs]

    f = Plt.Figure()
    Plt.Axis(
        f[1, 1],
        xscale = log,
        yscale = log,
        xticks = [10^2, 10^3, 10^4, 10^5, 10^6],
        yticks = [10^3, 1, 10^-3, 10^-6],
        xlabel = "λ",
        ylabel = "q",
        title = "q vs λ",
        height = 300,
        width = 400,
    )

    l1 = Plt.lines!(λs, qs, linewidth = 3, color = "Black", label = "q")
    l2 = Plt.lines!(
        λs,
        guesses,
        linewidth = 2,
        linestyle = :dash,
        color = "Red",
        label = "q_approximated",
    )

    Plt.axislegend("Legend", position = :lb)

    Plt.resize_to_layout!(f)
    Plt.save("SolverInitialGuess.svg", f)

end

lambda_guess_plot(FT(0.5), FT(200))
