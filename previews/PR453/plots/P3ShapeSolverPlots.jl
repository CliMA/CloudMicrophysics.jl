import CairoMakie as Plt
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

const PSP3 = CMP.ParametersP3
p3 = CMP.ParametersP3(FT)
F_liq = FT(0) # preserving original P3
# TODO: investigate for F_liq != 0

function guess_value(λ::FT, p1::FT, p2::FT, L1::FT, L2::FT)
    return L1 * (λ / p1)^((log(L1) - log(L2)) / (log(p1) - log(p2)))
end

function lambda_guess_plot()
    N = FT(1e6)

    F_rim_s = [FT(0.0), FT(0.5), FT(0.8)]
    ρ_r_s = [FT(200), FT(400), FT(800)]

    f = Plt.Figure()

    for i in 1:length(F_rim_s)
        for j in 1:length(ρ_r_s)
            F_rim = F_rim_s[i]
            ρ_r = ρ_r_s[j]

            Plt.Axis(
                f[i, j],
                xlabel = "log(L/N)",
                ylabel = "log(λ)",
                title = string(
                    "λ vs L/N for F_rim = ",
                    F_rim,
                    " and ρ_r = ",
                    ρ_r,
                ),
                height = 300,
                width = 400,
            )

            logλs = FT(1):FT(0.01):FT(6)
            λs = [10^logλ for logλ in logλs]
            th = P3.thresholds(p3, ρ_r, F_rim)
            Ls_over_N = [
                P3.L_over_N_gamma(
                    p3,
                    F_rim,
                    F_liq,
                    log(λ),
                    P3.DSD_μ(p3, λ),
                    th,
                ) for λ in λs
            ]
            guesses = [FT(0) for λ in λs]

            for i in 1:length(λs)
                (min,) = P3.get_bounds(
                    N,
                    Ls_over_N[i] * N,
                    P3.DSD_μ_approx(p3, Ls_over_N[i] * N, N, ρ_r, F_rim, F_liq),
                    F_rim,
                    F_liq,
                    p3,
                    th,
                )
                guesses[i] = exp(min)
            end


            Plt.lines!(
                log10.(Ls_over_N),
                log10.(λs),
                linewidth = 3,
                color = "Black",
                label = "true",
            )
            Plt.lines!(
                log10.(Ls_over_N),
                log10.(guesses),
                linewidth = 2,
                linestyle = :dash,
                color = "Red",
                label = "approximated",
            )

            Plt.axislegend("Legend", position = :lb)
        end
    end

    Plt.resize_to_layout!(f)
    Plt.save("SolverInitialGuess.svg", f)

end

lambda_guess_plot()
