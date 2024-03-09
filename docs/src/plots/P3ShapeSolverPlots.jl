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

function lambda_guess_plot()
    N = FT(1e6)

    F_r_s = [FT(0.0), FT(0.5), FT(0.8)]
    ρ_r_s = [FT(200), FT(400), FT(800)]

    f = Plt.Figure()

    for i in 1:length(F_r_s)
        for j in 1:length(ρ_r_s)
            F_r = F_r_s[i]
            ρ_r = ρ_r_s[j]

            Plt.Axis(
                f[i, j],
                xlabel = "log(q/N)",
                ylabel = "log(λ)",
                title = string("λ vs q/N for F_r = ", F_r, " and ρ_r = ", ρ_r),
                height = 300,
                width = 400,
            )


            logλs = FT(1):FT(0.01):FT(6)
            λs = [10^logλ for logλ in logλs]
            th = P3.thresholds(p3, ρ_r, F_r)
            qs = [
                P3.q_over_N_gamma(p3, F_r, log(λ), P3.DSD_μ(p3, λ), th) for
                λ in λs
            ]
            guesses = [FT(0) for λ in λs]

            for i in 1:length(λs)
                (min,) = P3.get_bounds(
                    N,
                    qs[i] * N,
                    P3.DSD_μ_approx(p3, qs[i] * N, N, ρ_r, F_r),
                    F_r,
                    p3,
                    th,
                )
                guesses[i] = exp(min)
            end


            Plt.lines!(
                log10.(qs),
                log10.(λs),
                linewidth = 3,
                color = "Black",
                label = "true",
            )
            Plt.lines!(
                log10.(qs),
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
