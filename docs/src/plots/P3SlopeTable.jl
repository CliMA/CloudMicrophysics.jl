import CairoMakie: CairoMakie, Makie
CairoMakie.activate!(type = "svg")

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3
FT = Float64

params = CMP.ParametersP3(FT)
ρ_rim = FT(400)
logλ_lo, logλ_hi = P3.logλ_bounds(params)
logλs = range(logλ_lo, logλ_hi; length = 400)

fig = Makie.with_theme(Makie.theme_minimal()) do
    fig = Makie.Figure(size = (900, 340), figure_padding = 20)
    ax_mass = Makie.Axis(fig[1, 1];
        title = "Slope table", xscale = log10, yscale = log10,
        xlabel = "λ (m⁻¹)", ylabel = "mean particle mass L/N (kg)",
    )
    ax_μ = Makie.Axis(fig[1, 2];
        title = "Shape parameter", xscale = log10,
        xlabel = "mean particle mass L/N (kg)", ylabel = "μ",
    )
    for (F_rim, color) in ((0, :black), (0.5, :royalblue), (0.9, :darkorange))
        particles = P3.P3Particles(params, FT(F_rim), ρ_rim)
        x̄ = exp.(P3.logLdivN.(Ref(particles), logλs))
        μ = P3.get_μ.(Ref(params.slope), logλs)
        Makie.lines!(ax_mass, exp.(logλs), x̄; color, label = "F_rim = $F_rim")
        Makie.lines!(ax_μ, x̄, μ; color, label = "F_rim = $F_rim")
    end
    Makie.vspan!(ax_μ, exp(-20.2), exp(-19.8); color = (:gray, 0.2))
    Makie.axislegend(ax_mass; position = :rt)
    Makie.axislegend(ax_μ; position = :rt)
    fig
end
Makie.save("P3SlopeTable.svg", fig)
