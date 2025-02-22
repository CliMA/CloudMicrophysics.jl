import CairoMakie: Makie
import CloudMicrophysics as CM
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

FT = Float64

logλs = @. log(10.0^FT(3:0.01:5))

function make_slope_plot(slope_law, title)
    fig = Makie.Figure(size = (400, 300), figure_padding = 20)

    ax = Makie.Axis(fig[1, 1]; title, xscale = log10, xlabel = "λ", ylabel = "μ")

    Makie.lines!(ax, exp.(logλs), P3.get_μ.(slope_law, logλs))

    return fig
end


# Power law parameterization
slope_power_law = CMP.ParametersP3(FT).slope

fig = Makie.with_theme(Makie.theme_minimal()) do
    make_slope_plot(slope_power_law, "μ as a function of λ (power law)")
end
Makie.save("P3SlopeParameterizations_power_law.svg", fig)

# Constant parameterization
# params = CMP.ParametersP3(FT; slope_law = :constant)
# slope_constant = params.slope

# fig = Makie.with_theme(Makie.theme_minimal()) do
#     make_slope_plot(slope_constant, "μ as a constant")
# end
# save("P3SlopeParameterizations_constant.svg", fig)
