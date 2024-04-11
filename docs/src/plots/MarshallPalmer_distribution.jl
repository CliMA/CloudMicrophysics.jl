import CairoMakie as MK

import SpecialFunctions as SF

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics1M as CM1

FT = Float64

function Marshall_Palmer_distribution(
    (; pdf, mass)::CMP.Rain{FT},
    q::FT,
    ρ::FT,
    r::FT,
)where{FT}

    n₀::FT = CM1.get_n0(pdf)
    λ = CM1.lambda(pdf, mass, q, ρ)
    # distribution
    n_r = n₀ * exp(-λ*r)

    return n_r
end

ρ = FT(1.01)
qᵣ_0 = FT(1.0e-4)
qᵣ_1 = FT(1.96e-4)
qᵣ_2 = FT(1.3e-4)
qᵣ_3 = FT(1.6e-4)

r_range = range(25e-6, stop = 1e-2, length = 1000)

n_r_0 = [Marshall_Palmer_distribution(rain, qᵣ_0, ρ, r) for r in r_range]
n_r_1 = [Marshall_Palmer_distribution(rain, qᵣ_1, ρ, r) for r in r_range]
n_r_2 = [Marshall_Palmer_distribution(rain, qᵣ_2, ρ, r) for r in r_range]
n_r_3 = [Marshall_Palmer_distribution(rain, qᵣ_3, ρ, r) for r in r_range]

fig = MK.Figure(resolution = (800, 600))
ax1 = MK.Axis(fig[1, 1], title = "Marshall-Palmer distribution", ylabel = "r [μm]", xlabel = "n(r) [mm/m]", limits = ((25e-6, 1e-2), nothing))
MK.lines!(ax1, n_r_0, r_range .* -6, label = "CM default a_w", color = :blue)
MK.lines!(ax1, n_r_1, r_range .* -6, label = "a_w using p(0,T)", color = :red)
MK.lines!(ax1, n_r_2, r_range .* -6, label = "CM default a_w_ice", color = :green)
MK.lines!(ax1, n_r_3, r_range .* -6, label = "a_w_ice using p(0,T)", color = :orange)
MK.axislegend()

MK.save("MarshallPalmer_distribution.svg", fig)
#! format: on