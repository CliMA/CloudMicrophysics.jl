import CairoMakie as MK

import SpecialFunctions as SF

import Thermodynamics as TD
import CloudMicrophysics as CM
import ClimaParams as CP

import Thermodynamics.Parameters as TDP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.Microphysics1M as CM1

FT = Float64

function Marshall_Palmer_distribution(
    (; pdf, mass)::CMP.Rain{FT},
    q::FT,
    ρ::FT,
    r::FT,
) where {FT}

    n₀::FT = CM1.get_n0(pdf)
    λ = CM1.lambda(pdf, mass, q, ρ)
    # distribution
    n_r = n₀ * exp(-λ * r)

    return n_r
end

ρ = FT(1.01)
qᵣ_0 = FT(1.0e-5)
qᵣ_1 = FT(2.0e-4)
qᵣ_2 = FT(1.0e-4)
qᵣ_3 = FT(1.0e-3)

r_range = range(25e-6, stop = 1e-2, length = 1000)

n_r_0 =
    [Marshall_Palmer_distribution(CMP.Rain(FT), qᵣ_0, ρ, r) for r in r_range]
n_r_1 =
    [Marshall_Palmer_distribution(CMP.Rain(FT), qᵣ_1, ρ, r) for r in r_range]
n_r_2 =
    [Marshall_Palmer_distribution(CMP.Rain(FT), qᵣ_2, ρ, r) for r in r_range]
n_r_3 =
    [Marshall_Palmer_distribution(CMP.Rain(FT), qᵣ_3, ρ, r) for r in r_range]

fig = MK.Figure(resolution = (1100, 600))
ax1 = MK.Axis(
    fig[1, 1],
    title = "Marshall-Palmer distribution",
    xlabel = "r [m]",
    ylabel = "n(r) [mm/m]",
    limits = ((0, 0.0015), nothing),
)
MK.lines!(ax1, r_range, n_r_0, label = "q⁰ᵣ = 10⁻⁵", color = :blue)
MK.lines!(ax1, r_range, n_r_1, label = "q¹ᵣ = 2x10⁻⁴", color = :red)
MK.lines!(ax1, r_range, n_r_2, label = "q²ᵣ = 10⁻⁴", color = :green)
MK.lines!(ax1, r_range, n_r_3, label = "q³ᵣ = 10⁻³", color = :orange)
ax2 = MK.Axis(
    fig[1, 2],
    title = "Marshall-Palmer distribution - semilog scale",
    ylabel = "n(r) [mm/m]",
    xlabel = "r [m]",
    xscale = log10,
)
MK.lines!(ax2, r_range, n_r_0, label = "q⁰ᵣ = 10⁻⁵", color = :blue)
MK.lines!(ax2, r_range, n_r_1, label = "q¹ᵣ = 2x10⁻⁴", color = :red)
MK.lines!(ax2, r_range, n_r_2, label = "q²ᵣ = 10⁻⁴", color = :green)
MK.lines!(ax2, r_range, n_r_3, label = "q³ᵣ = 10⁻³", color = :orange)
MK.axislegend()

MK.save("MarshallPalmer_distribution.svg", fig)
#! format: on
