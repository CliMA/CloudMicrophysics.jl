import CairoMakie: Makie
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.P3Scheme as P3

# Direct evaluation of the analytical solution for ρ_d, which loses accuracy as F_rim → 0.
function ρ_d_direct(β_va, F_rim, ρ_rim)
    p = 1 / (3 - β_va)
    Fᵤ = 1 - F_rim
    k = Fᵤ^(-p)
    den = (β_va - 2) * (k - 1) / (Fᵤ * k - 1) - Fᵤ
    return ρ_rim * F_rim / den
end

mass = CMP.ParametersP3(Float32).mass
β_va = mass.β_va
ρ_rim = 400

ρ_g(ρ_d, F_rim) = F_rim * ρ_rim + (1 - F_rim) * ρ_d
F_rims = [10.0^e for e in -7:0.05:-0.005]
ρ_g_ref(F_rim) = Float64(ρ_g(ρ_d_direct(big(β_va), big(F_rim), big(ρ_rim)), big(F_rim)))
rel_err(ρ_g_f32, F_rim) = abs(Float64(ρ_g_f32) - ρ_g_ref(F_rim)) / abs(ρ_g_ref(F_rim))

err_stable = [rel_err(ρ_g(P3.get_ρ_d(mass, Float32(F), Float32(ρ_rim)), Float32(F)), F) for F in F_rims]
err_direct = [rel_err(ρ_g(ρ_d_direct(β_va, Float32(F), Float32(ρ_rim)), Float32(F)), F) for F in F_rims]

fig = Makie.Figure(size = (760, 460))
ax = Makie.Axis(
    fig[1, 1];
    xscale = log10,
    yscale = log10,
    xlabel = "F_rim",
    ylabel = "relative error of ρ_g (Float32)",
    title = "Direct evaluation vs. numerically stable rewrite",
)
Makie.lines!(ax, F_rims, max.(err_direct, 1e-16), label = "direct evaluation", linewidth = 2)
Makie.lines!(ax, F_rims, max.(err_stable, 1e-16), label = "stable rewrite", linewidth = 2)
Makie.hlines!(ax, [eps(Float32)], color = :gray, linestyle = :dash, label = "eps(Float32)")
Makie.axislegend(ax; position = :rt)
Makie.save("P3RhoDStability.svg", fig)
