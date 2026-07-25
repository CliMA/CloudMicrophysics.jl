import CairoMakie as MK
import CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams

FT = Float64

tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
mp  = CMP.Microphysics1MParams(FT)

q_min, q_max = 1e-8, 1e-2          # wider range to capture threshold sweep
q_icl_range  = range(q_min, stop = q_max, length = 500)
ρ_air = 1.2
q_lcl, q_tot = 5e-4, 20e-3
q_rai = 1e-3
q_sno = 1e-4
T     = 255.0
limits = (0, q_max * 1e3, 0, nothing)

MK.set_theme!(MK.theme_minimal())

acnv0 = mp.options.snow_autoconversion.acnv1M
#τ₀, q_thr₀, k₀ = acnv0.τ, acnv0.q_threshold, acnv0.k
τ₀, k₀ = acnv0.τ, acnv0.k
q_thr₀ = 1e-3
@info "Default snow autoconversion parameters" τ₀ q_thr₀ k₀

acnv_rate(acnv) = map(q_icl_range) do q
    CM1.conv_q_icl_to_q_sno(
        CMP.NoSupersaturation(acnv), mp, tps,
        (; q_tot, q_lcl, q_icl = q, q_rai, q_sno),
        (; ρ = ρ_air, T),
    )
end

# much wider sweeps so the effect is actually visible
mult_τ    = FT(5)   .^ range(-1, 1, length = 5)   # timescale: 0.2x – 5x
mult_qthr = FT(10)  .^ range(-1, 1, length = 5)   # threshold: 0.1x – 10x
mult_k    = FT(100) .^ range(-1, 1, length = 5)   # steepness: 0.01x – 100x

τ_values     = τ₀     .* mult_τ
q_thr_values = q_thr₀ .* mult_qthr
k_values     = k₀     .* mult_k

fig = MK.Figure(size = (1350, 430))
axs = [MK.Axis(fig[1, i]; xlabel = "q_icl [g/kg]", limits) for i in 1:3]
axs[1].ylabel = "autoconversion rate [kg/kg/s]"
axs[1].title  = "τ varied (q_thr = $(round(q_thr₀ * 1e3, sigdigits=3)) g/kg, k = $(round(k₀, sigdigits=3)))"
axs[2].title  = "q_thr varied (τ = $(round(τ₀, sigdigits=3)) s, k = $(round(k₀, sigdigits=3)))"
axs[3].title  = "k varied (τ = $(round(τ₀, sigdigits=3)) s, q_thr = $(round(q_thr₀ * 1e3, sigdigits=3)) g/kg)"

for (τ, m) in zip(τ_values, mult_τ)
    MK.lines!(
        axs[1], q_icl_range * 1e3, acnv_rate(CMP.Acnv1M{FT}(τ, q_thr₀, k₀));
        label = "τ = $(round(τ, sigdigits=3)) s ($(round(m, digits=2))x)",
    )
end
for (q_thr, m) in zip(q_thr_values, mult_qthr)
    MK.lines!(
        axs[2], q_icl_range * 1e3, acnv_rate(CMP.Acnv1M{FT}(τ₀, q_thr, k₀));
        label = "q_thr = $(round(q_thr * 1e3, sigdigits=3)) g/kg ($(round(m, digits=2))x)",
    )
end
for (k, m) in zip(k_values, mult_k)
    MK.lines!(
        axs[3], q_icl_range * 1e3, acnv_rate(CMP.Acnv1M{FT}(τ₀, q_thr₀, k));
        label = "k = $(round(k, sigdigits=3)) ($(round(m, digits=2))x)",
    )
end
foreach(ax -> MK.axislegend(ax; position = :lt), axs)

MK.save("snow_autoconversion_sensitivity.svg", fig)
@info "Saved snow_autoconversion_sensitivity.svg"