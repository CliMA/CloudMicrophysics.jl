import CairoMakie as MK

import CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.ThermodynamicsInterface as TDI
import ClimaParams

FT = Float64

tps = TDI.TD.Parameters.ThermodynamicsParameters(FT)
mp = CMP.Microphysics1MParams(FT)

# example values
q_min, q_max = 1e-8, 5e-3
q_lcl_range = range(q_min, stop = q_max, length = 100)
ρ_air = 1.2
q_icl, q_tot = 5e-4, 20e-3
q_rai = 1e-3
q_sno = 1e-4
T = 273.15
limits = (0, q_max * 1e3, 0, nothing)

MK.set_theme!(MK.theme_minimal())

# default Kessler1M autoconversion parameters (populated from ClimaParams TOML)
acnv0 = mp.options.rain_autoconversion.acnv1M
τ₀, q_thr₀, k₀ = acnv0.τ, acnv0.q_threshold, FT(1)  # k₀ overridden from ClimaParams default (10) to 1
@info "Default rain autoconversion parameters" τ₀ q_thr₀ k₀

# rate curve for a given Acnv1M parameter struct.
# conv_q_lcl_to_q_rai dispatches on (and reads its parameters from) the option
# in the first argument, so we can pass a modified Kessler1M and reuse `mp`.
acnv_rate(acnv) = map(q_lcl_range) do q
    CM1.conv_q_lcl_to_q_rai(
        CMP.Kessler1M(acnv), mp, tps,
        (; q_tot, q_lcl = q, q_icl, q_rai, q_sno),
        (; ρ = ρ_air, T),
    )
end

# one-at-a-time parameter sweeps (keep the other two at their defaults)
# 5 points geometrically spaced from 1/5x to 5x the default, i.e. multipliers
# 0.2, ~0.447, 1, ~2.236, 5 (equal log-spacing, symmetric around the default)
mult = FT(5) .^ range(-1, 1, length = 5)
τ_values     = τ₀     .* mult   # autoconversion timescale [s]
q_thr_values = q_thr₀ .* mult   # threshold [kg/kg]
k_values     = k₀     .* mult   # smooth-transition steepness [-]

fig = MK.Figure(size = (1350, 430))
axs = [MK.Axis(fig[1, i]; xlabel = "q_lcl [g/kg]", limits) for i in 1:3]
axs[1].ylabel = "autoconversion rate [1/s]"
axs[1].title = "τ varied (q_thr = $(q_thr₀ * 1e3) g/kg, k = $(k₀))"
axs[2].title = "q_thr varied (τ = $(τ₀) s, k = $(k₀))"
axs[3].title = "k varied (τ = $(τ₀) s, q_thr = $(q_thr₀ * 1e3) g/kg)"

for (τ, m) in zip(τ_values, mult)
    MK.lines!(
        axs[1], q_lcl_range * 1e3, acnv_rate(CMP.Acnv1M{FT}(τ, q_thr₀, k₀));
        label = "τ = $(round(τ, sigdigits=3)) s ($(round(m, digits=2))x)",
    )
end
for (q_thr, m) in zip(q_thr_values, mult)
    MK.lines!(
        axs[2], q_lcl_range * 1e3, acnv_rate(CMP.Acnv1M{FT}(τ₀, q_thr, k₀));
        label = "q_thr = $(round(q_thr * 1e3, sigdigits=3)) g/kg ($(round(m, digits=2))x)",
    )
end
for (k, m) in zip(k_values, mult)
    MK.lines!(
        axs[3], q_lcl_range * 1e3, acnv_rate(CMP.Acnv1M{FT}(τ₀, q_thr₀, k));
        label = "k = $(round(k, sigdigits=3)) ($(round(m, digits=2))x)",
    )
end
foreach(ax -> MK.axislegend(ax; position = :lt), axs)

MK.save("autoconversion_sensitivity.svg", fig)