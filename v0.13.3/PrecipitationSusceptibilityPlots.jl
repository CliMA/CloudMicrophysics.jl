using CairoMakie
CairoMakie.activate!(type = "svg")

import CloudMicrophysics as CM
import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

const FT = Float64

toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
prs = cloud_microphysics_parameters(toml_dict)

scheme = CMT.SB2006Type()

q_tot = FT(0.5e-3)

# Start from small nonzero value to more clearly see the asymptotes
q_rai = range(0.0001 * q_tot, q_tot, 1000)
q_liq = q_tot .- q_rai

N_liq = FT(1e8)
ρ = FT(1)

τ = q_rai ./ q_tot

aut_rates =
    CMPS.precipitation_susceptibility_autoconversion.(
        Ref(prs),
        Ref(scheme),
        q_liq,
        q_rai,
        Ref(ρ),
        Ref(N_liq),
    )

acc_rates =
    CMPS.precipitation_susceptibility_accretion.(
        Ref(prs),
        Ref(scheme),
        q_liq,
        q_rai,
        Ref(ρ),
        Ref(N_liq),
    )

fig = Figure()

ax = Axis(fig[1, 1])

ax.xlabel = "q_rai / (q_liq + q_rai)"
ax.ylabel = "Precipitation susceptibility"

l1 = lines!(ax, τ, [r.d_ln_pp_d_ln_q_liq for r in aut_rates], color = :red)
l2 = lines!(ax, τ, [r.d_ln_pp_d_ln_q_rai for r in aut_rates], color = :brown)
l3 = lines!(ax, τ, [r.d_ln_pp_d_ln_q_liq for r in acc_rates], color = :blue)
l4 = lines!(ax, τ, [r.d_ln_pp_d_ln_q_rai for r in acc_rates], color = :green)

Legend(
    fig[1, 2],
    [l1, l2, l3, l4],
    [
        "autoconversion, q_liq",
        "autoconversion, q_rai",
        "accretion, q_liq",
        "accretion, q_rai",
    ],
)

save("Glassmeier-Lohmann_Fig2.svg", fig)
