using CairoMakie
CairoMakie.activate!(type = "svg")

import CLIMAParameters
import Thermodynamics as TD
import CloudMicrophysics
import CloudMicrophysics.Microphysics1M as CM1
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Parameters as CMP

FT = Float64

const tps = TD.Parameters.ThermodynamicsParameters(FT)
const aps = CMP.AirProperties(FT)

const KK2000 = CMP.KK2000(FT)
const B1994 = CMP.B1994(FT)
const TC1980 = CMP.TC1980(FT)
const LD2004 = CMP.LD2004(FT)
const VarTSc = CMP.VarTimescaleAcnv(FT)
const SB2006 = CMP.SB2006(FT)

const ce = CMP.CollisionEff(FT)

const liquid = CMP.CloudLiquid(FT)
const rain = CMP.Rain(FT)
const blk1mvel = CMP.Blk1MVelType(FT)

include(
    joinpath(pkgdir(CloudMicrophysics), "docs", "src", "plots", "Wooddata.jl"),
)

# Example values
q_liq_range = range(1e-8, stop = 1e-3, length = 1000)
q_rai_range = range(1e-8, stop = 1e-3, length = 1000)
N_d_range = range(1e7, stop = 1e9, length = 1000)
q_liq = 5e-4
q_rai = 5e-4
ρ_air = 1.0 # kg m^-3

q_liq_KK2000 = [
    CM2.conv_q_liq_to_q_rai(KK2000, q_liq, ρ_air, N_d = 1e8) for
    q_liq in q_liq_range
]
q_liq_B1994 = [
    CM2.conv_q_liq_to_q_rai(B1994, q_liq, ρ_air, N_d = 1e8) for
    q_liq in q_liq_range
]
q_liq_TC1980 = [
    CM2.conv_q_liq_to_q_rai(TC1980, q_liq, ρ_air, N_d = 1e8) for
    q_liq in q_liq_range
]
q_liq_LD2004 = [
    CM2.conv_q_liq_to_q_rai(LD2004, q_liq, ρ_air, N_d = 1e8) for
    q_liq in q_liq_range
]
q_liq_VarTimeScaleAcnv = [
    CM2.conv_q_liq_to_q_rai(VarTSc, q_liq, ρ_air, N_d = 1e8) for
    q_liq in q_liq_range
]
q_liq_SB2006 = [
    CM2.autoconversion(SB2006.acnv, q_liq, q_rai, ρ_air, 1e8).dq_rai_dt for
    q_liq in q_liq_range
]
q_liq_K1969 =
    [CM1.conv_q_liq_to_q_rai(rain.acnv1M, q_liq) for q_liq in q_liq_range]

N_d_KK2000 = [
    CM2.conv_q_liq_to_q_rai(KK2000, 5e-4, ρ_air, N_d = N_d) for N_d in N_d_range
]
N_d_B1994 = [
    CM2.conv_q_liq_to_q_rai(B1994, 5e-4, ρ_air, N_d = N_d) for N_d in N_d_range
]
N_d_TC1980 = [
    CM2.conv_q_liq_to_q_rai(TC1980, 5e-4, ρ_air, N_d = N_d) for N_d in N_d_range
]
N_d_LD2004 = [
    CM2.conv_q_liq_to_q_rai(LD2004, 5e-4, ρ_air, N_d = N_d) for N_d in N_d_range
]
N_d_VarTimeScaleAcnv = [
    CM2.conv_q_liq_to_q_rai(VarTSc, 5e-4, ρ_air, N_d = N_d) for N_d in N_d_range
]
N_d_SB2006 = [
    CM2.autoconversion(SB2006.acnv, q_liq, q_rai, ρ_air, N_d).dq_rai_dt for
    N_d in N_d_range
]

accKK2000_q_liq =
    [CM2.accretion(KK2000, q_liq, q_rai, ρ_air) for q_liq in q_liq_range]
accB1994_q_liq =
    [CM2.accretion(B1994, q_liq, q_rai, ρ_air) for q_liq in q_liq_range]
accTC1980_q_liq = [CM2.accretion(TC1980, q_liq, q_rai) for q_liq in q_liq_range]
accSB2006_q_liq = [
    CM2.accretion(SB2006, q_liq, q_rai, ρ_air, 1e8).dq_rai_dt for
    q_liq in q_liq_range
]
accK1969_q_liq = [
    CM1.accretion(liquid, rain, blk1mvel.rain, ce, q_liq, q_rai, ρ_air) for
    q_liq in q_liq_range
]

accKK2000_q_rai =
    [CM2.accretion(KK2000, q_liq, q_rai, ρ_air) for q_rai in q_rai_range]
accB1994_q_rai =
    [CM2.accretion(B1994, q_liq, q_rai, ρ_air) for q_rai in q_rai_range]
accTC1980_q_rai = [CM2.accretion(TC1980, q_liq, q_rai) for q_rai in q_rai_range]
accSB2006_q_rai = [
    CM2.accretion(SB2006, q_liq, q_rai, ρ_air, 1e8).dq_rai_dt for
    q_rai in q_rai_range
]
accK1969_q_rai = [
    CM1.accretion(liquid, rain, blk1mvel.rain, ce, q_liq, q_rai, ρ_air) for
    q_rai in q_rai_range
]

fig = Figure(resolution = (900, 600))

ax1 = Axis(fig[1, 1]; yscale = log10)
ax2 = Axis(fig[1, 2]; xscale = log10, yscale = log10)
ax3 = Axis(fig[2, 1]; yscale = log10)
ax4 = Axis(fig[2, 2]; yscale = log10)

ylims!(ax1, [1e-13, 1e-5])
ylims!(ax2, [1e-13, 1e-5])
ylims!(ax3, [1e-7, 5e-6])
ylims!(ax4, [1e-7, 5e-6])

l1 = lines!(ax1, q_liq_range * 1e3, q_liq_KK2000, color = :red)
l2 = lines!(ax1, q_liq_range * 1e3, q_liq_B1994, color = :green)
l3 = lines!(ax1, q_liq_range * 1e3, q_liq_TC1980, color = :blue)
l4 = lines!(ax1, q_liq_range * 1e3, q_liq_LD2004, color = :purple)
l5 = lines!(ax1, q_liq_range * 1e3, q_liq_K1969, color = :black)
l6 =
    lines!(ax1, KK2000_x_q_liq, KK2000_y_q_liq, color = :red, linestyle = :dash)
l7 =
    lines!(ax1, B1994_x_q_liq, B1994_y_q_liq, color = :green, linestyle = :dash)
l8 = lines!(
    ax1,
    TC1980_x_q_liq,
    TC1980_y_q_liq,
    color = :blue,
    linestyle = :dash,
)
l9 = lines!(
    ax1,
    LD2004_x_q_liq,
    LD2004_y_q_liq,
    color = :purple,
    linestyle = :dash,
)

l10 = lines!(ax2, N_d_range * 1e-6, N_d_KK2000, color = :red)
l11 = lines!(ax2, N_d_range * 1e-6, N_d_B1994, color = :green)
l12 = lines!(ax2, N_d_range * 1e-6, N_d_TC1980, color = :blue)
l13 = lines!(ax2, N_d_range * 1e-6, N_d_LD2004, color = :purple)
l14 = lines!(ax2, KK2000_x_N_d, KK2000_y_N_d, color = :red, linestyle = :dash)
l15 = lines!(ax2, B1994_x_N_d, B1994_y_N_d, color = :green, linestyle = :dash)
l16 = lines!(ax2, TC1980_x_N_d, TC1980_y_N_d, color = :blue, linestyle = :dash)
l17 =
    lines!(ax2, LD2004_x_N_d, LD2004_y_N_d, color = :purple, linestyle = :dash)

l18 = lines!(ax3, q_liq_range * 1e3, accKK2000_q_liq, color = :red)
l19 = lines!(ax3, q_liq_range * 1e3, accB1994_q_liq, color = :green)
l20 = lines!(ax3, q_liq_range * 1e3, accTC1980_q_liq, color = :blue)
l21 = lines!(ax3, q_liq_range * 1e3, accK1969_q_liq, color = :black)

l22 = lines!(ax4, q_rai_range * 1e3, accKK2000_q_rai, color = :red)
l23 = lines!(ax4, q_rai_range * 1e3, accB1994_q_rai, color = :green)
l24 = lines!(ax4, q_rai_range * 1e3, accTC1980_q_rai, color = :blue)
l25 = lines!(ax4, q_rai_range * 1e3, accK1969_q_rai, color = :black)

l26 = lines!(ax1, q_liq_range * 1e3, q_liq_SB2006, color = :cyan)
l27 = lines!(ax2, N_d_range * 1e-6, N_d_SB2006, color = :cyan)
l28 = lines!(ax3, q_liq_range * 1e3, accSB2006_q_liq, color = :cyan)
l28 = lines!(ax4, q_rai_range * 1e3, accSB2006_q_rai, color = :cyan)

l29 = lines!(ax1, q_liq_range * 1e3, q_liq_VarTimeScaleAcnv, color = :orange)
l30 = lines!(ax2, N_d_range * 1e-6, N_d_VarTimeScaleAcnv, color = :orange)

ax1.xlabel = "q_liq [g/kg]"
ax1.ylabel = "autoconversion rate [1/s]"
ax2.xlabel = "N_d [1/cm3]"
ax2.ylabel = "autoconversion rate [1/s]"
ax3.xlabel = "q_liq [g/kg] (q_rai = 0.5 g/kg)"
ax3.ylabel = "accretion rate [1/s]"
ax4.xlabel = "q_rai [g/kg] (q_liq = 0.5 g/kg)"
ax4.ylabel = "accretion rate [1/s]"

Legend(
    fig[1, 3],
    [l1, l2, l3, l4, l26, l5, l6, l7, l8, l9, l29],
    [
        "KK2000",
        "B1994",
        "TC1980",
        "LD2004",
        "SB2006",
        "K1969",
        "Wood_KK2000",
        "Wood_B1994",
        "Wood_TC1980",
        "Wood_LD2004",
        "SA2023",
    ],
)
save("Autoconversion_accretion.svg", fig)
