using GLMakie

import CloudMicrophysics
import CLIMAParameters

const CMT = CloudMicrophysics.CommonTypes
const CM1 = CloudMicrophysics.Microphysics1M
const CM2 = CloudMicrophysics.Microphysics2M
const CP = CLIMAParameters
const CMP = CloudMicrophysics.Parameters

include(joinpath(pkgdir(CloudMicrophysics), "test", "create_parameters.jl"))
FT = Float64
toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
const param_set = cloud_microphysics_parameters(toml_dict)
thermo_params = CMP.thermodynamics_params(param_set)

const liquid = CMT.LiquidType()
const rain = CMT.RainType()

include(joinpath(pkgdir(CloudMicrophysics), "docs", "src", "Wooddata.jl"))

# Example values
q_liq_range = range(1e-8, stop=1e-3, length=1000)
q_rai_range = range(1e-8, stop=1e-3, length=1000)
N_d_range = range(1e7, stop=1e9, length=1000)
q_liq = 5e-4
q_rai = 5e-4
ρ_air = 1.2 # kg m^-3

q_liq_KK2000 = [CM2.conv_q_liq_to_q_rai_KK2000(param_set, q_liq, 1.0, N_d = 1e8) for q_liq in q_liq_range]
q_liq_B1994 = [CM2.conv_q_liq_to_q_rai_B1994(param_set, q_liq, N_d = 1e8) for q_liq in q_liq_range]
q_liq_TC1980 = [CM2.conv_q_liq_to_q_rai_TC1980(param_set, q_liq, N_d = 1e8) for q_liq in q_liq_range]
q_liq_LD2004 = [CM2.conv_q_liq_to_q_rai_LD2004(param_set, q_liq, N_d = 1e8) for q_liq in q_liq_range]
q_liq_K1969 = [CM1.conv_q_liq_to_q_rai(param_set, q_liq) for q_liq in q_liq_range]

N_d_KK2000 = [CM2.conv_q_liq_to_q_rai_KK2000(param_set, 5e-4, 1.0, N_d = N_d) for N_d in N_d_range]
N_d_B1994 = [CM2.conv_q_liq_to_q_rai_B1994(param_set, 5e-4, N_d = N_d) for N_d in N_d_range]
N_d_TC1980 = [CM2.conv_q_liq_to_q_rai_TC1980(param_set, 5e-4, N_d = N_d) for N_d in N_d_range]
N_d_LD2004 = [CM2.conv_q_liq_to_q_rai_LD2004(param_set, 5e-4, N_d = N_d) for N_d in N_d_range]

accKK2000_q_liq = [CM2.accretion_KK2000(param_set, q_liq, q_rai, 1.0) for q_liq in q_liq_range]
accB1994_q_liq = [CM2.accretion_B1994(param_set, q_liq, q_rai) for q_liq in q_liq_range]
accTC1980_q_liq = [CM2.accretion_TC1980(param_set, q_liq, q_rai) for q_liq in q_liq_range]
accK1969_q_liq = [CM1.accretion(param_set, liquid, rain, q_liq, q_rai, ρ_air) for q_liq in q_liq_range]

accKK2000_q_rai = [CM2.accretion_KK2000(param_set, q_liq, q_rai, 1.0) for q_rai in q_rai_range]
accB1994_q_rai = [CM2.accretion_B1994(param_set, q_liq, q_rai) for q_rai in q_rai_range]
accTC1980_q_rai = [CM2.accretion_TC1980(param_set, q_liq, q_rai) for q_rai in q_rai_range]
accK1969_q_rai = [CM1.accretion(param_set, liquid, rain, q_liq, q_rai, ρ_air) for q_rai in q_rai_range]

fig = Figure()

ax1 = Axis(fig[1, 1]; yscale = log10)
ax2 = Axis(fig[1, 2]; xscale = log10, yscale = log10)
ax3 = Axis(fig[2, 1])
ax4 = Axis(fig[2, 2])

ylims!(ax1, [1e-13, 1e-5])
ylims!(ax2, [1e-13, 1e-5])
ylims!(ax3, [1e-13, 5e-6])
ylims!(ax4, [1e-13, 5e-6])

l1 = lines!(ax1, q_liq_range * 1e3, q_liq_KK2000, color = :red)
l2 = lines!(ax1, q_liq_range * 1e3, q_liq_B1994, color = :green)
l3 = lines!(ax1, q_liq_range * 1e3, q_liq_TC1980, color = :blue)
l4 = lines!(ax1, q_liq_range * 1e3, q_liq_LD2004, color = :purple)
l5 = lines!(ax1, q_liq_range * 1e3, q_liq_K1969, color = :black)
l6 = lines!(ax1, KK2000_x_q_liq, KK2000_y_q_liq, color = :red, linestyle = :dash)
l7 = lines!(ax1, B1994_x_q_liq, B1994_y_q_liq, color = :green, linestyle = :dash)
l8 = lines!(ax1, TC1980_x_q_liq, TC1980_y_q_liq, color = :blue, linestyle = :dash)
l9 = lines!(ax1, LD2004_x_q_liq, LD2004_y_q_liq, color = :purple, linestyle = :dash)

l10 = lines!(ax2, N_d_range * 1e-6, N_d_KK2000, color = :red)
l11 = lines!(ax2, N_d_range * 1e-6, N_d_B1994, color = :green)
l12 = lines!(ax2, N_d_range * 1e-6, N_d_TC1980, color = :blue)
l13 = lines!(ax2, N_d_range * 1e-6, N_d_LD2004, color = :purple)
l14 = lines!(ax2, KK2000_x_N_d, KK2000_y_N_d, color = :red, linestyle = :dash)
l15 = lines!(ax2, B1994_x_N_d, B1994_y_N_d, color = :green, linestyle = :dash)
l16 = lines!(ax2, TC1980_x_N_d, TC1980_y_N_d, color = :blue, linestyle = :dash)
l17 = lines!(ax2, LD2004_x_N_d, LD2004_y_N_d, color = :purple, linestyle = :dash)

l18 = lines!(ax3, q_liq_range * 1e3, accKK2000_q_liq, color = :red)
l19 = lines!(ax3, q_liq_range * 1e3, accB1994_q_liq, color = :green)
l20 = lines!(ax3, q_liq_range * 1e3, accTC1980_q_liq, color = :blue)
l21 = lines!(ax3, q_liq_range * 1e3, accK1969_q_liq, color = :black)

l22 = lines!(ax4, q_rai_range * 1e3, accKK2000_q_rai, color = :red)
l23 = lines!(ax4, q_rai_range * 1e3, accB1994_q_rai, color = :green)
l24 = lines!(ax4, q_rai_range * 1e3, accTC1980_q_rai, color = :blue)
l25 = lines!(ax4, q_rai_range * 1e3, accK1969_q_rai, color = :black)

ax1.xlabel = "q_liq [g/kg]"
ax1.ylabel = "autoconversion rate [kg m^-3 s^-1]"
ax2.xlabel = "N_d [cm^-3]"
ax2.ylabel = "autoconversion rate [kg m^-3 s^-1]"
ax3.xlabel = "q_liq [g/kg] (q_rai = 0.5e-4)"
ax3.ylabel = "accretion rate [?]"
ax4.xlabel = "q_rai [g/kg] (q_liq = 0.5e-4)"
ax4.ylabel = "accretion rate [?]"

Legend(
    fig[1, 3],
    [l1, l2, l3, l4, l5, l6, l7, l8, l9],
    ["KK2000", "B1994", "TC1980", "LD2004 (2.25 not 7.5)", "K1969", "Wood_KK2000", "Wood_B1994", "Wood_TC1980", "Wood_LD2004"]
)

#display(fig)
save("tmp.png", fig)
