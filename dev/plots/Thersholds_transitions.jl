import Plots
import CloudMicrophysics
import ClimaParams

const PL = Plots
const CM1 = CloudMicrophysics.Microphysics1M
const CM2 = CloudMicrophysics.Microphysics2M
const CP = ClimaParams
const CMP = CloudMicrophysics.Parameters

FT = Float64

rain = []
B1994 = []
TC1980 = []
LD2004 = []

k_thrshld_stpnss_values = [5.0, 2.0, 12.0]
for i in 1:3
    override_file = joinpath("override_dict.toml")
    open(override_file, "w") do io
        println(io, "[threshold_smooth_transition_steepness]")
        println(io, "alias = \"k_thrshld_stpnss\"")
        println(io, "value = " * string(k_thrshld_stpnss_values[i]))
        println(io, "type = \"float\"")
    end
    toml_dict = CP.create_toml_dict(FT; override_file)
    isfile(override_file) && rm(override_file; force = true)

    push!(rain, CMP.Rain(toml_dict))
    push!(B1994, CMP.B1994(toml_dict))
    push!(TC1980, CMP.TC1980(toml_dict))
    push!(LD2004, CMP.LD2004(toml_dict))
end

# example values
q_liq_range = range(1e-8, stop = 1.5e-3, length = 1000)
N_d_range = range(1e7, stop = 1e9, length = 1000)
ρ_air = 1.0 # kg m^-3
N_d = 1e8

q_liq_K1969 =
    [CM1.conv_q_liq_to_q_rai(rain[1].acnv1M, q_liq) for q_liq in q_liq_range]
q_liq_K1969_s = [
    CM1.conv_q_liq_to_q_rai(rain[1].acnv1M, q_liq, true) for
    q_liq in q_liq_range
]

q_liq_TC1980 = [
    CM2.conv_q_liq_to_q_rai(TC1980[2], q_liq, ρ_air, N_d) for
    q_liq in q_liq_range
]
q_liq_TC1980_s = [
    CM2.conv_q_liq_to_q_rai(TC1980[2], q_liq, ρ_air, N_d, true) for
    q_liq in q_liq_range
]
q_liq_LD2004 = [
    CM2.conv_q_liq_to_q_rai(LD2004[2], q_liq, ρ_air, N_d) for
    q_liq in q_liq_range
]
q_liq_LD2004_s = [
    CM2.conv_q_liq_to_q_rai(LD2004[2], q_liq, ρ_air, N_d, true) for
    q_liq in q_liq_range
]

N_d_B1994 =
    [CM2.conv_q_liq_to_q_rai(B1994[3], 5e-4, ρ_air, N_d) for N_d in N_d_range]
N_d_B1994_s = [
    CM2.conv_q_liq_to_q_rai(B1994[3], 5e-4, ρ_air, N_d, true) for
    N_d in N_d_range
]

PL.plot(
    q_liq_range * 1e3,
    q_liq_K1969,
    linewidth = 2,
    xlabel = "q_liq [g/kg]",
    ylabel = "autoconversion rate [1/s]",
    label = "K1969 without smoothing",
)
PL.plot!(
    q_liq_range * 1e3,
    q_liq_K1969_s,
    linewidth = 2,
    label = "K1969 with smoothing",
)
PL.savefig("q_liq_K1969.svg") # hide

PL.plot(
    q_liq_range * 1e3,
    q_liq_TC1980,
    linewidth = 2,
    xlabel = "q_liq [g/kg]",
    ylabel = "autoconversion rate [1/s]",
    label = "TC1980 without smoothing",
    yaxis = :log,
    ylim = (1e-10, 1e-5),
)
PL.plot!(
    q_liq_range * 1e3,
    q_liq_TC1980_s,
    linewidth = 2,
    label = "TC1980 with smoothing",
)
PL.plot!(
    q_liq_range * 1e3,
    q_liq_LD2004,
    linewidth = 2,
    label = "LD2004 without smoothing",
)
PL.plot!(
    q_liq_range * 1e3,
    q_liq_LD2004_s,
    linewidth = 2,
    label = "LD2004 with smoothing",
)
PL.savefig("q_liq_TC1980_LD2004.svg") # hide

PL.plot(
    N_d_range * 1e-6,
    N_d_B1994,
    linewidth = 2,
    xlabel = "N_d [1/cm3]",
    ylabel = "autoconversion rate [1/s]",
    label = "B1994 without smoothing",
    xaxis = :log,
    yaxis = :log,
    ylim = (1e-13, 1e-5),
)
PL.plot!(
    N_d_range * 1e-6,
    N_d_B1994_s,
    linewidth = 2,
    label = "B1994 with smoothing",
)
PL.savefig("N_d_B1994.svg") # hide
