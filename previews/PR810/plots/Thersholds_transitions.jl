import CairoMakie as MK
import CloudMicrophysics
import ClimaParams

const CM1 = CloudMicrophysics.Microphysics1M
const CM2 = CloudMicrophysics.Microphysics2M
const CO = CloudMicrophysics.Common
const CP = ClimaParams
const CMP = CloudMicrophysics.Parameters

FT = Float64

rain = []
B1994 = []
TC1980 = []
LD2004 = []

k_thrshld_stpnss_values = [5.0, 2.0, 12.0]
for i in 1:3
    override_file = Dict(
        "threshold_smooth_transition_steepness" =>
            Dict("value" => k_thrshld_stpnss_values[i], "type" => "float"),
    )
    toml_dict = CP.create_toml_dict(FT; override_file)

    push!(rain, CMP.Rain(toml_dict))
    push!(B1994, CMP.B1994(toml_dict))
    push!(TC1980, CMP.TC1980(toml_dict))
    push!(LD2004, CMP.LD2004(toml_dict))
end

# example values
q_lcl_range = range(1e-8, stop = 1.5e-3, length = 1000)
N_d_range = range(1e7, stop = 1e9, length = 1000)
ρ_air = 1.0 # kg m^-3
N_d = 1e8

mp = CMP.Microphysics1MParams(CP.create_toml_dict(FT))
# Kessler threshold and timescale for quiescent air (w = 0)
acnv_q_threshold = CM1.rain_autoconversion_threshold(mp.processes.rain_autoconversion, mp, FT(0))
acnv_τ = CM1.rain_autoconversion_timescale(mp.processes.rain_autoconversion, mp, FT(0))
acnv_k = mp.process_params.rain_autoconversion.k

q_lcl_K1969 = [
    max(0, q_lcl - acnv_q_threshold) / acnv_τ
    for q_lcl in q_lcl_range
]
q_lcl_K1969_s = [
    CO.logistic_function_integral(q_lcl, acnv_q_threshold, acnv_k) /
    acnv_τ
    for q_lcl in q_lcl_range
]

q_lcl_TC1980 = [
    CM2.conv_q_lcl_to_q_rai(TC1980[2], q_lcl, ρ_air, N_d) for
    q_lcl in q_lcl_range
]
q_lcl_TC1980_s = [
    CM2.conv_q_lcl_to_q_rai(TC1980[2], q_lcl, ρ_air, N_d, true) for
    q_lcl in q_lcl_range
]
q_lcl_LD2004 = [
    CM2.conv_q_lcl_to_q_rai(LD2004[2], q_lcl, ρ_air, N_d) for
    q_lcl in q_lcl_range
]
q_lcl_LD2004_s = [
    CM2.conv_q_lcl_to_q_rai(LD2004[2], q_lcl, ρ_air, N_d, true) for
    q_lcl in q_lcl_range
]

N_d_B1994 =
    [CM2.conv_q_lcl_to_q_rai(B1994[3], 5e-4, ρ_air, N_d) for N_d in N_d_range]
N_d_B1994_s = [
    CM2.conv_q_lcl_to_q_rai(B1994[3], 5e-4, ρ_air, N_d, true) for
    N_d in N_d_range
]

# rates below the threshold are exactly zero, which a log axis cannot show
positive(x) = x > 0 ? x : NaN

fig = MK.Figure(size = (700, 450))
ax = MK.Axis(fig[1, 1]; xlabel = "q_lcl [g/kg]", ylabel = "autoconversion rate [1/s]")
MK.lines!(ax, q_lcl_range * 1e3, q_lcl_K1969; linewidth = 2, label = "K1969 without smoothing")
MK.lines!(ax, q_lcl_range * 1e3, q_lcl_K1969_s; linewidth = 2, label = "K1969 with smoothing")
MK.axislegend(ax; position = :lt)
MK.save("q_lcl_K1969.svg", fig)

fig = MK.Figure(size = (700, 450))
ax = MK.Axis(fig[1, 1]; xlabel = "q_lcl [g/kg]", ylabel = "autoconversion rate [1/s]", yscale = log10)
MK.lines!(ax, q_lcl_range * 1e3, positive.(q_lcl_TC1980); linewidth = 2, label = "TC1980 without smoothing")
MK.lines!(ax, q_lcl_range * 1e3, positive.(q_lcl_TC1980_s); linewidth = 2, label = "TC1980 with smoothing")
MK.lines!(ax, q_lcl_range * 1e3, positive.(q_lcl_LD2004); linewidth = 2, label = "LD2004 without smoothing")
MK.lines!(ax, q_lcl_range * 1e3, positive.(q_lcl_LD2004_s); linewidth = 2, label = "LD2004 with smoothing")
MK.ylims!(ax, 1e-10, 1e-5)
MK.axislegend(ax; position = :rb)
MK.save("q_lcl_TC1980_LD2004.svg", fig)

fig = MK.Figure(size = (700, 450))
ax = MK.Axis(
    fig[1, 1];
    xlabel = "N_d [1/cm3]",
    ylabel = "autoconversion rate [1/s]",
    xscale = log10,
    yscale = log10,
    xticks = [10, 100, 1000],
)
MK.lines!(ax, N_d_range * 1e-6, positive.(N_d_B1994); linewidth = 2, label = "B1994 without smoothing")
MK.lines!(ax, N_d_range * 1e-6, positive.(N_d_B1994_s); linewidth = 2, label = "B1994 with smoothing")
MK.ylims!(ax, 1e-13, 1e-5)
MK.axislegend(ax; position = :rt)
MK.save("N_d_B1994.svg", fig)
nothing
