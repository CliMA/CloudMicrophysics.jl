import Plots as PL
import ClimaParams as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.MicrophysicsNonEq as CMNe

FT = Float32

# set constants
ip = CMP.Frostenberg2023(FT)
aps = CMP.AirProperties(FT)
q_icl = FT(1e-6)

A = FT(1e-5)
override_file = Dict(
    "sublimation_deposition_timescale" =>
        Dict("value" => A, "type" => "float"),
)
override_toml_dict = CP.create_toml_dict(FT; override_file)
ice = CMP.CloudIce(override_toml_dict)

T_range = range(233, stop = 271, length = 500)

τᵢ_values = [CMNe.τ_Frostenberg(ice, aps, ip, q_icl, T) for T in T_range]

PL.plot(
    T_range,
    τᵢ_values,
    xlabel = "T [K]",
    ylabel = "τᵢ [s]",
    label = "τᵢ(T)",
    legend = :topright,
    yaxis = :log,
)

PL.savefig("τᵢ_Frostenberg.svg")
