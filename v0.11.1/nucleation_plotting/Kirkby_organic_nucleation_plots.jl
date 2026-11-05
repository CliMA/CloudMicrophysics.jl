using Plots
import CLIMAParameters as CP

include("../../src/Nucleation.jl")
using .Nucleation

# Testing for CLOUD-experiment based nucleation rates.

FT = Float64
toml_dict = CP.create_toml_dict(FT)
param_names = ("a_1", "a_2", "a_3", "a_4", "a_5")
params = CP.get_parameter_values!(toml_dict, param_names)
params = (; params...)

HOM_concentrations = 10 .^ (6:0.125:8.7)

rates = map(HOM_concentrations) do HOM_conc
    sum(Nucleation.organic_nucleation_rate(0, HOM_conc, params)) * 1e-6
end

Plots.plot()
Plots.plot!(
    # title = title,
    HOM_concentrations,
    rates,
    xaxis = :log,
    yaxis = :log,
    lw = 3,
    ylims = (1e-4, 1e3),
    ylabel = "Nucleation rate (cm⁻³ s⁻¹)",
    xlabel = "HOM (cm⁻³)",
    label = "Kirkby parameterization",
)

Kirkby_points = [
    (2288812.4132151115, 0.0007608176801027033),
    (4933062.77556296, 0.005456422936842963),
    (7518401.669201947, 0.012808563896703628),
    (10899312.238623586, 0.04030351703529813),
    (11080975.881537117, 0.04150759064830022),
    (11740563.42748787, 0.047381882130014596),
    (17739769.94293252, 0.1491314951513984),
    (20588244.583238784, 0.17285245766668833),
    (21109938.829561356, 0.14929680441169513),
    (37322975.847007304, 0.5861159835346361),
    (46638708.493595704, 1.152374306087539),
    (108522641.6563595, 1.4866541939329214),
    (239753316.3650208, 11.816579222488736),
    (229881600.35871544, 20.029717576261014),
    (304770868.9495659, 14.532070496789618),
]


Plots.plot!(
    Kirkby_points,
    xticks = 10 .^ (5:10),
    seriestype = :scatter,
    label = "CLOUD Points",
)
Plots.svg("Kirkby_organic_nucleation");
