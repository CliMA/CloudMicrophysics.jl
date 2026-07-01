using Plots

import CLIMAParameters as CP
import CloudMicrophysics.Nucleation as Nucleation
import CloudMicrophysics.Parameters as CMP

FT = Float64
params = CMP.MixedNucleationParameters(FT)

# Units: 1/m³
bioOxOrg_concentrations = 10 .^ (5.8:0.125:8.5)

nucleation_rates = map(bioOxOrg_concentrations) do bioOxOrg_conc
    Nucleation.organic_and_h2so4_nucleation_rate_bioOxOrg_prescribed(
        2.6e6,
        bioOxOrg_conc,
        params,
    ) * 1e6
end

Plots.plot(xaxis = :log, yaxis = :log, lw = 3)
Plots.plot!(
    bioOxOrg_concentrations,
    nucleation_rates,
    label = "Riccobono parameterization",
)
CLOUD_points = [
    (5360344.038850841, 0.09948622490099472),
    (8968652.585232794, 0.3074108277542233),
    (31220049.06625331, 0.41892723198582127),
    (189317581.7233987, 1.054767073233024),
    (1194554.6269400613, 0.02575861817013657),
    (1913819.0312476552, 0.01580159786357796),
    (2995132.266946077, 0.017881741897758208),
    (969552.1151059747, 0.01014906869624897),
    (855010.5080933397, 0.021947616507750668),
    (1084499.6602613404, 0.03132689161710575),
    (1241807.5569017627, 0.05371756362194626),
    (763822.0791489122, 0.012957872278423216),
    (1135505.5592274622, 0.019432271996165488),
    (903852.3947998419, 0.05636713954136084),
    (1795038.1716049775, 0.06228330897095073),
    (2410058.720090212, 0.07866396091996666),
    (3388495.4821890164, 0.05452490332715843),
    (3506842.011886563, 0.04268266991569145),
    (9761038.49207006, 0.4279903441430223),
    (9646144.458902434, 0.6411390209745732),
    (34160788.5981378, 0.7543305822696003),
    (33011936.663243562, 0.8736791862602847),
    (27894160.336738832, 0.2215198967777188),
    (24209707.65594941, 0.20076866070325255),
    (209532141.16201726, 1.7647289054732143),
    (196732207.6712905, 2.987642282932369),
]

Plots.plot!(
    CLOUD_points,
    seriestype = :scatter,
    label = "CLOUD data",
    ylabel = "Nucleation rate (cm⁻³ s⁻¹)",
    xlabel = "[BioOxOrg] (cm⁻³)",
    xlims = [10^5.4, 10^8.4],
    ylims = [0.001, 10],
    xticks = [1e6, 1e7, 1e8],
    yticks = [0.01, 0.1, 1, 10],
)
Plots.svg("Riccobono_nucleation");
