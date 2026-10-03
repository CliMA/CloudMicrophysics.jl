import CairoMakie as MK

import ClimaParams as CP
import CloudMicrophysics as CM
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.Nucleation as Nucleation

FT = Float64
params = CMP.H2S04NucleationParameters(FT)

"""
    vehkamaki_nucleation_timestep(rh, temp, so4)

 - `rh` - Relative humidity, expressed as a percentage between 0 and 1
 - `temp` - Temperature (K)
 - `so4` - Total gas phase concentration of sulfuric acid (m⁻³)

Calculates the rate of binary H2SO4-H2O nucleation for a single timestep under the given conditions.
The parameterization is valid for the temperature range 230.15–300.15 K,
relative humidities 0.01–100% and total sulfuric acid concentrations 1e10-1e17 m⁻³.

All newly formed particles will be placed into the Aitken mode.
The initial nucleation rate comes from Vehkamaki et al, 2002 doi:10.1029/2002JD002184
Adjustment factor comes from Kerminen and Kulmala, 2002 doi:10.1029/2002JD002184
"""
function vehkamaki_nucleation_timestep(rh, temp, so4)
    so4 = so4 / 1e6 # convert to cm⁻³
    ln_so4 = log(so4)
    log_rh = log(rh)
    x_crit = (
        0.740997 - 0.00266379 * temp - 0.00349998 * ln_so4 +
        0.0000504022 * temp * ln_so4 +
        0.00201048 * log_rh - 0.000183289 * temp * log_rh +
        0.00157407 * log_rh^2 - 0.0000179059 * temp * log_rh^2 +
        0.000184403 * log_rh^3 - 1.50345e-6 * temp * log_rh^3
    )

    a_T = (
        0.14309 + 2.21956 * temp - 0.0273911 * temp^2 +
        0.0000722811 * temp^3 +
        5.91822 / x_crit
    )
    b_T = (
        0.117489 + 0.462532 * temp - 0.0118059 * temp^2 +
        0.0000404196 * temp^3 +
        15.7963 / x_crit
    )
    c_T = (
        -0.215554 - 0.0810269 * temp + 0.00143581 * temp^2 -
        4.7758e-6 * temp^3 - 2.91297 / x_crit
    )
    d_T = (
        -3.58856 + 0.049508 * temp - 0.00021382 * temp^2 + 3.10801e-7 * temp^3 - 0.0293333 / x_crit
    )
    e_T = (
        1.14598 - 0.600796 * temp + 0.00864245 * temp^2 -
        0.0000228947 * temp^3 - 8.44985 / x_crit
    )
    f_T = (
        2.15855 + 0.0808121 * temp - 0.000407382 * temp^2 -
        4.01957e-7 * temp^3 + 0.721326 / x_crit
    )
    g_T = (
        1.6241 - 0.0160106 * temp +
        0.0000377124 * temp^2 +
        3.21794e-8 * temp^3 - 0.0113255 / x_crit
    )
    h_T = (
        9.71682 - 0.115048 * temp +
        0.000157098 * temp^2 +
        4.00914e-7 * temp^3 +
        0.71186 / x_crit
    )
    i_T = (
        -1.05611 + 0.00903378 * temp - 0.0000198417 * temp^2 +
        2.46048e-8 * temp^3 - 0.0579087 / x_crit
    )
    j_T = (
        -0.148712 + 0.00283508 * temp - 9.24619e-6 * temp^2 +
        5.00427e-9 * temp^3 - 0.0127081 / x_crit
    )

    # Nucleation rate 1/(m^3 s)
    J =
        exp(
            a_T +
            b_T * log_rh +
            c_T * log_rh^2 +
            d_T * log_rh^3 +
            e_T * ln_so4 +
            f_T * log_rh * ln_so4 +
            g_T * log_rh^2 * ln_so4 +
            h_T * ln_so4^2 +
            i_T * log_rh * ln_so4^2 +
            j_T * ln_so4^3,
        ) * 1e6
    return J
end


function plot_vehk_cloud_comparison(humidities, temp, h2so4_concentrations; ticks...)
    cloud_rates = map(h2so4_concentrations) do h2so4_concentration
        sum(
            Nucleation.h2so4_nucleation_rate(
                h2so4_concentration * 1e6,
                nh3_concentration,
                negative_ion_concentration,
                temp,
                params,
            ),
        ) * 1e-6
    end
    fig = MK.Figure(size = (700, 500))
    ax = MK.Axis(
        fig[1, 1];
        title = "$temp K",
        xlabel = "[H2SO4] (cm⁻³)",
        ylabel = "Nucleation rate (cm⁻³ s⁻¹)",
        xscale = log10,
        yscale = log10,
        ticks...,
    )
    MK.lines!(ax, h2so4_concentrations, cloud_rates; linewidth = 3, label = "CLOUD")
    for rh in humidities
        vehk_rates = map(h2so4_concentrations) do x
            vehkamaki_nucleation_timestep(rh, temp, x * 1e6) / 1e6
        end
        MK.lines!(
            ax,
            h2so4_concentrations,
            vehk_rates;
            linewidth = 3,
            label = "Vehkamaki $(round(rh*100;digits=3))% RH",
        )
    end
    MK.axislegend(ax; position = :lt)
    MK.save("CLOUD_Vehk_comparison_$temp.svg", fig)
end

humidities = (0.55, 0.382)
nh3_concentration = 0
negative_ion_concentration = 0
h2so4_concentrations = 10 .^ range(log10(8e9), log10(3.7e10), length = 200)
# at 298 K the automatic log ticks fall on fractional powers of ten
plot_vehk_cloud_comparison(
    humidities,
    298,
    h2so4_concentrations;
    xticks = ([1e10, 2e10, 3e10], [MK.rich("$(k)×10", MK.superscript("10")) for k in 1:3]),
    yticks = (10.0 .^ (0:2:8), [MK.rich("10", MK.superscript(string(k))) for k in 0:2:8]),
)
humidities = (0.523, 0.382)
h2so4_concentrations = 10 .^ range(6, 9, length = 200)
plot_vehk_cloud_comparison(humidities, 236, h2so4_concentrations)
nothing
