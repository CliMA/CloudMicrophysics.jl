module Nucleation

"""
    nucleation_timestep(rh, temp, so4)

Calculates the nucleation rate for a single timestep under the given conditions.
All newly formed particles will be placed into the Aitken mode.

The initial nucleation rate comes from Vehkamaki et al, 2002 doi:10.1029/2002JD002184
Adjustment factor comes from Kerminen and Kulmala, 2002 doi:10.1029/2002JD002184

# Inputs
rh: Relative humidity, expressed as a percentage between 0 and 1
temp: Temperature (K)
so4: Total gas phase concentration of sulfuric acid (cm^-3)

# Examples
```julia-repl
julia> nucleation_timestep(0.55, 236, 1e7)
162.96173344575243
```
"""

function nucleation_timestep(rh, temp, so4)
    # Vehkamaki et al 2002 for initial parameterization of nucleation rate J
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

    # Nucleation rate 1/(cm^3 s)
    J = exp(
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
    )
    return J
end

end
