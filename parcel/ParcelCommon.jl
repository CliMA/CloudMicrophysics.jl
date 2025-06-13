import Interpolations as Intp

# Saturation ratio over ice
ξ(tps, T) =
    TDI.saturation_vapor_pressure_over_liquid(tps, T) /
    TDI.saturation_vapor_pressure_over_ice(tps, T)

# Vapour partial pressure
eᵥ(qᵥ, p_air, R_air, Rᵥ) = qᵥ * p_air * Rᵥ / R_air

# Saturation pressure over ice given saturation pressure over liquid
S_i(tps, T, S_liq) = ξ(tps, T) * S_liq

# Interpolating for cooling/expansion rate
function AIDA_rate(model_t, data_t, data)
    data_at_t = Intp.linear_interpolation(data_t, data)

    if model_t < data_t[end]
        return (data_at_t(model_t + 1) - data_at_t(model_t))
    else
        return 0
    end
end
