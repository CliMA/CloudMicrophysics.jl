# Saturation ratio over ice
ξ(tps, T) =
    TD.saturation_vapor_pressure(tps, T, TD.Liquid()) /
    TD.saturation_vapor_pressure(tps, T, TD.Ice())

# Vapour partial pressure
eᵥ(qᵥ, p_air, R_air, Rᵥ) = qᵥ * p_air * Rᵥ / R_air

# Saturation pressure over ice given saturation pressure over liquid
S_i(tps, T, S_liq) = ξ(tps, T) * S_liq
