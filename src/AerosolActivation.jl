"""
    Aerosol activation scheme, which includes:

  - mean hygroscopicity for each mode of the aerosol size distribution
  - critical supersaturation for each mode of the aerosol size distribution
  - maximum supersaturation
  - total number of particles actived
  - total mass of particles actived
"""
module AerosolActivation

import SpecialFunctions as SF

import Thermodynamics as TD

import CloudMicrophysics.AerosolModel as AM
import CloudMicrophysics.Common as CO

import CLIMAParameters as CP
import CLIMAParameters.Planet as CP_planet

const APS = CP.AbstractParameterSet

export mean_hygroscopicity
export max_supersaturation
export N_activated_per_mode
export M_activated_per_mode
export total_N_activated
export total_M_activated

"""
    coeff_of_curvature(param_set, T)

  - `param_set` - abstract set with Earth's parameters
  - `T` - air temperature

Returns a curvature coefficient.
"""
function coeff_of_curvature(param_set::APS, T::FT) where {FT <: Real}

    _molmass_water::FT = CP_planet.molmass_water(param_set)
    _gas_constant::FT = CP.gas_constant()
    _ρ_cloud_liq::FT = CP_planet.ρ_cloud_liq(param_set)
    _surface_tension::FT = CP_planet.surface_tension_coeff(param_set)

    return 2 * _surface_tension * _molmass_water / _ρ_cloud_liq /
           _gas_constant / T
end

"""
    mean_hygroscopicity(param_set, ad)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct

Returns a tuple of mean hygroscopicities
(one tuple element for each aerosol size distribution mode).
"""
function mean_hygroscopicity(param_set::APS, ad::AM.AerosolDistribution)

    _molmass_water = CP_planet.molmass_water(param_set)
    _ρ_cloud_liq = CP_planet.ρ_cloud_liq(param_set)

    return ntuple(length(ad.Modes)) do i

        mode_i = ad.Modes[i]

        nom = sum(1:(mode_i.n_components)) do j
            mode_i.mass_mix_ratio[j] *
            mode_i.dissoc[j] *
            mode_i.osmotic_coeff[j] *
            mode_i.soluble_mass_frac[j] / mode_i.molar_mass[j]
        end

        den = sum(1:(mode_i.n_components)) do j
            mode_i.mass_mix_ratio[j] / mode_i.aerosol_density[j]
        end

        nom / den * _molmass_water / _ρ_cloud_liq
    end
end

"""
    critical_supersaturation(param_set, ad, T)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature

Returns a tuple of critical supersaturations
(one tuple element for each aerosol size distribution mode).
"""
function critical_supersaturation(
    param_set::APS,
    ad::AM.AerosolDistribution,
    T::FT,
) where {FT <: Real}

    A::FT = coeff_of_curvature(param_set, T)
    B = mean_hygroscopicity(param_set, ad)

    return ntuple(length(ad.Modes)) do i
        2 / sqrt(B[i]) * (A / 3 / ad.Modes[i].r_dry)^(3 / 2)
    end
end

"""
    max_supersaturation(param_set, ad, T, p, w)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity

Returns the maximum supersaturation.
"""
function max_supersaturation(
    param_set::APS,
    ad::AM.AerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    _grav::FT = CP_planet.grav(param_set)
    _molmass_water::FT = CP_planet.molmass_water(param_set)
    _molmass_dryair::FT = CP_planet.molmass_dryair(param_set)
    _gas_constant::FT = CP.gas_constant()
    _cp_d::FT = CP_planet.cp_d(param_set)
    _ρ_cloud_liq::FT = CP_planet.ρ_cloud_liq(param_set)

    L::FT = TD.latent_heat_vapor(param_set, T)
    p_vs::FT = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
    G::FT = CO.G_func(param_set, T, TD.Liquid()) / _ρ_cloud_liq

    # eq 11, 12 in Razzak et al 1998
    α::FT =
        _grav * _molmass_water * L / _cp_d / _gas_constant / T^2 -
        _grav * _molmass_dryair / _gas_constant / T
    γ::FT =
        _gas_constant * T / p_vs / _molmass_water +
        _molmass_water * L^2 / _cp_d / p / _molmass_dryair / T

    A::FT = coeff_of_curvature(param_set, T)
    ζ::FT = 2 * A / 3 * sqrt(α * w / G)

    Sm = critical_supersaturation(param_set, ad, T)

    tmp::FT = sum(1:length(ad.Modes)) do i

        mode_i = ad.Modes[i]

        f::FT = 0.5 * exp(2.5 * (log(mode_i.stdev))^2)
        g::FT = 1 + 0.25 * log(mode_i.stdev)
        η::FT = (α * w / G)^(3 / 2) / (2 * pi * _ρ_cloud_liq * γ * mode_i.N)

        1 / (Sm[i])^2 *
        (f * (ζ / η)^(3 / 2) + g * (Sm[i]^2 / (η + 3 * ζ))^(3 / 4))
    end

    return FT(1) / sqrt(tmp)
end

"""
    N_activated_per_mode(param_set, ad, T, p, w)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity

Returns the number of activated aerosol particles
in each aerosol size distribution mode.
"""
function N_activated_per_mode(
    param_set::APS,
    ad::AM.AerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    smax::FT = max_supersaturation(param_set, ad, T, p, w)
    sm = critical_supersaturation(param_set, ad, T)

    return ntuple(length(ad.Modes)) do i

        mode_i = ad.Modes[i]
        u_i::FT = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        mode_i.N * (1 / 2) * (1 - SF.erf(u_i))
    end
end

"""
    M_activated_per_mode(param_set, ad, T, p, w)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity

Returns the mass of activated aerosol particles
per mode of the aerosol size distribution.
"""
function M_activated_per_mode(
    param_set::APS,
    ad::AM.AerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    smax = max_supersaturation(param_set, ad, T, p, w)
    sm = critical_supersaturation(param_set, ad, T)

    return ntuple(length(ad.Modes)) do i

        mode_i = ad.Modes[i]

        avg_molar_mass_i = sum(1:(mode_i.n_components)) do j
            mode_i.molar_mass[j] * mode_i.mass_mix_ratio[j]
        end

        u_i = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        avg_molar_mass_i * 1 / 2 *
        (1 - SF.erf(u_i - 3 * sqrt(2) / 2 * log(mode_i.stdev)))
    end
end

"""
    total_N_activated(param_set, ad, T, p, w)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity

Returns the total number of activated aerosol particles.
"""
function total_N_activated(
    param_set::APS,
    ad::AM.AerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    return sum(N_activated_per_mode(param_set, ad, T, p, w))

end

"""
    total_M_activated(param_set, ad, T, p, w)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity

Returns the total mass of activated aerosol particles.
"""
function total_M_activated(
    param_set::APS,
    ad::AM.AerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    return sum(M_activated_per_mode(param_set, ad, T, p, w))

end

end # module AerosolActivation.jl
