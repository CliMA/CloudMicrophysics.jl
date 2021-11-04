"""
    Aerosol activation scheme, which includes:

  - mean hygroscopicity for each mode of the aerosol size distribution
  - critical supersaturation for each mode of the aerosol size distribution
  - maximum supersaturation
  - total number of particles actived
  - total mass of particles actived
"""
module AerosolActivation

import SpecialFunctions

import Thermodynamics
import CloudMicrophysics
import CLIMAParameters

const SF = SpecialFunctions
const TD = Thermodynamics
const AM = CloudMicrophysics.AerosolModel
const CO = CloudMicrophysics.Common
const CP = CLIMAParameters
const CP_planet = CLIMAParameters.Planet
const APS = CP.AbstractParameterSet

export mean_hygroscopicity_parameter
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
    mean_hygroscopicity_parameter(param_set, ad)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct

Returns a tuple of hygroscopicity parameters
(one tuple element for each aerosol size distribution mode).
The tuple is computed either as mass-weighted B parameters
(Abdul-Razzak and Ghan 2000)
or volume weighted kappa parameters (Petters and Kreidenweis 2007).
Implemented via a dispatch based on aerosol distribution mode type.
"""
function mean_hygroscopicity_parameter(
    param_set::APS,
    ad::AM.AerosolDistribution{NTuple{N, T}},
) where {N, T <: AM.Mode_B}

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
function mean_hygroscopicity_parameter(
    param_set::APS,
    ad::AM.AerosolDistribution{NTuple{N, T}},
) where {N, T <: AM.Mode_κ}

    return ntuple(length(ad.Modes)) do i

        modei = ad.Modes[i]
        sum(1:(modei.n_components)) do j
            modei.vol_mix_ratio[j] * modei.kappa[j]
        end
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
    ad::AM.AbstractAerosolDistribution,
    T::FT,
) where {FT <: Real}

    A::FT = coeff_of_curvature(param_set, T)
    hygro = mean_hygroscopicity_parameter(param_set, ad)

    return ntuple(length(ad.Modes)) do i
        2 / sqrt(hygro[i]) * (A / 3 / ad.Modes[i].r_dry)^(3 / 2)
    end
end

"""
    max_supersaturation(param_set, ad, T, p, w, q)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the maximum supersaturation.
"""
function max_supersaturation(
    param_set::APS,
    ad::AM.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    _grav::FT = CP_planet.grav(param_set)
    _ρ_cloud_liq::FT = CP_planet.ρ_cloud_liq(param_set)

    ϵ::FT = 1 / CP_planet.molmass_ratio(param_set)
    R_m::FT = TD.gas_constant_air(param_set, q)
    cp_m::FT = TD.cp_m(param_set, q)

    L::FT = CP_planet.LH_v0(param_set)

    p_vs::FT = TD.saturation_vapor_pressure(param_set, T, TD.Liquid())
    G::FT = CO.G_func(param_set, T, TD.Liquid()) / _ρ_cloud_liq

    # eq 11, 12 in Razzak et al 1998
    # but following eq 10 from Rogers 1975
    α::FT = L * _grav * ϵ / R_m / cp_m / T^2 - _grav / R_m / T
    γ::FT = R_m * T / ϵ / p_vs + ϵ * L^2 / cp_m / T / p

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
    N_activated_per_mode(param_set, ad, T, p, w, q)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the number of activated aerosol particles
in each aerosol size distribution mode.
"""
function N_activated_per_mode(
    param_set::APS,
    ad::AM.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    smax::FT = max_supersaturation(param_set, ad, T, p, w, q)
    sm = critical_supersaturation(param_set, ad, T)

    return ntuple(length(ad.Modes)) do i

        mode_i = ad.Modes[i]
        u_i::FT = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        mode_i.N * (1 / 2) * (1 - SF.erf(u_i))
    end
end

"""
    M_activated_per_mode(param_set, ad, T, p, w, q)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the mass of activated aerosol particles
per mode of the aerosol size distribution.
"""
function M_activated_per_mode(
    param_set::APS,
    ad::AM.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    smax = max_supersaturation(param_set, ad, T, p, w, q)
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
    total_N_activated(param_set, ad, T, p, w, q)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the total number of activated aerosol particles.
"""
function total_N_activated(
    param_set::APS,
    ad::AM.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    return sum(N_activated_per_mode(param_set, ad, T, p, w, q))

end

"""
    total_M_activated(param_set, ad, T, p, w, q)

  - `param_set` - abstract set with Earth's parameters
  - `ad` - aerosol distribution struct
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the total mass of activated aerosol particles.
"""
function total_M_activated(
    param_set::APS,
    ad::AM.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    return sum(M_activated_per_mode(param_set, ad, T, p, w, q))

end

end # module AerosolActivation.jl
