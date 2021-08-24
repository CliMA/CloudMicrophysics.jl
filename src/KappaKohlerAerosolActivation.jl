"""
    Kappa-Kohler based Aerosol activation scheme, which includes:

  - critical supersaturation for each mode of the aerosol size distribution
  - maximum supersaturation
  - total number of particles actived
  - total mass of particles actived
"""
# module KappaKohlerAerosolActivation

using SpecialFunctions

using Thermodynamics

include("/home/idularaz/CloudMicrophysics.jl/src/KappaKohlerAerosolModel.jl")
using CloudMicrophysics.Microphysics_1M: G_func

using CLIMAParameters
using CLIMAParameters: gas_constant
using CLIMAParameters.Planet:
    ρ_cloud_liq,
    R_v,
    grav,
    molmass_water,
    molmass_dryair,
    cp_d,
    surface_tension_coeff
using CLIMAParameters.Atmos.Microphysics: K_therm, D_vapor

const APS = AbstractParameterSet

#export KK_max_supersaturation
# KK_N_activated_per_mode
#export KK_M_activated_per_mode
#export KK_total_N_activated
#export KK_total_M_activated

"""
    coeff_of_curvature(param_set, T)

  - `param_set` - abstract set with Earth's parameters
  - `T` - air temperature

Returns a curvature coefficient.
"""
function KK_coeff_of_curvature(param_set::APS, T::FT) where {FT <: Real}

    _molmass_water::FT = molmass_water(param_set)
    _gas_constant::FT = gas_constant()
    _ρ_cloud_liq::FT = ρ_cloud_liq(param_set)
    _surface_tension::FT = surface_tension_coeff(param_set)

    return 2 * _surface_tension * _molmass_water / _ρ_cloud_liq /
           _gas_constant / T
end

"""
    kappa(param_set, ad)
 
 - 'param_set' - abstract set with Earth's parameters
 - 'ad' - aerosol distribution

Returns volume-weighted kappa parameter. 
"""
function kappa(param_set::APS, ad::KappaKohlerAerosolDistribution)
    return ntuple(length(ad.KK_Modes)) do i
        modei = ad.KK_Modes[i]
        sum(1:modei.n_components) do i
            modei.vol_mix_ratio[i] * modei.kappa[i]
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
function KK_critical_supersaturation(param_set::APS, ad::KappaKohlerAerosolDistribution, T::FT) where {FT <: Real}
    kappa_avg = kappa(param_set, ad)
    A = KK_coeff_of_curvature(param_set, T)
    return ntuple(length(ad.KK_Modes)) do i
        modei = ad.KK_Modes[i]
        2 / sqrt(kappa_avg[i]) * (A / 3 / ad.KK_Modes[i].r_dry)^(3 / 2)
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
function KK_max_supersaturation(
    param_set::APS,
    ad::KappaKohlerAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    _grav::FT = grav(param_set)
    _molmass_water::FT = molmass_water(param_set)
    _molmass_dryair::FT = molmass_dryair(param_set)
    _gas_constant::FT = gas_constant()
    _cp_d::FT = cp_d(param_set)
    _ρ_cloud_liq::FT = ρ_cloud_liq(param_set)

    L::FT = latent_heat_vapor(param_set, T)
    p_vs::FT = saturation_vapor_pressure(param_set, T, Liquid())
    G::FT = G_func(param_set, T, Liquid()) / _ρ_cloud_liq

    # eq 11, 12 in Razzak et al 1998
    α::FT =
        _grav * _molmass_water * L / _cp_d / _gas_constant / T^2 -
        _grav * _molmass_dryair / _gas_constant / T
    γ::FT =
        _gas_constant * T / p_vs / _molmass_water +
        _molmass_water * L^2 / _cp_d / p / _molmass_dryair / T

    A::FT = KK_coeff_of_curvature(param_set, T)
    ζ::FT = 2 * A / 3 * sqrt(α * w / G)

    # Sm = KK_critical_supersaturation(param_set, ad, T)
    Sm = KK_critical_supersaturation(param_set, ad, T)

    tmp::FT = sum(1:length(ad.KK_Modes)) do i

        mode_i = ad.KK_Modes[i]

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
function KK_N_activated_per_mode(
    param_set::APS,
    ad::KappaKohlerAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    smax::FT = KK_max_supersaturation(param_set, ad, T, p, w)
    
    sm = KK_critical_supersaturation(param_set, ad, T)
    return ntuple(length(ad.KK_Modes)) do i

        mode_i = ad.KK_Modes[i]
        u_i::FT = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        mode_i.N * (1 / 2) * (1 - erf(u_i))
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
function KK_M_activated_per_mode(
    param_set::APS,
    ad::KappaKohlerAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    smax = KK_max_supersaturation(param_set, ad, T, p, w)
    sm = KK_critical_supersaturation(param_set, ad, T)
    return ntuple(length(ad.KK_Modes)) do i

        mode_i = ad.KK_Modes[i]

        avg_molar_mass_i = sum(1:(mode_i.n_components)) do j
            mode_i.molar_mass[j] * mode_i.mass_mix_ratio[j]
        end

        u_i = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        avg_molar_mass_i * 1 / 2 *
        (1 - erf(u_i - 3 * sqrt(2) / 2 * log(mode_i.stdev)))
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
function KK_total_N_activated(
    param_set::APS,
    ad::KappaKohlerAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    return sum(KK_N_activated_per_mode(param_set, ad, T, p, w))

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
function KK_total_M_activated(
    param_set::APS,
    ad::KappaKohlerAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
) where {FT <: Real}

    return sum(KK_M_activated_per_mode(param_set, ad, T, p, w))

end

# end # module AerosolActivation.jl
