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
import Thermodynamics.Parameters as TDP

import ..Common as CO
import ..AerosolModel as AM
import ..Parameters as CMP

export mean_hygroscopicity_parameter,
    max_supersaturation,
    N_activated_per_mode,
    M_activated_per_mode,
    total_N_activated,
    total_M_activated

"""
    coeff_of_curvature(ap, T)

  - `ap` - a struct with aerosol activation parameters
  - `T` - air temperature

Returns a curvature coefficient.
"""
function coeff_of_curvature(
    ap::CMP.AerosolActivationParameters,
    T::FT,
) where {FT}
    return FT(2) * ap.σ * ap.M_w / ap.ρ_w / ap.R / T
end

"""
    mean_hygroscopicity_parameter(ap, ad)

  - `ap` - a struct with aerosol activation parameters
  - `ad` - a struct with aerosol distribution (B or κ based)

Returns a tuple of hygroscopicity parameters
(one tuple element for each aerosol size distribution mode).
The tuple is computed either as mass-weighted B parameters
(Abdul-Razzak and Ghan 2000)
or volume weighted kappa parameters (Petters and Kreidenweis 2007).
Implemented via a dispatch based on aerosol distribution mode type.
"""
function mean_hygroscopicity_parameter(
    ap::CMP.AerosolActivationParameters,
    ad::AM.AerosolDistribution{NTuple{N, T}},
) where {N, T <: AM.Mode_B}
    return ntuple(Val(AM.n_modes(ad))) do i
        FT = eltype(ap)
        mode_i = ad.modes[i]

        nom = FT(0)
        @inbounds for j in 1:(AM.n_components(mode_i))
            nom +=
                mode_i.mass_mix_ratio[j] *
                mode_i.dissoc[j] *
                mode_i.osmotic_coeff[j] *
                mode_i.soluble_mass_frac[j] / mode_i.molar_mass[j]
        end

        den = FT(0)
        @inbounds for j in 1:(AM.n_components(mode_i))
            den += mode_i.mass_mix_ratio[j] / mode_i.aerosol_density[j]
        end

        nom / den * ap.M_w / ap.ρ_w
    end
end
function mean_hygroscopicity_parameter(
    ap::CMP.AerosolActivationParameters,
    ad::AM.AerosolDistribution{NTuple{N, T}},
) where {N, T <: AM.Mode_κ}

    return ntuple(Val(AM.n_modes(ad))) do i
        FT = eltype(ap)
        mode_i = ad.modes[i]

        result = FT(0)
        @inbounds for j in 1:(AM.n_components(mode_i))
            result += mode_i.vol_mix_ratio[j] * mode_i.kappa[j]
        end
        result
    end
end

"""
    critical_supersaturation(ap, ad, T)

  - `ap` - a set with aerosol activation parameters
  - `ad` - a struct with aerosol distribution
  - `T` - air temperature

Returns a tuple of critical supersaturations
(one tuple element for each aerosol size distribution mode).
"""
function critical_supersaturation(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    T::FT,
) where {FT}
    A::FT = coeff_of_curvature(ap, T)
    hygro = mean_hygroscopicity_parameter(ap, ad)

    return ntuple(Val(AM.n_modes(ad))) do i
        2 / sqrt(hygro[i]) * (A / 3 / ad.modes[i].r_dry)^FT(3 / 2)
    end
end

"""
    max_supersaturation(ap, ad, aip, tps, T, p, w, q)

  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - a struct with aerosol distribution
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T`   - air temperature
  - `p`   - air pressure
  - `w`   - vertical velocity
  - `q`   - phase partition

Returns the maximum supersaturation.
"""
function max_supersaturation(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT}
    ϵ::FT = 1 / TD.Parameters.molmass_ratio(tps)
    R_m::FT = TD.gas_constant_air(tps, q)
    cp_m::FT = TD.cp_m(tps, q)

    L::FT = TD.latent_heat_vapor(tps, T)
    p_vs::FT = TD.saturation_vapor_pressure(tps, T, TD.Liquid())
    G::FT = CO.G_func(aip, tps, T, TD.Liquid()) / ap.ρ_w

    # eq 11, 12 in Razzak et al 1998
    # but following eq 10 from Rogers 1975
    α::FT = L * ap.g * ϵ / R_m / cp_m / T^2 - ap.g / R_m / T
    γ::FT = R_m * T / ϵ / p_vs + ϵ * L^2 / cp_m / T / p

    A::FT = coeff_of_curvature(ap, T)
    ζ::FT = 2 * A / 3 * sqrt(α * w / G)

    Sm = critical_supersaturation(ap, ad, T)

    tmp::FT = FT(0)
    @inbounds for i in 1:AM.n_modes(ad)

        mode_i = ad.modes[i]

        f::FT = ap.f1 * exp(ap.f2 * (log(mode_i.stdev))^2)
        g::FT = ap.g1 + ap.g2 * log(mode_i.stdev)
        η::FT = (α * w / G)^FT(3 / 2) / (FT(2 * pi) * ap.ρ_w * γ * mode_i.N)

        tmp +=
            1 / (Sm[i])^2 *
            (f * (ζ / η)^ap.p1 + g * (Sm[i]^2 / (η + 3 * ζ))^ap.p2)
    end

    return FT(1) / sqrt(tmp)
end

"""
    N_activated_per_mode(ap, ad, aip, tps, T, p, w, q)

  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - aerosol distribution struct
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T`   - air temperature
  - `p`   - air pressure
  - `w`   - vertical velocity
  - `q`   - phase partition

Returns the number of activated aerosol particles
in each aerosol size distribution mode.
"""
function N_activated_per_mode(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT}
    smax::FT = max_supersaturation(ap, ad, aip, tps, T, p, w, q)
    sm = critical_supersaturation(ap, ad, T)

    return ntuple(Val(AM.n_modes(ad))) do i

        mode_i = ad.modes[i]
        u_i::FT = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        mode_i.N * FT(0.5) * (1 - SF.erf(u_i))
    end
end

"""
    M_activated_per_mode(ap, ad, aip, tps, T, p, w, q)

  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - a struct with aerosol distribution parameters
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the mass of activated aerosol particles
per mode of the aerosol size distribution.
"""
function M_activated_per_mode(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT}
    smax::FT = max_supersaturation(ap, ad, aip, tps, T, p, w, q)
    sm = critical_supersaturation(ap, ad, T)

    return ntuple(Val(AM.n_modes(ad))) do i

        mode_i = ad.modes[i]

        avg_molar_mass_i = FT(0)
        @inbounds for j in 1:(AM.n_components(mode_i))
            avg_molar_mass_i += mode_i.molar_mass[j] * mode_i.mass_mix_ratio[j]
        end

        u_i = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        avg_molar_mass_i * 1 / 2 *
        (1 - SF.erf(u_i - 3 * sqrt(2) / 2 * log(mode_i.stdev)))
    end
end

"""
    total_N_activated(ap, ad, aip, tps, T, p, w, q)

  - `ap` - a struct with aerosol activation parameters
  - `ad` - aerosol distribution struct
  - `aip` - a struct with air properties
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the total number of activated aerosol particles.
"""
function total_N_activated(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT}
    return sum(N_activated_per_mode(ap, ad, aip, tps, T, p, w, q))
end

"""
    total_M_activated(ap, ad, aip, tps, T, p, w, q)

  - `ap` - a struct with aerosol activation parameters
  - `ad` - aerosol distribution struct
  - `aip` - a struct with air properties
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q` - phase partition

Returns the total mass of activated aerosol particles.
"""
function total_M_activated(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDP.ThermodynamicsParameters,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT}
    return sum(M_activated_per_mode(ap, ad, aip, tps, T, p, w, q))
end

end # module AerosolActivation.jl
