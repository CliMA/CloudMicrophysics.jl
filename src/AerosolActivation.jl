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
import UnrolledUtilities as UU

import ..ThermodynamicsInterface as TDI
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
    max_supersaturation(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)

  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - a struct with aerosol distribution
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T`   - air temperature
  - `p`   - air pressure
  - `w`   - vertical velocity
  - `q_tot` - total water specific content
  - `q_liq` - liquid water specific content
  - `q_ice` - ice water specific content
  - `N_liq` - liquid water number concentration
  - `N_ice` - ice water number concentration

Returns the maximum supersaturation.
"""
function max_supersaturation(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    N_liq::FT,
    N_ice::FT,
) where {FT}
    R_v::FT = TDI.Rᵥ(tps)
    R_m::FT = TDI.Rₘ(tps, q_tot, q_liq, q_ice)
    cp_m::FT = TDI.cpₘ(tps, q_tot, q_liq, q_ice)

    Lᵥ::FT = TDI.Lᵥ(tps, T)
    ρ_air = TDI.air_density(tps, T, p, q_tot, q_liq, q_ice)
    p_v::FT = (q_tot - q_liq - q_ice) * ρ_air * R_v * T
    p_vs::FT = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    G::FT = CO.G_func_liquid(aip, tps, T) / ap.ρ_w

    # eq 11, 12 in Razzak et al 1998
    # but following eq A11 from Korolev and Mazin 2003
    α::FT = p_v / p_vs * (Lᵥ * ap.g / R_v / cp_m / T^2 - ap.g / R_m / T)
    γ::FT = R_v * T / p_vs + p_v / p_vs * R_m * Lᵥ^2 / R_v / cp_m / T / p

    A::FT = coeff_of_curvature(ap, T)
    ζ::FT = 2 * A / 3 * sqrt(α * w / G)

    Sm = critical_supersaturation(ap, ad, T)

    tmp::FT = FT(0)
    @inbounds for i in 1:AM.n_modes(ad)

        mode_i = ad.modes[i]

        f::FT = ap.f1 * exp(ap.f2 * (log(mode_i.stdev))^2)
        g::FT = ap.g1 + ap.g2 * log(mode_i.stdev)
        η::FT = sqrt(α * w / G)^3 / (FT(2 * pi) * ap.ρ_w * γ * mode_i.N)

        tmp +=
            1 / (Sm[i])^2 *
            (f * (ζ / η)^ap.p1 + g * (Sm[i]^2 / (η + 3 * ζ))^ap.p2)
    end
    S_max_ARG::FT = FT(1) / sqrt(tmp)

    r_liq::FT = N_liq < eps(FT) ? FT(0) : cbrt(ρ_air * q_liq / N_liq / ap.ρ_w / FT(4 / 3 * π))
    K_liq::FT = FT(4 * π) * ap.ρ_w * N_liq * r_liq * G * γ

    Lₛ::FT = TDI.Lₛ(tps, T)
    γᵢ::FT = R_v * T / p_vs + p_v / p_vs * R_m * Lᵥ * Lₛ / R_v / cp_m / T / p
    r_ice::FT = N_ice < eps(FT) ? FT(0) : cbrt(ρ_air * q_ice / N_ice / ap.ρ_i / FT(4 / 3 * π))
    ρᵢGᵢ::FT = CO.G_func_ice(aip, tps, T)
    ξ::FT = TDI.saturation_vapor_pressure_over_liquid(tps, T) / TDI.saturation_vapor_pressure_over_ice(tps, T)
    K_ice::FT = FT(4 * π) * N_ice * r_ice * ρᵢGᵢ * γᵢ

    S_max::FT = S_max_ARG * (α * w - K_ice * (ξ - FT(1))) / (α * w + (K_liq + K_ice * ξ) * S_max_ARG)

    return max(FT(0), S_max)
end
function max_supersaturation(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return max_supersaturation(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, FT(0), FT(0))
end

"""
    N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)

  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - aerosol distribution struct
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T`   - air temperature
  - `p`   - air pressure
  - `w`   - vertical velocity
  - `q_tot` - total water specific content
  - `q_liq` - liquid water specific content
  - `q_ice` - ice water specific content
  - `N_liq` - liquid water number concentration
  - `N_ice` - ice water number concentration

Returns the number of activated aerosol particles
in each aerosol size distribution mode.
"""
function N_activated_per_mode(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    N_liq::FT,
    N_ice::FT,
) where {FT}
    smax::FT = max_supersaturation(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)
    sm = critical_supersaturation(ap, ad, T)

    return ntuple(Val(AM.n_modes(ad))) do i

        mode_i = ad.modes[i]
        u_i = 2 * log(sm[i] / smax) / 3 / sqrt(FT(2)) / log(mode_i.stdev)

        mode_i.N * FT(0.5) * (1 - SF.erf(u_i))
    end
end
function N_activated_per_mode(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, FT(0), FT(0))
end

"""
    M_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)

  - `ap`  - a struct with aerosol activation parameters
  - `ad`  - a struct with aerosol distribution parameters
  - `aip` - a struct with air parameters
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q_tot` - total water specific content
  - `q_liq` - liquid water specific content
  - `q_ice` - ice water specific content
  - `N_liq` - liquid water number concentration
  - `N_ice` - ice water number concentration

Returns the mass of activated aerosol particles
per mode of the aerosol size distribution.
"""
function M_activated_per_mode(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    N_liq::FT,
    N_ice::FT,
) where {FT}
    smax::FT = max_supersaturation(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)
    sm = critical_supersaturation(ap, ad, T)

    return ntuple(Val(AM.n_modes(ad))) do i
        mode_i = ad.modes[i]
        Mᵢ = UU.unrolled_sum(mode_i.molar_mass .* mode_i.mass_mix_ratio)
        σᵢ = mode_i.stdev
        fac = 3log(σᵢ) * √(FT(2)) / 2  # 3√2/2 log(σᵢ), shared factor in `erf`
        u_i = log(sm[i] / smax) / fac

        # erfc(x) ≡ 1 - erf(x), but more accurate for large x
        Mᵢ / 2 * SF.erfc(u_i - fac)
    end
end
function M_activated_per_mode(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return M_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, FT(0), FT(0))
end

"""
    total_N_activated(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)

  - `ap` - a struct with aerosol activation parameters
  - `ad` - aerosol distribution struct
  - `aip` - a struct with air properties
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q_tot` - total water specific content
  - `q_liq` - liquid water specific content
  - `q_ice` - ice water specific content
  - `N_liq` - liquid water number concentration
  - `N_ice` - ice water number concentration

Returns the total number of activated aerosol particles.
"""
function total_N_activated(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    N_liq::FT,
    N_ice::FT,
) where {FT}
    return sum(N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice))
end
function total_N_activated(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return sum(N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice))
end

"""
    total_M_activated(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice)

  - `ap` - a struct with aerosol activation parameters
  - `ad` - aerosol distribution struct
  - `aip` - a struct with air properties
  - `tps` - a struct with thermodynamics parameters
  - `T` - air temperature
  - `p` - air pressure
  - `w` - vertical velocity
  - `q_tot` - total water specific content
  - `q_liq` - liquid water specific content
  - `q_ice` - ice water specific content
  - `N_liq` - liquid water number concentration
  - `N_ice` - ice water number concentration

Returns the total mass of activated aerosol particles.
"""
function total_M_activated(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
    N_liq::FT,
    N_ice::FT,
) where {FT}
    return sum(M_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice, N_liq, N_ice))
end
function total_M_activated(
    ap::CMP.AerosolActivationParameters,
    ad::CMP.AerosolDistributionType,
    aip::CMP.AirProperties,
    tps::TDI.PS,
    T::FT,
    p::FT,
    w::FT,
    q_tot::FT,
    q_liq::FT,
    q_ice::FT,
) where {FT}
    return sum(M_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice))
end

end # module AerosolActivation.jl
