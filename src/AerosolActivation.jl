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

import DataFrames as DF
import MLJ

import ..CommonTypes as CT
import ..Common as CO
import ..AerosolModel as AM
import ..Parameters as CMP
const APS = CMP.AbstractCloudMicrophysicsParameters

export mean_hygroscopicity_parameter
export max_supersaturation

export N_activated_per_mode
export M_activated_per_mode

export total_N_activated
export total_M_activated

"""
    MLEmulatedAerosolActivation

The type for aerosol activation schemes that are emulated with an ML model
"""
struct MLEmulatedAerosolActivation <: CT.AbstractAerosolActivation
    machine::MLJ.Machine
end

function MLEmulatedAerosolActivation(emulator_filepath::String)
    machine = MLJ.machine(emulator_filepath)
    return MLEmulatedAerosolActivation(machine)
end

"""
    coeff_of_curvature(param_set, T)

  - `param_set` - abstract set with Earth's parameters
  - `T` - air temperature

Returns a curvature coefficient.
"""
function coeff_of_curvature(param_set::APS, T::FT) where {FT <: Real}

    _molmass_water::FT = CMP.molmass_water(param_set)
    _gas_constant::FT = CMP.gas_constant(param_set)
    _ρ_cloud_liq::FT = CMP.ρ_cloud_liq(param_set)
    _surface_tension::FT = CMP.surface_tension_coeff(param_set)

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
    return ntuple(Val(AM.n_modes(ad))) do i
        _molmass_water = CMP.molmass_water(param_set)
        _ρ_cloud_liq = CMP.ρ_cloud_liq(param_set)

        FT = eltype(param_set)

        mode_i = ad.Modes[i]

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

        nom / den * _molmass_water / _ρ_cloud_liq
    end
end
function mean_hygroscopicity_parameter(
    param_set::APS,
    ad::AM.AerosolDistribution{NTuple{N, T}},
) where {N, T <: AM.Mode_κ}

    return ntuple(Val(AM.n_modes(ad))) do i
        FT = eltype(param_set)

        mode_i = ad.Modes[i]
        _result = FT(0)
        @inbounds for j in 1:(AM.n_components(mode_i))
            _result += mode_i.vol_mix_ratio[j] * mode_i.kappa[j]
        end
        _result
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
    ad::CT.AbstractAerosolDistribution,
    T::FT,
) where {FT <: Real}
    A::FT = coeff_of_curvature(param_set, T)
    hygro = mean_hygroscopicity_parameter(param_set, ad)
    critical_supersaturation(param_set, ad, T, A, hygro)
end
function critical_supersaturation(
    param_set::APS,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    A,
    hygro,
) where {FT <: Real}
    return ntuple(Val(AM.n_modes(ad))) do i
        2 / sqrt(hygro[i]) * (A / 3 / ad.Modes[i].r_dry)^FT(3 / 2)
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
    scheme::CT.ARG2000Type,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    thermo_params = CMP.thermodynamics_params(param_set)
    _grav::FT = CMP.grav(param_set)
    _ρ_cloud_liq::FT = CMP.ρ_cloud_liq(param_set)

    _ϵ::FT = 1 / CMP.molmass_ratio(param_set)
    R_m::FT = TD.gas_constant_air(thermo_params, q)
    cp_m::FT = TD.cp_m(thermo_params, q)

    L::FT = TD.latent_heat_vapor(thermo_params, T)
    p_vs::FT = TD.saturation_vapor_pressure(thermo_params, T, TD.Liquid())
    G::FT = CO.G_func(param_set, T, TD.Liquid()) / _ρ_cloud_liq

    f_coeff_1::FT = CMP.f_coeff_1_ARG2000(param_set)
    f_coeff_2::FT = CMP.f_coeff_2_ARG2000(param_set)
    g_coeff::FT = CMP.g_coeff_ARG2000(param_set)

    # eq 11, 12 in Razzak et al 1998
    # but following eq 10 from Rogers 1975
    α::FT = L * _grav * _ϵ / R_m / cp_m / T^2 - _grav / R_m / T
    γ::FT = R_m * T / _ϵ / p_vs + _ϵ * L^2 / cp_m / T / p

    A::FT = coeff_of_curvature(param_set, T)
    ζ::FT = 2 * A / 3 * sqrt(α * w / G)

    Sm = critical_supersaturation(param_set, ad, T)

    tmp::FT = FT(0)
    @inbounds for i in 1:AM.n_modes(ad)

        mode_i = ad.Modes[i]

        f::FT = f_coeff_1 * exp(f_coeff_2 * (log(mode_i.stdev))^2)
        g::FT = 1 + g_coeff * log(mode_i.stdev)
        η::FT =
            (α * w / G)^FT(3 / 2) / (FT(2 * pi) * _ρ_cloud_liq * γ * mode_i.N)

        tmp +=
            1 / (Sm[i])^2 *
            (f * (ζ / η)^FT(3 / 2) + g * (Sm[i]^2 / (η + 3 * ζ))^FT(3 / 4))
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
    scheme::CT.AbstractParameterizedAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}
    smax::FT = max_supersaturation(param_set, scheme, ad, T, p, w, q)
    sm = critical_supersaturation(param_set, ad, T)
    N_activated_per_mode(param_set, scheme, ad, T, p, w, q, smax, sm)
end

function N_activated_per_mode(
    param_set::APS,
    scheme::CT.AbstractParameterizedAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
    smax,
    sm,
) where {FT <: Real}
    return ntuple(Val(AM.n_modes(ad))) do i

        mode_i = ad.Modes[i]
        u_i::FT = 2 * log(sm[i] / smax) / 3 / sqrt(2) / log(mode_i.stdev)

        mode_i.N * (1 / 2) * (1 - SF.erf(u_i))
    end
end

function N_activated_per_mode(
    param_set::APS,
    scheme::MLEmulatedAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}
    hygro = mean_hygroscopicity_parameter(param_set, ad)
    return ntuple(Val(AM.n_modes(ad))) do i
        # Model predicts activation of the first mode. So, swap each mode
        # with the first mode repeatedly to predict all activations.
        modes_perm = collect(1:AM.n_modes(ad))
        modes_perm[[1, i]] = modes_perm[[i, 1]]
        per_mode_data = [
            (;
                Symbol("mode_$(j)_N") => ad.Modes[modes_perm[j]].N,
                Symbol("mode_$(j)_mean") => ad.Modes[modes_perm[j]].r_dry,
                Symbol("mode_$(j)_stdev") => ad.Modes[modes_perm[j]].stdev,
                Symbol("mode_$(j)_kappa") => hygro[modes_perm[j]],
            ) for j in 1:AM.n_modes(ad)
        ]
        additional_data = (;
            :velocity => w,
            :initial_temperature => T,
            :initial_pressure => p,
        )
        X = DF.DataFrame([merge(reduce(merge, per_mode_data), additional_data)])
        max(FT(0), min(FT(1), MLJ.predict(scheme.machine, X)[1])) *
        ad.Modes[i].N
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
    scheme::CT.AbstractParameterizedAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}
    smax = max_supersaturation(param_set, scheme, ad, T, p, w, q)
    sm = critical_supersaturation(param_set, ad, T)
    M_activated_per_mode(param_set, scheme, ad, T, p, w, q, smax, sm)
end

function M_activated_per_mode(
    param_set::APS,
    scheme::CT.AbstractParameterizedAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
    smax,
    sm,
) where {FT <: Real}
    return ntuple(Val(AM.n_modes(ad))) do i

        mode_i = ad.Modes[i]

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
    scheme::CT.AbstractAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    return sum(N_activated_per_mode(param_set, scheme, ad, T, p, w, q))

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
    scheme::CT.AbstractAerosolActivation,
    ad::CT.AbstractAerosolDistribution,
    T::FT,
    p::FT,
    w::FT,
    q::TD.PhasePartition{FT},
) where {FT <: Real}

    return sum(M_activated_per_mode(param_set, scheme, ad, T, p, w, q))

end

end # module AerosolActivation.jl
