"""
    AerosolActivation

Aerosol activation scheme, which includes:

- mean hygroscopicity for each mode of the aerosol size distribution
- critical supersaturation for each mode of the aerosol size distribution
- maximum supersaturation
- total number of particles activated
- total mass of particles activated
"""
module AerosolActivation

import SpecialFunctions as SF
import UnrolledUtilities as UU

import ..ThermodynamicsInterface as TDI
import ..Common as CO
import ..AerosolModel as AM
import ..Parameters as CMP
import ..Utilities as UT

export mean_hygroscopicity_parameter,
    max_supersaturation,
    cloud_droplet_activation_rate,
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

  - `ap` - a struct with aerosol activation parameters
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
    max_supersaturation(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice)

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
    N_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice)

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
    M_activated_per_mode(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice)

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
    total_N_activated(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice)

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
    total_M_activated(ap, ad, aip, tps, T, p, w, q_tot, q_liq, q_ice)

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

"""
    _activation_lognormal_argument(S_crit, S, stdev)

The lognormal tail argument `u = 2 log(S_crit / S) / (3 sqrt(2) log(stdev))` whose complementary
error function gives the activated fraction of a mode.

Written once because it is shared between the parcel-maximum entry points above and the
supersaturation-taking ones below, and guarded because both are reached from inside a GPU kernel
where a `DomainError` aborts the whole kernel: `S` is floored at `eps(FT)` so the logarithm stays
in domain, and `log(stdev)` likewise, so a monodisperse mode gives `u = ±Inf` and an activated
fraction of exactly zero or one instead of dividing by zero. Below the `S` floor the argument is
large and positive, `erf` is at or within rounding of one, and both the activated number and its
derivative fall to a negligible fraction of the aerosol budget, which is the correct answer there.
"""
@inline function _activation_lognormal_argument(S_crit, S, stdev)
    FT = UT.promote_typeof(S_crit, S, stdev)
    S_pos = max(S, eps(FT))
    S_crit_pos = max(S_crit, eps(FT))
    lnσ = max(log(stdev), eps(FT))
    return 2 * log(S_crit_pos / S_pos) / (3 * sqrt(FT(2)) * lnσ)
end

"""
    N_activated_per_mode(ap, ad, T, S)
    total_N_activated(ap, ad, T, S)

Activated number concentration [1/m³], per mode and summed, at a GIVEN supersaturation `S`
rather than at the maximum an adiabatic parcel rising at some updraft speed would reach.

The parcel-maximum entry points above take the updraft speed and solve for the supersaturation
themselves, and `max_supersaturation` is not analytically invertible in the updraft speed, so a
caller that already knows the supersaturation it wants to activate at has no way to say so. That
caller exists: a cell whose AMBIENT supersaturation already exceeds the parcel maximum has more
aerosol above their critical supersaturation than the parcel calculation admits, and evaluating
the activated number at the parcel maximum understates it. These methods are the tail of
[`N_activated_per_mode`](@ref) with the parcel solve replaced by the argument.

`S` is a fractional supersaturation, not a percentage, and not a saturation ratio: `0.002` is
0.2 percent, the order a cloud base reaches.
"""
@inline function N_activated_per_mode(
    ap::CMP.AerosolActivationParameters, ad::CMP.AerosolDistributionType, T, S,
)
    sm = critical_supersaturation(ap, ad, T)
    return ntuple(Val(AM.n_modes(ad))) do i
        mode_i = ad.modes[i]
        u_i = _activation_lognormal_argument(sm[i], S, mode_i.stdev)
        mode_i.N * (1 - SF.erf(u_i)) / 2
    end
end

@inline total_N_activated(
    ap::CMP.AerosolActivationParameters, ad::CMP.AerosolDistributionType, T, S,
) = sum(N_activated_per_mode(ap, ad, T, S))

"""
    ∂N_activated_∂S(ap, ad, T, S)

Derivative of [`total_N_activated`](@ref) with respect to the supersaturation, [1/m³] per unit
`S`, in closed form:

```math
∂N/∂S = Σᵢ \\frac{2 Nᵢ}{3 \\sqrt{2π} \\log(σᵢ) S} e^{-uᵢ²}
```

Positive everywhere: more supersaturation activates more aerosol. Supplied analytically rather
than left to the linearization because it is the only route by which the substep's vapor budget
reaches the activation source, and the exponential makes it vanish smoothly on both tails.
"""
@inline function ∂N_activated_∂S(
    ap::CMP.AerosolActivationParameters, ad::CMP.AerosolDistributionType, T, S,
)
    sm = critical_supersaturation(ap, ad, T)
    FT = UT.promote_typeof(S, sm[1])
    S_pos = max(S, eps(FT))
    return sum(
        ntuple(Val(AM.n_modes(ad))) do i
            mode_i = ad.modes[i]
            u_i = _activation_lognormal_argument(sm[i], S, mode_i.stdev)
            lnσ = max(log(mode_i.stdev), eps(FT))
            2 * mode_i.N * exp(-u_i^2) / (3 * sqrt(2 * FT(π)) * lnσ * S_pos)
        end,
    )
end

"""
    ACTIVATION_MIN_UPDRAFT

Updraft speed [m/s] below which the adiabatic-parcel branch of the activation supersaturation is
switched off.

`max_supersaturation` forms `sqrt(α w / G)`, so it throws a `DomainError` at `w < 0` and returns
`NaN` at `w == 0` exactly, where `ζ/η` is `0/0`. Inside a GPU kernel a `DomainError` aborts the
whole kernel and hides the state that caused it, so the parcel branch is selected rather than
entered: below this speed the ambient supersaturation is the only term, which is the physically
right answer in still or subsiding air.
"""
const ACTIVATION_MIN_UPDRAFT = 1e-4

"""
    cloud_droplet_activation_rate(ap, pa, aip, tps, T, p, w, ρₐ, q_tot, q_liq, q_ice, n_lcl, x_seed)

Cloud droplet activation as a nucleation-class source: droplet number and the mass those droplets
carry, together.

```math
∂ₜn = \\max(0, N_{act}(S_{eff})/ρₐ - n) / τ_{act}(S_{eff}), \\qquad ∂ₜq = x_{seed} ∂ₜn
```

# The supersaturation it activates at, and the regime that is valid in

`S_eff = max(S_ambient, S_max(w))`, the larger of the supersaturation the cell already carries and
the maximum an adiabatic parcel rising at `w` could itself generate, with the parcel term admitted
only where the cell is at or above liquid saturation, since that is the state the ARG parcel
calculation starts from and no droplet passes its critical radius below it.

**The parcel logic is valid for `S_ambient <= S_max`.** There, the resolved updraft is what
generates the supersaturation, the ARG closure solves the parcel balance for its peak, and the
activated number is that peak's Koehler cut through the aerosol spectrum. That is the regime the
scheme was built for and it is where `S_eff = S_max`.

**Above it the scheme is outside that regime and says so.** A cell whose ambient supersaturation
already exceeds anything its own updraft could produce is not a rising parcel; its supersaturation
was put there by something else, radiative cooling being the common case. What is returned there is
the KOEHLER/TWOMEY LOWER BOUND: every aerosol particle whose critical supersaturation is below the
ambient `S` is above its barrier and activates, which is a bound rather than a solution because it
neglects the depletion those droplets would themselves cause. It is capped by the aerosol budget
`N_a`, since the sum over modes cannot exceed the particles that exist. A lower bound is the right
thing to return when the closure's own assumptions do not hold: it is defensible in the direction
it errs.

Returning zero there instead would reverse the physics: if the ambient supersaturation is already
above the parcel maximum then MORE aerosol are above their critical supersaturation than the parcel
calculation admits, so the activated number is at least the parcel number and never zero. The
reversal is also self-reinforcing, since the condition is met more firmly the more supersaturated
the air becomes, and it can switch activation off across a whole persistently supersaturated
domain.

The returned `outside_parcel_regime` is that distinction, per cell, and it is the DIAGNOSABLE
SIGNAL for it: a host or an analysis that wants to know whether a run is spending its time outside
the parcel-valid regime counts the fraction of cells where it is true, rather than assuming one way
or the other.

# Why it is a relaxation, and not an instantaneous adjustment

An aerosol particle does not become a cloud droplet the moment it is above its critical
supersaturation; it has to grow by diffusion through the Köhler barrier to the activation radius,
and that takes `τ_act = r_act² / (2 G S)` with `G` the volumetric diffusional growth coefficient
the condensation closure already uses. Both factors are state functions, so no timestep appears
anywhere in the rate: at 0.2 percent supersaturation this is seconds, at 20 percent it is tens of
milliseconds, and the substep's implicit update recovers the instantaneous limit by itself
wherever `h ≫ τ_act`. The reciprocal is carried rather than the timescale so that a vanishing
supersaturation gives a rate of exactly zero instead of a division by one.

# Why it carries mass

The droplets it creates are droplets. A source that supplies number alone leaves the category in
a state the scheme has no size for, and the number adjustment then correctly drains it, so the
population can never establish itself. Each new droplet arrives at `x_seed`, the distribution's
own minimum droplet mass, which is exactly a 1 μm droplet and exactly the size the activation
radius describes. This is the same pairing rule the deposition-nucleation source follows on the
ice side, one phase over.

# Totality

Every argument that reaches `max_supersaturation` is sanitized first, and the parcel branch is
SELECTED rather than entered, because that function throws a `DomainError` on a negative updraft
or on a state whose condensate exceeds its total water, and returns `NaN` at zero updraft. The
pre-existing-droplet sink is deliberately not passed to it: that term takes the cube root of the
mean droplet volume, whose derivative is infinite at the mass-free state this source exists to
rescue, and the depletion it represents is already carried by `S_ambient` through the vapor
budget.

With `w` or `p` at zero, which is what the defaulted arguments of the tendency entry points
supply, the parcel branch contributes nothing and the AMBIENT branch is evaluated as usual. The
source is therefore NOT switched off by a zero updraft, and deliberately so: radiatively driven
supersaturation persists at negligible resolved vertical velocity, and an updraft gate is what
keeps a cloud from ever forming there. The rate is inert only where the air is subsaturated over
liquid or the prescribed aerosol has no particles.

# Arguments
 - `ap`: [`CMP.AerosolActivationParameters`](@ref), the ARG2000 fit
 - `pa`: [`CMP.PrescribedAerosol`](@ref), the aerosol modes to activate, or `nothing`
 - `aip`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `T`: temperature [K]
 - `p`: air pressure [Pa]; non-positive switches the parcel branch off
 - `w`: vertical velocity [m/s]; at or below [`ACTIVATION_MIN_UPDRAFT`](@ref) the parcel branch is off
 - `ρₐ`: air density [kg/m³]
 - `q_tot`: total water specific content [kg/kg]
 - `q_liq`: liquid water specific content [kg/kg]
 - `q_ice`: ice water specific content [kg/kg]
 - `n_lcl`: cloud droplet specific number [1/kg]
 - `x_seed`: mass of a newly activated droplet [kg]

# Returns
A `NamedTuple` `(; ∂ₜn_lcl, ∂ₜq_lcl, inv_τ_act, ∂ₜn_∂S, outside_parcel_regime, qᵥ_sat)`.

`inv_τ_act` and `∂ₜn_∂S` are for the Jacobian: the relaxation rate that goes on the number
diagonal, and the closed-form derivative of the number source with respect to the supersaturation.

`outside_parcel_regime` is true when the ambient supersaturation is at or above the parcel maximum,
i.e. when the cell is outside the regime the ARG parcel closure is valid in and the returned number
is the Koehler/Twomey lower bound rather than the parcel solution. It serves twice: the Jacobian
uses it to gate the vapor coupling, which exists only on the ambient branch, and it is the cheap
per-cell signal from which a fraction-of-domain diagnostic is counted.
"""
@inline function cloud_droplet_activation_rate(
    ap::CMP.AerosolActivationParameters, pa::CMP.PrescribedAerosol,
    aip::CMP.AirProperties, tps::TDI.PS,
    T, p, w, ρₐ, q_tot, q_liq, q_ice, n_lcl, x_seed,
)
    FT = UT.promote_typeof(T, p, w, ρₐ, q_tot, q_liq, q_ice, n_lcl)
    ad = AM.aerosol_distribution(pa)

    q_liq_s = UT.clamp_to_nonneg(q_liq)
    q_ice_s = UT.clamp_to_nonneg(q_ice)
    qᵥ = UT.clamp_to_nonneg(q_tot - q_liq_s - q_ice_s)
    qᵥ_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρₐ)
    S_ambient = qᵥ / max(qᵥ_sat, eps(FT)) - 1

    # The parcel branch, selected rather than entered, and evaluated on a physically enveloped
    # copy of the state rather than on the state itself. `max_supersaturation` is not a total
    # function, and inside a GPU kernel a `DomainError` aborts the whole kernel and hides the
    # state that caused it, so the guarding is done here. Four ways it fails:
    #
    #   w < 0                 sqrt of a negative, throws
    #   w = 0 exactly         ζ/η is 0/0, returns NaN, which `max(0, ·)` does not sanitize
    #   no aerosol at all     the mode sum is zero, `1/sqrt(0)` is Inf, `Inf/Inf` is NaN
    #   condensate of order 1 α ∝ (Lᵥ R_m / (R_v cp_m T) − 1) turns negative, because cp_m is
    #                         then the heat capacity of water rather than of air; sqrt throws
    #
    # The envelope is admissible because the parcel maximum is a secondary term: it is selected
    # only where it EXCEEDS the ambient supersaturation, which is near-saturated ascending air
    # carrying little condensate, and it is exactly there that the envelope is inert.
    #
    # The parcel branch also requires the cell to be at or above liquid saturation. The ARG
    # calculation asks what peak supersaturation a parcel rising through CLOUD BASE reaches, and
    # it assumes the parcel is at saturation when it starts; applied to subsaturated air it
    # returns a positive peak for a parcel that would never reach saturation at all, and no
    # droplet can pass its critical radius there whatever the updraft does. Without this the
    # source would fire over most of a domain, since most of a domain is subsaturated and half of
    # it is rising.
    n_aer = pa.N_accum + pa.N_coarse
    parcel_ok = (w > oftype(w, ACTIVATION_MIN_UPDRAFT)) & (p > 0) & (qᵥ > 0) &
                (n_aer > 0) & (S_ambient >= 0)
    q_liq_p = min(q_liq_s, FT(0.05))
    q_ice_p = min(q_ice_s, FT(0.05))
    q_tot_p = clamp(FT(q_tot), q_liq_p + q_ice_p + FT(1e-6), FT(0.2))
    T_p = FT(clamp(T, oftype(T, 150), oftype(T, 350)))
    p_p = FT(clamp(p, oftype(p, 1e3), oftype(p, 1.2e5)))
    w_p = ifelse(parcel_ok, clamp(FT(w), FT(ACTIVATION_MIN_UPDRAFT), FT(100)), one(FT))
    S_parcel = max_supersaturation(
        ap, ad, aip, tps,
        T_p, p_p, w_p, q_tot_p, q_liq_p, q_ice_p, zero(FT), zero(FT),
    )
    S_eff = max(S_ambient, ifelse(parcel_ok, S_parcel, zero(FT)))
    outside_parcel_regime = S_ambient >= ifelse(parcel_ok, S_parcel, zero(FT))

    # Diffusional growth through the Köhler barrier to the activation radius.
    ρ_w = ap.ρ_w
    r_seed = cbrt(3 * x_seed / (4 * FT(π) * ρ_w))
    G_volumetric = CO.G_func_liquid(aip, tps, T) / ρ_w
    inv_τ_act = 2 * G_volumetric * max(S_eff, zero(FT)) / r_seed^2

    n_act = total_N_activated(ap, ad, T, S_eff) / ρₐ
    deficit = max(n_act - n_lcl, zero(FT))
    active = (S_eff > 0) & (deficit > 0)

    ∂ₜn_lcl = ifelse(active, deficit * inv_τ_act, zero(FT))
    ∂ₜq_lcl = x_seed * ∂ₜn_lcl

    # ∂(∂ₜn)/∂S = (∂n_act/∂S)/τ + deficit ∂(1/τ)/∂S, and ∂(1/τ)/∂S = (1/τ)/S.
    ∂n_act_∂S = ∂N_activated_∂S(ap, ad, T, S_eff) / ρₐ
    ∂ₜn_∂S = ifelse(
        active,
        ∂n_act_∂S * inv_τ_act + deficit * inv_τ_act / max(S_eff, eps(FT)),
        zero(FT),
    )

    return (; ∂ₜn_lcl, ∂ₜq_lcl, inv_τ_act = ifelse(active, inv_τ_act, zero(FT)),
        ∂ₜn_∂S, outside_parcel_regime, qᵥ_sat)
end

"""
    cloud_droplet_activation_rate(ap, ::Nothing, aip, tps, T, p, w, ρₐ, q_tot, q_liq, q_ice, n_lcl, x_seed)

Droplet activation for a parameter set carrying no aerosol population: the same `NamedTuple`, with
every rate and every derivative exactly zero.

`nothing` is the default aerosol of [`CMP.PrescribedAerosol`](@ref)'s consumers, so this is the
method a configuration that has not stated an aerosol population reaches. It is a computed zero
rather than an error because an aerosol-free configuration is a legitimate one: the droplet number
is then supplied by whatever else the host prescribes, and every other microphysical process still
runs.

`qᵥ_sat` is the saturation vapor specific content over liquid, which does not depend on the
aerosol, so it is returned as it is on the populated branch and stays usable by a caller that
reads it.
"""
@inline function cloud_droplet_activation_rate(
    ap::CMP.AerosolActivationParameters, ::Nothing,
    aip::CMP.AirProperties, tps::TDI.PS,
    T, p, w, ρₐ, q_tot, q_liq, q_ice, n_lcl, x_seed,
)
    FT = UT.promote_typeof(T, p, w, ρₐ, q_tot, q_liq, q_ice, n_lcl)
    o = zero(FT)
    qᵥ_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, T, ρₐ)
    return (; ∂ₜn_lcl = o, ∂ₜq_lcl = o, inv_τ_act = o, ∂ₜn_∂S = o,
        outside_parcel_regime = false, qᵥ_sat)
end

end # module AerosolActivation.jl
