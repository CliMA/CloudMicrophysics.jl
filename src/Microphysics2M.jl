"""
Double-moment bulk microphysics parametrizations including:
 - autoconversion, accretion, self-collection, breakup, mean terminal velocity of raindrops,
    and rain evaporation rates from Seifert and Beheng 2006.
 - number concentration adjustment from Horn 2012.
 - additional double-moment bulk microphysics autoconversion and accretion rates
   from: Khairoutdinov and Kogan 2000, Beheng 1994, Tripoli and Cotton 1980, and
   Liu and Daum 2004.
"""
module Microphysics2M

import SpecialFunctions as SF
import RootSolvers as RS
import ForwardDiff as FD

import ..ThermodynamicsInterface as TDI
import ..Common as CO
import ..Parameters as CMP
import ..DistributionTools as DT
import ..Utilities as UT

import ..DistributionTools: size_distribution

export activation_droplet_mass,
    autoconversion,
    accretion,
    cloud_liquid_self_collection,
    autoconversion_and_cloud_liquid_self_collection,
    rain_terminal_velocity,
    conv_q_lcl_to_q_rai,
    rain_evaporation,
    ∂rain_evaporation_∂N_rai_∂q_rai,
    rain_self_collection,
    rain_breakup,
    rain_self_collection_and_breakup,
    rain_equilibrium_number,
    rain_number_relaxation,
    size_distribution,
    get_size_distribution_bounds,
    number_tendency_from_mass_limits,
    orphan_mass_drain,
    orphan_mass_drain_ice,
    orphan_mass_inv_timescale,
    orphan_mass_inv_timescale_ice

"""
    pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)

Return the parameters of the rain drop diameter distribution

    n_r(D) = N_0 * exp(- D / Dr_mean)

 where
 - `D` is the diameter of the raindrop,
 - `N_0` [1/m⁴] is the intercept parameter of the distribution,
 - `Dr_mean` [m] is the mean diameter of the raindrops.

 Note: in SB2006, Eq. (83) the distribution is given as:

    f(D) = N_0 * exp(- λ_r D)

 where `λ_r ≡ 1 / Dr_mean` [1/m] is the inverse of the mean diameter of the raindrops.

# Arguments
 - `pdf_r`: struct containing size distribution parameters for rain, one of
        [`CMP.RainParticlePDF_SB2006_notlimited`](@ref) (the honest inversion) or
        [`CMP.RainParticlePDF_SB2006_windowed`](@ref) (one bound, on the mean drop mass).
 - `qᵣ`: rain water specific content [kg/kg]
 - `ρₐ`: air density [kg/m³]
 - `Nᵣ`: number of rain drops [1/m³]

# Returns
 - A `NamedTuple` with the fields `(; N₀r, Dr_mean, xr_mean)`
"""
function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_notlimited, qᵣ, ρₐ, Nᵣ)
    FT = UT.promote_typeof(qᵣ, ρₐ, Nᵣ)
    (; ρw) = pdf_r
    safe_qᵣ = max(qᵣ, UT.ϵ_numerics_2M_M(FT))
    safe_Nᵣ = max(Nᵣ, UT.ϵ_numerics_2M_N(FT))
    Lᵣ = ρₐ * safe_qᵣ

    xr_mean = Lᵣ / safe_Nᵣ
    λr = cbrt(π * ρw / xr_mean)
    N₀r = λr * safe_Nᵣ

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    # The predicate is evaluated once into `cond`; the three `ifelse`es are predicated
    # selects reusing it, not three re-evaluations. This replaces an early-return branch,
    # so all warp lanes stay on one instruction stream (no GPU warp divergence) - the
    # three extra selects are far cheaper than a data-dependent branch (cf. PR #749).
    cond = Nᵣ < UT.ϵ_numerics_2M_N(FT) || qᵣ < UT.ϵ_numerics_2M_M(FT)
    return (;
        N₀r = ifelse(cond, zero(N₀r), N₀r),
        Dr_mean = ifelse(cond, zero(Dr_mean), Dr_mean),
        xr_mean = ifelse(cond, zero(xr_mean), xr_mean),
    )
end

function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_limited, qᵣ, ρₐ, Nᵣ)
    FT = UT.promote_typeof(qᵣ, ρₐ, Nᵣ)
    (; xr_min, xr_max, N0_min, N0_max, λ_min, λ_max, ρw) = pdf_r
    safe_qᵣ = max(qᵣ, UT.ϵ_numerics_2M_M(FT))
    safe_Nᵣ = max(Nᵣ, UT.ϵ_numerics_2M_N(FT))
    Lᵣ = ρₐ * safe_qᵣ

    # Sequence of limiting steps in Seifert and Beheng 2006:
    x̃r = clamp(Lᵣ / safe_Nᵣ, xr_min, xr_max)  # Eq. (94)
    N₀r = clamp(safe_Nᵣ * cbrt(π * ρw / x̃r), N0_min, N0_max)  # Eq. (95)
    λr = clamp(sqrt(sqrt(π * ρw * N₀r / Lᵣ)), λ_min, λ_max)  # Eq. (96)
    xr_mean = clamp(Lᵣ * λr / N₀r, xr_min, xr_max)  # Eq. (97)

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    cond = Nᵣ < UT.ϵ_numerics_2M_N(FT) && qᵣ < UT.ϵ_numerics_2M_M(FT)
    return (;
        N₀r = ifelse(cond, zero(N₀r), N₀r),
        Dr_mean = ifelse(cond, zero(Dr_mean), Dr_mean),
        xr_mean = ifelse(cond, zero(xr_mean), xr_mean),
    )
end
function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_windowed, qᵣ, ρₐ, Nᵣ)
    FT = UT.promote_typeof(qᵣ, ρₐ, Nᵣ)
    (; xr_min, xr_max, ρw) = pdf_r
    safe_qᵣ = max(qᵣ, UT.ϵ_numerics_2M_M(FT))
    safe_Nᵣ = max(Nᵣ, UT.ϵ_numerics_2M_N(FT))
    Lᵣ = ρₐ * safe_qᵣ

    xr_mean = clamp(Lᵣ / safe_Nᵣ, xr_min, xr_max)
    Nᵣ_bounded = Lᵣ / xr_mean
    λr = cbrt(π * ρw / xr_mean)
    N₀r = λr * Nᵣ_bounded

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    cond = Nᵣ < UT.ϵ_numerics_2M_N(FT) || qᵣ < UT.ϵ_numerics_2M_M(FT)
    return (;
        N₀r = ifelse(cond, zero(N₀r), N₀r),
        Dr_mean = ifelse(cond, zero(Dr_mean), Dr_mean),
        xr_mean = ifelse(cond, zero(xr_mean), xr_mean),
    )
end

"""
    pdf_rain_parameters_mass(pdf_r::CMP.RainParticlePDF_SB2006, qᵣ, ρₐ, Nᵣ)

Return the parameters of the rain drop diameter distribution in terms of mass.

As a function of diameter, the size distribution is given by:

    n(D) = N₀r * exp(-D / Dr_mean)

In terms of mass (`x`), the size distribution is given by:

    f(x) = n(D(x)) * ∂D∂x(x)
         = N₀ * exp(-D(x) / Dr_mean) * (6 / (π * ρw))^(1/3) / 3 * x^(-2/3)
         = N₀ / 3 * (6 / (π * ρw))^(1/3) * x^(-2/3) * exp(- (6 / (π * ρw))^(1/3) / Dr_mean * x^(1/3))

 where
 - `D(x) = (6x / (π * ρw))^(1/3)` is the diameter of a raindrop of mass `x`.
 - `∂D∂x(x) = (6 / (π * ρw))^(1/3) / 3 * x^(-2/3)` is the derivative of the diameter with respect to the mass.

If we write the general form of the size distribution as:

    f(x) = A * x^ν * exp(-B * x^μ)

 then we have that:
 - `A = N₀ / 3 * (6 / (π * ρw))^(1/3)`
 - `B = (6 / (π * ρw))^(1/3) / Dr_mean`
 - `ν = -2/3`
 - `μ = 1/3`

# Returns
 - A `NamedTuple` with the fields `(; Ar, Br)`, where `Ar = A` and `Br = B`
   expressed via the mean raindrop mass:
   `Ar = Nᵣ * (6 / xr_mean)^(1/3) / 3` and `Br = (6 / xr_mean)^(1/3)`.
"""
function pdf_rain_parameters_mass(pdf_r::CMP.RainParticlePDF_SB2006, qᵣ, ρₐ, Nᵣ)
    (; xr_mean) = pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)
    Br = cbrt(6 / xr_mean)
    Ar = Nᵣ * Br / 3
    return (; Ar, Br)
end

"""
    activation_droplet_mass(pdf_c)

Mass of a freshly activated cloud droplet [kg], the smallest droplet the size distribution
resolves. Doubles as the floor on the mean droplet mass in
[`log_pdf_cloud_parameters_mass`](@ref) and as the mass a droplet activation source pairs with
its number tendency.

# Arguments
- `pdf_c`: [`CMP.CloudParticlePDF_SB2006`](@ref)

# Returns
- Activation droplet mass [kg]
"""
@inline activation_droplet_mass(pdf_c) = pdf_c.xc_min

"""
    cloud_mean_droplet_mass_and_number(pdf_c, q, ρₐ, N)

Return `(x̄, N_eff)`: the cloud mean droplet mass bounded to
`[activation_droplet_mass(pdf_c), xc_max]`, and the droplet number consistent with it. At the
upper bound the number is rescaled to `ρₐ q / xc_max`, the same choice
[`pdf_rain_parameters`](@ref) makes for rain; at the lower bound `N` is returned unchanged, since
a mass-preserving rescale there would return zero droplets for a mass-free, number-carrying
population. `N_eff` equals `max(N, ϵ_numerics_2M_N)` wherever the upper bound does not bind.

# Arguments
 - `pdf_c`: Size distribution parameters for cloud droplets, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: Liquid mass content [kg/kg]
 - `ρₐ`: Air density [kg/m³]
 - `N`: Number concentration of the particle [1/m³]

# Returns
 - `(x̄, N_eff)`: mean droplet mass [kg] within `[xc_min, xc_max]`, and the number [1/m³]
   consistent with it
"""
function cloud_mean_droplet_mass_and_number(pdf_c, q, ρₐ, N)
    FT = UT.promote_typeof(q, ρₐ, N)
    (; xc_max) = pdf_c
    safe_N = max(N, UT.ϵ_numerics_2M_N(FT))
    L = ρₐ * UT.clamp_to_nonneg(q)
    x_raw = L / safe_N
    x̄ = clamp(x_raw, FT(activation_droplet_mass(pdf_c)), FT(xc_max))
    N_eff = ifelse(x_raw > FT(xc_max), L / FT(xc_max), safe_N)
    return (x̄, N_eff)
end

"""
    log_pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)

Return the log of the parameters of the generalized gamma distribution of the form

    f(x) = A * x^ν * exp(-B * x^μ),  [Eq. (79) in Seifert and Beheng 2006, but using the symbol `B` instead of `λ`]

where

    B = [  x̄ Γ(z₁) / Γ(z₂) ]^(-μ)
    A = μ N B^(z₁) / Γ(z₁)
    z₁ = (ν + 1) / μ
    z₂ = (ν + 2) / μ

That is,

    log(B) = - μ [ log(x̄) + logΓ(z₁) - logΓ(z₂) ]
    log(A) = log(μ) + log(N) + z₁ * log(B) - logΓ(z₁)

The distribution is degenerate when `N` falls below [`UT.ϵ_numerics_2M_N`](@ref); `logA = -Inf`
and `logB = Inf` then send every moment to zero. Presence is decided by `N` alone rather than by
`q`, since a two-moment category with `N` droplets and negligible mass content is a real
population still acquiring its mass, not a degenerate one.

The mean droplet mass is bounded to `[activation_droplet_mass(pdf_c), xc_max]` through
[`cloud_mean_droplet_mass_and_number`](@ref). `xc_max` = 2.6e-10 kg is, at liquid-water density,
a sphere of 79 μm, the raindrop minimum and the diameter at which
[`CMP.StokesRegimeVelType`](@ref)'s creeping-flow fall-speed law stops being valid.

# Arguments
 - `pdf_c`: Size distribution parameters for cloud droplets, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: Liquid mass content [kg/kg]
 - `ρₐ`: Air density [kg/m³]
 - `N`: Number concentration of the particle [1/m³]

# Returns
 - `(logA, logB)`: Log of the parameters of the generalized gamma distribution
"""
function log_pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)
    FT = UT.promote_typeof(q, ρₐ, N)
    (; νc, μc, loggamma_z1, loggamma_z2) = pdf_c
    x̄, N_eff = cloud_mean_droplet_mass_and_number(pdf_c, q, ρₐ, N)
    logx̄ = log(x̄)
    z1 = (νc + 1) / μc
    # loggamma_z1 = SF.loggamma(z1) (pre-computed in pdf_c)
    # loggamma_z2 = SF.loggamma(z2) (pre-computed in pdf_c)
    logB = -μc * (logx̄ + loggamma_z1 - loggamma_z2)
    logA = log(μc) + log(N_eff) + z1 * logB - loggamma_z1

    cond = N < UT.ϵ_numerics_2M_N(FT)
    return (ifelse(cond, oftype(logA, -Inf), logA), ifelse(cond, oftype(logB, Inf), logB))
end

"""
    pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)

Return the parameters of the size distribution of cloud particles in terms of mass.

See [`log_pdf_cloud_parameters_mass`](@ref) for more details.
"""
function pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)
    logA, logB = log_pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)
    return (; Ac = exp(logA), Bc = exp(logB))
end

"""
    pdf_cloud_parameters(pdf_c, q, ρₐ, N)

Return the parameters of the size distribution of cloud particles in terms of diameter.

The size distribution is given by:

    n(D) = N₀c * D^νcD * exp(-λc * D^μcD)

where
- `νcD = 3νc + 2`
- `μcD = 3μc`

# Arguments
 - `pdf_c`: Size distribution parameters for cloud droplets, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: Liquid mass content [kg/kg]
 - `ρₐ`: Air density [kg/m³]
 - `N`: Number concentration of the particle [1/m³]

# Returns
 - A `NamedTuple` with the fields `(; logN₀c, λc, νcD, μcD)`:
   Parameters of the generalized gamma distribution in terms of diameter
"""
function pdf_cloud_parameters(pdf_c, q, ρₐ, N)
    FT = UT.promote_typeof(q, ρₐ, N)
    logAc, logBc = log_pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)
    (; νc, μc, ρw) = pdf_c
    k_m = ρw * π / 6
    # Convert from mass-based to diameter-based distribution
    logN₀c = logAc + log(FT(3)) + (νc + 1) * log(k_m)
    λc = exp(logBc) * k_m^μc
    return (; logN₀c, λc, νcD = 3νc + 2, μcD = 3μc)
end

"""
    cloud_condensation_timescale(pdf_c, aps, tps, Tₐ, ρₐ, q_lcl, N_lcl)

Compute the condensation relaxation timescale of the cloud droplet population
from its capacitance integral,

```math
τ_{cond} = \\frac{ρₐ q_{v,sl}}{2π G_l ∫ D n(D) dD},
```

with spherical capacitance `C = D/2` and unit ventilation. The diameter moment
of the generalized gamma distribution is closed form,
`∫ D n(D) dD = N₀/μ λ^{-(ν+2)/μ} Γ((ν+2)/μ)`. The timescale diverges as the
population vanishes and shrinks as the integrated droplet surface grows.

# Arguments
 - `pdf_c`: cloud droplet size distribution parameters, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `aps`: [`CMP.AirProperties`](@ref)
 - `tps`: thermodynamics parameters
 - `Tₐ`: temperature (K)
 - `ρₐ`: air density
 - `q_lcl`: cloud liquid mass content [kg/kg]
 - `N_lcl`: cloud droplet number concentration [1/m³]

# Returns
- Condensation timescale [s], bounded above at [`CLOUD_COND_TIMESCALE_MAX`](@ref). A value at
  the bound means the capacitance integral underflowed and there is no population to relax; test
  with [`cloud_condensation_is_degenerate`](@ref) rather than comparing to a literal.
"""
@inline function cloud_condensation_timescale(
    pdf_c::CMP.CloudParticlePDF_SB2006, aps::CMP.AirProperties, tps::TDI.PS,
    Tₐ, ρₐ, q_lcl, N_lcl,
)
    FT = UT.promote_typeof(q_lcl, ρₐ, N_lcl, Tₐ)
    G = CO.G_func_liquid(aps, tps, Tₐ)
    qᵥ_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, Tₐ, ρₐ)
    (; logN₀c, λc, νcD, μcD) = pdf_cloud_parameters(pdf_c, q_lcl, ρₐ, N_lcl)
    z = (νcD + 2) / μcD
    log_moment = logN₀c - log(μcD) - z * log(λc) + SF.loggamma(z)
    denom = 2 * FT(π) * G * exp(log_moment)
    return min(ρₐ * qᵥ_sat_liq / max(denom, floatmin(FT)), CLOUD_COND_TIMESCALE_MAX(FT))
end

"""
    CLOUD_COND_TIMESCALE_MAX(FT)

Upper bound on the cloud condensation/evaporation relaxation timescale [s], applied by
[`cloud_condensation_timescale`](@ref) as the droplet population vanishes. A returned value at
this bound signals a degenerate population, not a slow relaxation; use
[`cloud_condensation_is_degenerate`](@ref) to detect it.
"""
@inline CLOUD_COND_TIMESCALE_MAX(::Type{FT}) where {FT} = FT(1e10)

"""
    cloud_condensation_is_degenerate(τ_cond)

`true` when `τ_cond` from [`cloud_condensation_timescale`](@ref) sits at
[`CLOUD_COND_TIMESCALE_MAX`](@ref): the capacitance integral underflowed and the
condensation/evaporation rate is exactly zero.
"""
@inline cloud_condensation_is_degenerate(τ_cond::FT) where {FT} =
    τ_cond >= CLOUD_COND_TIMESCALE_MAX(FT)

"""
    log_size_distribution_mass(pdf::CMP.CloudParticlePDF_SB2006, q_c, ρₐ, N_c)

Return the log of the size distribution, as a function of mass, of the form

    f(x) = A * x^ν * exp(-B * x^μ)

that is, the function

    log(f(x)) = log(A) + ν * log(x) - B * x^μ

"""
function log_size_distribution_mass(pdf::CMP.CloudParticlePDF_SB2006, q_c, ρₐ, N_c)
    (; νc, μc) = pdf
    logA, logB = log_pdf_cloud_parameters_mass(pdf, q_c, ρₐ, N_c)
    B = exp(logB)
    logpsd(x) = logA + νc * log(x) - B * x^μc
    return logpsd
end
size_distribution_mass(pdf, q_c, ρₐ, N_c) = exp ∘ log_size_distribution_mass(pdf, q_c, ρₐ, N_c)

"""
    size_distribution(pdf::CMP.RainParticlePDF_SB2006, q, ρₐ, N)

Return `n(D)`, a function that computes the size distribution for rain particles at diameter `D`

# Arguments
- `pdf`: Rain size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
- `q`: Rain water specific content [kg/kg]
- `ρₐ`: Density of air [kg/m³]
- `N`: Rain water number concentration [1/m³]
"""
function DT.size_distribution(pdf::CMP.RainParticlePDF_SB2006, q, ρₐ, N)
    (; N₀r, Dr_mean) = pdf_rain_parameters(pdf, q, ρₐ, N)
    function n(D)
        v = N₀r * exp(-D / Dr_mean)
        return ifelse(iszero(N₀r), zero(v), v)
    end
    return n
end

"""
    size_distribution(pdf::CMP.CloudParticlePDF_SB2006, q, ρₐ, N)

Return `n(D)`, a function that computes the size distribution for cloud particles at diameter `D`

The size distribution is given by:

    n(D) = N₀c * D^(3νc + 2) * exp(-λc * D^(3μc))

# Arguments
 - `pdf`: Cloud size distribution parameters, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: Cloud water specific content [kg/kg]
 - `ρₐ`: Density of air [kg/m³]
 - `N`: Cloud water number concentration [1/m³]

"""
function DT.size_distribution(pdf::CMP.CloudParticlePDF_SB2006, q, ρₐ, N)
    (; logN₀c, λc, νcD, μcD) = pdf_cloud_parameters(pdf, q, ρₐ, N)
    # zero the active value (not logN₀c) so the closure is type-concrete under mixed Dual/plain D
    function n(D)
        v = exp(logN₀c + νcD * log(D) - λc * D^μcD)
        return ifelse(logN₀c == -Inf, zero(v), v)
    end
    return n
end

"""
    size_distribution_value(pdf, q, ρₐ, N, D)

Return the size distribution value for a cloud or rain particle of diameter `D`.

See [`size_distribution`](@ref) for more details.
"""
function size_distribution_value(pdf, q, ρₐ, N, D)
    n = size_distribution(pdf, q, ρₐ, N)
    return n(D)
end

"""
    get_size_distribution_bounds(pdf, q, ρₐ, N, p = eps(eltype(q)))

Return the minimum and maximum diameters of a cloud or rain particle,
set at the `p`-th and `(1 - p)`-th quantiles of the size distribution.

# Arguments
 - `pdf`: Size distribution parameters for cloud or rain,
    [`CMP.RainParticlePDF_SB2006`](@ref) or [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: specific content of cloud or rain water [kg/kg]
 - `ρₐ`: density of air
 - `N`: number concentration of cloud or rain drops [1/m³]
 - `p`: probability level (0 ≤ p ≤ 1), default is `eps(eltype(q))`

# Returns
 - `D_min, D_max`: minimum and maximum diameters of a cloud or rain particle,
    at the `p`-th and `(1 - p)`-th quantiles of the size distribution.
    All inputs and output diameters are in base SI units.
    The bounds are calculated through quantile functions of the size distribution.
"""
function get_size_distribution_bounds(
    pdf::CMP.RainParticlePDF_SB2006, q, ρₐ, N, p = eps(eltype(q)),
)
    FT = UT.promote_typeof(q, ρₐ, N, p)
    (; Dr_mean) = pdf_rain_parameters(pdf, q, ρₐ, N)
    (isfinite(Dr_mean) && Dr_mean > 0) || return (FT(0), FT(0))
    D_min = DT.exponential_quantile(Dr_mean, p)
    D_max = DT.exponential_quantile(Dr_mean, 1 - p)
    return D_min, D_max
end
function get_size_distribution_bounds(
    pdf::CMP.CloudParticlePDF_SB2006, q, ρₐ, N, p = eps(eltype(q)),
)
    FT = UT.promote_typeof(q, ρₐ, N, p)
    (; λc, νcD, μcD) = pdf_cloud_parameters(pdf, q, ρₐ, N)
    (isfinite(λc) && λc > 0 && μcD > 0) || return (FT(0), FT(0))
    D_min = FT(DT.generalized_gamma_quantile(νcD, μcD, λc, p))
    D_max = FT(DT.generalized_gamma_quantile(νcD, μcD, λc, 1 - p))
    return D_min, D_max
end


### ----- ###
### RATES ###
### ----- ###

"""
A structure containing the rates of change of the specific contents and number
densities of cloud liquid water and rain water.
"""
@kwdef struct LclRaiRates{FT}
    "Rate of change of the cloud liquid water specific content"
    dq_lcl_dt::FT = FT(0)
    "Rate of change of the cloud liquid water number density"
    dN_lcl_dt::FT = FT(0)
    "Rate of change of the rain water specific content"
    dq_rai_dt::FT = FT(0)
    "Rate of change of the rain water number density"
    dN_rai_dt::FT = FT(0)
end
LclRaiRates(dq_lcl_dt, dN_lcl_dt, dq_rai_dt, dN_rai_dt) =
    LclRaiRates(promote(dq_lcl_dt, dN_lcl_dt, dq_rai_dt, dN_rai_dt)...)

"""
    mean_mass_bound_factor(x, x_max; onset = 1 // 2)

Factor in `[0, 1]` for a mean-particle-mass-dependent rate: `1` for
`x ≤ onset * x_max`, decreasing smoothly (continuous value and slope) to `0`
as `x` increases from `onset * x_max` to `x_max`, and identically `0` for
`x ≥ x_max`.
"""
@inline function mean_mass_bound_factor(x, x_max; onset = 1 // 2)
    FT = UT.promote_typeof(x, x_max)
    r = x / x_max
    s = clamp((r - FT(onset)) / (1 - FT(onset)), zero(FT), one(FT))
    return 1 - s^2 * (3 - 2 * s)
end

"""
    autoconversion(acnv, pdf_c, q_lcl, q_rai, ρ, N_lcl)

Compute autoconversion rates

# Arguments
 - `acnv`: Autoconversion parameters, [`CMP.AcnvSB2006`](@ref)
 - `pdf_c`: Cloud size distribution parameters, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q_lcl`: Cloud liquid water specific content [kg/kg]
 - `q_rai`: Rain water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `N_lcl`: Cloud droplet number density [1/m³]

# Returns
 - [`LclRaiRates`](@ref) with `q_lcl`, `N_lcl`, `q_rai`, `N_rai` tendencies due to
    collisions between cloud droplets (autoconversion)
"""
function autoconversion(
    acnv::CMP.AcnvSB2006, pdf_c::CMP.CloudParticlePDF_SB2006, q_lcl, q_rai, ρ, N_lcl,
)
    FT = UT.promote_typeof(q_lcl, q_rai, ρ, N_lcl)
    (; kcc, x_star, ρ0, A, a, b) = acnv
    (; νc) = pdf_c

    safe_q_lcl = max(q_lcl, UT.ϵ_numerics_2M_M(FT))
    safe_N_lcl = max(N_lcl, UT.ϵ_numerics_2M_N(FT))
    L_lcl = ρ * safe_q_lcl
    x_lcl = min(x_star, L_lcl / safe_N_lcl)
    bound_factor = mean_mass_bound_factor(L_lcl / safe_N_lcl, x_star)
    safe_q_rai = max(0, q_rai)
    τ = 1 - safe_q_lcl / (safe_q_lcl + safe_q_rai)  # Eq. (5) from SB2006
    # τ^a has a vertical tangent at τ = 0; the ifelse keeps the ForwardDiff
    # derivative w.r.t. q_rai finite at q_rai = 0 (and the code branch-free)
    ϕ_au = ifelse(q_rai < UT.ϵ_numerics_2M_M(FT), zero(τ), A * τ^a * (1 - τ^a)^b)

    # Eq. (4) from SB2006, scaled by `bound_factor` so the whole event rate
    # (mass and number together) vanishes continuously as the mean droplet
    # mass approaches `x_star` from below, instead of saturating at a fixed
    # value once `x_lcl` reaches the `min(x_star, ...)` clamp above.
    dL_rai_dt =
        kcc / 20 / x_star * (νc + 2) * (νc + 4) / (νc + 1)^2 *
        L_lcl^2 * x_lcl^2 * (1 + ϕ_au / (1 - τ)^2) * ρ0 / ρ * bound_factor
    dN_rai_dt = dL_rai_dt / x_star
    dL_lcl_dt = -dL_rai_dt
    dN_lcl_dt = -2 * dN_rai_dt

    cond = q_lcl < UT.ϵ_numerics_2M_M(FT) || N_lcl < UT.ϵ_numerics_2M_N(FT)
    return LclRaiRates(
        dq_lcl_dt = ifelse(cond, zero(FT), dL_lcl_dt / ρ),
        dN_lcl_dt = ifelse(cond, zero(FT), dN_lcl_dt),
        dq_rai_dt = ifelse(cond, zero(FT), dL_rai_dt / ρ),
        dN_rai_dt = ifelse(cond, zero(FT), dN_rai_dt),
    )
end

"""
    accretion(scheme::CMP.SB2006, q_lcl, q_rai, ρ, N_lcl)

Compute accretion rate

# Arguments
 - `scheme`: [`CMP.SB2006`](@ref) 2-moment scheme parameters
   (the accretion parameters [`CMP.AccrSB2006`](@ref) are taken from its `accr` field)
 - `q_lcl`: Cloud liquid water specific content [kg/kg]
 - `q_rai`: Rain water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `N_lcl`: Cloud droplet number density [1/m³]

# Returns
 - [`LclRaiRates`](@ref) with `q_lcl`, `N_lcl`, `q_rai`, `N_rai` tendencies due to
    collisions between raindrops and cloud droplets (accretion)
"""
function accretion((; accr)::CMP.SB2006, q_lcl, q_rai, ρ, N_lcl)
    FT = UT.promote_typeof(q_lcl, q_rai, ρ, N_lcl)
    (; kcr, τ0, ρ0, c) = accr
    safe_q_lcl = max(q_lcl, UT.ϵ_numerics_2M_M(FT))
    safe_q_rai = max(q_rai, UT.ϵ_numerics_2M_M(FT))
    safe_N_lcl = max(N_lcl, UT.ϵ_numerics_2M_N(FT))

    L_lcl = ρ * safe_q_lcl
    L_rai = ρ * safe_q_rai
    x_lcl = L_lcl / safe_N_lcl
    τ = 1 - safe_q_lcl / (safe_q_lcl + safe_q_rai)  # Eq. (5) from SB2006
    ϕ_ac = (τ / (τ + τ0))^c          # Eq. (8) from SB2006

    dL_rai_dt = kcr * L_lcl * L_rai * ϕ_ac * sqrt(ρ0 / ρ)  # Eq. (7) from SB2006
    dN_rai_dt = zero(FT)
    dL_lcl_dt = -dL_rai_dt
    dN_lcl_dt = dL_lcl_dt / x_lcl

    cond = q_lcl < UT.ϵ_numerics_2M_M(FT) || q_rai < UT.ϵ_numerics_2M_M(FT) || N_lcl < UT.ϵ_numerics_2M_N(FT)
    return LclRaiRates(
        dq_lcl_dt = ifelse(cond, zero(FT), dL_lcl_dt / ρ),
        dN_lcl_dt = ifelse(cond, zero(FT), dN_lcl_dt),
        dq_rai_dt = ifelse(cond, zero(FT), dL_rai_dt / ρ),
        dN_rai_dt = ifelse(cond, zero(FT), dN_rai_dt),
    )
end

"""
    cloud_liquid_self_collection(acnv, pdf_c, q_lcl, ρ, N_lcl, dN_lcl_dt_au)

Compute cloud liquid self-collection rate

# Arguments
 - `acnv`: 2-moment autoconversion parameterization, [`CMP.AcnvSB2006`](@ref)
 - `pdf_c`: Cloud size distribution parameters, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q_lcl`: Cloud liquid water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `N_lcl`: Cloud droplet number density [1/m³]
 - `dN_lcl_dt_au`: Rate of change of cloud droplets number density due to autoconversion [1/m³/s]

# Returns
 - The cloud droplets number density tendency due to collisions of cloud droplets
    that produce larger cloud droplets (self-collection)
"""
function cloud_liquid_self_collection(
    acnv::CMP.AcnvSB2006, pdf_c::CMP.CloudParticlePDF_SB2006, q_lcl, ρ, N_lcl, dN_lcl_dt_au,
)
    FT = UT.promote_typeof(q_lcl, ρ, N_lcl, dN_lcl_dt_au)
    (; kcc, ρ0, x_star) = acnv
    (; νc) = pdf_c

    L_lcl = ρ * q_lcl
    safe_N_lcl = max(N_lcl, UT.ϵ_numerics_2M_N(FT))
    bound_factor = mean_mass_bound_factor(L_lcl / safe_N_lcl, x_star)
    # Eq. (9) from SB2006, scaled by `bound_factor` so the sink vanishes
    # continuously as the mean droplet mass approaches `x_star` from below.
    dN_lcl_dt_sc = -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * L_lcl^2 * bound_factor - dN_lcl_dt_au

    cond = q_lcl < UT.ϵ_numerics_2M_M(FT)
    return ifelse(cond, FT(0), dN_lcl_dt_sc)
end

"""
    autoconversion_and_cloud_liquid_self_collection(scheme, q_lcl, q_rai, ρ, N_lcl)

Compute autoconversion and cloud liquid self-collection rates

# Arguments
 - `scheme`: 2-moment rain autoconversion parameterization, [`CMP.SB2006`](@ref)
 - `q_lcl`: Cloud liquid water specific content [kg/kg]
 - `q_rai`: Rain water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `N_lcl`: Cloud droplet number density [1/m³]

# Returns
 - `(au, sc)`: A `NamedTuple` containing the autoconversion rate and the
    cloud liquid self-collection rate.
"""
function autoconversion_and_cloud_liquid_self_collection(
    (; acnv, pdf_c)::CMP.SB2006, q_lcl, q_rai, ρ, N_lcl,
)

    au = autoconversion(acnv, pdf_c, q_lcl, q_rai, ρ, N_lcl)
    sc = cloud_liquid_self_collection(acnv, pdf_c, q_lcl, ρ, N_lcl, au.dN_lcl_dt)

    return (; au, sc)
end

"""
    rain_self_collection(pdf, self, q_rai, ρ, N_rai)

Compute the rain self-collection rate

# Arguments
 - `pdf`: Rain size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
 - `self`: Rain self-collection parameters, [`CMP.SelfColSB2006`](@ref)
 - `q_rai`: Rain water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `N_rai`: Raindrops number density [1/m³]

# Returns
 - The raindrops number density tendency due to collisions of raindrops that
    produce larger raindrops (self-collection).
"""
function rain_self_collection(
    pdf::CMP.RainParticlePDF_SB2006, self::CMP.SelfColSB2006, q_rai, ρ, N_rai,
)
    FT = UT.promote_typeof(q_rai, ρ, N_rai)
    (; krr, κrr, d) = self
    (; ρ0) = pdf

    safe_q_rai = max(q_rai, UT.ϵ_numerics_2M_M(FT))
    safe_N_rai = max(N_rai, UT.ϵ_numerics_2M_N(FT))
    L_rai = ρ * safe_q_rai
    (; Br) = pdf_rain_parameters_mass(pdf, safe_q_rai, ρ, safe_N_rai)
    dN_rai_dt_sc = -krr * N_rai * L_rai * √(ρ0 / ρ) * (1 + κrr / Br)^d  # Eq. (11) from SB2006

    cond = q_rai < UT.ϵ_numerics_2M_M(FT) || N_rai < UT.ϵ_numerics_2M_N(FT)
    return ifelse(cond, FT(0), dN_rai_dt_sc)
end

"""
    rain_breakup(pdf, brek, q_rai, ρ, N_rai, dN_rai_dt_sc)

Compute the raindrops number density tendency due to breakup of raindrops

# Arguments
 - `pdf`: Rain size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
 - `brek`: Rain breakup parameters, [`CMP.BreakupSB2006`](@ref)
 - `q_rai`: Rain water specific content
 - `ρ`: Air density
 - `N_rai`: Raindrops number density
 - `dN_rai_dt_sc`: Rate of change of raindrops number density due to self-collection

# Returns
 - The raindrops number density tendency due to breakup of raindrops that produce
    smaller raindrops
"""
function rain_breakup(
    pdf::CMP.RainParticlePDF_SB2006, brek::CMP.BreakupSB2006, q_rai, ρ, N_rai, dN_rai_dt_sc,
)
    FT = UT.promote_typeof(q_rai, ρ, N_rai, dN_rai_dt_sc)
    (; Deq, Dr_th, kbr, κbr) = brek
    (; ρw) = pdf

    safe_q_rai = max(q_rai, UT.ϵ_numerics_2M_M(FT))
    safe_N_rai = max(N_rai, UT.ϵ_numerics_2M_N(FT))

    (; xr_mean) = pdf_rain_parameters(pdf, safe_q_rai, ρ, safe_N_rai)
    Dr = cbrt(xr_mean * 6 / (π * ρw))  # mean volume raindrop diameter
    ΔD = Dr - Deq

    # Dr < Dr_th:  below the threshold diameter, breakup is neglected
    # Dr ≤ Deq:    below the equilibrium diameter, breakup is a linear function
    # Dr > Deq:    above the equilibrium diameter, breakup is an exponential function
    # TODO: the exponential branch disagrees with the docs. SB2006 Eq. (13)
    # prints 2 exp(κbr ΔD) - 1, which is discontinuous at ΔD = 0;
    # docs/src/Microphysics2M.md reads it as 2 (exp(κbr ΔD) - 1), while the
    # code uses exp(κbr ΔD) - 1. Both amended forms are continuous but
    # differ by a factor of 2. Decide which is intended and align the docs.
    Φ_br = ifelse(Dr < Dr_th, FT(-1), ifelse(Dr ≤ Deq, kbr * ΔD, exp(κbr * ΔD) - 1))
    dN_rai_dt_br = -(Φ_br + 1) * dN_rai_dt_sc  # Eq. (13) from SB2006

    cond = q_rai < UT.ϵ_numerics_2M_M(FT) || N_rai < UT.ϵ_numerics_2M_N(FT)
    return ifelse(cond, FT(0), dN_rai_dt_br)
end

"""
    rain_self_collection_and_breakup(params, q_rai, ρ, N_rai)

Compute the raindrops self-collection and breakup rates.

# Arguments
 - `params`: 2-moment rain size distribution parameters, [`CMP.SB2006`](@ref)
    including raindrop size distribution, self collection, and breakup parameters
 - `q_rai`: Rain water specific content
 - `ρ`: Air density
 - `N_rai`: Raindrops number density

# Returns
- `(sc, br)`: A `NamedTuple` containing the raindrops self-collection and breakup rates, respectively.
"""
function rain_self_collection_and_breakup(
    (; pdf_r, self, brek)::CMP.SB2006, q_rai, ρ, N_rai,
)

    sc = rain_self_collection(pdf_r, self, q_rai, ρ, N_rai)
    br = rain_breakup(pdf_r, brek, q_rai, ρ, N_rai, sc)

    return (; sc, br)
end

"""
    rain_equilibrium_number(brek, pdf_r, q_rai, ρ)

The rain number density `N_eq` [1/m³] at which self-collection and breakup balance:
`L_rai / x_eq`, with `x_eq = (π/6) ρw Deq³` the mass of a drop at the SB2006 fitted
collisional-equilibrium mean-volume diameter `Deq`.

`N_eq` is defined for every `L_rai ≥ 0` and is exactly zero on an empty state, so a relaxation
toward it cannot manufacture number without mass and needs no special empty-state arm.
"""
@inline function rain_equilibrium_number(
    brek::CMP.BreakupSB2006, pdf_r::CMP.RainParticlePDF_SB2006, q_rai, ρ,
)
    FT = UT.promote_typeof(q_rai, ρ)
    return ρ * q_rai / _rain_equilibrium_mass(FT, brek, pdf_r)
end

"""Mass of a drop at the collisional-equilibrium mean-volume diameter, `x_eq = (π/6) ρw Deq³`."""
@inline _rain_equilibrium_mass(::Type{FT}, brek, pdf_r) where {FT} =
    FT(π) / 6 * FT(pdf_r.ρw) * FT(brek.Deq)^3

"""
    rain_number_relaxation(pdf_r, self, brek, q_rai, ρ, N_rai)

Write the SB2006 self-collection/breakup pair as a relaxation toward the equilibrium number
density,

    ∂ₜN_rai = -(N_rai - N_eq(L)) / τ_eff,    τ_eff = (N_eq - N_rai) / f,

where `f = sc + br` is the pair's net number tendency and `N_eq` is
[`rain_equilibrium_number`](@ref). `∂ₜN_rai` equals `sc + br` exactly; `1/τ_eff` supplies the
Jacobian's diagonal entry `∂(∂ₜN_rai)/∂N_rai = -1/τ_eff`.

At `N_rai = N_eq` the quotient is `0/0`; the removable limit is the linearization of the pair at
equilibrium,

    1/τ_relax(L) = (κ Deq / 3) k_rr √(ρ0/ρ) (1 + κrr/Br(x_eq))^d L,

with `κ = κbr` above `Deq` and `κ = kbr` below, matching the breakup fit's two branches (C¹ kink
at `Deq`).

Requires a rain PSD whose parameters are honest functions of `(L, N)`, i.e. `λ` and `N₀` not
clamped independently of the state's own moments; both concrete `RainParticlePDF_SB2006`
variants satisfy this by construction.

# Arguments
 - `pdf_r`: rain size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
 - `self`: rain self-collection parameters, [`CMP.SelfColSB2006`](@ref)
 - `brek`: rain breakup parameters, [`CMP.BreakupSB2006`](@ref)
 - `q_rai`: rain water specific content [kg/kg]
 - `ρ`: air density [kg/m³]
 - `N_rai`: raindrop number density [1/m³]

# Returns
 - `(; ∂ₜN_rai, N_eq, inv_τ_eff)`, the pair's net number tendency [1/(m³ s)], the equilibrium
   number density [1/m³], and the relaxation rate `1/τ_eff` [1/s].
"""
@inline function rain_number_relaxation(
    pdf_r::CMP.RainParticlePDF_SB2006, self::CMP.SelfColSB2006, brek::CMP.BreakupSB2006,
    q_rai, ρ, N_rai,
)
    FT = UT.promote_typeof(q_rai, ρ, N_rai)
    (; krr, κrr, d) = self
    (; Deq, kbr, κbr) = brek
    (; ρ0) = pdf_r

    sc = rain_self_collection(pdf_r, self, q_rai, ρ, N_rai)
    br = rain_breakup(pdf_r, brek, q_rai, ρ, N_rai, sc)
    ∂ₜN_rai = sc + br

    N_eq = rain_equilibrium_number(brek, pdf_r, q_rai, ρ)
    Δ = N_eq - N_rai

    Br_eq = cbrt(6 / _rain_equilibrium_mass(FT, brek, pdf_r))
    κ = ifelse(Δ ≥ 0, κbr, kbr)
    inv_τ_lin =
        κ * Deq / 3 * krr * sqrt(ρ0 / ρ) * (1 + κrr / Br_eq)^d * (ρ * q_rai)

    near_eq = abs(Δ) ≤ sqrt(eps(FT)) * N_eq
    inv_τ_eff = ifelse(near_eq, inv_τ_lin, ∂ₜN_rai / ifelse(near_eq, one(FT), Δ))
    # Nonnegative for an honest PSD (see the docstring); floored for robustness against any
    # future `RainParticlePDF_SB2006` variant. `∂ₜN_rai` is unaffected either way.
    inv_τ_eff = max(zero(FT), inv_τ_eff)

    inv_τ_eff = ifelse(
        q_rai < UT.ϵ_numerics_2M_M(FT) || N_rai < UT.ϵ_numerics_2M_N(FT),
        zero(FT), inv_τ_eff,
    )
    return (; ∂ₜN_rai, N_eq, inv_τ_eff)
end

"""
    cloud_terminal_velocity(pdf_c, vel_params, q_liq, ρₐ, N_liq)

Compute the number-averaged and mass-averaged terminal velocities of cloud droplets
assuming a gamma size distribution for droplet mass and the analytical Stokes-regime terminal
velocity of spherical particles.

# Arguments
- `pdf_c`: Cloud droplet size distribution parameters, [`CMP.CloudParticlePDF_SB2006`](@ref).
- `vel_params`: Terminal velocity parameters, [`CMP.StokesRegimeVelType`](@ref).
- `q_liq`: Cloud liquid water specific content [kg kg⁻¹].
- `ρₐ`: Air density [kg m⁻³].
- `N_liq`: Cloud droplet number concentration [m⁻³].

# Returns
A tuple containing the number- and mass-weighted mean fall velocities of cloud droplets in [m/s].
Individual droplet terminal velocities follow v_{term}(D) = (1/18) (ρw - ρₐ) g D^2 / μ_air with
μ_air = ρₐ * ν_air and assuming constant ν_air.
"""
function cloud_terminal_velocity(
    pdf_c::CMP.CloudParticlePDF_SB2006,
    (; ρw, grav, ν_air)::CMP.StokesRegimeVelType,
    q_liq, ρₐ, N_liq,
)
    FT = UT.promote_typeof(q_liq, ρₐ, N_liq)
    (; νc, μc) = pdf_c
    safe_q_liq = max(q_liq, UT.ϵ_numerics_2M_M(FT))
    safe_N_liq = max(N_liq, UT.ϵ_numerics_2M_N(FT))
    (; Bc) = pdf_cloud_parameters_mass(pdf_c, safe_q_liq, ρₐ, safe_N_liq)
    # BOTH halves of the canonicalized pair, because the distribution is built from both.
    # `log_pdf_cloud_parameters_mass` forms `logB` from `x̄` and `logA` from `N_eff`, so the
    # distribution these moments integrate carries a mass concentration of `x̄ * N_eff` and NOT
    # `ρₐ q`. The two agree wherever the upper bound binds, since there `N_eff = L / xc_max` and
    # `x̄ = xc_max`. They do not agree where the LOWER bound binds: `x̄` is raised to the
    # activation droplet mass while the number is left alone, so `x̄ * N_eff` exceeds `ρₐ q` and
    # normalising the mass-weighted moment by the latter divides a bounded numerator by an
    # unbounded denominator. Measured on a state with 1.66e-6 kg/kg of cloud liquid, that returned
    # a fall speed of 21884 m/s where rain was 6.10 and ice 2.08.
    x̄, N_eff = cloud_mean_droplet_mass_and_number(pdf_c, safe_q_liq, ρₐ, safe_N_liq)
    L_eff = x̄ * N_eff

    terminal_velocity_prefactor = FT(1 / 18) * cbrt((FT(6) / ρw / FT(π))^2) * (ρw / ρₐ - 1) * grav / ν_air
    vt0 = terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_eff, FT(2 / 3)) / N_eff
    vt1 = terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_eff, FT(5 / 3)) / L_eff

    cond = N_liq < UT.ϵ_numerics_2M_N(FT) || q_liq < UT.ϵ_numerics_2M_M(FT)
    return (ifelse(cond, FT(0), vt0), ifelse(cond, FT(0), vt1))
end

"""
    rain_terminal_velocity(scheme::CMP.SB2006, vel, q_rai, ρ, N_rai)

Compute the raindrops terminal velocity.

# Arguments
 - `scheme`: [`CMP.SB2006`](@ref) 2-moment scheme parameters
   (the rain size distribution parameters are taken from its `pdf_r` field)
 - `vel`: Terminal velocity parameters,
   [`CMP.SB2006VelType`](@ref) or [`CMP.Chen2022VelTypeRain`](@ref)
 - `q_rai`: Rain water specific content
 - `ρ`: Air density
 - `N_rai`: Raindrops number density

# Returns
A tuple containing the number and mass weighted mean fall velocities of raindrops in [m/s],
assuming an exponential size distribution from Seifert and Beheng 2006.
The fall velocity of individual rain drops is parameterized:
 - assuming an empirical relation similar to Rogers (1993) for `vel::CMP.SB2006VelType`,
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171 for `vel::CMP.Chen2022VelTypeRain`.
"""
function rain_terminal_velocity(
    (; pdf_r)::CMP.SB2006, (; ρ0, aR, bR, cR)::CMP.SB2006VelType, q_rai, ρ, N_rai,
)
    FT = UT.promote_typeof(q_rai, ρ, N_rai)
    safe_q_rai = max(q_rai, UT.ϵ_numerics_2M_M(FT))
    safe_N_rai = max(N_rai, UT.ϵ_numerics_2M_N(FT))

    (; Dr_mean) = pdf_rain_parameters(pdf_r, safe_q_rai, ρ, safe_N_rai)
    _pa0, _pb0, _pa1, _pb1 =
        _sb_rain_terminal_velocity_helper(pdf_r, 1 / Dr_mean, aR, bR, cR)

    vt0 = max(0, sqrt(ρ0 / ρ) * (aR * _pa0 - bR * _pb0 / (1 + cR * Dr_mean)))
    vt1 = max(0, sqrt(ρ0 / ρ) * (aR * _pa1 - bR * _pb1 / (1 + cR * Dr_mean)^4))

    cond_N = N_rai < UT.ϵ_numerics_2M_N(FT)
    cond_q = q_rai < UT.ϵ_numerics_2M_M(FT)
    return (ifelse(cond_N, FT(0), vt0), ifelse(cond_q, FT(0), vt1))
end
function rain_terminal_velocity(
    (; pdf_r)::CMP.SB2006, vel::CMP.Chen2022VelTypeRain, q_rai, ρ, N_rai,
)
    FT = UT.promote_typeof(q_rai, ρ, N_rai)
    aiu, bi, ciu = CO.Chen2022_vel_coeffs(vel, ρ)
    safe_q_rai = max(q_rai, UT.ϵ_numerics_2M_M(FT))
    safe_N_rai = max(N_rai, UT.ϵ_numerics_2M_N(FT))
    (; Dr_mean) = pdf_rain_parameters(pdf_r, safe_q_rai, ρ, safe_N_rai)

    # It should be (ϕ^κ * vt0, ϕ^κ * vt3), but for rain drops ϕ = 1 and κ = 0
    vt0 = sum(map((a, b, c) -> CO.Chen2022_exponential_pdf(a, b, c, Dr_mean, 0), aiu, bi, ciu))
    vt3 = sum(map((a, b, c) -> CO.Chen2022_exponential_pdf(a, b, c, Dr_mean, 3), aiu, bi, ciu))

    cond_N = N_rai < UT.ϵ_numerics_2M_N(FT)
    cond_q = q_rai < UT.ϵ_numerics_2M_M(FT)
    return (ifelse(cond_N, FT(0), max(0, vt0)), ifelse(cond_q, FT(0), max(0, vt3)))
end
# The `SB2006VelType` moment factors. The individual-drop fit `v = aR - bR exp(-cR D)` is
# NEGATIVE below the diameter where it crosses zero, so the honest bulk moments integrate only
# over the range where it is positive; that is the second method below, shared by the windowed
# and unbounded variants. The cascade variant returns the UNTRUNCATED factors instead, which is
# not a variant of the limiting at all but a second, undocumented difference riding along with
# it.
#
# The truncated form is the implementation. The untruncated method exists only to keep the
# cascade reproducing what it has always produced, and it retires when the cascade does. It is
# not a supported alternative and nothing new should dispatch to it. Production sedimentation
# uses `Chen2022VelTypeRain`, so the stakes are documentation rather than results, but the
# asymmetry must not persist silently: any golden that moves on a `SB2006VelType` path when the
# default distribution changes then has TWO causes, the limiting and this truncation, not one.
#
# Dispatching the truncated form on the ABSTRACT `RainParticlePDF_SB2006` is what collapses
# them, because the cascade type is a subtype of it and is silently caught. The methods are
# therefore written against the concrete types.
function _sb_rain_terminal_velocity_helper(
    ::CMP.RainParticlePDF_SB2006_limited, λr, aR, bR, cR,
)
    FT = eltype(λr)
    return (FT(1), FT(1), FT(1), FT(1))
end
function _sb_rain_terminal_velocity_helper(
    ::Union{
        CMP.RainParticlePDF_SB2006_notlimited,
        CMP.RainParticlePDF_SB2006_windowed,
    },
    λr, aR, bR, cR,
)
    # Integrate velocity of particles over a range of r with
    # positive terminal velocity (v = aR - bR exp(-lambda D))
    _rc = -1 / (2 * cR) * log(aR / bR)
    _Γ_1(t) = exp(-t)
    _Γ_4(t) = (t^3 + 3 * t^2 + 6 * t + 6) * exp(-t)
    _pa0 = _Γ_1(2 * _rc * λr)
    _pb0 = _Γ_1(2 * _rc * (λr + cR))
    _pa1 = _Γ_4(2 * _rc * λr) / 6
    _pb1 = _Γ_4(2 * _rc * (λr + cR)) / 6
    return (_pa0, _pb0, _pa1, _pb1)
end

"""
    Γ_incl(a, x)

Returns the approximation of an incomplete gamma function for a ∈ {-1.0, -0.101}, and x in [0.067 1.82]
"""
function Γ_incl(a, x)
    FT = UT.promote_typeof(a, x)
    #return exp(-x) / ((FT(1.5) - FT(0.54) * a) * x^(FT(0.46) - FT(0.75) * a))
    return exp(-x) / (
        (FT(0.33) - FT(0.7) * a) * x^(FT(0.08) - FT(0.93) * a) +
        (FT(1.34) - FT(0.1) * a) * x^(FT(0.8) - a)
    )
end

"""
    rain_evaporation(scheme, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)

Compute the evaporation of raindrop number and mass.

# Arguments
  - `scheme`: precipitation formation scheme, [`CMP.SB2006`](@ref). Notably, need the fields:
    + `pdf_r`: Raindrop size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
    + `evap`: evaporation parameterization scheme, [`CMP.EvaporationSB2006`](@ref)
  - `aps`: air properties, [`CMP.AirProperties`](@ref)
  - `tps`: thermodynamics parameters, [`ThermodynamicsParameters`](@extref Thermodynamics.Parameters.ThermodynamicsParameters)
  - `q_tot`, `q_lcl`, `q_icl`, `q_rai`, `q_sno`: total water,
     cloud liquid water, cloud ice, rain and snow specific contents, [kg kg⁻¹]
  - `ρ`: air density [kg m⁻³]
  - `N_rai`: raindrops number density [m⁻³]
  - `T`: air temperature [K]

# Returns
  - A NamedTuple `(; ∂ₜρn_rai, ∂ₜq_rai)` with
    + `∂ₜρn_rai`: tendency of raindrops number density [m⁻³ s⁻¹]
    + `∂ₜq_rai`: tendency of rain water specific content [kg kg⁻¹ s⁻¹]

These are computed assuming a power law velocity relation for the
fall velocity of individual drops and an exponential drop size distribution.
"""
function rain_evaporation(
    (; pdf_r, evap)::CMP.SB2006, aps::CMP.AirProperties, tps::TDI.PS,
    q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T,
)
    # the early return below must match the main path's type for any mix of
    # plain-float and Dual arguments
    FT = UT.promote_typeof(q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)
    ϵₘ = UT.ϵ_numerics_2M_M(FT)
    ϵₙ = UT.ϵ_numerics_2M_N(FT)

    S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)

    (; ν_air, D_vapor) = aps
    (; α, β, ρ0) = evap
    ρw = pdf_r.ρw
    x_star = pdf_r.xr_min
    G = CO.G_func_liquid(aps, tps, T)

    safe_q_rai = max(q_rai, ϵₘ)
    safe_N_rai = max(N_rai, ϵₙ)

    (; xr_mean) = pdf_rain_parameters(pdf_r, safe_q_rai, ρ, safe_N_rai)
    Dr = cbrt(6 * xr_mean / (π * ρw))

    t_star = cbrt(FT(6) * x_star / xr_mean)
    a_vent_0 = evap.a_vent_0_coeff * Γ_incl(FT(-1), t_star)
    b_vent_0 = evap.b_vent_0_coeff * Γ_incl(evap.β_vent_0, t_star)

    a_vent_1 = evap.a_vent_1
    b_vent_1 = evap.b_vent_1

    N_Re = α * xr_mean^β * sqrt(ρ0 / ρ) * Dr / ν_air
    # cbrt(Schmidt) and sqrt(Reynolds) are identical for both ventilation moments, so compute once
    # (libm cbrt is not reliably CSE'd across the two uses).
    cbrt_Sc = cbrt(ν_air / max(D_vapor, UT.ϵ_numerics(FT)))
    sqrt_N_Re = sqrt(N_Re)
    Fv0 = a_vent_0 + b_vent_0 * cbrt_Sc * sqrt_N_Re
    Fv1 = a_vent_1 + b_vent_1 * cbrt_Sc * sqrt_N_Re

    ∂ₜρn_rai = min(zero(FT), 2 * FT(π) * G * S * N_rai * Dr * Fv0 / xr_mean)
    ∂ₜq_rai = min(zero(FT), 2 * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)

    # When xr = 0, ∂ₜρn_rai becomes NaN. We replace NaN with 0 which is the limit of
    # ∂ₜρn_rai for xr -> 0.
    ∂ₜρn_rai = ifelse(q_rai < ϵₘ || xr_mean / x_star < eps(FT) || N_rai ≤ ϵₙ || S ≥ 0, FT(0), ∂ₜρn_rai)
    ∂ₜq_rai = ifelse(q_rai < ϵₘ || N_rai ≤ ϵₙ || S ≥ 0, FT(0), ∂ₜq_rai)

    return (; ∂ₜρn_rai, ∂ₜq_rai)
end

"""
    ∂rain_evaporation_∂N_rai_∂q_rai(sb, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)

Returns the leading-order derivatives of the rain evaporation tendencies with
respect to rain specific content `q_rai` and rain number concentration N_rai.

Uses a donor-based leading-order approximation:
- ∂(∂ₜρn_rai/ρ)/∂N_rai ≈ ∂ₜρn_rai / N_rai  (number tendency, first)
- ∂(∂ₜq_rai)/∂q_rai ≈ ∂ₜq_rai / q_rai  (mass tendency, second)

# Returns
`NamedTuple` with fields `(; ∂N_rai, ∂q_rai)`.
"""
@inline function ∂rain_evaporation_∂N_rai_∂q_rai(
    sb::CMP.SB2006, aps::CMP.AirProperties, tps::TDI.PS,
    q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T,
)
    FT = eltype(q_tot)
    result = rain_evaporation(sb, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)
    ∂N_rai = ifelse(N_rai > UT.ϵ_numerics_2M_N(FT), result.∂ₜρn_rai / N_rai, zero(result.∂ₜρn_rai))
    ∂q_rai = ifelse(q_rai > UT.ϵ_numerics_2M_M(FT), result.∂ₜq_rai / q_rai, zero(result.∂ₜq_rai))
    return (; ∂N_rai, ∂q_rai)
end

"""
    number_tendency_from_mass_limits(params, q, n, sat_excess = 0)

Compute the specific number tendency (rate of change) to relax the mean
particle mass, `x = q / n` [kg], towards the physical bounds `[x_min, x_max]`
[kg].

The relaxation tendency is given by

    ∂n/∂t = (n_target - n) / τ

where `n_target` is the specific number that corresponds to the nearest
valid mean particle mass,

    n_target = q / clamp(x, x_min, x_max)

for `q > 0`. At `q ≤ 0` the mean mass does not exist and the target is decided by `sat_excess`
instead: `n_target = n` where the vapor is in excess, and `n_target = 0` where it is not. This
follows the number-concentration exchange with a background reservoir of cloud condensation
nuclei, `∂N_CCN/∂t = -∂N/∂t` (see the Number concentration adjustment section of the
`Microphysics2M` documentation, and [Horn2012](@cite)): draining number at zero mass returns
droplets to that reservoir, which happens only once the air can no longer sustain them
(subsaturation), not merely because a clamp zeroed their mass.

Only the sign of `sat_excess` is read, through `FD.value`, since the excess is state-dependent
and differentiating the branch itself would put a spurious derivative of a switch into the
Jacobian. The default `sat_excess = 0` drains at zero mass; the retention arm is currently
supplied for cloud droplets only.

# Arguments
  - `params`: Number concentration adjustment parameters, a `NamedTuple` with fields:
    + `x_min`: Minimum allowed mean particle mass [kg]
    + `x_max`: Maximum allowed mean particle mass [kg]
    + `τ`: Relaxation timescale [s]
  - `q`: Specific mass (mass mixing ratio) [kg/kg]
  - `n`: Specific number (number mixing ratio) [1/kg]
  - `sat_excess`: Vapor specific content in excess of saturation over the
    species' own phase [kg/kg]. Only its sign is used.

# Returns
- The rate of change of specific number [1/(kg·s)] needed to bring the mean mass within the valid bounds.
"""
function number_tendency_from_mass_limits(
    (; x_min, x_max, τ), q, n, sat_excess = zero(q); invent_from_zero = true,
)
    # `q > 0` is a presence test, not a smallness threshold: `q/x_max` is a positive number at
    # every positive mass, so a smallness threshold on that arm would relax `n` up from zero and
    # manufacture particles at the largest mass the window allows. The orphan mass at `q ≤ 0` is
    # drained separately, as mass, by `orphan_mass_drain` and `orphan_mass_drain_ice`.
    # `invent_from_zero` keeps the previous behavior (relax `n` up from zero at any positive mass)
    # for callers whose corner doctrine is undecided; every production caller passes `false`.
    FT = UT.promote_typeof(q, n)
    orphan = !invent_from_zero & !(FD.value(n) > zero(FT))
    n_target = ifelse(
        orphan,
        zero(FT),
        ifelse(
            q > zero(FT),
            clamp(FT(n), q / x_max, q / x_min),
            ifelse(FD.value(sat_excess) > 0, FT(n), zero(FT)),
        ),
    )
    return (n_target - n) / τ
end

"""
    orphan_mass_drain(aps, tps, Tₐ, ρₐ, q, x_min, ρ_w, sat_excess)

The mass tendency [kg/kg/s] that removes orphan condensate - mass whose number
concentration is absent, so no particle carries it - by evaporating it to vapor at the rate a
population of minimum-mass particles would evaporate at:

    ∂ₜq = -q / τ_orphan,     1/τ_orphan = 2π G D_min (q_sat - qᵥ) / (x_min q_sat)

`D_min` is the diameter of a particle of mass `x_min`. This is the condensation closure's own
`∂ₜq = (qᵥ - q_sat)/τ` evaluated at the number `N = ρₐ q / x_min` that makes the mean mass
exactly `x_min`; `N` cancels out of the result, so no invented population number reaches any
rate. Minimum mass gives the largest surface-to-mass ratio and thus the fastest possible
evaporation; ventilation is omitted, which can only make the drain slower than the true rate.

The drain is zero at or above saturation: droplet activation supplies a real number there, and
the orphan mass is adopted by the population that arrives within an activation timescale.
"""
@inline function orphan_mass_drain(
    aps::CMP.AirProperties, tps::TDI.PS, Tₐ, ρₐ, q, x_min, ρ_w, sat_excess,
)
    inv_τ = orphan_mass_inv_timescale(aps, tps, Tₐ, ρₐ, x_min, ρ_w, sat_excess)
    return -UT.clamp_to_nonneg(q) * inv_τ
end

"""
    orphan_mass_drain_ice(aps, tps, Tₐ, ρₐ, q, x_min, ρ_i, sat_excess)

The ice-phase [`orphan_mass_drain`](@ref): the same minimum-mass relaxation with the vapor
diffusivity and saturation taken over ice, so orphan ice mass sublimates to vapor at the rate
a population of nucleation-mass crystals would. `sat_excess` is the vapor excess over ice
saturation; at or above ice saturation the drain is zero and the mass is left in place until
a number source, transport, or drying air resolves it.
"""
@inline function orphan_mass_drain_ice(
    aps::CMP.AirProperties, tps::TDI.PS, Tₐ, ρₐ, q, x_min, ρ_i, sat_excess,
)
    inv_τ = orphan_mass_inv_timescale_ice(aps, tps, Tₐ, ρₐ, x_min, ρ_i, sat_excess)
    return -UT.clamp_to_nonneg(q) * inv_τ
end

"""
    orphan_mass_inv_timescale(aps, tps, Tₐ, ρₐ, x_min, ρ_w, sat_excess)

`1/τ_orphan` [1/s] of [`orphan_mass_drain`](@ref), which is the whole of that
rate's mass dependence-free part. Split out because the substep's manual Jacobian
carries `-1/τ_orphan` as an exact diagonal and has to read the same number the
primal rate was built from; a linearization recomputing its own copy is how f and
J drift apart. The ice phase shares the expression through
[`orphan_mass_inv_timescale_ice`](@ref), differing only in the diffusional growth
factor and the saturation reference.
"""
@inline function orphan_mass_inv_timescale(
    aps::CMP.AirProperties, tps::TDI.PS, Tₐ, ρₐ, x_min, ρ_w, sat_excess,
)
    G = CO.G_func_liquid(aps, tps, Tₐ)
    qᵥ_sat = TDI.saturation_vapor_specific_content_over_liquid(tps, Tₐ, ρₐ)
    return _orphan_mass_inv_timescale(G, qᵥ_sat, x_min, ρ_w, sat_excess)
end

"""
    orphan_mass_inv_timescale_ice(aps, tps, Tₐ, ρₐ, x_min, ρ_i, sat_excess)

`1/τ_orphan` [1/s] of [`orphan_mass_drain_ice`](@ref), split out for the manual
Jacobian exactly as [`orphan_mass_inv_timescale`](@ref) is for the warm phase.
"""
@inline function orphan_mass_inv_timescale_ice(
    aps::CMP.AirProperties, tps::TDI.PS, Tₐ, ρₐ, x_min, ρ_i, sat_excess,
)
    G = CO.G_func_ice(aps, tps, Tₐ)
    qᵥ_sat = TDI.saturation_vapor_specific_content_over_ice(tps, Tₐ, ρₐ)
    return _orphan_mass_inv_timescale(G, qᵥ_sat, x_min, ρ_i, sat_excess)
end

@inline function _orphan_mass_inv_timescale(G, qᵥ_sat, x_min, ρ_x, sat_excess)
    FT = UT.promote_typeof(G, qᵥ_sat, sat_excess)
    D_min = cbrt(6 * FT(x_min) / (FT(π) * FT(ρ_x)))
    subsaturation = max(-sat_excess, zero(FT)) / max(qᵥ_sat, floatmin(FT))
    return 2 * FT(π) * G * D_min * subsaturation / FT(x_min)
end

"""
    number_bounded_by_mass_limits((; x_min, x_max), q, n, sat_excess = 0;
        invent_from_zero = true)

Specific number bounded by the mean-particle-mass limits `[x_min, x_max]` [kg]:
`clamp(n, q / x_max, q / x_min)` for `q > 0`, and at `q ≤ 0` the number itself where
`sat_excess` is positive, zero where it is not. This is the `n_target` of
[`number_tendency_from_mass_limits`](@ref) exactly, so process rates evaluated at this number
are consistent with the adjusted mean mass:

    number_tendency_from_mass_limits(p, q, n, s) ==
        (number_bounded_by_mass_limits(p, q, n, s) - n) / τ
"""
function number_bounded_by_mass_limits(
    (; x_min, x_max), q, n, sat_excess = zero(q); invent_from_zero = true,
)
    FT = UT.promote_typeof(q, n)
    orphan = !invent_from_zero & !(FD.value(n) > zero(FT))
    return ifelse(
        orphan,
        zero(FT),
        ifelse(
            q > zero(FT),
            clamp(FT(n), q / x_max, q / x_min),
            ifelse(FD.value(sat_excess) > 0, FT(n), zero(FT)),
        ),
    )
end

# Additional double moment autoconversion and accretion parametrizations:
# - Khairoutdinov and Kogan (2000)
# - Beheng (1994)
# - Tripoli and Cotton (1980)
# - Liu and Daum (2004)

"""
    conv_q_lcl_to_q_rai(scheme, q_lcl, ρ, N_d, smooth_transition = false)

Compute the `q_rai` tendency due to collisions between cloud droplets
(autoconversion).

# Arguments
 - `scheme` - autoconversion scheme parameters
   (the autoconversion parameters are taken from its `acnv` field):
   - Khairoutdinov and Kogan (2000) for `scheme::CMP.KK2000`
   - Beheng (1994) for `scheme::CMP.B1994`
   - Tripoli and Cotton (1980) for `scheme::CMP.TC1980`
   - Liu and Daum (2004) for `scheme::CMP.LD2004`
 - `q_lcl` - cloud liquid water specific content
 - `ρ` - air density
 - `N_d` - prescribed cloud droplet number concentration
 - `smooth_transition` - for the `CMP.B1994`, `CMP.TC1980` and `CMP.LD2004`
   schemes, an optional flag that smoothes their threshold behavior if set
   to `true`. The default value is `false`.
"""
function conv_q_lcl_to_q_rai((; acnv)::CMP.KK2000, q_lcl, ρ, N_d)
    q_lcl = max(0, q_lcl)
    (; A, a, b, c) = acnv
    return A * q_lcl^a * N_d^b * ρ^c
end
function conv_q_lcl_to_q_rai((; acnv)::CMP.B1994, q_lcl, ρ, N_d, smooth_transition = false)
    q_lcl = max(0, q_lcl)
    (; C, a, b, c, N_0, k, d_low, d_high) = acnv
    d = zero(q_lcl)
    if smooth_transition
        d_low_acnv_fraction = CO.logistic_function(N_d, N_0, k)
        d_high_acnv_fraction = 1 - d_low_acnv_fraction
        d = d_low_acnv_fraction * d_low + d_high_acnv_fraction * d_high
    else
        d = N_d >= N_0 ? d_low : d_high
    end
    return C * d^a * (q_lcl * ρ)^b * N_d^c / ρ
end
function conv_q_lcl_to_q_rai((; acnv)::CMP.TC1980, q_lcl, ρ, N_d, smooth_transition = false)
    #TODO - The original paper is actually formulated for mixing ratios, not specific contents
    q_lcl = max(0, q_lcl)
    (; m0_liq_coeff, me_liq, D, a, b, r_0, k) = acnv
    q_liq_threshold = m0_liq_coeff * N_d / ρ * r_0^me_liq
    output =
        smooth_transition ? CO.logistic_function(q_lcl, q_liq_threshold, k) :
        CO.heaviside(q_lcl - q_liq_threshold)
    return D * q_lcl^a * N_d^b * output
end
function conv_q_lcl_to_q_rai(
    (; ρ_w, R_6C_0, E_0, k)::CMP.LD2004, q_lcl, ρ, N_d, smooth_transition = false,
)
    if q_lcl <= UT.ϵ_numerics_2M_M(eltype(q_lcl))
        return zero(UT.promote_typeof(q_lcl, ρ, N_d))
    else
        # Mean volume radius in microns (assuming spherical cloud droplets)
        r_vol = cbrt(3 * q_lcl * ρ / 4 / π / ρ_w / N_d) * 1_000_000

        # Assumed size distribution: modified gamma distribution
        β_6 = cbrt((r_vol + 3) / r_vol)
        E = E_0 * β_6^6
        R_6 = β_6 * r_vol
        R_6C = R_6C_0 / cbrt(sqrt(q_lcl * ρ)) / sqrt(R_6)  # cbrt(sqrt(x)) = x^(1/6)

        output =
            smooth_transition ? CO.logistic_function(R_6, R_6C, k) :
            CO.heaviside(R_6 - R_6C)
        return E * (q_lcl * ρ)^3 / N_d / ρ * output
    end
end
"""
    accretion(scheme::CMP.KK2000, q_lcl, q_rai, ρ)
    accretion(scheme::CMP.B1994, q_lcl, q_rai, ρ)
    accretion(scheme::CMP.TC1980, q_lcl, q_rai)

Compute the accretion rate of rain.

# Arguments
 - `scheme` - accretion scheme parameters
   (the accretion parameters are taken from its `accr` field):
   - Khairoutdinov and Kogan (2000) for `scheme::CMP.KK2000`
   - Beheng (1994) for `scheme::CMP.B1994`
   - Tripoli and Cotton (1980) for `scheme::CMP.TC1980`
 - `q_lcl` - cloud liquid water specific content
 - `q_rai` - rain water specific content
 - `ρ` - air density (for `CMP.KK2000` and `CMP.B1994` only)
"""
function accretion((; accr)::CMP.KK2000, q_lcl, q_rai, ρ)
    q_lcl = max(0, q_lcl)
    q_rai = max(0, q_rai)
    (; A, a, b) = accr
    return A * (q_lcl * q_rai)^a * ρ^b
end

function accretion((; accr)::CMP.B1994, q_lcl, q_rai, ρ)
    q_lcl = max(0, q_lcl)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_lcl * ρ * q_rai
end

function accretion((; accr)::CMP.TC1980, q_lcl, q_rai)
    #TODO - The original paper is actually formulated for mixing ratios, not specific contents
    q_lcl = max(0, q_lcl)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_lcl * q_rai
end

end # module
