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

import ..ThermodynamicsInterface as TDI
import ..Common as CO
import ..Parameters as CMP
import ..DistributionTools as DT
import ..Utilities as UT

import ..DistributionTools: size_distribution

export autoconversion,
    accretion,
    cloud_liquid_self_collection,
    autoconversion_and_cloud_liquid_self_collection,
    rain_terminal_velocity,
    conv_q_lcl_to_q_rai,
    rain_evaporation,
    rain_self_collection,
    rain_breakup,
    rain_self_collection_and_breakup,
    size_distribution,
    get_size_distribution_bounds,
    number_increase_for_mass_limit,
    number_decrease_for_mass_limit

"""
    pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)

Return the parameters of the rain drop diameter distribution

    n_r(D) = N_0 * exp(- D / Dr_mean)

 where
 - `D` is the diameter of the raindrop,
 - `N_0` [1/m³] is the number concentration of raindrops,
 - `Dr_mean` [m] is the mean diameter of the raindrops.

 Note: in SB2006, Eq. (83) the distribution is given as:

    f(D) = N_0 * exp(- λ_r D)

 where `λ_r ≡ 1 / Dr_mean` [1/m] is the inverse of the mean diameter of the raindrops.

# Arguments
 - `pdf_r`: struct containing size distribution parameters for rain.
        Can either be [`CMP.RainParticlePDF_SB2006_notlimited`](@ref) or [`CMP.RainParticlePDF_SB2006_limited`](@ref).
        For the latter, the values for `N_0`, `Dr_mean`, and `xr_mean` are limited to be within provided ranges.
 - `qᵣ`: mass of rain water [kg]
 - `ρₐ`: air density [kg/m³]
 - `Nᵣ`: number of rain drops [1/m³]

# Returns
 - A `NamedTuple` with the fields `(; N₀r, Dr_mean, xr_mean)`
"""
function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_notlimited, qᵣ, ρₐ, Nᵣ)
    (iszero(Nᵣ) && iszero(qᵣ)) && return (; N₀r = zero(Nᵣ), Dr_mean = zero(qᵣ), xr_mean = zero(qᵣ))
    (; ρw) = pdf_r
    Lᵣ = ρₐ * qᵣ

    xr_mean = Lᵣ / Nᵣ
    λr = cbrt(π * ρw / xr_mean)
    N₀r = λr * Nᵣ

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    return (; N₀r, Dr_mean, xr_mean)
end
function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_limited, qᵣ, ρₐ, Nᵣ)
    FT = eltype(qᵣ)
    (Nᵣ < UT.ϵ_numerics_2M_N(FT) && qᵣ < UT.ϵ_numerics_2M_M(FT)) &&
        return (; N₀r = zero(Nᵣ), Dr_mean = zero(qᵣ), xr_mean = zero(qᵣ))
    (; xr_min, xr_max, N0_min, N0_max, λ_min, λ_max, ρw) = pdf_r
    Lᵣ = ρₐ * max(0, qᵣ)

    # Sequence of limiting steps in Seifert and Beheng 2006:
    x̃r = clamp(Lᵣ / Nᵣ, xr_min, xr_max)  # Eq. (94)
    N₀r = clamp(Nᵣ * cbrt(π * ρw / x̃r), N0_min, N0_max)  # Eq. (95)
    λr = clamp(sqrt(sqrt(π * ρw * N₀r / Lᵣ)), λ_min, λ_max)  # Eq. (96)
    xr_mean = clamp(Lᵣ * λr / N₀r, xr_min, xr_max)  # Eq. (97)

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    return (; N₀r, Dr_mean, xr_mean)
end

"""
    pdf_rain_parameters_mass(pdf_r::CMP.RainParticlePDF_SB2006, qᵣ, ρₐ, Nᵣ)

Return the parameters of the rain drop diameter distribution in terms of mass.

As a function of diameter, the size distribution is given by:

    n(D) = N₀r * exp(-D / Dr_mean)

In terms of mass (`x`), the size distribution is given by:

    f(x) = n(D(x)) * ∂D∂x(x)
         = N₀ * exp(-D(x) / Dr_mean) * 2 / (π * ρw) * x^(-2/3)
         = N₀ * 2 / (π * ρw) * x^(-2/3) * exp(- (6 / (π * ρw))^(1/3) / Dr_mean * x^(1/3))

 where
 - `D(x) = (6x / (π * ρw))^(1/3)` is the diameter of a raindrop of mass `x`.
 - `∂D∂x(x) = (6 / (π * ρw))^(1/3) * x^(-2/3)` is the derivative of the diameter with respect to the mass.

If we write the general form of the size distribution as:

    f(x) = A * x^ν * exp(-B * x^μ)

 then we have that:
 - `A = N₀ * 2 / (π * ρw)`
 - `B = (6 / (π * ρw))^(1/3) / Dr_mean`
 - `ν = -2/3`
 - `μ = 1/3`
"""
function pdf_rain_parameters_mass(pdf_r::CMP.RainParticlePDF_SB2006, qᵣ, ρₐ, Nᵣ)
    (; xr_mean) = pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)
    Br = cbrt(6 / xr_mean)
    Ar = Nᵣ * Br / 3
    return (; Ar, Br)
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

# Arguments
 - `pdf_c`: Size distribution parameters for cloud droplets, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: Liquid mass content [kg/kg]
 - `ρₐ`: Air density [kg/m³]
 - `N`: Number concentration of the particle [1/m³]

# Returns
 - `(logA, logB)`: Log of the parameters of the generalized gamma distribution
"""
function log_pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)
    FT = eltype(q)
    L = ρₐ * q
    # If L or N are (essentially) zero, return `A=0` (no number per mass), `B=∞` (zero mass "length" scale)
    (N < UT.ϵ_numerics_2M_N(FT) || L < UT.ϵ_numerics_2M_M(FT)) && return (log(zero(N)), log(1 / zero(q)))
    (; νc, μc) = pdf_c
    logx̄ = log(L / N)
    z1 = (νc + 1) / μc
    z2 = (νc + 2) / μc
    logB = -μc * (logx̄ + SF.loggamma(z1) - SF.loggamma(z2))
    logA = log(μc) + log(N) + z1 * logB - SF.loggamma(z1)
    return (logA, logB)
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
 - `(logN₀c, λc, νcD, μcD)`: Parameters of the generalized gamma distribution in terms of diameter
"""
function pdf_cloud_parameters(pdf_c, q, ρₐ, N)
    FT = eltype(pdf_c)
    logAc, logBc = log_pdf_cloud_parameters_mass(pdf_c, q, ρₐ, N)
    (; νc, μc, ρw) = pdf_c
    k_m = ρw * π / 6
    # Convert from mass-based to diameter-based distribution
    logN₀c = logAc + log(FT(3)) + (νc + 1) * log(k_m)
    λc = exp(logBc) * k_m^μc
    return (; logN₀c, λc, νcD = 3νc + 2, μcD = 3μc)
end

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
    return n(D) = iszero(N₀r) ? zero(D) : N₀r * exp(-D / Dr_mean)
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
    return n(D) = logN₀c == -Inf ? zero(D) : exp(logN₀c + νcD * log(D) - λc * D^μcD)
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
    get_size_distribution_bounds(pdf, q, ρₐ, N, p)

Return the minimum and maximum diameters of a cloud or rain particle such that
the size distribution is within (1 - p) to (p) probability of the true size distribution.

# Arguments
 - `pdf`: Size distribution parameters for cloud or rain,
    [`CMP.RainParticlePDF_SB2006`](@ref) or [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q`: mass mixing ratio of cloud or rain water
 - `ρₐ`: density of air
 - `N`: number mixing ratio of cloud or rain
 - `p`: probability level (0 ≤ p ≤ 1)

# Returns
 - `D_min, D_max`: minimum and maximum diameters of a cloud or rain particle such that
    the size distribution is within (1 - p) to (p) probability of the true size distribution.
    All inputs and output diameters are in base SI units.
    The bounds are calculated through quantile functions of the size distribution.
"""
function get_size_distribution_bounds(
    pdf::CMP.RainParticlePDF_SB2006{FT},
    q, ρₐ, N, p = eps(FT),
) where {FT}
    (; Dr_mean) = pdf_rain_parameters(pdf, q, ρₐ, N)
    iszero(Dr_mean) && return (FT(0), FT(0))
    D_min = DT.exponential_quantile(Dr_mean, p)
    D_max = DT.exponential_quantile(Dr_mean, 1 - p)
    return D_min, D_max
end
function get_size_distribution_bounds(
    pdf::CMP.CloudParticlePDF_SB2006{FT},
    q, ρₐ, N, p = eps(FT),
) where {FT}
    (; λc, νcD, μcD) = pdf_cloud_parameters(pdf, q, ρₐ, N)

    D_min = DT.generalized_gamma_quantile(νcD, μcD, λc, p)
    D_max = DT.generalized_gamma_quantile(νcD, μcD, λc, 1 - p)
    return D_min, D_max
end


### ----- ###
### RATES ###
### ----- ###

"""
A structure containing the rates of change of the specific contents and number
densities of cloud liquid water and rain water.
"""
Base.@kwdef struct LclRaiRates{FT}
    "Rate of change of the cloud liquid water specific content"
    dq_lcl_dt::FT = FT(0)
    "Rate of change of the cloud liquid water number density"
    dN_lcl_dt::FT = FT(0)
    "Rate of change of the rain water specific content"
    dq_rai_dt::FT = FT(0)
    "Rate of change of the rain water number density"
    dN_rai_dt::FT = FT(0)
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
    acnv::CMP.AcnvSB2006{FT}, pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    q_lcl, q_rai, ρ, N_lcl,
) where {FT}

    if q_lcl < UT.ϵ_numerics_2M_M(FT) || N_lcl < UT.ϵ_numerics_2M_N(FT)
        return LclRaiRates{FT}()
    end

    (; kcc, x_star, ρ0, A, a, b) = acnv
    (; νc) = pdf_c

    L_lcl = ρ * q_lcl
    x_lcl = min(x_star, L_lcl / N_lcl)
    q_rai = max(FT(0), q_rai)
    τ = 1 - q_lcl / (q_lcl + q_rai)  # Eq. (5) from SB2006
    ϕ_au = A * τ^a * (1 - τ^a)^b

    dL_rai_dt =
        kcc / 20 / x_star * (νc + 2) * (νc + 4) / (νc + 1)^2 *
        L_lcl^2 * x_lcl^2 * (1 + ϕ_au / (1 - τ)^2) * ρ0 / ρ  # Eq. (4) from SB2006
    dN_rai_dt = dL_rai_dt / x_star
    dL_lcl_dt = -dL_rai_dt
    dN_lcl_dt = -2 * dN_rai_dt

    return LclRaiRates(
        dq_lcl_dt = dL_lcl_dt / ρ,
        dN_lcl_dt = dN_lcl_dt,
        dq_rai_dt = dL_rai_dt / ρ,
        dN_rai_dt = dN_rai_dt,
    )
end

"""
    accretion(accr, q_lcl, q_rai, ρ, N_lcl)

Compute accretion rate

# Arguments
 - `accr`: Accretion parameters, [`CMP.AccrSB2006`](@ref)
 - `q_lcl`: Cloud liquid water specific content [kg/kg]
 - `q_rai`: Rain water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `N_lcl`: Cloud droplet number density [1/m³]

# Returns
 - [`LclRaiRates`](@ref) with `q_lcl`, `N_lcl`, `q_rai`, `N_rai` tendencies due to
    collisions between raindrops and cloud droplets (accretion)
"""
function accretion((; accr)::CMP.SB2006{FT}, q_lcl, q_rai, ρ, N_lcl) where {FT}

    if q_lcl < UT.ϵ_numerics_2M_M(FT) || q_rai < UT.ϵ_numerics_2M_M(FT) || N_lcl < UT.ϵ_numerics_2M_N(FT)
        return LclRaiRates{FT}()
    end

    (; kcr, τ0, ρ0, c) = accr
    L_lcl = ρ * q_lcl
    L_rai = ρ * q_rai
    x_lcl = L_lcl / N_lcl
    τ = 1 - q_lcl / (q_lcl + q_rai)  # Eq. (5) from SB2006
    ϕ_ac = (τ / (τ + τ0))^c          # Eq. (8) from SB2006

    dL_rai_dt = kcr * L_lcl * L_rai * ϕ_ac * sqrt(ρ0 / ρ)  # Eq. (7) from SB2006
    dN_rai_dt = zero(N_lcl)
    dL_lcl_dt = -dL_rai_dt
    dN_lcl_dt = dL_lcl_dt / x_lcl

    return LclRaiRates(
        dq_lcl_dt = dL_lcl_dt / ρ,
        dN_lcl_dt = dN_lcl_dt,
        dq_rai_dt = dL_rai_dt / ρ,
        dN_rai_dt = dN_rai_dt,
    )
end

"""
    cloud_liquid_self_collection(acnv, pdf_c, q_lcl, ρ, dN_lcl_dt_au)

Compute cloud liquid self-collection rate

# Arguments
 - `acnv`: 2-moment autoconversion parameterization, [`CMP.AcnvSB2006`](@ref)
 - `pdf_c`: Cloud size distribution parameters, [`CMP.CloudParticlePDF_SB2006`](@ref)
 - `q_lcl`: Cloud liquid water specific content [kg/kg]
 - `ρ`: Air density [kg/m³]
 - `dN_lcl_dt_au`: Rate of change of cloud droplets number density due to autoconversion [1/m³/s]

# Returns
 - The cloud droplets number density tendency due to collisions of cloud droplets
    that produce larger cloud droplets (self-collection)
"""
function cloud_liquid_self_collection(
    acnv::CMP.AcnvSB2006{FT}, pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    q_lcl, ρ, dN_lcl_dt_au,
) where {FT}

    if q_lcl < UT.ϵ_numerics_2M_M(FT)
        return FT(0)
    end
    (; kcc, ρ0) = acnv
    (; νc) = pdf_c

    L_lcl = ρ * q_lcl

    # Eq. (9) from SB2006
    dN_lcl_dt_sc = -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * L_lcl^2 - dN_lcl_dt_au

    return dN_lcl_dt_sc
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
    (; acnv, pdf_c)::CMP.SB2006{FT},
    q_lcl, q_rai, ρ, N_lcl,
) where {FT}

    au = autoconversion(acnv, pdf_c, q_lcl, q_rai, ρ, N_lcl)
    sc = cloud_liquid_self_collection(acnv, pdf_c, q_lcl, ρ, au.dN_lcl_dt)

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
    pdf::CMP.RainParticlePDF_SB2006{FT}, self::CMP.SelfColSB2006{FT},
    q_rai, ρ, N_rai,
) where {FT}

    if q_rai < UT.ϵ_numerics_2M_M(FT) || N_rai < UT.ϵ_numerics_2M_N(FT)
        return FT(0)
    end

    (; krr, κrr, d) = self
    (; ρ0) = pdf

    L_rai = ρ * q_rai
    (; Br) = pdf_rain_parameters_mass(pdf, q_rai, ρ, N_rai)
    dN_rai_dt_sc = -krr * N_rai * L_rai * √(ρ0 / ρ) * (1 + κrr / Br)^d  # Eq. (11) from SB2006

    return dN_rai_dt_sc
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
    pdf::CMP.RainParticlePDF_SB2006{FT}, brek::CMP.BreakupSB2006{FT},
    q_rai, ρ, N_rai, dN_rai_dt_sc,
) where {FT}

    if q_rai < UT.ϵ_numerics_2M_M(FT) || N_rai < UT.ϵ_numerics_2M_N(FT)
        return FT(0)
    end
    (; Deq, Dr_th, kbr, κbr) = brek
    (; ρw) = pdf
    (; xr_mean) = pdf_rain_parameters(pdf, q_rai, ρ, N_rai)
    Dr = cbrt(xr_mean * 6 / (π * ρw))  # mean volume raindrop diameter
    ΔD = Dr - Deq
    Φ_br = if Dr < Dr_th  # Below the threshold diameter, breakup is neglected
        FT(-1)
    elseif Dr ≤ Deq  # Below the equilibrium diameter, breakup is parameterized as a linear function
        kbr * ΔD
    else
        exp(κbr * ΔD) - 1  # Above the equilibrium diameter, breakup is parameterized as an exponential function
    end
    dN_rai_dt_br = -(Φ_br + 1) * dN_rai_dt_sc  # Eq. (13) from SB2006

    return dN_rai_dt_br
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
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    (; ρw, grav, ν_air)::CMP.StokesRegimeVelType{FT},
    q_liq, ρₐ, N_liq,
) where {FT}

    if N_liq < UT.ϵ_numerics_2M_N(FT) || q_liq < UT.ϵ_numerics_2M_M(FT)
        return (FT(0), FT(0))
    end

    (; νc, μc) = pdf_c
    (; Bc) = pdf_cloud_parameters_mass(pdf_c, q_liq, ρₐ, N_liq)

    terminal_velocity_prefactor = FT(1 / 18) * (6 / ρw / π)^(2 // 3) * (ρw / ρₐ - 1) * grav / ν_air
    vt0 =
        N_liq < UT.ϵ_numerics_2M_N(FT) ? FT(0) :
        terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_liq, FT(2 / 3)) / N_liq
    vt1 =
        q_liq < UT.ϵ_numerics_2M_M(FT) ? FT(0) :
        terminal_velocity_prefactor * DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_liq, FT(5 / 3)) / ρₐ / q_liq

    return (vt0, vt1)

end

"""
    rain_terminal_velocity(SB2006, vel, q_rai, ρ, N_rai)

Compute the raindrops terminal velocity.

# Arguments
 - `pdf_r`: Rain size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
 - `vel`: Terminal velocity parameters, [`CMP.Chen2022VelTypeRain`](@ref)
 - `q_rai`: Rain water specific content
 - `ρ`: Air density
 - `N_rai`: Raindrops number density

# Returns
A tuple containing the number and mass weigthed mean fall velocities of raindrops in [m/s].
Assuming an exponential size distribution from Seifert and Beheng 2006 for `scheme == SB2006Type`
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical relation similar to Rogers (1993) for `velo_scheme == SB2006VelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171 for `velo_scheme == Chen2022Type`
"""
function rain_terminal_velocity(
    (; pdf_r)::CMP.SB2006{FT}, (; ρ0, aR, bR, cR)::CMP.SB2006VelType{FT},
    q_rai, ρ, N_rai,
) where {FT}
    # TODO: Input argument list needs to be redesigned

    (; Dr_mean) = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai)
    _pa0, _pb0, _pa1, _pb1 =
        _sb_rain_terminal_velocity_helper(pdf_r, 1 / Dr_mean, aR, bR, cR)

    vt0 =
        N_rai < UT.ϵ_numerics_2M_N(FT) ? FT(0) :
        max(FT(0), sqrt(ρ0 / ρ) * (aR * _pa0 - bR * _pb0 / (1 + cR * Dr_mean)))
    vt1 =
        q_rai < UT.ϵ_numerics_2M_M(FT) ? FT(0) :
        max(FT(0), sqrt(ρ0 / ρ) * (aR * _pa1 - bR * _pb1 / (1 + cR * Dr_mean)^4))
    return (vt0, vt1)
end
function rain_terminal_velocity(
    (; pdf_r)::CMP.SB2006{FT}, vel::CMP.Chen2022VelTypeRain{FT},
    q_rai, ρ, N_rai,
) where {FT}
    # coefficients from Table B1 from Chen et. al. 2022
    aiu, bi, ciu = CO.Chen2022_vel_coeffs(vel, ρ)
    # size distribution parameter
    (; Dr_mean) = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai)

    # eq 20 from Chen et al 2022
    vt0 = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, Dr_mean, 0))
    vt3 = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, Dr_mean, 3))

    vt0 = N_rai < UT.ϵ_numerics_2M_N(FT) ? FT(0) : max(FT(0), vt0)
    vt3 = q_rai < UT.ϵ_numerics_2M_M(FT) ? FT(0) : max(FT(0), vt3)
    # It should be (ϕ^κ * vt0, ϕ^κ * vt3), but for rain drops ϕ = 1 and κ = 0
    return (vt0, vt3)
end
function _sb_rain_terminal_velocity_helper(
    pdf_r::CMP.RainParticlePDF_SB2006_limited{FT}, λr, aR, bR, cR,
) where {FT}
    return (FT(1), FT(1), FT(1), FT(1))
end
function _sb_rain_terminal_velocity_helper(
    pdf_r::CMP.RainParticlePDF_SB2006_notlimited{FT}, λr, aR, bR, cR,
) where {FT}
    # Integrate velocity of particles over a range of r with
    # positive terminal velocity (v = aR - bR exp(-lambda D))
    _rc = -1 / (2 * cR) * log(aR / bR)
    _Γ_1(t) = exp(-t)
    _Γ_4(t) = (t^3 + 3 * t^2 + 6 * t + 6) * exp(-t)
    _pa0::FT = _Γ_1(2 * _rc * λr)
    _pb0::FT = _Γ_1(2 * _rc * (λr + cR))
    _pa1::FT = _Γ_4(2 * _rc * λr) / 6
    _pb1::FT = _Γ_4(2 * _rc * (λr + cR)) / 6
    return (_pa0, _pb0, _pa1, _pb1)
end

"""
    Γ_incl(a, x)

Returns the approximation of an incomplete gamma function for a ∈ {-1.0, -0.101}, and x in [0.067 1.82]
"""
function Γ_incl(a::FT, x::FT) where {FT}
    #return exp(-x) / ((FT(1.5) - FT(0.54) * a) * x^(FT(0.46) - FT(0.75) * a))
    return exp(-x) / (
        (FT(0.33) - FT(0.7) * a) * x^(FT(0.08) - FT(0.93) * a) +
        (FT(1.34) - FT(0.1) * a) * x^(FT(0.8) - a)
    )
end

"""
    rain_evaporation(evap, aps, tps, q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T)

 - `evap` - evaporation parameterization scheme
 - `aps` - air properties
 - `tps` - thermodynamics parameters
 - `q_tot`, `q_lcl`, `q_icl`, `q_rai`, `q_sno` - total water,
    cloud liquid water, cloud ice, rain and snow specific contents,
 - `ρ` - air density
 - `N_rai` - raindrops number density
 - `T` - air temperature

Returns a named tuple containing the tendency of raindrops number density and rain water
specific content due to rain rain_evaporation, assuming a power law velocity relation for
fall velocity of individual drops and an exponential size distribution, for `scheme == SB2006Type`
"""
function rain_evaporation(
    (; pdf_r, evap)::CMP.SB2006{FT}, aps::CMP.AirProperties,
    tps::TDI.PS,
    q_tot, q_lcl, q_icl, q_rai, q_sno, ρ, N_rai, T,
) where {FT}

    evap_rate_0 = FT(0)
    evap_rate_1 = FT(0)
    S = TDI.supersaturation_over_liquid(tps, q_tot, q_lcl + q_rai, q_icl + q_sno, ρ, T)

    if (N_rai > UT.ϵ_numerics_2M_N(FT) && S < FT(0))

        (; ν_air, D_vapor) = aps
        (; av, bv, α, β, ρ0) = evap
        x_star = pdf_r.xr_min
        ρw = pdf_r.ρw
        G = CO.G_func_liquid(aps, tps, T)

        (; xr_mean) = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai)
        Dr = cbrt(6 * xr_mean / (π * ρw))

        t_star = cbrt(6 * x_star / xr_mean)
        a_vent_0 = av * Γ_incl(FT(-1), t_star) / FT(6)^(-2 // 3)
        b_vent_0 = bv * Γ_incl(-1 // 2 + 3 // 2 * β, t_star) / FT(6)^(β / 2 - 1 // 2)

        a_vent_1 = av * SF.gamma(FT(2)) / cbrt(FT(6))
        b_vent_1 = bv * SF.gamma(5 // 2 + 3 // 2 * β) / FT(6)^(β / 2 + 1 // 2)

        N_Re = α * xr_mean^β * sqrt(ρ0 / ρ) * Dr / ν_air
        Fv0 = a_vent_0 + b_vent_0 * cbrt(ν_air / D_vapor) * sqrt(N_Re)
        Fv1 = a_vent_1 + b_vent_1 * cbrt(ν_air / D_vapor) * sqrt(N_Re)

        evap_rate_0 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv0 / xr_mean)
        evap_rate_1 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)

        # When xr = 0 evap_rate_0 becomes NaN. We replace NaN with 0 which is the limit of
        # evap_rate_0 for xr -> 0.
        evap_rate_0 = xr_mean / x_star < eps(FT) ? FT(0) : evap_rate_0
        evap_rate_1 = q_rai < UT.ϵ_numerics_2M_M(FT) ? FT(0) : evap_rate_1

    end

    return (; evap_rate_0, evap_rate_1)
end

"""
    number_increase_for_mass_limit(numadj, x_max, q, ρ, N)

Compute the tendency (rate of change) of number concentration `N` required to ensure that
the mean particle mass `x = ρq / N` does not exceed the upper limit `x_max`. Returns a positive
tendency when the mean mass is too high (`x > x_max`), and zero otherwise.
The method is based on Horn (2012, DOI: [10.5194/gmd-5-345-2012](https://doi.org/10.5194/gmd-5-345-2012)).

# Arguments
- `numadj`: Number concentration adjustment parameters ([CMP.NumberAdjustmentHorn2012](@ref))
- `q`: Mass mixing ratio [kg/kg]
- `ρ`: Air density [kg/m³]
- `N`: Number concentration [1/m³]
- `x_max`: Maximum allowed mean particle mass [kg]

# Returns
- The rate of change of number concentration [1/(m³·s)] needed to bring the mean mass within the upper bound.
"""
function number_increase_for_mass_limit(
    (; τ)::CMP.NumberAdjustmentHorn2012{FT}, x_max, q, ρ, N,
) where {FT}
    return max(FT(0), ρ * q / x_max - N) / τ
end

"""
    number_decrease_for_mass_limit(numadj, x_min, q, ρ, N)

Compute the tendency (rate of change) of number concentration `N` required to ensure that
the mean particle mass `x = ρq / N` does not fall below the lower limit `x_min`. Returns a negative
tendency when the mean mass is too low (`x < x_min`), and zero otherwise.
The method is based on Horn (2012, DOI: [10.5194/gmd-5-345-2012](https://doi.org/10.5194/gmd-5-345-2012)).

# Arguments
- `numadj`: Number concentration adjustment parameters ([CMP.NumberAdjustmentHorn2012](@ref))
- `q`: Mass mixing ratio [kg/kg]
- `ρ`: Air density [kg/m³]
- `N`: Number concentration [1/m³]
- `x_min`: Minimum allowed mean particle mass [kg]

# Returns
- The rate of change of number concentration [1/(m³·s)] needed to bring the mean mass within the lower bound.
"""
function number_decrease_for_mass_limit(
    (; τ)::CMP.NumberAdjustmentHorn2012{FT}, x_min, q, ρ, N,
) where {FT}

    # Avoid NaN when both q and x_min are 0; use typed Inf to avoid type promotion
    N_max = iszero(x_min) ? FT(Inf) : ρ * q / x_min

    return min(FT(0), N_max - N) / τ
end

# Additional double moment autoconversion and accretion parametrizations:
# - Khairoutdinov and Kogan (2000)
# - Beheng (1994)
# - Tripoli and Cotton (1980)
# - Liu and Daum (2004)
# - variable timescale autoconversion Azimi (2023)

"""
    conv_q_lcl_to_q_rai(acnv, q_lcl, ρ, N_d; smooth_transition)

 - `acnv` - 2-moment rain autoconversion parameterization
 - `q_lcl` - cloud liquid water specific content
 - `ρ` - air density
 - `N_d` - prescribed cloud droplet number concentration

Returns the q_rai tendency due to collisions between cloud droplets
(autoconversion), parametrized following:
 - Khairoutdinov and Kogan (2000) for `scheme == KK2000Type`
 - Beheng (1994) for `scheme == B1994Type`
 - Tripoli and Cotton (1980) for `scheme == TC1980Type`
 - Liu and Daum (2004) for `scheme ==LD2004Type`

The `Beheng1994Type`, `TC1980Type` and `LD2004Type` of schemes
additionally accept `smooth_transition` flag that
smoothes their thershold behaviour if set to `true`.
The default value is `false`.
"""
function conv_q_lcl_to_q_rai((; acnv)::CMP.KK2000{FT}, q_lcl, ρ, N_d) where {FT}
    q_lcl = max(0, q_lcl)
    (; A, a, b, c) = acnv
    return A * q_lcl^a * N_d^b * ρ^c
end
function conv_q_lcl_to_q_rai(
    (; acnv)::CMP.B1994{FT},
    q_lcl,
    ρ,
    N_d,
    smooth_transition = false,
) where {FT}
    q_lcl = max(0, q_lcl)
    (; C, a, b, c, N_0, k, d_low, d_high) = acnv
    d = FT(0)
    if smooth_transition
        d_low_acnv_fraction = CO.logistic_function(N_d, N_0, k)
        d_high_acnv_fraction = FT(1) - d_low_acnv_fraction
        d = d_low_acnv_fraction * d_low + d_high_acnv_fraction * d_high
    else
        d = N_d >= N_0 ? d_low : d_high
    end
    return C * d^a * (q_lcl * ρ)^b * N_d^c / ρ
end
function conv_q_lcl_to_q_rai(
    (; acnv)::CMP.TC1980{FT},
    q_lcl,
    ρ,
    N_d,
    smooth_transition = false,
) where {FT}
    #TODO - The original paper is actually formulated for mixing ratios, not specific contents
    q_lcl = max(0, q_lcl)
    (; m0_liq_coeff, me_liq, D, a, b, r_0, k) = acnv
    q_liq_threshold::FT = m0_liq_coeff * N_d / ρ * r_0^me_liq
    output =
        smooth_transition ? CO.logistic_function(q_lcl, q_liq_threshold, k) :
        CO.heaviside(q_lcl - q_liq_threshold)
    return D * q_lcl^a * N_d^b * output
end
function conv_q_lcl_to_q_rai(
    (; ρ_w, R_6C_0, E_0, k)::CMP.LD2004{FT},
    q_lcl,
    ρ,
    N_d,
    smooth_transition = false,
) where {FT}
    if q_lcl <= UT.ϵ_numerics_2M_M(FT)
        return FT(0)
    else
        # Mean volume radius in microns (assuming spherical cloud droplets)
        r_vol =
            (FT(3) * (q_lcl * ρ) / FT(4) / FT(π) / ρ_w / N_d)^FT(1 / 3) *
            FT(1e6)

        # Assumed size distribution: modified gamma distribution
        β_6 = ((r_vol + FT(3)) / r_vol)^FT(1 / 3)
        E = E_0 * β_6^6
        R_6 = β_6 * r_vol
        R_6C = R_6C_0 / (q_lcl * ρ)^FT(1 / 6) / R_6^FT(1 / 2)

        output =
            smooth_transition ? CO.logistic_function(R_6, R_6C, k) :
            CO.heaviside(R_6 - R_6C)
        return E * (q_lcl * ρ)^3 / N_d / ρ * output
    end
end
function conv_q_lcl_to_q_rai(
    (; τ, α)::CMP.VarTimescaleAcnv{FT},
    q_lcl::FT,
    ρ::FT,
    N_d::FT,
) where {FT}
    return max(0, q_lcl) / (τ * (N_d / 1e8)^α)
end

"""
    accretion(accretion_scheme, q_lcl, q_rai, ρ)

 - `accretion_scheme` - type for 2-moment rain accretion parameterization
 - `q_lcl` - cloud liquid water specific content
 - `q_rai` - rain water specific content
 - `ρ` - air density (for `KK2000Type` and `Beheng1994Type`)

 Returns the accretion rate of rain, parametrized following
 - Khairoutdinov and Kogan (2000) for `scheme == KK2000Type`
 - Beheng (1994) for `scheme == B1994Type`
 - Tripoli and Cotton (1980) for `scheme == TC1980Type`
"""
function accretion((; accr)::CMP.KK2000{FT}, q_lcl, q_rai, ρ) where {FT}
    q_lcl = max(0, q_lcl)
    q_rai = max(0, q_rai)
    (; A, a, b) = accr
    return A * (q_lcl * q_rai)^a * ρ^b
end

function accretion((; accr)::CMP.B1994{FT}, q_lcl, q_rai, ρ) where {FT}
    q_lcl = max(0, q_lcl)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_lcl * ρ * q_rai
end

function accretion((; accr)::CMP.TC1980{FT}, q_lcl, q_rai) where {FT}
    #TODO - The original paper is actually formulated for mixing ratios, not specific contents
    q_lcl = max(0, q_lcl)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_lcl * q_rai
end

end # module
