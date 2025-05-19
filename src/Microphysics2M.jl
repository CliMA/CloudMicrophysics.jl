"""
Double-moment bulk microphysics parametrizations including:
 - autoconversion, accretion, self-collection, breakup, mean terminal velocity of raindrops and rain
    evaporation rates from Seifert and Beheng 2006,
 - additional double-moment bulk microphysics autoconversion and accretion rates
   from: Khairoutdinov and Kogan 2000, Beheng 1994, Tripoli and Cotton 1980, and
   Liu and Daum 2004.
"""
module Microphysics2M

import SpecialFunctions as SF
import RootSolvers as RS

import Thermodynamics as TD
import Thermodynamics.Parameters as TDP

import ..Common as CO
import ..Parameters as CMP

export autoconversion,
    accretion,
    liquid_self_collection,
    autoconversion_and_liquid_self_collection,
    rain_terminal_velocity,
    conv_q_liq_to_q_rai,
    rain_evaporation,
    rain_self_collection,
    rain_breakup,
    rain_self_collection_and_breakup,
    size_distribution

"""
A structure containing the rates of change of the specific contents and number
densities of liquid and rain water.
"""
Base.@kwdef struct LiqRaiRates{FT}
    "Rate of change of the liquid water specific content"
    dq_liq_dt::FT = FT(0)
    "Rate of change of the liquid water number density"
    dN_liq_dt::FT = FT(0)
    "Rate of change of the rain water specific content"
    dq_rai_dt::FT = FT(0)
    "Rate of change of the rain water number density"
    dN_rai_dt::FT = FT(0)
end

# Double-moment bulk microphysics autoconversion, accretion, self-collection, breakup,
# mean terminal velocity of raindrops, and rain evaporation rates from Seifert and Beheng 2001

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
        Can either be [`RainParticlePDF_SB2006_notlimited`](@ref) or [`RainParticlePDF_SB2006_limited`](@ref).
        For the latter, the values for ``N_0``, ``Dr_mean``, and ``xr_mean`` are limited to be within provided ranges.
 - `qᵣ`: mass of rain water [kg]
 - `ρₐ`: air density [kg/m³]
 - `Nᵣ`: number of rain drops [1/m³]

# Returns
 - A `NamedTuple` with the fields `(; N₀r, Dr_mean, xr_mean)`
"""
function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_notlimited, qᵣ, ρₐ, Nᵣ)
    (; ρw) = pdf_r
    Lᵣ = ρₐ * qᵣ

    xr_mean = Lᵣ / Nᵣ 
    λr = ∛(π * ρw / xr_mean)
    N₀r = λr * Nᵣ

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    return (; N₀r, Dr_mean, xr_mean)
end
function pdf_rain_parameters(pdf_r::CMP.RainParticlePDF_SB2006_limited, qᵣ, ρₐ, Nᵣ)
    FT = eltype(pdf_r)
    (; xr_min, xr_max, N0_min, N0_max, λ_min, λ_max, ρw) = pdf_r
    Lᵣ = ρₐ * max(0, qᵣ)

    # Sequence of limiting steps in Seifert and Beheng 2006:
    x̃r      = clamp(Lᵣ / Nᵣ,                 xr_min, xr_max)  # Eq. (94)  # TODO: Ill-defined for 0 / 0
    N₀r     = clamp(Nᵣ * ∛(π * ρw / x̃r),     N0_min, N0_max)  # Eq. (95)
    λr      = clamp((π * ρw * N₀r)^FT(1 / 4), λ_min,  λ_max)  # Eq. (96)
    xr_mean = clamp(Lᵣ * λr / N₀r,           xr_min, xr_max)  # Eq. (97)

    Dr_mean = 1 / λr  # The inverse of λr is the mean diameter of the raindrops (units: `m`)
    return (; N₀r, Dr_mean, xr_mean)
end

"""
    pdf_rain_parameters_mass(pdf_r::CMP.RainParticlePDF_SB2006_notlimited, qᵣ, ρₐ, Nᵣ)

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
    (; N₀r, Dr_mean) = pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)
    Ar = N₀r * 2 / (π * ρw)
    Br = (6 / (π * ρw))^(1/3) / Dr_mean
    return (; Ar, Br)
end

"""
    log_pdf_cloud_parameters_mass(L, N, ν, μ)

Return the log of the parameters of the generalized gamma distribution of the form

    f(x) = A * x^ν * exp(-B * x^μ),  [Eq. (79) in Seifert and Beheng 2006, but using the symbol `B` instead of `λ`]

where

    B_c = [  x̄_c Γ((ν_c + 1) / μ_c) / Γ((ν_c + 2) / μ_c) ]^(-μ_c)
    A_c = μ_c N_c B_c^((ν_c + 1) / μ_c) / Γ((ν_c + 1) / μ_c)

That is,

    log(B_c) = - μ_c [ log(x̄_c) + logΓ((ν_c + 1) / μ_c) - logΓ((ν_c + 2) / μ_c)) ]
    log(A_c) = log(μ_c) + log(N_c) + (ν_c + 1) / μ_c * log(B_c) - logΓ((ν_c + 1) / μ_c)

# Arguments
 - `L`: Liquid mass content [kg/m³]
 - `N`: Number concentration of the particle [1/m³]
 - `ν`: Exponent of the mass distribution
 - `μ`: Exponent of the size distribution

# Returns
 - `(logA, logB)`: Log of the parameters of the generalized gamma distribution
"""
function log_pdf_cloud_parameters_mass(L, N, ν, μ)
    logx̄ = log(L / N)
    logB = -μ * (logx̄ + SF.loggamma((ν + 1) / μ) - SF.loggamma((ν + 2) / μ))
    logA = log(μ) + log(N) + (ν + 1) / μ * logB - SF.loggamma((ν + 1) / μ)
    return (logA, logB)
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
    L_c = ρₐ * q_c
    logA, logB = log_pdf_cloud_parameters_mass(L_c, N_c, νc, μc)
    B = exp(logB)
    logpsd(x) = logA + νc * log(x) - B * x^μc
    return logpsd
end

"""
    size_distribution(pdf, q, ρₐ, N)

Return a function in diameter `D` that computes the size distribution value for rain or cloud particles.

# Arguments
- `pdf`: Struct containing size distribution parameters of cloud or rain
- `q`: Cloud or rain water specific content [kg/kg]
- `ρₐ`: Density of air [kg/m³]
- `N`: Cloud or rain water number concentration [1/m³]
- `D`: Particle size (i.e. maximum dimension of particle) [m]
"""
function size_distribution(pdf::CMP.RainParticlePDF_SB2006, q, ρₐ, N)
    (; N₀r, D_mean) = pdf_rain_parameters(pdf, q, ρₐ, N)
    return rain_psd(D) = N₀r * exp(-D / D_mean)
end

"""
    size_distribution(pdf::CMP.CloudParticlePDF_SB2006, q, ρₐ, N)

Return the size distribution, as a function of diameter, of the form

    n(D) = f_c(x(D)) * ∂x∂D(D)

where 
- `f_c(x)` is the size distribution in terms of mass
- `x(D)` is the mass of the particle as a function of diameter
- `∂x∂D(D)` is the derivative of the mass of a spherical particle with respect to diameter

"""
function size_distribution(pdf::CMP.CloudParticlePDF_SB2006, q, ρₐ, N)
    (; ρw) = pdf
    logpsd_mass = log_size_distribution_mass(pdf, q, ρₐ, N)
    ∂mass_∂D(D) = ρw * π / 2 * D^2
    mass(D) = ρw * CO.volume_sphere_D(D)
    return psd(D) = exp(logpsd_mass(mass(D))) * ∂mass_∂D(D)
end

"""
    size_distribution_value(pdf, q, ρₐ, N, D)

Return the size distribution value for a cloud or rain particle of diameter `D`.

# Arguments
- `pdf`: struct containing size distribution parameters for cloud or rain
- `q`: mass mixing ratio of cloud or rain water
- `ρₐ`: density of air
- `N`: number mixing ratio of cloud or rain
- `D`: diameter of the particle

See [`size_distribution`](@ref) for more details.
"""
function size_distribution_value(pdf, q, ρₐ, N, D)
    psd = size_distribution(pdf, q, ρₐ, N)
    return psd(D)
end

"""
    get_size_distribution_bound(pdf, q, N, ρₐ, tolerance)

 - pdf_r - struct containing size distribution parameters for cloud or rain
 - q - mass mixing ratio of cloud or rain water
 - N - number mixing ratio of cloud or rain
 - ρₐ - density of air
 - tolerance - tolerance for integration error

    Returns D_max value such that (1 - tolerance) = 1/N * ∫ N'(D) dD from 0 to D_max.
    All inputs and output D_max are in base SI units.
    For rain size distribution D_max is obtained analytically.
    For cloud size distribution D_max is calculated through a linear approximation
    of the bounds from numerical solutions.
"""
function get_size_distribution_bound(
    pdf::CMP.RainParticlePDF_SB2006,
    q, N, ρₐ, tolerance,
)
    (; Dr_mean) = pdf_rain_parameters(pdf, q, ρₐ, N)
    return -Dr_mean * log(tolerance)
end
function get_size_distribution_bound(
    pdf::CMP.CloudParticlePDF_SB2006{FT},
    q, N, ρₐ, tolerance,
) where {FT}
    L_c = ρₐ * q
    (; νc, μc) = pdf
    _, logB = log_pdf_cloud_parameters_mass(L_c, N, νc, μc)

    # This should give the mean diameter. Then use essentially same UB as for rain
    BD = exp(logB) * (ρ_w * π / 6)^μc  # constants in exponential of Eq. (6) in 2M docs
    z1 = (3νc + 4) / (3μc)
    z2 = (3νc + 3) / (3μc)
    mean_D = BD^(-1/(3μc)) * SF.gamma(z1) / SF.gamma(z2)
    return -mean_D * log(tolerance)

    ###

    # ψc = pdf_cloud_parameters(pdf, q, ρₐ, N).ψc
    # cloud_λ = pdf_cloud_parameters(pdf, q, ρₐ, N).Ec^(1 / ψc) * 1e3 # converting to m
    # cloud_problem(x) =
    #     tolerance -
    #     exp(-exp(x)^ψc * cloud_λ^ψc) * (
    #         1 +
    #         exp(x)^ψc * cloud_λ^ψc +
    #         1 / 2 * exp(x)^(2 * ψc) * cloud_λ^(2 * ψc)
    #     )
    # guess =
    #     log(0.5) +
    #     (log(0.00025) - log(0.5)) / (log(1e12) - log(1e2)) *
    #     (log(cloud_λ^3) - log(10^2))
    # log_cloud_x =
    #     RS.find_zero(
    #         cloud_problem,
    #         RS.NewtonsMethodAD(guess),
    #         RS.CompactSolution(),
    #         RS.RelativeSolutionTolerance(eps(FT)),
    #         5,
    #     ).root
    # return exp(log_cloud_x)
end

"""
    autoconversion(scheme, q_liq, q_rai, ρ, N_liq)

 - `acnv`, `pdf_c` - structs with autoconversion and cloud size distribution parameters
 - `q_liq` - cloud water specific content
 - `q_rai` - rain water specific content
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between cloud droplets (autoconversion) for `scheme == SB2006Type`
"""
function autoconversion(
    acnv::CMP.AcnvSB2006{FT},
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    q_liq,
    q_rai,
    ρ,
    N_liq,
) where {FT}

    if q_liq < eps(FT) || N_liq < eps(FT)
        return LiqRaiRates{FT}()
    end

    (; kcc, x_star, ρ0, A, a, b) = acnv
    (; νc) = pdf_c

    L_liq = ρ * q_liq
    x_liq = min(x_star, L_liq / N_liq)
    q_rai = max(FT(0), q_rai)
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_au = A * τ^a * (FT(1) - τ^a)^b

    dL_rai_dt =
        kcc / 20 / x_star * (νc + 2) * (νc + 4) / (νc + 1)^2 *
        L_liq^2 *
        x_liq^2 *
        (1 + ϕ_au / (1 - τ)^2) *
        ρ0 / ρ
    dN_rai_dt = dL_rai_dt / x_star
    dL_liq_dt = -dL_rai_dt
    dN_liq_dt = -2 * dN_rai_dt

    return LiqRaiRates(
        dq_liq_dt = dL_liq_dt / ρ,
        dN_liq_dt = dN_liq_dt,
        dq_rai_dt = dL_rai_dt / ρ,
        dN_rai_dt = dN_rai_dt,
    )
end

"""
    accretion(scheme, q_liq, q_rai, ρ, N_liq)

 - `scheme` - type for 2-moment accretion parameterization
 - `q_liq` - cloud water specific content
 - `q_rai` - rain water specific content
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a LiqRaiRates object containing `q_liq`, `N_liq`, `q_rai`, `N_rai` tendencies due to
collisions between raindrops and cloud droplets (accretion) for `scheme == SB2006Type`
"""
function accretion((; accr)::CMP.SB2006{FT}, q_liq, q_rai, ρ, N_liq) where {FT}

    if q_liq < eps(FT) || q_rai < eps(FT) || N_liq < eps(FT)
        return LiqRaiRates{FT}()
    end

    (; kcr, τ0, ρ0, c) = accr
    L_liq = ρ * q_liq
    L_rai = ρ * q_rai
    x_liq = L_liq / N_liq
    τ = FT(1) - q_liq / (q_liq + q_rai)
    ϕ_ac = (τ / (τ + τ0))^c

    dL_rai_dt = kcr * L_liq * L_rai * ϕ_ac * sqrt(ρ0 / ρ)
    dN_rai_dt = zero(N_liq)
    dL_liq_dt = -dL_rai_dt
    dN_liq_dt = dL_liq_dt / x_liq

    return LiqRaiRates(
        dq_liq_dt = dL_liq_dt / ρ,
        dN_liq_dt = dN_liq_dt,
        dq_rai_dt = dL_rai_dt / ρ,
        dN_rai_dt = dN_rai_dt,
    )
end

"""
    liquid_self_collection(scheme, q_liq, ρ, dN_liq_dt_au)

 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_liq` - cloud water specific content
 - `ρ` - air density
 - `dN_liq_dt_au` - rate of change of cloud droplets number density due to autoconversion

Returns the cloud droplets number density tendency due to collisions of cloud droplets
that produce larger cloud droplets (self-collection) for `scheme == SB2006Type`
"""
function liquid_self_collection(
    acnv::CMP.AcnvSB2006{FT},
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    q_liq::FT,
    ρ::FT,
    dN_liq_dt_au::FT,
) where {FT}

    if q_liq < eps(FT)
        return FT(0)
    end
    (; kcc, ρ0) = acnv
    (; νc) = pdf_c

    L_liq = ρ * q_liq

    # Eq. (9) from Seifert and Beheng 2006
    dN_liq_dt_sc =
        -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * L_liq^2 - dN_liq_dt_au

    return dN_liq_dt_sc
end

"""
    autoconversion_and_liquid_self_collection(scheme, q_liq, q_rai, ρ, N_liq)

Compute autoconversion and liquid self-collection rates for the [`CMP.SB2006`](@ref) scheme.

# Arguments
 - `scheme`: type for 2-moment rain autoconversion parameterization
 - `q_liq`: cloud water specific content
 - `q_rai`: rain water specific content
 - `ρ`: air density
 - `N_liq`: cloud droplet number density

# Returns
`NamedTuple` containing a `LiqRaiRates` object for the autoconversion rate and the liquid self-collection rate.
"""
function autoconversion_and_liquid_self_collection(
    (; acnv, pdf_c)::CMP.SB2006{FT},
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT}

    au = autoconversion(acnv, pdf_c, q_liq, q_rai, ρ, N_liq)
    sc = liquid_self_collection(acnv, pdf_c, q_liq, ρ, au.dN_liq_dt)

    return (; au, sc)
end

"""
    rain_self_collection(scheme, q_rai, ρ, N_rai)

# Arguments
 - `scheme`: type for 2-moment rain self-collection parameterization
 - `q_rai`: rain water specific content
 - `ρ`: air density
 - `N_rai`: raindrops number density

# Returns
The raindrops number density tendency due to collisions of raindrops that produce larger raindrops (self-collection).
"""
function rain_self_collection(
    pdf::CMP.RainParticlePDF_SB2006{FT},
    self::CMP.SelfColSB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}

    if q_rai < eps(FT) || N_rai < eps(FT)
        return FT(0)
    end

    (; krr, κrr, d) = self
    (; ρ0) = pdf

    L_rai = ρ * q_rai
    (; Dr_mean) = pdf_rain_parameters(pdf, q_rai, ρ, N_rai)
    # Note: `λᵣ` from Eq. (11) in the paper is related to `Dr_mean` by `λᵣ ≡ 1/Dr_mean`
    dN_rai_dt_sc = -krr * N_rai * L_rai * sqrt(ρ0 / ρ) * (1 + κrr * Dr_mean)^d  # Eq. (11)

    return dN_rai_dt_sc
end

"""
    rain_breakup(scheme, q_rai, ρ, dN_rai_dt_sc)

 - `scheme` - type for 2-moment liquid self-collection parameterization
 - `q_rai` - rain water specific content
 - `ρ` - air density
 - `N_rai` - raindrops number density
 - `dN_rai_dt_sc` - rate of change of raindrops number density due to self-collection

Returns the raindrops number density tendency due to breakup of raindrops
that produce smaller raindrops for `scheme == SB2006Type`
"""
function rain_breakup(
    pdf::CMP.RainParticlePDF_SB2006{FT},
    brek::CMP.BreakupSB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
    dN_rai_dt_sc::FT,
) where {FT}

    if q_rai < eps(FT) || N_rai < eps(FT)
        return FT(0)
    end
    (; Deq, Dr_th, kbr, κbr) = brek
    ρw = pdf.ρw
    (; xr_mean) = pdf_rain_parameters(pdf, q_rai, ρ, N_rai)
    Dr = (xr_mean * 6 / FT(π) / ρw)^FT(1 / 3)
    ΔD = Dr - Deq
    phi_br =
        (Dr < Dr_th) ? FT(-1) : ((ΔD <= 0) ? kbr * ΔD : 2 * (exp(κbr * ΔD) - 1))
    dN_rai_dt_br = -(phi_br + 1) * dN_rai_dt_sc

    return dN_rai_dt_br
end

"""
    rain_self_collection_and_breakup(SB2006, q_rai, ρ, N_rai)

 - `SB2006` - a struct with SB2006 parameters for raindrops size
    distribution, self collection, and breakup
 - `q_rai` - rain water specific content
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns a named tupple containing the raindrops self-collection and breakup rates
for `scheme == SB2006Type`
"""
function rain_self_collection_and_breakup(
    (; pdf_r, self, brek)::CMP.SB2006, q_rai, ρ, N_rai,
)

    sc = rain_self_collection(pdf_r, self, q_rai, ρ, N_rai)
    br = rain_breakup(pdf_r, brek, q_rai, ρ, N_rai, sc)

    return (; sc, br)
end

"""
    rain_terminal_velocity(SB2006, vel, q_rai, ρ, N_rai)

 - `SB2006` - a struct with SB2006 rain size distribution parameters
 - `vel` - a struct with terminal velocity parameters
 - `q_rai` - rain water specific content [kg/kg]
 - `ρ` - air density [kg/m³]
 - `N_rai` - raindrops number density [1/m³]

Returns a tuple containing the number and mass weigthed mean fall velocities of raindrops in [m/s].
Assuming an exponential size distribution from Seifert and Beheng 2006 for `scheme == SB2006Type`
Fall velocity of individual rain drops is parameterized:
 - assuming an empirical relation similar to Rogers (1993) for `velo_scheme == SB2006VelType`
 - following Chen et. al 2022, DOI: 10.1016/j.atmosres.2022.106171 for `velo_scheme == Chen2022Type`
"""
function rain_terminal_velocity(
    (; pdf_r)::CMP.SB2006{FT},
    (; ρ0, aR, bR, cR)::CMP.SB2006VelType{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}
    # TODO: Input argument list needs to be redesigned

    (; Dr_mean) = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai)
    _pa0, _pb0, _pa1, _pb1 =
        _sb_rain_terminal_velocity_helper(pdf_r, 1/Dr_mean, aR, bR, cR)

    vt0 =
        N_rai < eps(FT) ? FT(0) :
        max(FT(0), sqrt(ρ0 / ρ) * (aR * _pa0 - bR * _pb0 / (1 + cR * Dr_mean)))
    vt1 =
        q_rai < eps(FT) ? FT(0) :
        max(FT(0), sqrt(ρ0 / ρ) * (aR * _pa1 - bR * _pb1 / (1 + cR * Dr_mean)^FT(4)))
    return (vt0, vt1)
end
function rain_terminal_velocity(
    (; pdf_r)::CMP.SB2006{FT},
    vel::CMP.Chen2022VelTypeRain{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}
    # coefficients from Table B1 from Chen et. al. 2022
    aiu, bi, ciu = CO.Chen2022_vel_coeffs_B1(vel, ρ)
    # size distribution parameter
    (; Dr_mean) = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai)

    # eq 20 from Chen et al 2022
    vt0 = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, Dr_mean, 0))
    vt3 = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, Dr_mean, 3))

    vt0 = N_rai < eps(FT) ? FT(0) : max(FT(0), vt0)
    vt3 = q_rai < eps(FT) ? FT(0) : max(FT(0), vt3)
    # It should be (ϕ^κ * vt0, ϕ^κ * vt3), but for rain drops ϕ = 1 and κ = 0
    return (vt0, vt3)
end
function _sb_rain_terminal_velocity_helper(
    pdf_r::CMP.RainParticlePDF_SB2006_limited{FT},
    λr,
    aR,
    bR,
    cR,
) where {FT}
    return (FT(1), FT(1), FT(1), FT(1))
end
function _sb_rain_terminal_velocity_helper(
    pdf_r::CMP.RainParticlePDF_SB2006_notlimited{FT},
    λr,
    aR,
    bR,
    cR,
) where {FT}
    # Integrate velocity of particles over a range of r with
    # positive terminal velocity (v = aR - bR exp(-lambda D))
    _rc = -1 / (2 * cR) * log(aR / bR)
    _Γ_1(t) = exp(-t)
    _Γ_4(t) = (t^3 + 3 * t^2 + 6 * t + 6) * exp(-t)
    _pa0::FT = _Γ_1(2 * _rc * λr)
    _pb0::FT = _Γ_1(2 * _rc * (λr + cR))
    _pa1::FT = _Γ_4(2 * _rc * λr) / FT(6)
    _pb1::FT = _Γ_4(2 * _rc * (λr + cR)) / FT(6)
    return (_pa0, _pb0, _pa1, _pb1)
end

"""
    Returns the approximation of an incomplete gamma function
    for a ∈ {-1.0, -0.101}, and x in [0.067 1.82]
"""
function Γ_incl(a::FT, x::FT) where {FT}
    #return exp(-x) / ((FT(1.5) - FT(0.54) * a) * x^(FT(0.46) - FT(0.75) * a))
    return exp(-x) / (
        (FT(0.33) - FT(0.7) * a) * x^(FT(0.08) - FT(0.93) * a) +
        (FT(1.34) - FT(0.1) * a) * x^(FT(0.8) - a)
    )
end

"""
    rain_evaporation(evap, aps, tps, q, q_rai, ρ, N_rai, T)

 - `evap` - evaporation parameterization scheme
 - `aps` - air properties
 - `tps` - thermodynamics parameters
 - `q` - phase partition
 - `q_rai` - rain specific content
 - `ρ` - air density
 - `N_rai` - raindrops number density
 - `T` - air temperature

Returns a named tuple containing the tendency of raindrops number density and rain water
specific content due to rain rain_evaporation, assuming a power law velocity relation for
fall velocity of individual drops and an exponential size distribution, for `scheme == SB2006Type`
"""
function rain_evaporation(
    (; pdf_r, evap)::CMP.SB2006{FT},
    aps::CMP.AirProperties{FT},
    tps::TDP.ThermodynamicsParameters{FT},
    q::TD.PhasePartition{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
    T::FT,
) where {FT}

    evap_rate_0 = FT(0)
    evap_rate_1 = FT(0)
    S = TD.supersaturation(tps, q, ρ, T, TD.Liquid())

    if ((q_rai > eps(FT) || N_rai > eps(FT)) && S < FT(0))

        (; ν_air, D_vapor) = aps
        (; av, bv, α, β, ρ0) = evap
        x_star = pdf_r.xr_min
        ρw = pdf_r.ρw
        G = CO.G_func(aps, tps, T, TD.Liquid())

        (; xr_mean) = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai)
        Dr = (FT(6) / FT(π) / ρw)^FT(1 / 3) * xr_mean^FT(1 / 3)

        t_star = (FT(6) * x_star / xr_mean)^FT(1 / 3)
        a_vent_0 = av * Γ_incl(FT(-1), t_star) / FT(6)^FT(-2 / 3)
        b_vent_0 =
            bv * Γ_incl(-FT(0.5) + FT(1.5) * β, t_star) /
            FT(6)^FT(β / 2 - FT(0.5))

        a_vent_1 = av * SF.gamma(FT(2)) / FT(6)^FT(1 / 3)
        b_vent_1 =
            bv * SF.gamma(FT(5 / 2) + FT(3 / 2) * β) / FT(6)^FT(β / 2 + 1 / 2)

        N_Re = α * xr_mean^β * sqrt(ρ0 / ρ) * Dr / ν_air
        Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
        Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

        evap_rate_0 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv0 / xr_mean)
        evap_rate_1 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)

        # When xr = 0 evap_rate_0 becomes NaN. We replace NaN with 0 which is the limit of
        # evap_rate_0 for xr -> 0.
        evap_rate_0 =
            N_rai < eps(FT) || xr_mean / x_star < eps(FT) ? FT(0) : evap_rate_0
        evap_rate_1 = q_rai < eps(FT) ? FT(0) : evap_rate_1
    end

    return (; evap_rate_0, evap_rate_1)
end

# Additional double moment autoconversion and accretion parametrizations:
# - Khairoutdinov and Kogan (2000)
# - Beheng (1994)
# - Tripoli and Cotton (1980)
# - Liu and Daum (2004)
# - variable timescale autoconversion Azimi (2023)

"""
    conv_q_liq_to_q_rai(acnv, q_liq, ρ, N_d; smooth_transition)

 - `acnv` - 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific content
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
function conv_q_liq_to_q_rai((; acnv)::CMP.KK2000{FT}, q_liq, ρ, N_d) where {FT}
    q_liq = max(0, q_liq)
    (; A, a, b, c) = acnv
    return A * q_liq^a * N_d^b * ρ^c
end
function conv_q_liq_to_q_rai(
    (; acnv)::CMP.B1994{FT},
    q_liq,
    ρ,
    N_d,
    smooth_transition = false,
) where {FT}
    q_liq = max(0, q_liq)
    (; C, a, b, c, N_0, k, d_low, d_high) = acnv
    d = FT(0)
    if smooth_transition
        d_low_acnv_fraction = CO.logistic_function(N_d, N_0, k)
        d_high_acnv_fraction = FT(1) - d_low_acnv_fraction
        d = d_low_acnv_fraction * d_low + d_high_acnv_fraction * d_high
    else
        d = N_d >= N_0 ? d_low : d_high
    end
    return C * d^a * (q_liq * ρ)^b * N_d^c / ρ
end
function conv_q_liq_to_q_rai(
    (; acnv)::CMP.TC1980{FT},
    q_liq,
    ρ,
    N_d,
    smooth_transition = false,
) where {FT}
    #TODO - The original paper is actually formulated for mixing ratios, not specific contents
    q_liq = max(0, q_liq)
    (; m0_liq_coeff, me_liq, D, a, b, r_0, k) = acnv
    q_liq_threshold::FT = m0_liq_coeff * N_d / ρ * r_0^me_liq
    output =
        smooth_transition ? CO.logistic_function(q_liq, q_liq_threshold, k) :
        CO.heaviside(q_liq - q_liq_threshold)
    return D * q_liq^a * N_d^b * output
end
function conv_q_liq_to_q_rai(
    (; ρ_w, R_6C_0, E_0, k)::CMP.LD2004{FT},
    q_liq,
    ρ,
    N_d,
    smooth_transition = false,
) where {FT}
    if q_liq <= eps(FT)
        return FT(0)
    else
        # Mean volume radius in microns (assuming spherical cloud droplets)
        r_vol =
            (FT(3) * (q_liq * ρ) / FT(4) / FT(π) / ρ_w / N_d)^FT(1 / 3) *
            FT(1e6)

        # Assumed size distribution: modified gamma distribution
        β_6 = ((r_vol + FT(3)) / r_vol)^FT(1 / 3)
        E = E_0 * β_6^6
        R_6 = β_6 * r_vol
        R_6C = R_6C_0 / (q_liq * ρ)^FT(1 / 6) / R_6^FT(1 / 2)

        output =
            smooth_transition ? CO.logistic_function(R_6, R_6C, k) :
            CO.heaviside(R_6 - R_6C)
        return E * (q_liq * ρ)^3 / N_d / ρ * output
    end
end
function conv_q_liq_to_q_rai(
    (; τ, α)::CMP.VarTimescaleAcnv{FT},
    q_liq::FT,
    ρ::FT,
    N_d::FT,
) where {FT}
    return max(0, q_liq) / (τ * (N_d / 1e8)^α)
end

"""
    accretion(accretion_scheme, q_liq, q_rai, ρ)

 - `accretion_scheme` - type for 2-moment rain accretion parameterization
 - `q_liq` - cloud water specific content
 - `q_rai` - rain water specific content
 - `ρ` - air density (for `KK2000Type` and `Beheng1994Type`)

 Returns the accretion rate of rain, parametrized following
 - Khairoutdinov and Kogan (2000) for `scheme == KK2000Type`
 - Beheng (1994) for `scheme == B1994Type`
 - Tripoli and Cotton (1980) for `scheme == TC1980Type`
"""
function accretion((; accr)::CMP.KK2000{FT}, q_liq, q_rai, ρ) where {FT}
    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)
    (; A, a, b) = accr
    return A * (q_liq * q_rai)^a * ρ^b
end

function accretion((; accr)::CMP.B1994{FT}, q_liq, q_rai, ρ) where {FT}
    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_liq * ρ * q_rai
end

function accretion((; accr)::CMP.TC1980{FT}, q_liq, q_rai) where {FT}
    #TODO - The original paper is actually formulated for mixing ratios, not specific contents
    q_liq = max(0, q_liq)
    q_rai = max(0, q_rai)
    (; A) = accr
    return A * q_liq * q_rai
end

end
