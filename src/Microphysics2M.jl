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
    size_distribution,
    get_size_distribution_bound

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
    pdf_cloud_parameters(pdf_c, qₗ, ρₐ, Nₗ)

 - `pdf_c` - a struct with SB2006 cloud droplets size distribution parameters
 - `qₗ` - cloud water specific content
 - `ρₐ` - air density
 - `Nₗ` cloud droplet number concentration

    Returns the mean mass of cloud droplets xc [μg],
    a multiplier needed to convert xc to base SI units χ [-],
    and the two cloud droplet mass distribution parameters
    Ac [m^-3, μg^-3] and Bc [μg^-1],
    It also returns four cloud droplet size distribution parameters
    Cc [(1/mm3)^4], Ec [1/mm3], ϕc [-], ψc [-] that are consistent with the
    assumed cloud droplet mass distribution.
"""
function pdf_cloud_parameters(
    pdf_c::CMP.CloudParticlePDF_SB2006{FT},
    qₗ::FT,
    ρₐ::FT,
    Nₗ::FT,
) where {FT}

    (; νc, μc, ρw) = pdf_c
    χ = 9 # convert mean droplet mass from kg to μg

    if qₗ < eps(FT) || Nₗ < eps(FT)
        return (
            xc = FT(0),
            χ,
            Ac = FT(0),
            Bc = FT(0),
            Cc = FT(0),
            Ec = FT(0),
            ϕc = FT(0),
            ψc = FT(0),
        )
    else
        xc = qₗ * ρₐ / Nₗ * FT(10)^χ                                             # μg
        Bc = (xc * SF.gamma((νc + 1) / μc) / SF.gamma((νc + 2) / μc))^(-μc)      # 1/μg
        Ac = μc * Nₗ * Bc^((νc + 1) / μc) / SF.gamma((νc + 1) / μc)              # 1/m3 1/μg3

        _Ac = Ac * FT(1e-9)                                                      # 1/mm3 1/μg3
        Cc = _Ac * ((FT(π) * ρw * FT(10)^(χ - 9)) / FT(2))^(νc + 1) / FT(3)^νc   # (1/mm3)^4
        Ec = Bc * (FT(π / 6) * ρw * FT(10)^(χ - 9))^μc                           # 1/mm3
        ϕc = 3 * νc + 2
        ψc = 3 * μc

        return (; xc, χ, Ac, Bc, Cc, Ec, ϕc, ψc)
    end
end

"""
    pdf_rain_parameters(pdf_r, qᵣ, ρₐ, Nᵣ)

 - `pdf_r` - a struct with SB2006 raindrops size distribution parameters
 - `qᵣ` - rain water specific content
 - `ρₐ` - air density
 - `Nᵣ` rain drop number concentration

    Returns the parameters of the rain drop diameter distribution λr [1/m], N₀r [1/m4],
    optionally limited within prescribed ranges.
    Also returns the mean mass of rain drops xr [kg],
    and the two rain drop mass distribution parameters Ar [1/m3 * (1/kg)^(νr+1)], Br [(1/kg)^μr].
"""
function pdf_rain_parameters(
    pdf_r::CMP.RainParticlePDF_SB2006{FT},
    qᵣ::FT,
    ρₐ::FT,
    Nᵣ::FT,
) where {FT}
    (; νr, μr, ρw) = pdf_r

    if qᵣ < eps(FT) || Nᵣ < eps(FT)
        return (λr = FT(0), xr = FT(0), Ar = FT(0), Br = FT(0))
    else
        λr = (Nᵣ * FT(π) * ρw / ρₐ / qᵣ)^FT(1 / 3)   #1/m
        N₀r = Nᵣ * λr                                #1/m4

        xr = qᵣ * ρₐ / Nᵣ
        Br = (xr * SF.gamma((νr + 1) / μr) / SF.gamma((νr + 2) / μr))^(-μr)
        Ar = μr * Nᵣ * Br^((νr + 1) / μr) / SF.gamma(FT(νr + 1) / μr)

        return (; λr, N₀r, xr, Ar, Br)
    end
end
function pdf_rain_parameters(
    pdf_r::CMP.RainParticlePDF_SB2006_limited{FT},
    qᵣ::FT,
    ρₐ::FT,
    Nᵣ::FT,
) where {FT}
    (; νr, μr, xr_min, xr_max, N0_min, N0_max, λ_min, λ_max, ρw) = pdf_r

    qᵣ = qᵣ > FT(0) ? qᵣ : FT(0)

    Lᵣ = ρₐ * qᵣ
    xr_0 = Nᵣ > eps(FT) ? Lᵣ / Nᵣ : FT(0)
    xr_hat = max(xr_min, min(xr_max, xr_0))
    N0 = max(N0_min, min(N0_max, Nᵣ * (FT(π) * ρw / xr_hat)^FT(1 / 3)))
    λr = max(λ_min, min(λ_max, (FT(π) * ρw * N0 / Lᵣ)^FT(1 / 4)))
    xr = max(xr_min, min(xr_max, Lᵣ * λr / N0))
    N₀r = Nᵣ * λr

    Br = (SF.gamma(FT(νr + 1) / μr) * xr / SF.gamma(FT(νr + 2) / μr))^(-μr)
    Ar = μr * Nᵣ * Br^(FT(νr + 1) / μr) / SF.gamma(FT(νr + 1) / μr)

    return (; λr, N₀r, xr, Ar, Br)
end

"""
    size_distribution(pdf, D, q, ρ, N)

 - pdf - struct containing size distribution parameters of cloud or rain
 - D - particle size (i.e. maximum dimension of particle)
 - q - cloud or rain water specific content
 - ρₐ - density of air
 - N - cloud or rain water number concentration

    Returns the size distribution value for rain or cloud particles for a given
    particle size
"""
function size_distribution(
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
    D::FT,
    q::FT,
    ρₐ::FT,
    N::FT,
) where {FT}
    (; N₀r, λr) = pdf_rain_parameters(pdf, q, ρₐ, N)
    return N₀r * exp(-λr * D)
end
function size_distribution(
    pdf::CMP.CloudParticlePDF_SB2006{FT},
    D::FT,
    q::FT,
    ρₐ::FT,
    N::FT,
) where {FT}
    (; Cc, Ec, ϕc, ψc) = pdf_cloud_parameters(pdf, q, ρₐ, N)
    # Convert values from mm to m assuming
    Ec *= FT(10^(3 * ψc))
    # TODO - overflow for Float32. We only need it in P3 scheme and it only
    # works for Float64 right now. The solution is to non-dimensionalize
    # all distributions by a reference mass/size, similar to the 1M scheme.
    Cc *= Float64(10^((ϕc + 4) * 3))
    return Cc * D^ϕc * exp(-Ec * D^ψc)
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
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
    q::FT,
    N::FT,
    ρₐ::FT,
    tolerance::FT,
) where {FT}
    (; λr) = pdf_rain_parameters(pdf, q, ρₐ, N)
    return -1 / λr * log(tolerance)
end
function get_size_distribution_bound(
    pdf::CMP.CloudParticlePDF_SB2006{FT},
    q::FT,
    N::FT,
    ρₐ::FT,
    tolerance::FT,
) where {FT}
    ψc = pdf_cloud_parameters(pdf, q, ρₐ, N).ψc
    cloud_λ = pdf_cloud_parameters(pdf, q, ρₐ, N).Ec^(1 / ψc) * 1e3 # converting to m
    cloud_problem(x) =
        tolerance -
        exp(-exp(x)^ψc * cloud_λ^ψc) * (
            1 +
            exp(x)^ψc * cloud_λ^ψc +
            1 / 2 * exp(x)^(2 * ψc) * cloud_λ^(2 * ψc)
        )
    guess =
        log(0.5) +
        (log(0.00025) - log(0.5)) / (log(1e12) - log(1e2)) *
        (log(cloud_λ^3) - log(10^2))
    log_cloud_x =
        RS.find_zero(
            cloud_problem,
            RS.NewtonsMethodAD(guess),
            RS.CompactSolution(),
            RS.RelativeSolutionTolerance(eps(FT)),
            5,
        ).root
    return exp(log_cloud_x)
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

    dN_liq_dt_sc =
        -kcc * (νc + 2) / (νc + 1) * (ρ0 / ρ) * L_liq^2 - dN_liq_dt_au

    return dN_liq_dt_sc
end

"""
    autoconversion_and_liquid_self_collection(scheme, q_liq, q_rai, ρ, N_liq)

 - `scheme` - type for 2-moment rain autoconversion parameterization
 - `q_liq` - cloud water specific content
 - `q_rai` - rain water specific content
 - `ρ` - air density
 - `N_liq` - cloud droplet number density

Returns a named tupple containing a LiqRaiRates object for the autoconversion rate and
the liquid self-collection rate for `scheme == SB2006Type`
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

 - `scheme` - type for 2-moment rain self-collection parameterization
 - `q_rai` - rain water specific content
 - `ρ` - air density
 - `N_rai` - raindrops number density

Returns the raindrops number density tendency due to collisions of raindrops
that produce larger raindrops (self-collection) for `scheme == SB2006Type`
"""
function rain_self_collection(
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
    self::CMP.SelfColSB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}

    if q_rai < eps(FT) || N_rai < eps(FT)
        return FT(0)
    end

    (; krr, κrr, d) = self
    (; ρ0, ρw) = pdf

    L_rai = ρ * q_rai
    Br = pdf_rain_parameters(pdf, q_rai, ρ, N_rai).Br
    dN_rai_dt_sc = -krr * N_rai * L_rai * sqrt(ρ0 / ρ) * (1 + κrr / Br)^d

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
    pdf::Union{
        CMP.RainParticlePDF_SB2006{FT},
        CMP.RainParticlePDF_SB2006_limited{FT},
    },
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
    xr = pdf_rain_parameters(pdf, q_rai, ρ, N_rai).xr
    Dr = (xr * 6 / FT(π) / ρw)^FT(1 / 3)
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
    (; pdf_r, self, brek)::CMP.SB2006{FT},
    q_rai::FT,
    ρ::FT,
    N_rai::FT,
) where {FT}

    sc = rain_self_collection(pdf_r, self, q_rai, ρ, N_rai)
    br = rain_breakup(pdf_r, brek, q_rai, ρ, N_rai, sc)

    return (; sc, br)
end

"""
    rain_particle_terminal_velocity(D, Chen2022, ρₐ)

 - D - maximum particle dimension
 - Chen2022 - a struct with terminal velocity parameters from Chen 2022
 - ρₐ - air density

Returns the terminal velocity of an individual rain drop as a function
of its size (maximum dimension) following Chen 2022 velocity parametrization.
Needed for numerical integrals in the P3 scheme.
"""
function rain_particle_terminal_velocity(
    D::FT,
    Chen2022::CMP.Chen2022VelTypeRain,
    ρₐ::FT,
) where {FT}
    (ai, bi, ci) = CO.Chen2022_vel_coeffs_B1(Chen2022, ρₐ)

    v = sum(@. sum(ai * D^bi * exp(-ci * D)))
    return v
end

"""
    rain_terminal_velocity(SB2006, vel, q_rai, ρ, N_rai)

 - `SB2006` - a struct with SB2006 rain size distribution parameters
 - `vel` - a struct with terminal velocity parameters
 - `q_rai` - rain water specific content [kg/kg]
 - `ρ` - air density [kg/m^3]
 - `N_rai` - raindrops number density [1/m^3]

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

    λr = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai).λr
    _pa0, _pb0, _pa1, _pb1 =
        _sb_rain_terminal_velocity_helper(pdf_r, λr, aR, bR, cR)

    vt0 =
        N_rai < eps(FT) ? FT(0) :
        max(FT(0), sqrt(ρ0 / ρ) * (aR * _pa0 - bR * _pb0 / (1 + cR / λr)))
    vt1 =
        q_rai < eps(FT) ? FT(0) :
        max(FT(0), sqrt(ρ0 / ρ) * (aR * _pa1 - bR * _pb1 / (1 + cR / λr)^FT(4)))
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
    λ = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai).λr

    # eq 20 from Chen et al 2022
    vt0 = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, λ, 0))
    vt3 = sum(CO.Chen2022_exponential_pdf.(aiu, bi, ciu, λ, 3))

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
    pdf_r::CMP.RainParticlePDF_SB2006{FT},
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

        xr = pdf_rain_parameters(pdf_r, q_rai, ρ, N_rai).xr
        Dr = (FT(6) / FT(π) / ρw)^FT(1 / 3) * xr^FT(1 / 3)

        t_star = (FT(6) * x_star / xr)^FT(1 / 3)
        a_vent_0 = av * Γ_incl(FT(-1), t_star) / FT(6)^FT(-2 / 3)
        b_vent_0 =
            bv * Γ_incl(-FT(0.5) + FT(1.5) * β, t_star) /
            FT(6)^FT(β / 2 - FT(0.5))

        a_vent_1 = av * SF.gamma(FT(2)) / FT(6)^FT(1 / 3)
        b_vent_1 =
            bv * SF.gamma(FT(5 / 2) + FT(3 / 2) * β) / FT(6)^FT(β / 2 + 1 / 2)

        N_Re = α * xr^β * sqrt(ρ0 / ρ) * Dr / ν_air
        Fv0 = a_vent_0 + b_vent_0 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)
        Fv1 = a_vent_1 + b_vent_1 * (ν_air / D_vapor)^FT(1 / 3) * sqrt(N_Re)

        evap_rate_0 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv0 / xr)
        evap_rate_1 = min(FT(0), FT(2) * FT(π) * G * S * N_rai * Dr * Fv1 / ρ)

        # When xr = 0 evap_rate_0 becomes NaN. We replace NaN with 0 which is the limit of
        # evap_rate_0 for xr -> 0.
        evap_rate_0 =
            N_rai < eps(FT) || xr / x_star < eps(FT) ? FT(0) : evap_rate_0
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
