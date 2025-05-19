"""
    Cloud diagnostics

 - radar reflectivity (1-moment and 2-moment)
 - effective radius  (1-moment and 2-moment)
"""
module CloudDiagnostics

import SpecialFunctions as SF

import ..Parameters as CMP
import ..Microphysics1M as CM1
import ..Microphysics2M as CM2

"""
    radar_reflectivity_1M(precip, q, ρ)

  - `precip` - struct with 1-moment microphysics rain free parameters
  - `q` - specific content of rain
  - `ρ` - air density

Returns logarithmic radar reflectivity for the 1-moment microphysics
based on the assumed rain particle size distribution.
Normalized by the reflectivty of 1 millimiter drop in a volume of 1m3.
The values are clipped at -150 dBZ.
"""
function radar_reflectivity_1M(
    (; pdf, mass)::CMP.Rain{FT},
    q::FT,
    ρ::FT,
) where {FT}

    # change units for accuracy
    n0 = CM1.get_n0(pdf) * FT(1e-12)
    λ_inv = CM1.lambda_inverse(pdf, mass, q, ρ) / FT(1e-3)

    Z = 720 * n0 * λ_inv^7
    log_10_Z₀ = FT(-18)
    log_Z = FT(10) * (log10(Z) - log_10_Z₀ - FT(9))

    return max(FT(-150), log_Z)
end

"""
    radar_reflectivity_2M(structs, q_liq, q_rai, N_liq, N_rai, ρ_air)

    - `structs` - structs microphysics 2-moment with SB2006 cloud droplets
                  and raindrops size distributions parameters
    - `q_liq` - cloud water specific content
    - `q_rai` - rain water specific content
    - `N_liq` - cloud droplet number density
    - `N_rai` - rain droplet number density
    - `ρ_air` - air density

Returns logarithmic radar reflectivity for the 2-moment microphysics SB2006
based on the assumed cloud and rain particle size distribuions.
Normalized by the reflectivty of 1 millimiter drop in a volume of 1m3.
The values are clipped at -150 dBZ.
"""
function radar_reflectivity_2M(
    (; pdf_c, pdf_r)::CMP.SB2006{FT},
    q_liq::FT,
    q_rai::FT,
    N_liq::FT,
    N_rai::FT,
    ρ_air::FT,
) where {FT}
    # free parameters
    (; νc, μc) = pdf_c
    (; νr, μr, ρw) = pdf_r
    C = FT(4 / 3 * π * ρw)
    log_10_Z₀ = -18

    # rain size distribution parameters
    (; Ar, Br) = CM2.pdf_rain_parameters_mass(pdf_r, q_rai, ρ_air, N_rai)

    # cloud size distribution parameters
    L_c = ρ_air * q_liq
    logAc, logBc = CM2.log_pdf_cloud_parameters_mass(L_c, N_liq, νc, μc)
    Ac = exp(logAc)
    Bc = exp(logBc)

    function radar_reflectivity(A, B, ν, μ)
        n = 2  # 2nd moment in mass
        z = (ν + 1 + n) / μ
        M_n = A * C^(ν+1) * (B * C^μ)^(-z) * SF.gamma(z) / μ
        return M_n
        # TODO: Consider calculating log(M_n)
        # logM_n = logA + (ν+1) * logC - z * (logB + μ * logC) + SF.loggamma(z) - log(μ)
    end

    Zc = radar_reflectivity(Ac, Bc, νc, μc)
    # TODO: Rewrite in terms of diameter
    Zr = radar_reflectivity(Ar, Br, νr, μr)

    return max(FT(-150), 10 * (log10(max(FT(0), Zc + Zr)) - log_10_Z₀))
end

"""
    effective_radius_2M(structs, q_liq, q_rai, N_liq, N_rai, ρ_air)

    - `structs` - structs with SB2006 cloud droplets and raindrops
                size distribution parameters
    - `q_liq` - cloud water specific content
    - `q_rai` - rain water specific content
    - `N_liq` - cloud droplet number density
    - `N_rai` - rain droplet number density
    - `ρ_air` - air density

Returns effective radius for the 2-moment microphysics scheme.
Computed based on the assumed cloud and rain particle size distributions.
"""
function effective_radius_2M(
    (; pdf_c, pdf_r)::CMP.SB2006{FT},
    q_liq::FT,
    q_rai::FT,
    N_liq::FT,
    N_rai::FT,
    ρ_air::FT,
) where {FT}
    # free parameters
    (; νc, μc) = pdf_c
    (; νr, μr, ρw) = pdf_r
    C = FT(4 / 3 * π * ρw)
    # rain size distribution parameters
    (; Ar, Br) = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ_air, N_rai)

    # cloud size distribution parameters (χ converts from μg to base SI)
    (; Ac, Bc, χ) = CM2.pdf_cloud_parameters(pdf_c, q_liq, ρ_air, N_liq)

    M3_c =
        Bc == 0 ? FT(0) :
        Ac *
        C^(νc + 1) *
        (Bc * C^μc)^(-(2 + νc) / μc) *
        SF.gamma((2 + νc) / μc) / μc / FT(10)^χ
    M3_r =
        Br == 0 ? FT(0) :
        Ar *
        C^(νr + 1) *
        (Br * C^μr)^(-(2 + νr) / μr) *
        SF.gamma((2 + νr) / μr) / μr

    M2_c =
        Bc == 0 ? FT(0) :
        Ac *
        C^(νc + 1) *
        (Bc * C^μc)^(-(5 + 3 * νc) / (3 * μc)) *
        SF.gamma((5 + 3 * νc) / (3 * μc)) / μc / 10^(χ * FT(2) / 3)
    M2_r =
        Br == 0 ? FT(0) :
        Ar *
        C^(νr + 1) *
        (Br * C^μr)^(-(5 + 3 * νr) / (3 * μr)) *
        SF.gamma((5 + 3 * νr) / (3 * μr)) / μr

    return M2_c + M2_r <= eps(FT) ? FT(0) : (M3_c + M3_r) / (M2_c + M2_r)
end

"""
    effective_radius_Liu_Hallet_97(q_liq, q_rai, N_liq, N_rai, ρ_air, ρ_w)

    - `q_liq` - cloud water specific content
    - `q_rai` - rain water specific content
    - `N_liq` - cloud droplet number density
    - `N_rai` - rain droplet number density
    - `ρ_air` - air density
    - `ρ_w` - water density

Returns effective radius using the "1/3" power law from Liu and Hallett (1997).
If not provided by the user, it is assumed that there is no rain present and that
the cloud droplet number concentration is 100 1/cm3.
"""
function effective_radius_Liu_Hallet_97(
    (; ρw)::Union{CMP.WaterProperties{FT}, CMP.CloudLiquid{FT}},
    ρ_air::FT,
    q_liq::FT,
    N_liq::FT,
    q_rai::FT,
    N_rai::FT,
) where {FT}

    k = FT(0.8)
    r_vol =
        ((N_liq + N_rai) < eps(FT)) ? FT(0) :
        (
            (FT(3) * (q_liq + q_rai) * ρ_air) /
            (FT(4) * π * ρw * (N_liq + N_rai))
        )^FT(1 / 3)

    return r_vol / k^FT(1 / 3)
end
function effective_radius_Liu_Hallet_97(
    wtr::Union{CMP.WaterProperties{FT}, CMP.CloudLiquid{FT}},
    ρ_air::FT,
    q_liq::FT,
) where {FT}
    return effective_radius_Liu_Hallet_97(
        wtr::Union{CMP.WaterProperties{FT}, CMP.CloudLiquid{FT}},
        ρ_air::FT,
        q_liq::FT,
        FT(100),
        FT(0),
        FT(0),
    )
end

"""
    effective_radius_const(cloud_params)

  - cloud_params - a struct with cloud liquid or cloud ice parameters

Returns a constant assumed effective radius for clouds
"""
function effective_radius_const(cloud_params::CMP.CloudLiquid{FT}) where {FT}
    return cloud_params.r_eff
end
function effective_radius_const(cloud_params::CMP.CloudIce{FT}) where {FT}
    return cloud_params.r_eff
end

end # end module
