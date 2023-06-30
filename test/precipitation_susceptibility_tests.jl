import Test as TT

import CloudMicrophysics as CM
import CLIMAParameters as CP
import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.CommonTypes as CMT
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

include(joinpath(pkgdir(CM), "test", "create_parameters.jl"))

const APS = CMP.AbstractCloudMicrophysicsParameters

FT = Float64

"""
Logarithmic derivative of universal function for autoconversion as described
    in Glassmeier & Lohmann, which is (1 + Phi(τ)/(1 - τ^2)) for Phi(τ)
    from SB2006.
"""
function d_ln_phi_au_d_ln_τ(
    param_set::APS,
    scheme::CMT.SB2006Type,
    τ::FT,
) where {FT <: Real}
    A::FT = CMP.A_phi_au_SB2006(param_set)
    a::FT = CMP.a_phi_au_SB2006(param_set)
    b::FT = CMP.b_phi_au_SB2006(param_set)
    return -(
        A *
        τ^a *
        (1 - τ^a)^(b - 1) *
        (a * (τ - 1) * ((b + 1) * τ^a - 1) - 2 * τ * (τ^a - 1))
    ) / (A * (τ - 1) * τ^a * (1 - τ^a)^b + (τ - 1)^3)
end

"""
Logarithmic derivative of universal function for accretion as described
    in Glassmeier & Lohmann
"""
function d_ln_phi_acc_d_ln_τ(
    param_set::APS,
    scheme::CMT.SB2006Type,
    τ::FT,
) where {FT <: Real}
    τ0::FT = CMP.τ0_phi_ac_SB2006(param_set)
    c::FT = CMP.c_phi_ac_SB2006(param_set)
    return (c * τ0) / (τ + τ0)
end


toml_dict = CP.create_toml_dict(FT; dict_type = "alias")
prs = cloud_microphysics_parameters(toml_dict)

TT.@testset "precipitation_susceptibility_SB2006" begin
    scheme = CMT.SB2006Type()

    q_liq = FT(0.5e-3)
    N_liq = FT(1e8)
    q_rai = FT(1e-5)
    ρ = FT(1)

    τ = FT(1) - q_liq / (q_liq + q_rai)

    aut_rates = CMPS.precipitation_susceptibility_autoconversion(
        prs,
        scheme,
        q_liq,
        q_rai,
        ρ,
        N_liq,
    )

    acc_rates = CMPS.precipitation_susceptibility_accretion(
        prs,
        scheme,
        q_liq,
        q_rai,
        ρ,
        N_liq,
    )

    TT.@test aut_rates.d_ln_pp_d_ln_N_liq ≈ -2
    TT.@test aut_rates.d_ln_pp_d_ln_q_liq ≈
             4 - (1 - τ) * d_ln_phi_au_d_ln_τ(prs, scheme, τ)
    TT.@test aut_rates.d_ln_pp_d_ln_q_rai ≈
             (1 - τ) * d_ln_phi_au_d_ln_τ(prs, scheme, τ)

    TT.@test acc_rates.d_ln_pp_d_ln_q_liq ≈
             1 - (1 - τ) * d_ln_phi_acc_d_ln_τ(prs, scheme, τ)
    TT.@test acc_rates.d_ln_pp_d_ln_q_rai ≈
             1 + (1 - τ) * d_ln_phi_acc_d_ln_τ(prs, scheme, τ)
end
