import Test as TT

import CloudMicrophysics.Parameters as CMP
import CloudMicrophysics.PrecipitationSusceptibility as CMPS

FT = Float64

"""
Logarithmic derivative of universal function for autoconversion as described
    in Glassmeier & Lohmann, which is (1 + Phi(τ)/(1 - τ^2)) for Phi(τ)
    from SB2006.
"""
d_ln_phi_au_d_ln_τ(
    (; A, a, b)::CMP.AcnvSB2006{FT},
    τ::FT,
) where {FT <: AbstractFloat} =
    -(
        A *
        τ^a *
        (1 - τ^a)^(b - 1) *
        (a * (τ - 1) * ((b + 1) * τ^a - 1) - 2 * τ * (τ^a - 1))
    ) / (A * (τ - 1) * τ^a * (1 - τ^a)^b + (τ - 1)^3)

"""
Logarithmic derivative of universal function for accretion as described
    in Glassmeier & Lohmann
"""
d_ln_phi_acc_d_ln_τ(
    (; τ0, c)::CMP.AccrSB2006{FT},
    τ::FT,
) where {FT <: AbstractFloat} = (c * τ0) / (τ + τ0)


TT.@testset "precipitation_susceptibility_SB2006" begin
    sb2006 = CMP.SB2006(FT)

    q_liq = FT(0.5e-3)
    N_liq = FT(1e8)
    q_rai = FT(1e-5)
    ρ = FT(1)

    τ = FT(1) - q_liq / (q_liq + q_rai)

    aut_rates = CMPS.precipitation_susceptibility_autoconversion(
        sb2006,
        q_liq,
        q_rai,
        ρ,
        N_liq,
    )

    acc_rates = CMPS.precipitation_susceptibility_accretion(
        sb2006,
        q_liq,
        q_rai,
        ρ,
        N_liq,
    )

    TT.@test aut_rates.d_ln_pp_d_ln_N_liq ≈ -2
    TT.@test aut_rates.d_ln_pp_d_ln_q_liq ≈
             4 - (1 - τ) * d_ln_phi_au_d_ln_τ(sb2006.acnv, τ)
    TT.@test aut_rates.d_ln_pp_d_ln_q_rai ≈
             (1 - τ) * d_ln_phi_au_d_ln_τ(sb2006.acnv, τ)

    TT.@test acc_rates.d_ln_pp_d_ln_q_liq ≈
             1 - (1 - τ) * d_ln_phi_acc_d_ln_τ(sb2006.accr, τ)
    TT.@test acc_rates.d_ln_pp_d_ln_q_rai ≈
             1 + (1 - τ) * d_ln_phi_acc_d_ln_τ(sb2006.accr, τ)
end
