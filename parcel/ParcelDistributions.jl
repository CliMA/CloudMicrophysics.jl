import CloudMicrophysics.Parameters as CMP

struct Monodisperse{FT} <: CMP.ParametersType{FT} end

struct Gamma{FT} <: CMP.ParametersType{FT} end

struct Lognormal{FT} <: CMP.ParametersType{FT}
    σ_g::FT
end

#! format: off
# Size distributiom moments
function distribution_moments(::Monodisperse, q, N, ρ, ρ_air)
    FT = typeof(q)
    if N == FT(0) || q == FT(0)
        r = FT(0)
        A = FT(0)
        V = FT(0)
    else
        r = cbrt(q / N / FT(4 / 3 * π) / ρ * ρ_air)
        A = 4 * FT(π) * r^2
        V = FT(4 / 3 * π) * r^3
    end
    return (; r, A, V)
end

function distribution_moments(::Gamma, q, N, ρ, ρ_air)
    FT = typeof(q)
    #integral(n(r)) = 1/λ^2
    #integral(r*n(r)) = 2/λ^3
    #<r> = 2/λ
    #integral(r^2*n(r)) = 6/λ^4
    #<r^2> = 6/λ^2
    #integral(r^3*n(r)) = 24/λ^5
    #<r^3> = 24/λ^3
    if N == FT(0)|| q == FT(0)
        r = FT(0)
        A = FT(0)
        V = FT(0)
    else
        λ = cbrt(32 * FT(π) * N / q * ρ / ρ_air)
        r = 2 / λ
        A = 4 * FT(π) * 6 / λ^2
        V = 4 / 3 * FT(π) * 24 / λ^3
    end
    return (; r, A, V)
end

function distribution_moments(params::Lognormal, q, N, ρ, ρ_air)
    FT = typeof(q)
    # M_0 = N_0
    # M_1 = N_0 * r_m * exp(1/2 * ln^2(σ_g))
    # M_2 = N_0 * r_m^2 * exp(2 * ln^2(σ_g))
    # M_3 = N_0 * r_m^3 * exp(9/2 * ln^2(σ_g))
    (; σ_g) = params
    if N == FT(0) || q == FT(0)
        r = FT(0)
        A = FT(0)
        V = FT(0)
    else
        r_median = cbrt(3 / 4 / FT(π) / N * q / ρ * ρ_air / exp(9 / 2 * (log(σ_g))^2))
        r = r_median * exp(0.5 * (log(σ_g))^2)
        A = 4 * FT(π) * r_median^2 * exp(2 * (log(σ_g))^2)
        V = 4 / 3 * FT(π) * r_median^3 * exp(9 / 2 * (log(σ_g))^2)
    end
    return (; r, A, V)
end
#! format: on
