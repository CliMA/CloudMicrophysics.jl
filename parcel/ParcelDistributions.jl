import CloudMicrophysics.Parameters as CMP

struct Monodisperse{FT} <: CMP.ParametersType{FT} end

struct MonodisperseMix{FT} <: CMP.ParametersType{FT} end

struct Gamma{FT} <: CMP.ParametersType{FT} end

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

function distribution_moments(::MonodisperseMix, q, N, ρ, ρ_air, q_mode1, N_mode1)
    FT = typeof(q)

    if N == 0 || q == 0
        r = FT(0)
        A = FT(0)
        V = FT(0)
    else
        r_mode1 = N_mode1 == FT(0) ? FT(0) : cbrt(q_mode1 / N_mode1 / FT(4 / 3 * π) / ρ * ρ_air)
        A_mode1 = 4 * FT(π) * r_mode1^2
        V_mode1 = FT(4 / 3 * π) * r_mode1^3
    
        if (N - N_mode1) <= FT(0) || (q - q_mode1) <= FT(0)
            r_mode2 = FT(0)
            A_mode2 = FT(0)
            V_mode2 = FT(0)
        else
            r_mode2 = cbrt((q - q_mode1) / (N - N_mode1) / FT(4 / 3 * π) / ρ * ρ_air)
            A_mode2 = 4 * FT(π) * r_mode2^2
            V_mode2 = FT(4 / 3 * π) * r_mode2^3
        end

        r = (N_mode1 * r_mode1 + (N - N_mode1) * r_mode2) / N
        A = (N_mode1 * A_mode1 + (N - N_mode1) * A_mode2) / N
        V = (N_mode1 * V_mode1 + (N - N_mode1) * V_mode2) / N
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

#! format: on
