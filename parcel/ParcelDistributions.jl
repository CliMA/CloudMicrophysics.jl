import CloudMicrophysics.Parameters as CMP

struct Monodisperse{FT} <: CMP.ParametersType{FT} end

struct Gamma{FT} <: CMP.ParametersType{FT} end

struct Lognormal{FT} <: CMP.ParametersType{FT} end

#! format: off
# Size distributiom moments
function distribution_moments(::Monodisperse, qₗ, Nₗ, ρₗ, ρ_air, qᵢ, Nᵢ, ρᵢ, distr_params)
    FT = typeof(qₗ)
    # liquid droplet
    if Nₗ == FT(0)
        rₗ = FT(0)
        Aₗ = FT(0)
        Vₗ = FT(0)
    else
        rₗ = cbrt(qₗ / Nₗ / FT(4 / 3 * π) / ρₗ * ρ_air)
        Aₗ = 4 * FT(π) * rₗ^2
        Vₗ = FT(4 / 3 * π) * rₗ^3
    end
    # ice crystals
    if Nᵢ == FT(0)
        rᵢ = FT(0)
    else
        rᵢ = cbrt(qᵢ / Nᵢ / FT(4 / 3 * π) / ρᵢ * ρ_air)
    end
    return (; rₗ, Aₗ, Vₗ, rᵢ)
end

function distribution_moments(::Gamma, qₗ, Nₗ, ρₗ, ρ_air, qᵢ, Nᵢ, ρᵢ, distr_params)
    #integral(n(r)) = 1/λ^2
    #integral(r*n(r)) = 2/λ^3
    #<r> = 2/λ
    #integral(r^2*n(r)) = 6/λ^4
    #<r^2> = 6/λ^2
    #integral(r^3*n(r)) = 24/λ^5
    #<r^3> = 24/λ^3
    # liquid droplets
    if Nₗ == FT(0)
        rₗ = FT(0)
        Aₗ = FT(0)
        Vₗ = FT(0)
    else
        λₗ = cbrt(32 * FT(π) * Nₗ / qₗ * ρₗ / ρ_air)
        rₗ = 2 / λₗ
        Aₗ = 4 * FT(π) * 6 / λₗ^2
        Vₗ = 4 / 3 * FT(π) * 24 / λₗ^3
    end
    # ice crystals
    if Nᵢ == FT(0)
        rᵢ = FT(0)
    else
        λᵢ = cbrt(32 * FT(π) * Nᵢ / qᵢ * ρᵢ / ρ_air)
        rᵢ = 2 / λᵢ
    end
    return (; rₗ, Aₗ, Vₗ, rᵢ)
end

function distribution_moments(::Lognormal, qₗ, Nₗ, ρₗ, ρ_air, qᵢ, Nᵢ, ρᵢ, distr_params)
    (; σ_g) = distr_params
    # liquid droplets
    if Nₗ == FT(0)
        rₗ = FT(0)
        Aₗ = FT(0)
        Vₗ = FT(0)
    else
        r_median = cbrt(3 / 4 / FT(π) / Nₗ * qₗ / ρₗ * ρ_air / exp(9 / 2 * (log(σ_g))^2))
        rₗ = r_median * exp(0.5 * (log(σ_g))^2)
        Aₗ = 4 * FT(π) * r_median^2 * exp(2 * (log(σ_g))^2)
        Vₗ = 4 / 3 * FT(π) * r_median^3 * exp(9 / 2 * (log(σ_g))^2)
    end
    # ice crystals; note this is gamma distribution
    if Nᵢ == FT(0)
        rᵢ = FT(0)
    else
        λᵢ = cbrt(32 * FT(π) * Nᵢ / qᵢ * ρᵢ / ρ_air)
        rᵢ = 2 / λᵢ
    end
    return (; rₗ, Aₗ, Vₗ, rᵢ)
end
#! format: on
