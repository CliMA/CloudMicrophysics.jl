# Some wrappers to cast types from SF.gamma
# (which returns Float64 even when the input is Float32)
Γ(a::FT, z::FT) where {FT <: Real} = FT(SF.gamma(a, z))
Γ(a::FT) where {FT <: Real} = FT(SF.gamma(a))

"""
    ∫_Γ(x₀, x_end, c1, c2, c3)

 - x₀ - lower bound
 - x_end - upper bound
 - c1, c2, c3 - respective constants

f(D, c1, c2, c3) = c1 * D ^ (c2) * exp(-c3 * D)

Integrates f(D, c1, c2, c3) dD from x₀ to x_end
"""
function ∫_Γ(x₀::FT, x_end::FT, c1::FT, c2::FT, c3::FT) where {FT}
    if x_end == Inf
        return c1 * c3^(-c2 - 1) * (Γ(1 + c2, x₀ * c3))
    elseif x₀ == 0
        return c1 * c3^(-c2 - 1) * (Γ(1 + c2) - Γ(1 + c2, x_end * c3))
    else
        return c1 * c3^(-c2 - 1) * (Γ(1 + c2, x₀ * c3) - Γ(1 + c2, x_end * c3))
    end
end

"""
    ∫_Γ(x₀, xₘ, x_end, c1, c2, c3, c4, c5, c6)

 - x₀ - lower bound
 - xₘ - switch point
 - x_end - upper bound
 - c1, c2, c3 - respective constants for the first part of the integral
 - c4, c5, c6 - respective constants for the second part of the integral

f(D, c1, c2, c3) = c1 * D ^ (c2) * exp(-c3 * D)

Integrates f(D, c1, c2, c3) dD from x₀ to xₘ and f(D, c4, c5, c6) dD from xₘ to x_end
"""
function ∫_Γ(
    x₀::FT,
    xₘ::FT,
    x_end::FT,
    c1::FT,
    c2::FT,
    c3::FT,
    c4::FT,
    c5::FT,
    c6::FT,
) where {FT}
    return ∫_Γ(x₀, xₘ, c1, c2, c3) + ∫_Γ(xₘ, x_end, c4, c5, c6)
end