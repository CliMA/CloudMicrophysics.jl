# helper functions for Γ_approx
function c₁(a::FT) where {FT <: Real}
    p1 = FT(9.4368392235e-03)
    p2 = -FT(1.0782666481e-04)
    p3 = -FT(5.8969657295e-06)
    p4 = FT(2.8939523781e-07)
    p5 = FT(1.0043326298e-01)
    p6 = FT(5.5637848465e-01)
    p = [p1, p2, p3, p4, p5, p6]

    ret = 1 + p[5] * (exp(-p[6] * a) - 1)
    for it in range(1, 4, step = 1)
        ret += p[it] * a^it
    end
    return ret
end
function c₂(a::FT) where {FT <: Real}
    q1 = FT(1.1464706419e-01)
    q2 = FT(2.6963429121)
    q3 = -FT(2.9647038257)
    q4 = FT(2.1080724954)
    q = [q1, q2, q3, q4]

    ret = q[1]
    for it in range(2, 4, step = 1)
        ret += q[it] / a^(it - 1)
    end
    return ret
end
function c₃(a::FT) where {FT <: Real}
    r1 = FT(0)
    r2 = FT(1.1428716184)
    r3 = -FT(6.6981186438e-03)
    r4 = FT(1.0480765092e-04)
    r = [r1, r2, r3, r4]

    ret = 0
    for it in range(1, 4, step = 1)
        ret += r[it] * a^(it - 1)
    end
    return ret
end
function c₄(a::FT) where {FT <: Real}
    s1 = FT(1.0356711153)
    s2 = FT(2.3423452308)
    s3 = -FT(3.6174503174e-01)
    s4 = -FT(3.1376557650)
    s5 = FT(2.9092306039)
    s = [s1, s2, s3, s4, s5]

    ret = 0
    for it in range(1, 5, step = 1)
        ret += s[it] / a^(it - 1)
    end
    return ret
end

"""
    Γ_lower(a, z)

An approximated lower incomplete Gamma function based on Blahak 2010:
doi:10.5194/gmd-3-329-2010
https://gmd.copernicus.org/articles/3/329/2010/gmd-3-329-2010.pdf
"""
function Γ_lower(a::FT, z::FT) where {FT <: Real}

    if isnan(z) || z <= FT(0)
        return FT(0)
    else
        W(a, z) = FT(0.5) + FT(0.5) * tanh(c₂(a) * (z - c₃(a)))

        return exp(-z) *
               z^a *
               (
                   1 / a +
                   c₁(a) * z / a / (a + 1) +
                   (c₁(a) * z)^2 / a / (a + 1) / (a + 2)
               ) *
               (1 - W(a, z)) + CO.Γ(a) * W(a, z) * (1 - c₄(a)^(-z))
    end
end

"""
    Γ_upper(a, z)

An approximated upper incomplete Gamma function computed as Γ(a) - Γ_lower(a, z)
"""
function Γ_upper(a::FT, z::FT) where {FT <: Real}
    return CO.Γ(a) - Γ_lower(a, z)
end

"""
    ∫_Γ(x₀, x_end, c1, c2, c3)

 - x₀ - lower bound
 - x_end - upper bound
 - c1, c2, c3 - respective constants

f(D, c1, c2, c3) = c1 * D ^ (c2) * exp(-c3 * D)

Integrates f(D, c1, c2, c3) dD from x₀ to x_end
"""
function ∫_Γ(x₀::FT, x_end::FT, c1::FT, c2::FT, c3::FT) where {FT}
    if x_end == FT(Inf)
        return c1 * c3^(-c2 - 1) * (Γ_upper(1 + c2, x₀ * c3))
    elseif x₀ == 0
        return c1 * c3^(-c2 - 1) * (Γ_lower(1 + c2, x_end * c3))
    else
        return c1 *
               c3^(-c2 - 1) *
               (Γ_upper(1 + c2, x₀ * c3) - Γ_upper(1 + c2, x_end * c3))
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
