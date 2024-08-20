export Γ_approx

# Some wrappers to cast types from SF.gamma
# (which returns Float64 even when the input is Float32)
Γ(a::FT, z::FT) where {FT <: Real} = FT(SF.gamma(a, z))
Γ(a::FT) where {FT <: Real} = FT(SF.gamma(a))

function c₁(a::FT) where {FT <: Real}
    # TODO - move to parameters
    p = [9.4368392235e-03, −1.0782666481e-04, −5.8969657295e-06, 2.8939523781e-07, 1.0043326298e-01, 5.5637848465e-01]

    ret = 1 + p[5] * (exp(-p[6] * a) - 1)
    for it in range(1, 4, step=1)
        ret += p[it] * a^it
    end
    return ret
end
function c₂(a::FT) where {FT <: Real}
    q = [1.1464706419e-01, 2.6963429121, −2.9647038257, 2.1080724954]
    #ret = 0
    #for it in range(1,4,step=1)
    #    ret += q[it] / a^(it-1)
    #end
    ret = q[1] + (q[2]/a) + (q[3]/a^2) + (q[4]/a^3)
    @info(" ", a, q[1] + q[2]/a + q[3]/a^2 + q[4]/a^3, ret)
    return ret
end
function c₃(a::FT) where {FT <: Real}
    r = [0.0, 1.1428716184, −6.6981186438e-03, 1.0480765092e-04]
    ret = 0
    for it in range(1,4,step=1)
       ret +=  r[it] * a^(it -1)
    end
    return ret
end
function c₄(a::FT) where {FT <: Real}
    s = [1.0356711153, 2.3423452308, −3.6174503174e-01, −3.1376557650, 2.9092306039]
    ret = 0
    for it in range(1,5,step=1)
        ret +=  s[it] / a^(it-1)
    end
    return ret
end

"""
    Γ_approx(a, z)

An approximated lower incomplete Gamma function based on Blahak 2010:
doi:10.5194/gmd-3-329-2010
https://gmd.copernicus.org/articles/3/329/2010/gmd-3-329-2010.pdf
"""
function Γ_approx(a::FT, z::FT) where {FT <: Real}

    W(a, z) = 0.5 + 0.5 * tanh(c₂(a) * (z - c₃(a)))

    @info(" ", c₁(a), c₂(a), c₃(a), c₄(a), W(a, z))

    tmp = exp(-z) * z^a
    tmp2 = (1/a + c₁(a) * z / a / (a+1) + (c₁(a) * z)^2 / a / (a+1) / (a+2))
    tmp3 = (1 - W(a, z)) + Γ(a) * W(a, z) * (1 - c₄(a)^(-z))

    @info(" ", tmp, tmp2, tmp3)

    return exp(-z) * z^a *
        (1/a + c₁(a) * z / a / (a+1) + (c₁(a) * z)^2 / a / (a+1) / (a+2)) *
        (1 - W(a, z)) + Γ(a) * W(a, z) * (1 - c₄(a)^(-z))
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
