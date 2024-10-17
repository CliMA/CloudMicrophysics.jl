"""
    α_va_si(p3)

 - p3 - a struct with P3 scheme parameters

Returns `α_va` coefficient for the assumed particle mass(size) relation for
large unrimed ice and dense nonspherical ice, in base SI units: kg m^(-β_va).
`β_va` is another coefficient of the mass(size) relation.
From measurements of mass grown by vapor diffusion and aggregation
in midlatitude cirrus by Brown and Francis (1995)
doi: 10.1175/1520-0426(1995)012<0410:IMOTIW>2.0.CO;2
"""
α_va_si(p3::PSP3{FT}) where {FT} = p3.α_va * 10^(6 * p3.β_va - 3)

"""
    D_th_helper(p3)

 - p3 - a struct with P3 scheme parameters

Returns the critical size separating spherical and nonspherical ice, in meters.
Eq. 8 in Morrison and Milbrandt (2015).
"""
D_th_helper(p3::PSP3{FT}) where {FT} =
    (FT(π) * p3.ρ_i / 6 / α_va_si(p3))^(1 / (p3.β_va - 3))

"""
    D_cr_helper(p3, F_rim, ρ_g)

 - p3 - a struct with P3 scheme parameters
 - F_rim - rime mass fraction L_rim / L_ice) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the size of equal mass for graupel and partially rimed ice, in meters.
Eq. 14 in Morrison and Milbrandt (2015).
"""
function D_cr_helper(p3::PSP3{FT}, F_rim::FT, ρ_g::FT) where {FT}
    α_va = α_va_si(p3)
    return (1 / (1 - F_rim) * 6 * α_va / FT(π) / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    D_gr_helper(p3, ρ_g)

 - p3 - a struct with P3 scheme parameters
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the size of equal mass for graupel and unrimed ice, in meters.
Eq. 15 in Morrison and Milbrandt (2015).
"""
function D_gr_helper(p3::PSP3{FT}, ρ_g::FT) where {FT}
    α_va = α_va_si(p3)
    return (6 * α_va / FT(π) / ρ_g)^(1 / (3 - p3.β_va))
end

"""
    ρ_g_helper(ρ_r, F_rim, ρ_d)

 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim / L_ice) [-]
 - ρ_g - is the effective density of a spherical graupel particle [kg/m^3]

Returns the density of total (deposition + rime) ice mass for graupel, in kg/m3
Eq. 16 in Morrison and Milbrandt (2015).
"""
ρ_g_helper(ρ_r::FT, F_rim::FT, ρ_d::FT) where {FT} =
    F_rim * ρ_r + (1 - F_rim) * ρ_d

"""
    ρ_d_helper(p3, D_cr, D_gr)

 - p3 - a struct with P3 scheme parameters
 - D_cr - is the size of equal mass for graupel and partially rimed ice, in meters
 - D_gr - the size of equal mass for graupel and unrimed ice, in meters

Returns the density of unrimed ice mass, in kg/m3
Eq. 17 in Morrison and Milbrandt (2015).
"""
function ρ_d_helper(p3::PSP3{FT}, D_cr::FT, D_gr::FT) where {FT}
    α_va = α_va_si(p3)
    β_m2 = p3.β_va - 2
    return 6 * α_va * (D_cr^β_m2 - D_gr^β_m2) / FT(π) / β_m2 /
           max(D_cr - D_gr, eps(FT))
end

"""
    thresholds(p3, ρ_r, F_rim)

 - p3 - a struct with P3 scheme parameters
 - ρ_r - rime density (L_rim/B_rim) [kg/m^3]
 - F_rim - rime mass fraction (L_rim / L_ice) [-]

Solves the nonlinear system consisting of D_cr, D_gr, ρ_g, ρ_d
for a given rime density and rime mass fraction.
Returns a named tuple containing:
 - D_cr - is the threshold size separating partially rimed ice and graupel [m],
 - D_gr - is the threshold size separating graupel and dense nonspherical ice [m],
 - ρ_g - is the effective density of a spherical graupel particle [kg/m3],
 - ρ_d - is the density of the unrimed portion of the particle [kg/m3],
"""
function thresholds(p3::PSP3{FT}, ρ_r::FT, F_rim::FT) where {FT}

    @assert F_rim >= FT(0)   # rime mass fraction must be positive ...
    @assert F_rim < FT(1)    # ... and there must always be some unrimed part

    if F_rim == FT(0)
        return (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0))
    else
        @assert ρ_r > FT(0)   # rime density must be positive ...
        @assert ρ_r <= p3.ρ_l # ... and as a bulk ice density can't exceed the density of water

        P3_problem(ρ_d) =
            ρ_d - ρ_d_helper(
                p3,
                D_cr_helper(p3, F_rim, ρ_g_helper(ρ_r, F_rim, ρ_d)),
                D_gr_helper(p3, ρ_g_helper(ρ_r, F_rim, ρ_d)),
            )

        ρ_d =
            RS.find_zero(
                P3_problem,
                RS.SecantMethod(FT(0), FT(1000)),
                RS.CompactSolution(),
            ).root
        ρ_g = ρ_g_helper(ρ_r, F_rim, ρ_d)

        return (;
            D_cr = D_cr_helper(p3, F_rim, ρ_g),
            D_gr = D_gr_helper(p3, ρ_g),
            ρ_g,
            ρ_d,
        )
    end
end

"""
    p3_F_liq_average(F_liq, X_ice, X_liq)

 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - X_ice - ice core parameterization (i.e. mass, etc)
 - X_liq - liquid part parameterization

Returns the liquid fraction weighted average of X_ice and X_liq.
"""
function p3_F_liq_average(F_liq::FT, X_ice::FT, X_liq::FT) where {FT}
    return (1 - F_liq) * X_ice + F_liq * X_liq
end

"""
    p3_density(p3, D, F_rim, th)

- p3 - a struct with P3 parameters
- D - maximum particle dimension [m]
- F_rim - rime mass fraction (L_rim / L_ice) [-]
- th - P3 scheme thresholds() output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns the density of a particle based on where it falls in the particle-size-based properties
regime. Following Morrison and Milbrandt (2015), the density of nonspherical particles is assumed to
be the particle mass divided by the volume of a sphere with the same D.
Needed for aspect ratio calculation, so we assume zero liquid fraction.
"""
function p3_density(p3::PSP3, D::FT, F_rim::FT, th) where {FT}
    if D_th_helper(p3) > D
        # small spherical ice
        return p3.ρ_i
    elseif F_rim == 0
        # large nonspherical unrimed ice
        return (6 * α_va_si(p3)) / FT(π) * D^(p3.β_va - 3)
    elseif th.D_gr > D >= D_th_helper(p3)
        # dense nonspherical ice
        return (6 * α_va_si(p3)) / FT(π) * D^(p3.β_va - 3)
    elseif th.D_cr > D >= th.D_gr
        # graupel
        return th.ρ_g
    else #elseif D >= th.D_cr
        # partially rimed ice
        return (6 * α_va_si(p3)) / (FT(π) * (1 - F_rim)) * D^(p3.β_va - 3)
    end
end

"""
    mass_(p3, D, ρ, F_rim)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension [m]
 - ρ - bulk ice density (ρ_i for small ice, ρ_g for graupel) [kg/m3]
 - F_rim - rime mass fraction [L_rim/L_ice]

Returns mass as a function of size for differen particle regimes [kg]
"""
# for spherical particles (small ice, completely rimed ice, or liquid on ice)
mass_s(D::FT, ρ::FT) where {FT <: Real} = FT(π) / 6 * ρ * D^3
# for large nonspherical ice (used for unrimed and dense types)
mass_nl(p3::PSP3, D::FT) where {FT <: Real} = α_va_si(p3) * D^p3.β_va
# for partially rimed ice
mass_r(p3::PSP3, D::FT, F_rim::FT) where {FT <: Real} =
    α_va_si(p3) / (1 - F_rim) * D^p3.β_va

"""
    p3_mass(p3, D, F_rim, F_liq, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_rim - rime mass fraction (L_rim / L_ice)
 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - th - P3 scheme thresholds() output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns mass(D) regime, used to create figures for the docs page.
"""
function p3_mass(
    p3::PSP3,
    D::FT,
    F_rim::FT,
    F_liq::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    if D_th_helper(p3) > D
        return p3_F_liq_average(F_liq, mass_s(D, p3.ρ_i), mass_s(D, p3.ρ_l))         # small spherical ice
    elseif F_rim == 0
        return p3_F_liq_average(F_liq, mass_nl(p3, D), mass_s(D, p3.ρ_l))            # large nonspherical unrimed ice
    elseif th.D_gr > D >= D_th_helper(p3)
        return p3_F_liq_average(F_liq, mass_nl(p3, D), mass_s(D, p3.ρ_l))            # dense nonspherical ice
    elseif th.D_cr > D >= th.D_gr
        return p3_F_liq_average(F_liq, mass_s(D, th.ρ_g), mass_s(D, p3.ρ_l))         # graupel
    else #elseif D >= th.D_cr
        return p3_F_liq_average(F_liq, mass_r(p3, D, F_rim), mass_s(D, p3.ρ_l))      # partially rimed ice
    end
end

"""
    p3_dmdD(p3, D, F_r, th)

 - p3 - a struct containing p3 parameters
 - D - maximum dimension of the particle
 - F_r - rime mass fraction (q_rim/ q_i)
 - th - thresholds as calculated by thresholds()

Returns dm(D)/dD for each particle regime, used in snow melting
"""
function p3_dmdD(
    p3::PSP3,
    D::FT,
    F_r::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT <: Real}
    D_th = D_th_helper(p3)
    if D_th > D
        return FT(π) / 2 * p3.ρ_i * D^2
    elseif F_r == 0
        return α_va_si(p3) * p3.β_va * D^(p3.β_va - 1)
    elseif th.D_gr > D >= D_th
        return α_va_si(p3) * p3.β_va * D^(p3.β_va - 1)
    elseif th.D_cr > D >= th.D_gr
        return FT(π) / 2 * th.ρ_g * D^2
    else # D >= th.D_cr
        return α_va_si(p3) / (1 - F_r) * p3.β_va * D^(p3.β_va - 1)
    end
end

"""
    A_(p3, D)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension

Returns particle projected area as a function of size for different particle regimes
"""
# for spherical particles
A_s(D::FT) where {FT} = FT(π) / 4 * D^2
# for nonspherical particles
A_ns(p3::PSP3, D::FT) where {FT} = p3.γ * D^p3.σ
# partially rimed ice
A_r(p3::PSP3, F_rim::FT, D::FT) where {FT} =
    F_rim * A_s(D) + (1 - F_rim) * A_ns(p3, D)

"""
    p3_area(p3, D, F_rim, F_liq, th)

 - p3 - a struct with P3 scheme parameters
 - D - maximum particle dimension
 - F_rim - rime mass fraction (L_rim / L_ice)
 - F_liq - liquid fraction (L_liq / L_p3_tot)
 - th - P3 scheme thresholds() output tuple (D_cr, D_gr, ρ_g, ρ_d)

Returns area(D), used to create figures for the documentation.
"""
function p3_area(
    p3::PSP3,
    D::FT,
    F_rim::FT,
    F_liq::FT,
    th = (; D_cr = FT(0), D_gr = FT(0), ρ_g = FT(0), ρ_d = FT(0)),
) where {FT}
    # Area regime:
    if D_th_helper(p3) > D
        return A_s(D)                                              # small spherical ice
    elseif F_rim == 0
        return p3_F_liq_average(F_liq, A_ns(p3, D), A_s(D))        # large nonspherical unrimed ice
    elseif th.D_gr > D >= D_th_helper(p3)
        return p3_F_liq_average(F_liq, A_ns(p3, D), A_s(D))        # dense nonspherical ice
    elseif th.D_cr > D >= th.D_gr
        return A_s(D)                                              # graupel
    else # D >= th.D_cr
        return p3_F_liq_average(F_liq, A_r(p3, F_rim, D), A_s(D))  # partially rimed ice
    end
end

"""
    ϕᵢ(p3, D, F_rim, th)

 - p3 - a struct containing P3 parameters
 - D - maximum dimension of ice particle [m]
 - F_rim - rime mass fraction (L_rim/ L_ice) [-]
 - th - P3 particle properties thresholds

Returns the aspect ratio (ϕ) for an ice particle with mass, cross-sectional area,
and ice density determined using the size-dependent particle property regimes
following Morrison and Milbrandt (2015). The density of nonspherical
particles is assumed to be equal to the particle mass divided by the volume of a
spherical particle with the same D_max.
Assuming zero liquid fraction and oblate shape.
"""
function ϕᵢ(p3::PSP3, D::FT, F_rim::FT, th) where {FT}
    F_liq = FT(0)
    mᵢ = p3_mass(p3, D, F_rim, F_liq, th)
    aᵢ = p3_area(p3, D, F_rim, F_liq, th)
    ρᵢ = p3_density(p3, D, F_rim, th)

    # TODO - prolate or oblate?
    ϕ_ob = min(1, 3 * sqrt(FT(π)) * mᵢ / (4 * ρᵢ * aᵢ^FT(1.5))) # κ =  1/3
    #ϕ_pr = max(1, 16 * ρᵢ^2 * aᵢ^3 / (9 * FT(π) * mᵢ^2))       # κ = -1/6

    return ifelse(D == 0, FT(0), ϕ_ob)
end
