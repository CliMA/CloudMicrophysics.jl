"""
    CloudDiagnostics

 - radar reflectivity (1-moment and 2-moment)
 - effective radius  (1-moment, 2-moment and P3)
"""
module CloudDiagnostics

import SpecialFunctions as SF

import ..Parameters as CMP
import ..Microphysics1M as CM1
import ..Microphysics2M as CM2
import ..P3Scheme as P3
import ..DistributionTools as DT
import ..Common as CO
import ..Quadrature as QD
import ..Utilities as UT

"""
    radar_reflectivity_1M(precip, q, ρ)

  - `precip` - struct with 1-moment microphysics rain free parameters
  - `q` - specific content of rain
  - `ρ` - air density

Returns logarithmic radar reflectivity for the 1-moment microphysics
based on the assumed rain particle size distribution.
Normalized by the reflectivity of a 1 millimeter drop in a volume of 1 m³.
The values are clipped at -150 dBZ.
"""
function radar_reflectivity_1M(
    (; pdf, mass)::CMP.Rain,
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
    radar_reflectivity_2M(structs, q_lcl, q_rai, N_lcl, N_rai, ρ_air)

 - `structs` - structs microphysics 2-moment with SB2006 cloud droplets
   and raindrops size distributions parameters
 - `q_lcl` - cloud liquid water specific content
 - `q_rai` - rain water specific content
 - `N_lcl` - cloud droplet number density
 - `N_rai` - rain droplet number density
 - `ρ_air` - air density

Returns logarithmic radar reflectivity for the 2-moment microphysics SB2006
based on the assumed cloud and rain particle size distributions.
Normalized by the reflectivity of a 1 millimeter drop in a volume of 1 m³.
The values are clipped at -150 dBZ.
"""
function radar_reflectivity_2M((; pdf_c, pdf_r)::CMP.SB2006, q_lcl, q_rai, N_lcl, N_rai, ρ_air)
    FT = eltype(q_lcl)
    # free parameters
    (; νc, μc) = pdf_c
    (; νr, μr, ρw) = pdf_r
    C = FT(4 / 3 * π * ρw)
    log_10_Z₀ = -18

    notvalid(B) = iszero(B) || !isfinite(B)  # TODO: Verify that this is the right limit

    # Rain and cloud size distribution parameters
    (; Br) = CM2.pdf_rain_parameters_mass(pdf_r, q_rai, ρ_air, N_rai)
    (; Bc) = CM2.pdf_cloud_parameters_mass(pdf_c, q_lcl, ρ_air, N_lcl)

    # 2nd moment in mass = 6th moment in radius
    n_mass = 2
    Zc = notvalid(Bc) ? FT(0) : DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_lcl, n_mass) / C^n_mass
    Zr = notvalid(Br) ? FT(0) : DT.generalized_gamma_Mⁿ(νr, μr, Br, N_rai, n_mass) / C^n_mass

    return max(FT(-150), 10 * (log10(max(FT(0), Zc + Zr)) - log_10_Z₀))
end

"""
    effective_radius_2M(structs, q_lcl, q_rai, N_lcl, N_rai, ρ_air)

 - `structs` - structs with SB2006 cloud droplets and raindrops
   size distribution parameters
 - `q_lcl` - cloud liquid water specific content
 - `q_rai` - rain water specific content
 - `N_lcl` - cloud droplet number density
 - `N_rai` - rain droplet number density
 - `ρ_air` - air density

Returns effective radius for the 2-moment microphysics scheme.
Computed based on the assumed cloud and rain particle size distributions.
"""
function effective_radius_2M((; pdf_c, pdf_r)::CMP.SB2006, q_lcl, q_rai, N_lcl, N_rai, ρ_air)
    FT = eltype(q_lcl)
    # free parameters
    (; νc, μc) = pdf_c
    (; νr, μr, ρw) = pdf_r
    C = FT(4 / 3 * π * ρw)
    # Rain and cloud size distribution parameters
    (; Br) = CM2.pdf_rain_parameters_mass(pdf_r, q_rai, ρ_air, N_rai)
    (; Bc) = CM2.pdf_cloud_parameters_mass(pdf_c, q_lcl, ρ_air, N_lcl)

    notvalid(B) = iszero(B) || !isfinite(B)
    # 3rd moment in radius = 1st moment in mass
    n_mass = 1
    M3_c = notvalid(Bc) ? FT(0) : DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_lcl, n_mass) / C
    M3_r = notvalid(Br) ? FT(0) : DT.generalized_gamma_Mⁿ(νr, μr, Br, N_rai, n_mass) / C

    # 2nd moment in radius = (2/3)rd moment in mass
    # (use a Float, not a Rational: a rational exponent forces a runtime
    #  Rational{Int64}/gcd construction in GPU kernels)
    n_mass = FT(2) / 3
    M2_c = notvalid(Bc) ? FT(0) : DT.generalized_gamma_Mⁿ(νc, μc, Bc, N_lcl, n_mass) / C^(n_mass)
    M2_r = notvalid(Br) ? FT(0) : DT.generalized_gamma_Mⁿ(νr, μr, Br, N_rai, n_mass) / C^(n_mass)

    return M2_c + M2_r <= UT.ϵ_numerics(FT) ? FT(0) : (M3_c + M3_r) / (M2_c + M2_r)
end

"""
    effective_radius_P3(state, logλ; quad)

Compute the effective radius of the P3 ice size distribution,

```math
r_e = \\frac{3}{4} \\frac{ρq_{ice} / ρ_i}{\\int aᵢ(D) N'(D) dD},
```

the ratio of ice volume to the projected area of the population, where the ice volume is
the ice water content over the solid ice density `ρ_i` and `aᵢ = ice_area(state, D)`.

The numerator is the PROGNOSTIC ice mass rather than its quadrature reconstruction, so
only the area moment is integrated. The two are not interchangeable: the mass is a state
variable the scheme conserves, while its reconstruction carries the quadrature's own
error, and a radius formed as the ratio of two integrals would let that error move the
optical properties of a cell whose mass never changed.

The result has no consumer inside CloudMicrophysics. It exists for a host computing ice
optical properties from the particle size distribution rather than from a prescribed
radius, which is what makes the four-moment P3 ice structure visible to radiation at all.

# Arguments
 - `state`: a [`P3Scheme.P3State`](@ref) object
 - `logλ`: the log of the slope parameter [log(1/m)]

# Keyword Arguments
 - `quad`: quadrature rule (a `Quadrature.QuadratureRule`)

# Returns
 - Effective radius [m], or zero where the area integral vanishes, following the
   convention of [`effective_radius_2M`](@ref).
"""
@inline function effective_radius_P3(state::P3.P3State, logλ; quad)
    FT = eltype(state)
    N′ = DT.size_distribution(state, logλ)
    bnds = P3.integral_bounds(state, logλ; p = FT(1e-6))
    ∫aN = QD.integrate(D -> P3.ice_area(state, D) * N′(D), bnds, quad)
    V = state.ρq_ice / state.params.ρ_i
    return ∫aN > FT(0) ? FT(3) / 4 * V / ∫aN : FT(0)
end

"""
    effective_radius_Liu_Hallet_97(wtr, ρ_air, q_lcl, N_lcl, q_rai, N_rai)
    effective_radius_Liu_Hallet_97(wtr, ρ_air, q_lcl)

 - `wtr` - a struct with water properties (contains the water density `ρw`)
 - `ρ_air` - air density
 - `q_lcl` - cloud water specific content
 - `N_lcl` - cloud droplet number density
 - `q_rai` - rain water specific content
 - `N_rai` - rain droplet number density

Returns effective radius using the "1/3" power law from Liu and Hallett (1997).
If not provided by the user, it is assumed that there is no rain present and that
the cloud droplet number concentration is 1e8 1/m3 (100 1/cm3).
"""
function effective_radius_Liu_Hallet_97(
    (; ρw)::Union{CMP.WaterProperties{FT}, CMP.CloudLiquid{FT}},
    ρ_air::FT,
    q_lcl::FT,
    N_lcl::FT,
    q_rai::FT,
    N_rai::FT,
) where {FT}

    k = FT(0.8)
    r_vol =
        ((N_lcl + N_rai) < UT.ϵ_numerics(FT)) ? FT(0) :
        (
            (FT(3) * (q_lcl + q_rai) * ρ_air) /
            (FT(4) * π * ρw * (N_lcl + N_rai))
        )^FT(1 / 3)

    return r_vol / k^FT(1 / 3)
end
function effective_radius_Liu_Hallet_97(
    wtr::Union{CMP.WaterProperties{FT}, CMP.CloudLiquid{FT}},
    ρ_air::FT,
    q_lcl::FT,
) where {FT}
    return effective_radius_Liu_Hallet_97(
        wtr::Union{CMP.WaterProperties{FT}, CMP.CloudLiquid{FT}},
        ρ_air::FT,
        q_lcl::FT,
        FT(1e8),
        FT(0),
        FT(0),
    )
end

"""
    effective_radius_const(cloud_params)

  - `cloud_params` - a struct with cloud liquid or cloud ice parameters

Returns a constant assumed effective radius for clouds
"""
function effective_radius_const(cloud_params::CMP.CloudLiquid{FT}) where {FT}
    return cloud_params.r_eff
end
function effective_radius_const(cloud_params::CMP.CloudIce{FT}) where {FT}
    return cloud_params.r_eff
end

"""
    rain_intercept_plausibility(range, pdf_r, q_rai, ρ_air, N_rai)

Report whether the rain intercept `N₀` implied by the state falls outside its observational
plausibility range, WITHOUT touching the state or any rate.

[SeifertBeheng2006](@cite) applies this range as a clamp inside the PSD inversion. Under
[`CMP.RainParticlePDF_SB2006_windowed`](@ref) the only bound is on the mean drop mass, and the
intercept range keeps its observational content here instead: a state that leaves the range is
reported and integrated unchanged, rather than being silently rewritten into one that does not
describe it.

An out-of-range intercept is not by itself an error. `N₀ = λ N_r` grows with the drop number at
fixed mean size, so ordinary heavy rain with many drops leaves the upper end and sparse large-drop
populations leave the lower end; what the flag identifies is where the SB2006 cascade WOULD have
intervened, which makes it the natural diagnostic for auditing the difference between the two.

# Arguments
 - `range`: the plausibility range, [`CMP.RainInterceptRange`](@ref)
 - `pdf_r`: rain size distribution parameters, [`CMP.RainParticlePDF_SB2006`](@ref)
 - `q_rai`: rain water specific content [kg/kg]
 - `ρ_air`: air density [kg/m³]
 - `N_rai`: raindrop number density [1/m³]

# Returns
 - `(; N₀r, below, above)`: the implied intercept [1/m⁴] and the two out-of-range flags. Both
   flags are `false` on an empty population, where the inversion returns a zero intercept and
   there is no distribution to call implausible.
"""
function rain_intercept_plausibility(
    (; N0_min, N0_max)::CMP.RainInterceptRange,
    pdf_r::CMP.RainParticlePDF_SB2006, q_rai, ρ_air, N_rai,
)
    (; N₀r) = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ_air, N_rai)
    populated = N₀r > 0
    return (; N₀r, below = populated & (N₀r < N0_min), above = populated & (N₀r > N0_max))
end

end # end module
