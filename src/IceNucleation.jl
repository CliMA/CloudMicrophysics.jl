# This file contains the `HetIceNucleation` and `HomIceNucleation` modules.

"""
Parameterization for homogeneous cloud ice nucleation
"""
module HomIceNucleation

import ..Parameters as CMP
import CloudMicrophysics.Common as CO

export homogeneous_J_cubic
export homogeneous_J_linear

"""
    homogeneous_J_cubic(ip, Δa_w)

Calculate the homogeneous freezing nucleation rate coefficient, `J` [m⁻³ s⁻¹],
for sulphuric acid solutions.

# Arguments
  - `ip`: The [`CMP.Koop2000`](@ref) struct with ice nucleation parameters,
    + `c₁`, `c₂`, `c₃`, `c₄`: cubic fit coefficients [-]
    + `Δa_w_min`: minimum change in water activity [-]
    + `Δa_w_max`: maximum change in water activity [-]
  - `Δa_w`: change in water activity [-].

Returns the homogeneous freezing nucleation rate coefficient,
`J`, in m⁻³ s⁻¹ for sulphuric acid solutions.
Parameterization based on [Koop2000](@cite), see doi.org/10.1038/35020537.
"""
function homogeneous_J_cubic(ip::CMP.Koop2000, Δa_w::FT) where {FT}
    (; Δa_w_min, Δa_w_max) = ip
    Δa_w_min ≤ Δa_w ≤ Δa_w_max || throw(
        DomainError(Δa_w,
            lazy"Change in water activity must be within the valid range: Δa_w ∈ [$Δa_w_min, $Δa_w_max], but Δa_w = $Δa_w",
        ),
    )
    return _homogeneous_J_cubic_fit(ip, Δa_w)
end

# The cubic fit of `homogeneous_J_cubic` without the range check
@inline function _homogeneous_J_cubic_fit((; c₁, c₂, c₃, c₄)::CMP.Koop2000, Δa_w)
    logJ = c₁ + c₂ * Δa_w - c₃ * Δa_w^2 + c₄ * Δa_w^3
    return 10^(logJ + 6)
end

"""
    homogeneous_J_linear(ip, Δa_w)

Calculate the homogeneous freezing nucleation rate coefficient, `J` [m⁻³ s⁻¹],
for sulphuric acid solutions.

# Arguments
  - `ip`: The [`CMP.Koop2000`](@ref) struct with ice nucleation parameters,
    + `linear_c₁`, `linear_c₂`: linear fit coefficients [-]
  - `Δa_w`: change in water activity [-].

Model is a linear fit of the [Koop2000](@cite) parameterization.
See: doi.org/10.1038/35020537
"""
function homogeneous_J_linear((; linear_c₁, linear_c₂)::CMP.Koop2000, Δa_w)
    logJ = linear_c₂ * Δa_w + linear_c₁
    return 10^(logJ + 6)
end

"""
    homogeneous_freezing_rate_coefficient(hom, tps, T)

Return the Koop (2000) homogeneous freezing rate coefficient of pure water [m⁻³(water) s⁻¹].

The water activity difference is `Δa_w = 1 - e_sat,ice(T) / e_sat,liq(T)`.
The coefficient is zero warmer than the fitted range of `Δa_w`, and it is the value at the cold edge colder than it.

# Arguments
 - `hom`: the [`CMP.Koop2000`](@ref) parameterization.
 - `tps`: thermodynamics parameters.
 - `T`: air temperature [K].
"""
@inline function homogeneous_freezing_rate_coefficient(hom::CMP.Koop2000, tps, T)
    (; Δa_w_min, Δa_w_max) = hom
    a_w = CO.a_w_ice(tps, T)
    # Where the saturation vapor pressure over liquid underflows, `a_w` is NaN, and `Δa_w = 1` is at the cold edge
    Δa_w = 1 - ifelse(isnan(a_w), zero(a_w), a_w)
    J = _homogeneous_J_cubic_fit(hom, clamp(Δa_w, Δa_w_min, Δa_w_max))
    return ifelse(Δa_w < Δa_w_min, zero(J), J)
end

end # end module

"""
Parameterization for heterogenous cloud ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP
import ..HomIceNucleation as CM_HomIce
import CloudMicrophysics.Common as CO
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.Utilities as UT
import CloudMicrophysics.Quadrature as QUAD
import ForwardDiff as FD

export dust_activated_number_fraction
export MohlerDepositionRate
export deposition_J
export ABIFM_J
export P3_deposition_N_i
export P3_het_N_i
export INP_concentration_frequency
export INP_concentration_mean
export rain_freezing_rate
export cloud_freezing_rate
export immersion_limit_rate
export deposition_rate
export n_active

"""
    dust_activated_number_fraction(dust, ip, Si, T)

Calculate the number fraction of mineral dust particles acting as deposition nuclei,

```n ice nuclei / n dust particles```

# Arguments
  - `dust`: a struct with dust parameters
  - `ip`: a struct with cloud ice nucleation parameters
  - `Si`: ice saturation ratio
  - `T`: air temperature [K]


From [Mohler2006](@cite) Table 2 (averages from different measurements
excluding those where a was not measured), see doi.org/10.5194/acp-6-3007-2006
"""
function dust_activated_number_fraction(
    dust::Union{CMP.DesertDust, CMP.ArizonaTestDust}, ip::CMP.Mohler2006, Si, T,
)
    @assert Si < ip.Sᵢ_max

    S₀ = ifelse(T > ip.T_thr, dust.S₀_warm, dust.S₀_cold)
    a = ifelse(T > ip.T_thr, dust.a_warm, dust.a_cold)
    return max(0, exp(a * (Si - S₀)) - 1)
end

"""
    MohlerDepositionRate(dust, ip, Si, T, dSi_dt, N_aer)

Calculate the cloud ice nucleation rate from deposition.

# Arguments
  - `dust`: a struct with dust parameters
  - `ip`: a struct with cloud ice nucleation parameters
  - `Si`: ice saturation
  - `T`: ambient temperature
  - `dSi_dt`: change in ice saturation over time
  - `N_aer`: number of unactivated aerosols

See [Mohler2006](@cite) Equation 5; doi.org/10.5194/acp-6-3007-2006
"""
function MohlerDepositionRate(
    dust::Union{CMP.DesertDust, CMP.ArizonaTestDust}, ip::CMP.Mohler2006,
    Si, T, dSi_dt, N_aer,
)
    @assert Si < ip.Sᵢ_max

    a = ifelse(T > ip.T_thr, dust.a_warm, dust.a_cold)
    return max(0, N_aer * a * dSi_dt)
end

"""
    deposition_J(dust, Δa_w)

Calculate the deposition nucleation rate coefficient, `J` [m⁻² s⁻¹],
for water vapor deposition onto different dust and mineral aerosol types.

# Arguments
  - `dust`: a struct with dust parameters (supported types: feldspar,
    ferrihydrite, kaolinite, illite, Arizona Test Dust, Saharan dust,
    Asian dust, and generic dust)
  - `Δa_w`: change in water activity [unitless].

# Returns
 - `J` [m⁻² s⁻¹]; zero for unsupported aerosol types.

See [China2017](@cite) for details on the parameterization.
"""
function deposition_J(
    dust::Union{
        CMP.Ferrihydrite, CMP.Feldspar, CMP.Kaolinite, CMP.Illite, CMP.ArizonaTestDust,
        CMP.SaharanDust, CMP.AsianDust, CMP.Dust,
    },
    Δa_w,
)
    logJ = dust.deposition_m * Δa_w + dust.deposition_c
    return 10^(logJ + 4) # converts cm⁻² s⁻¹ to m⁻² s⁻¹
end
deposition_J(::CMP.AerosolType, Δa_w) = zero(eltype(Δa_w))

"""
    ABIFM_J(dust, Δa_w)

Compute the heterogeneous ice nucleation rate coefficient, `J` [m⁻² s⁻¹]
    for the given `dust` type and solution water activity, `Δa_w`,
    using the "a_w based immersion freezing model" (ABIFM)

# Arguments
 - `dust`: The given mineral in liquid solution; currently supports:
    + `DesertDust`, `Illite`, `Kaolinite`, `Dust`, `ArizonaTestDust`,
      `MiddleEasternDust`, `AsianDust`
    + all other `AerosolType`s are not supported and will return zero
 - `Δa_w`: change in water activity [unitless].

# Returns
 - `J`: heterogeneous ice nucleation rate coefficient [m⁻² s⁻¹]

The free parameters `m` and `c` are taken from Knopf & Alpert 2013
see: doi.org/10.1039/C3FD00035D
"""
function ABIFM_J(
    dust::Union{
        CMP.DesertDust, CMP.Illite, CMP.Kaolinite, CMP.Dust,
        CMP.ArizonaTestDust, CMP.MiddleEasternDust, CMP.AsianDust,
    },
    Δa_w,
)
    logJ = dust.ABIFM_m * Δa_w + dust.ABIFM_c
    return 10^(logJ + 4) # `+4` converts cm⁻² s⁻¹ to m⁻² s⁻¹
end
ABIFM_J(::CMP.AerosolType, Δa_w) = zero(eltype(Δa_w))

"""
    P3_deposition_N_i(ip, T)

Calculate the number of ice crystals nucleated via deposition nucleation with units of m⁻³.

# Arguments
 - `ip`: a struct with ice nucleation parameters:
    + `c₁`: constant [L⁻¹]
    + `c₂`: constant [K⁻¹]
    + `T₀`: freezing temperature [K]
    + `T_dep_thres`: lower cutoff temperature [K]
 - `T`: air temperature [K].

# Returns
 - `Nᵢ`: number of ice crystals nucleated via deposition nucleation with units of m⁻³.

From Thompson et al 2004 eqn 2 as used in Morrison & Milbrandt 2015,

```
Nᵢ = c₁ exp(c₂ (T₀ - T))
```

where, in Thompson et al 2004, `c₁ = 0.005`, `c₂ = 0.304`, `T₀ = 273.15 K`,
and `T` is the air temperature [K].
The nucleation number is at most the value at `T = T_dep_thres`, and is zero above `T₀ = 0°C`.
"""
function P3_deposition_N_i((; c₁, c₂, T₀, T_dep_thres)::CMP.MorrisonMilbrandt2014, T)
    T′ = max(T_dep_thres, T)  # clamp T to T_thres ≤ T
    Nᵢ = 1000 * c₁ * exp(c₂ * (T₀ - T′))  # 1000 converts L⁻¹ to m⁻³
    return ifelse(T < T₀, Nᵢ, zero(Nᵢ))  # only allow deposition nucleation below T₀ (0°C)
end

"""
    P3_het_N_i(ip, T, Nₗ, Vₗ, Δt)

Compute number of ice crystals formed from heterogeneous condensation freezing

# Arguments
 - `ip`: The [`CMP.MorrisonMilbrandt2014`](@ref) paramterization, where:
    + `het_a`: empirical parameter [K⁻¹]
    + `het_B`: water-type dependent parameter [m⁻³ s⁻¹]
    + `T₀`: freezing temperature [K]
 - `T`: air temperature [K],
 - `Nₗ`: number of droplets [m⁻³],
 - `Vₗ`: volume of droplets to be frozen [m³],
 - `Δt`: timestep [s].

# Returns
 - `Nᵢ`: number of ice crystals [m⁻³] heterogeneously nucleated
    from cloud droplets in `Δt` seconds.

From Pruppacher & Klett 1997 eqn (9-51) as used in [MorrisonMilbrandt2015](@cite):

```
ln N₀ / Nᵤ(t) = B Vₗ [exp(aTₛ)] t
```

where `N₀=Nᵤ(t=0)` is the initial number of cloud droplets, `a` and `B` are
empirical parameters, `Tₛ` is the temperature difference between the freezing
point and the air temperature, and `Vₗ` is the volume of cloud droplets to be
frozen. Rearranged in terms of `Nᵤ(t)`:

```
Nᵤ(t) = N₀ exp(-B Vₗ [exp(aTₛ)] t)
```
"""
function P3_het_N_i((; het_a, het_B, T₀)::CMP.MorrisonMilbrandt2014, T, Nₗ, Vₗ, Δt)
    Tₛ = T₀ - T
    return Nₗ * (1 - exp(-het_B * Vₗ * Δt * exp(het_a * Tₛ)))
end

"""
    INP_concentration_frequency(params, INPC, T)

Calculate the relative frequency of a given INP concentration as a function of temperature.

# Arguments
 - `params`: a struct with INPC(T) distribution parameters
 - `INPC`: concentration of ice nucleating particles [m^-3]
 - `T`: air temperature [K]

For details see: [Frostenberg2023](@cite), doi.org/10.5194/acp-23-10883-2023
"""
function INP_concentration_frequency(params::CMP.Frostenberg2023, INPC, T)
    (; T_freeze, σ) = params
    T ≥ T_freeze && return zero(INPC)
    μ = INP_concentration_mean(params, T)
    return exp(-(log(INPC) - μ)^2 / 2σ^2) / √(π * 2σ^2)
end

"""
    INP_concentration_mean(params, T)

Calculate the mean log(INPC) as a function of temperature.

# Arguments
  - `params`: The [`CMP.Frostenberg2023`](@ref) INPC(T) distribution parameters, including
    + `T_freeze`: freezing temperature [K]
    + `b`: temperature normalization coefficient [°C⁻¹]
  - `T`: air temperature [K]

Following Eq. (1) of [Frostenberg2023](@cite), the mean log(INPC) is
```
μ(T) = 9 log(-b · T_celsius / 10)
```
with the corresponding INPC obtained by exponentiating. At the default `b = 1`
this is the marine-dataset curve `log((-T_celsius / 10)^9)`.

For details see: [Frostenberg2023](@cite), doi.org/10.5194/acp-23-10883-2023
"""
function INP_concentration_mean((; T_freeze, b)::CMP.Frostenberg2023, T)
    T_celsius = min(T - T_freeze, 0)
    return 9log(-b * T_celsius / 10)  # = log((-b * T_celsius / 10)^9)
end

# Whether a liquid population can freeze: both moments are positive and it is supercooled
@inline function _liquid_freezing_is_active(q, n, T, T_freeze)
    # Compare values, since a ForwardDiff dual with zero value and nonzero partials exceeds zero
    return (FD.value(n) > zero(FD.value(n))) & (FD.value(q) > zero(FD.value(q))) & (T < T_freeze)
end

"""
    drop_freezing_heat_timescale(vent, aps, tps, Tₐ, ρₐ, qᵥ, ρw, v_term)

Return a function `τ_heat(D)` [s], the time a nucleated drop of diameter `D` [m] needs to shed its latent heat.

The drop sheds the rest of its latent heat of fusion by conduction and evaporation through its ventilated surface.
It is zero at and above freezing. See the Ice Nucleation docs page.

# Arguments
 - `vent`: ventilation coefficients, [`CMP.VentilationFactor`](@ref).
 - `aps`: air properties, [`CMP.AirProperties`](@ref).
 - `tps`: thermodynamics parameters.
 - `Tₐ`: air temperature [K].
 - `ρₐ`: air density [kg m⁻³].
 - `qᵥ`: water vapor specific humidity [kg kg⁻¹].
 - `ρw`: liquid water density [kg m⁻³].
 - `v_term`: function `v_term(D)` that returns the drop fall speed [m s⁻¹].
"""
@inline function drop_freezing_heat_timescale(vent, aps, tps, Tₐ, ρₐ, qᵥ, ρw, v_term::VT) where {VT}
    FT = UT.promote_typeof(Tₐ, ρₐ, qᵥ, ρw)
    (; K_therm, D_vapor) = aps
    T_freeze = TDI.T_freeze(tps)
    ΔT = T_freeze - Tₐ
    F_v = CO.ventilation_factor(vent, aps, v_term)

    # Vapor density deficit between the drop surface, saturated over liquid at T_freeze, and the air
    ρᵥ_sfc = ρₐ * TDI.p2q(tps, T_freeze, ρₐ, TDI.saturation_vapor_pressure_over_liquid(tps, T_freeze))
    Δρᵥ = max(ρᵥ_sfc - ρₐ * qᵥ, zero(FT))
    # Heat to shed per unit mass, from liquid at Tₐ to ice at T_freeze [J kg⁻¹]
    ΔL = max(TDI.Lf(tps, T_freeze) - TDI.cp_l(tps) * ΔT, zero(FT))

    Q = K_therm * ΔT + TDI.Lᵥ(tps, T_freeze) * D_vapor * Δρᵥ
    # x_drop ΔL / (2π D) with x_drop = (π/6) ρw D³
    c_heat = ρw * ΔL / 12
    function τ_heat(D)
        τ = c_heat * D^2 / (F_v(D) * Q)
        return ifelse(ΔT > 0, τ, zero(τ))
    end
    return τ_heat
end

"""
    liquid_freezing_moments(τ_heat, v_dend, ρw, J, n, D_at, quad)

Return the number and mass freezing rates `(; ∂ₜn_frz, ∂ₜq_frz)` of a drop population [kg⁻¹ s⁻¹, kg kg⁻¹ s⁻¹].

Each drop freezes at the rate `1 / (τ_nuc + τ_heat + τ_dend)`, the reciprocal of the times to nucleate, to shed
the latent heat, and to grow the freezing front across the drop. See the Ice Nucleation docs page.

# Arguments
 - `τ_heat`: function of the drop diameter, from [`drop_freezing_heat_timescale`](@ref).
 - `v_dend`: speed of the freezing front [m s⁻¹].
 - `ρw`: liquid water density [kg m⁻³].
 - `J`: the volumetric nucleation coefficient [m⁻³(water) s⁻¹], summed over the nucleation pathways.
 - `n`: specific number concentration of the population [kg⁻¹].
 - `D_at`: maps a quadrature node to a drop diameter [m].
 - `quad`: the [`Quadrature.GaussLaguerre`](@ref) rule whose weight function is the normalized size distribution.
"""
@inline function liquid_freezing_moments(
    τ_heat::TH, v_dend, ρw, J, n, D_at::DF, quad::QUAD.GaussLaguerre,
) where {TH, DF}
    FT = UT.promote_typeof(ρw, J, n)
    V_1 = FT(π) / 6
    ∂ₜn_frz = zero(FT)
    ∂ₜq_frz = zero(FT)
    for i in eachindex(quad.nodes)
        D = D_at(FT(quad.nodes[i]))
        w = FT(quad.weights[i])
        inv_τ_nuc = J * V_1 * D^3
        # Compare values, since a ForwardDiff dual with zero value and nonzero partials exceeds zero
        present = FD.value(inv_τ_nuc) > zero(FD.value(inv_τ_nuc))
        # Substitutes keep the discarded branch and its partials finite
        Ds = ifelse(present, D, one(D))
        inv_s = ifelse(present, inv_τ_nuc, one(inv_τ_nuc))
        τ_freeze = τ_heat(Ds) + Ds / v_dend
        r = ifelse(present, 1 / (1 / inv_s + τ_freeze), zero(inv_τ_nuc))
        ∂ₜn_frz += w * r
        ∂ₜq_frz += w * r * V_1 * ρw * Ds^3
    end
    return (; ∂ₜn_frz = n * ∂ₜn_frz, ∂ₜq_frz = n * ∂ₜq_frz)
end

"""
    rain_freezing_rate(mp, tps, micro, thermo)

Return the raindrop freezing rates from immersion (Bigg 1953) and homogeneous (Koop 2000) nucleation.

The freezing of each drop is limited by the time to shed its latent heat and to grow the freezing front across it;
see the Ice Nucleation docs page. The size distribution integral uses the rule `mp.ice.rain_freezing_quad`.
The number freezes at most at the fraction of the mass that freezes, `∂ₜn_frz ≤ n_rai ∂ₜq_frz / q_rai`, because
the two-moment size distribution has a fixed shape and cannot lose its small drops alone.

# Arguments
 - `mp`: the [`CMP.Microphysics2MParams`](@ref) parameters with P3 ice.
 - `tps`: thermodynamics parameters.
 - `micro`: the contents per kg of air; the rate reads `q_tot`, `q_lcl`, `q_rai`, `n_rai` and `q_ice`.
 - `thermo`: the air state; the rate reads `ρ` and `T`.

# Returns
 - `(; ∂ₜn_frz, ∂ₜq_frz)`, the specific number rate [kg⁻¹ s⁻¹] and mass rate [kg kg⁻¹ s⁻¹].
"""
@inline function rain_freezing_rate(mp, tps, micro, thermo)
    (; ρ, T) = thermo
    (; q_rai, n_rai) = micro
    FT = UT.promote_typeof(q_rai, n_rai, ρ, T)
    T_freeze = TDI.T_freeze(tps)
    pdf_r = mp.ice.rain_pdf
    aps = mp.warm_rain.air_properties
    (; α, β, ρ0) = mp.warm_rain.seifert_beheng.evap
    (; ρw) = pdf_r
    (; Dr_mean) = CM2.pdf_rain_parameters(pdf_r, q_rai, ρ, n_rai * ρ)
    qᵥ = TDI.q_vap(micro.q_tot, micro.q_lcl + q_rai, micro.q_ice)

    # Immersion and homogeneous nucleation act in parallel, so their coefficients add
    J = mp.ice.rain_freezing(T, T_freeze) + CM_HomIce.homogeneous_freezing_rate_coefficient(mp.ice.homogeneous, tps, T)

    # The SB2006 fall speed `α x^β √(ρ0/ρ)` of a drop of mass `x = (π/6) ρw D³`
    x_per_D³ = FT(π) / 6 * ρw
    v_term = D -> α * (x_per_D³ * D^3)^β * sqrt(ρ0 / ρ)
    τ_heat = drop_freezing_heat_timescale(mp.ice.scheme.vent, aps, tps, T, ρ, qᵥ, ρw, v_term)
    # Gauss-Laguerre in `u = D / Dr_mean`, whose weight is the size distribution
    (; ∂ₜn_frz, ∂ₜq_frz) = liquid_freezing_moments(
        τ_heat, mp.ice.dendritic_freezing.v_dend, ρw, J, n_rai, u -> Dr_mean * u, mp.ice.rain_freezing_quad,
    )
    ∂ₜn_frz = min(∂ₜn_frz, n_rai * UT.guarded_quotient(∂ₜq_frz, q_rai))

    cond = _liquid_freezing_is_active(q_rai, n_rai, T, T_freeze)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(FT))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(FT))
    return (; ∂ₜn_frz, ∂ₜq_frz)
end

"""
    cloud_freezing_rate(mp, tps, micro, thermo)

Return the cloud droplet freezing rates from immersion (Bigg 1953) and homogeneous (Koop 2000) nucleation.

The per-drop freezing rate of [`rain_freezing_rate`](@ref) is integrated over the SB2006 cloud droplet size
distribution, with the droplet fall speed in the Stokes regime, by the rule `mp.ice.cloud_freezing_quad`.
The INP budget of the target spectrum `mp.ice.ice_nucleation`, [`immersion_limit_rate`](@ref), bounds the
immersion coefficient so that the immersion number rate without the freezing time limits does not exceed it.
The homogeneous coefficient has no bound.

# Arguments
 - `mp`: the [`CMP.Microphysics2MParams`](@ref) parameters with P3 ice.
 - `tps`: thermodynamics parameters.
 - `micro`: the contents per kg of air; the rate reads `q_tot`, `q_lcl`, `n_lcl`, `q_rai`, `q_ice` and the
   entries the depletion model reads.
 - `thermo`: the air state; the rate reads `ρ` and `T`.

# Returns
 - `(; ∂ₜn_frz, ∂ₜq_frz)`, the specific number rate [kg⁻¹ s⁻¹] and mass rate [kg kg⁻¹ s⁻¹].
"""
@inline function cloud_freezing_rate(mp, tps, micro, thermo)
    (; ρ, T) = thermo
    (; q_lcl, n_lcl) = micro
    FT = UT.promote_typeof(q_lcl, n_lcl, ρ, T)
    T_freeze = TDI.T_freeze(tps)
    pdf_c = mp.ice.cloud_pdf
    aps = mp.warm_rain.air_properties
    (; ρw) = pdf_c
    (; λc, μcD) = CM2.pdf_cloud_parameters(pdf_c, q_lcl, ρ, n_lcl * ρ)
    qᵥ = TDI.q_vap(micro.q_tot, q_lcl + micro.q_rai, micro.q_ice)

    # The immersion number rate at coefficient `J` is `J q_lcl / ρw`, so the budget bounds the coefficient
    n_act = n_active(mp.ice.inp_depletion_model, micro)
    ∂ₜn_lim = immersion_limit_rate(mp.ice.ice_nucleation, T, ρ; n_active = n_act).∂ₜn_frz
    J_het_max = ifelse(q_lcl > 0, ∂ₜn_lim * ρw / q_lcl, oftype(∂ₜn_lim, Inf))
    J_het = min(mp.ice.rain_freezing(T, T_freeze), J_het_max)
    J = J_het + CM_HomIce.homogeneous_freezing_rate_coefficient(mp.ice.homogeneous, tps, T)

    stokes = CMP.StokesRegimeVelType(; ρw, ν_air = aps.ν_air, grav = TDI.grav(tps))
    v_stokes = CO.particle_terminal_velocity(stokes, ρ)
    τ_heat = drop_freezing_heat_timescale(mp.ice.scheme.vent, aps, tps, T, ρ, qᵥ, ρw, v_stokes)
    # Gauss-Laguerre in `t = λc D^μcD`, so `D = (t/λc)^(1/μcD)`
    (; ∂ₜn_frz, ∂ₜq_frz) = liquid_freezing_moments(
        τ_heat, mp.ice.dendritic_freezing.v_dend, ρw, J, n_lcl, t -> (t / λc)^(1 / μcD), mp.ice.cloud_freezing_quad,
    )

    cond = _liquid_freezing_is_active(q_lcl, n_lcl, T, T_freeze)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(FT))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(FT))
    return (; ∂ₜn_frz, ∂ₜq_frz)
end

# ---------------------------------------------------------------------------
# Deposition nucleation toward an INP target spectrum
# ---------------------------------------------------------------------------

"""
    is_active(inp, T, S_i)

Return whether deposition nucleation proceeds for the INP target spectrum `inp`.

`inp` is a [`CMP.INPTargetSpectrum`](@ref), at air temperature `T` [K] and ice supersaturation
`S_i = q_vap/q_sat_ice - 1`.
"""
function is_active end

"""
    delivery_rate(inp, mp, tps, T, S_i)

Return the rate [s⁻¹] at which the INP deficit of the target spectrum `inp` becomes ice crystals.

The deficit is that between the target INP concentration of `inp` and the activated concentration.
"""
function delivery_rate end

# Target concentrations [m⁻³]: Cooper (1986) with a ceiling, the median of the Frostenberg 2023 climatology, and
# the base target with its log shifted and with a ceiling
@inline ((; a, b, T₀, N_max)::CMP.ExponentialSupercoolingINP)(T) = min(a * exp(b * (T₀ - T)), N_max)
@inline (inp::CMP.Frostenberg2023)(T) = exp(INP_concentration_mean(inp, T))
@inline (inp::CMP.StochasticINPSpectrum)(T) = min(exp(log(inp.base(T)) + inp.log_shift), inp.N_max)

@inline is_active(inp::CMP.ExponentialSupercoolingINP, T, S_i) = (T < inp.T_thr) & (S_i >= inp.S_thr)
@inline is_active(inp::CMP.Frostenberg2023, T, S_i) = (T < inp.T_thr) & (S_i >= inp.S_thr)
@inline is_active(inp::CMP.StochasticINPSpectrum, T, S_i) = (T >= inp.T_min) & is_active(inp.base, T, S_i)

# The Cooper spectrum relaxes at the rate at which vapor deposition grows a crystal to the radius of a newly
# nucleated crystal, `r ∂ₜr = G_ice S_i / ρ_i`
@inline function delivery_rate(inp::CMP.ExponentialSupercoolingINP, mp, tps, T, S_i)
    (; ρ_i, D_nuc) = mp.ice.scheme
    r_nuc = D_nuc / 2
    return 2 * CO.G_func_ice(mp.warm_rain.air_properties, tps, T) * max(S_i, 0) / (ρ_i * r_nuc^2)
end
@inline delivery_rate(inp::CMP.Frostenberg2023, mp, tps, T, S_i) = inv(oftype(T, inp.τ))
@inline delivery_rate(inp::CMP.StochasticINPSpectrum, mp, tps, T, S_i) = delivery_rate(inp.base, mp, tps, T, S_i)

# The mass rate of nucleated crystals of mass `m_nuc` at the number rate `∂ₜn_frz`; the Frostenberg 2023 rate is
# at most half the vapor excess `q_excess` per relaxation time `1/k`
@inline _nucleated_mass_rate(::CMP.INPTargetSpectrum, m_nuc, ∂ₜn_frz, q_excess, k) = m_nuc * ∂ₜn_frz
@inline _nucleated_mass_rate(::CMP.Frostenberg2023, m_nuc, ∂ₜn_frz, q_excess, k) =
    min(m_nuc * ∂ₜn_frz, q_excess * k / 2)
@inline _nucleated_mass_rate(inp::CMP.StochasticINPSpectrum, m_nuc, ∂ₜn_frz, q_excess, k) =
    _nucleated_mass_rate(inp.base, m_nuc, ∂ₜn_frz, q_excess, k)

"""
    deposition_rate(mp, tps, micro, thermo)

Compute the deposition nucleation of pristine ice toward the INP target spectrum `mp.ice.ice_nucleation`.

The nucleation relaxes the activated concentration toward the target concentration `N_t(T) = inp(T)` of the
[`CMP.INPTargetSpectrum`](@ref) `inp`:

```
∂ₜn_frz = max(0, N_t(T)/ρ - n_active) ⋅ k   [kg⁻¹ s⁻¹]
∂ₜq_frz = m_nuc ⋅ ∂ₜn_frz                   [kg kg⁻¹ s⁻¹]
```

with `n_active` the activated concentration of the depletion model `mp.ice.inp_depletion_model` (see
[`n_active`](@ref)), `k` the relaxation rate of `inp`, and `m_nuc` the mass of a newly nucleated crystal, a
solid-ice sphere of the P3 diameter `D_nuc`. For [`CMP.Frostenberg2023`](@ref) and a spectrum based on it, the
mass rate is at most half the vapor excess over ice saturation per relaxation time. Both rates are zero outside
the activation conditions of `inp`. See the Ice Nucleation docs page.

# Arguments
 - `mp`: the [`CMP.Microphysics2MParams`](@ref) parameters with P3 ice.
 - `tps`: thermodynamics parameters.
 - `micro`: the contents per kg of air; the rate reads `q_tot`, `q_lcl`, `q_rai`, `q_ice` and the entries the
   depletion model reads.
 - `thermo`: the air state; the rate reads `ρ` and `T`.

# Returns
 - `(; ∂ₜn_frz, ∂ₜq_frz)`, the specific number rate [kg⁻¹ s⁻¹] and mass rate [kg kg⁻¹ s⁻¹].
"""
@inline function deposition_rate(mp, tps, micro, thermo)
    inp = mp.ice.ice_nucleation
    (; ρ, T) = thermo
    q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    q_vap = TDI.q_vap(micro.q_tot, micro.q_lcl + micro.q_rai, micro.q_ice)
    S_i = q_vap / q_sat_ice - 1
    n_act = n_active(mp.ice.inp_depletion_model, micro)
    m_nuc = CO.nucleated_ice_mass(mp.ice.scheme)
    k = delivery_rate(inp, mp, tps, T, S_i)
    ∂ₜn_frz = max(0, inp(T) / ρ - n_act) * k
    ∂ₜn_frz = ifelse(is_active(inp, T, S_i), ∂ₜn_frz, zero(∂ₜn_frz))
    ∂ₜq_frz = _nucleated_mass_rate(inp, m_nuc, ∂ₜn_frz, max(0, q_vap - q_sat_ice), k)
    return (; ∂ₜn_frz, ∂ₜq_frz)
end

"""
    immersion_limit_rate(inp, T, ρ; n_active = 0)

Return the upper limit on the cloud droplet immersion freezing number rate from the INP budget of `inp`.

The limit is `max(0, N_t(T)/ρ - n_active) / τ` [kg⁻¹ s⁻¹], as a `NamedTuple` `(; ∂ₜn_frz)`, with `N_t` the
target concentration and `τ` the relaxation timescale of the [`CMP.INPTargetSpectrum`](@ref) `inp`,
`ρ` the air density [kg m⁻³], and `n_active` the activated concentration [kg⁻¹].
For [`CMP.Frostenberg2023`](@ref), the limit is zero at and above `T_freeze`.
For a [`CMP.StochasticINPSpectrum`](@ref) on a Frostenberg 2023 base, `N_t` is the shifted target, and the limit
is also zero below `T_min`.
A spectrum without an INP budget, such as [`CMP.ExponentialSupercoolingINP`](@ref), imposes no limit,
`∂ₜn_frz = Inf`.
"""
@inline function immersion_limit_rate(
    inp::Union{CMP.Frostenberg2023, CMP.StochasticINPSpectrum{<:Any, <:CMP.Frostenberg2023}}, T, ρ;
    n_active = zero(T),
)
    FT = UT.promote_typeof(T, ρ, n_active)
    ∂ₜn_frz = max(zero(FT), inp(T) / ρ - n_active) / _budget_timescale(inp)
    return (; ∂ₜn_frz = ifelse(_immersion_budget_is_open(inp, T), ∂ₜn_frz, zero(FT)))
end

# A spectrum without an INP budget, such as the Cooper spectrum, imposes no limit
@inline immersion_limit_rate(::CMP.INPTargetSpectrum, T, ρ; n_active = zero(T)) =
    (; ∂ₜn_frz = UT.promote_typeof(T, ρ, n_active)(Inf))

@inline _budget_timescale(inp::CMP.Frostenberg2023) = inp.τ
@inline _budget_timescale(inp::CMP.StochasticINPSpectrum) = _budget_timescale(inp.base)

@inline _immersion_budget_is_open(inp::CMP.Frostenberg2023, T) = T < inp.T_freeze
@inline _immersion_budget_is_open(inp::CMP.StochasticINPSpectrum, T) =
    (T >= inp.T_min) & _immersion_budget_is_open(inp.base, T)

"""
    n_active(model, micro)

Return the specific concentration [kg⁻¹] of the ice nucleating particles already activated.

[`deposition_rate`](@ref) and the immersion budget subtract it from the target concentration.
`micro` holds the contents per kg of air.
For [`CMP.NIceProxyDepletion`](@ref), it is the specific ice number `micro.n_ice`.
"""
@inline n_active(::CMP.NIceProxyDepletion, micro) = micro.n_ice

end # end module
