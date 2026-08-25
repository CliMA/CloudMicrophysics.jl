# This file contains the `HetIceNucleation` and `HomIceNucleation` modules.

"""
Parameterization for homogeneous cloud ice nucleation
"""
module HomIceNucleation

import ..Parameters as CMP

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
function homogeneous_J_cubic((; c₁, c₂, c₃, c₄, Δa_w_min, Δa_w_max)::CMP.Koop2000, Δa_w::FT) where {FT}
    Δa_w_min ≤ Δa_w ≤ Δa_w_max || throw(
        DomainError(Δa_w,
            lazy"Change in water activity must be within the valid range: Δa_w ∈ [$Δa_w_min, $Δa_w_max], but Δa_w = $Δa_w",
        ),
    )
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

end # end module

"""
Parameterization for heterogenous cloud ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP
# `HomIceNucleation` is defined ABOVE this module in this file so that the Koop coefficient
# can be imported here by the ordinary relative path. The two modules are independent -
# `HomIceNucleation` imports only `..Parameters` - so the ordering is free to be chosen.
import ..HomIceNucleation as CM_HomIce
import CloudMicrophysics.Common as CO
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.Utilities as UT
import ForwardDiff as FD

export dust_activated_number_fraction
export MohlerDepositionRate
export deposition_J
export ABIFM_J
export P3_deposition_N_i
export P3_het_N_i
export INP_concentration_frequency
export INP_concentration_mean
export liquid_freezing_rate
export drop_freezing_heat_timescale
export drop_freezing_dendrite_timescale
export DENDRITE_GROWTH_VELOCITY
export PRANDTL_NUMBER_AIR
export homogeneous_freezing_rate_coefficient
export rain_freezing_rate
export rain_freezing_quadrature
export RAIN_FREEZING_QUADRATURE_ORDER
export cloud_freezing_rate
export cloud_freezing_quadrature
export CLOUD_FREEZING_QUADRATURE_ORDER
export immersion_limit_rate
export deposition_rate
export is_active
export delivery_rate
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
    + `log_a`: log of the INPC normalization coefficient `a` [m³]
  - `T`: air temperature [K]

Following Eq. (1) of [Frostenberg2023](@cite), `log(a · INPC)` is normally
distributed with mean `μ(T) = log(-(b · T_celsius / 10)^9)`, so the mean
`log(INPC)` returned here is
```
μ(T) - log(a) = 9 log(-b · T_celsius / 10) - log(a)
```
with the corresponding INPC obtained by exponentiating. The parameters `a` and
`b` are read from `ClimaParams`; at their defaults `a = b = 1` this reduces to
the marine-dataset curve `log((-T_celsius / 10)^9)`.

For details see: [Frostenberg2023](@cite), doi.org/10.5194/acp-23-10883-2023
"""
function INP_concentration_mean((; T_freeze, b, log_a)::CMP.Frostenberg2023, T)
    T_celsius = min(T - T_freeze, 0)
    return 9log(-b * T_celsius / 10) - log_a  # = log((-b * T_celsius / 10)^9) - log(a)
end

"""
    liquid_freezing_rate(opt, pdf, tps, q, ρ, N, T)

Compute the rate of liquid water freezing into ice.

# Arguments
 - `opt`: The [`CMP.RainFreezing`](@ref) parameterization.
 - `pdf`: The liquid water particle size distribution (PSD) PDF.
 - `tps`: Thermodynamics parameters.
 - `q`: Liquid water specific content [kg(water) kg⁻¹(air)].
 - `ρ`: Air density [kg(air) m⁻³(air)].
 - `N`: Liquid water number concentration [m⁻³(air)].
 - `T`: Air temperature [K].

# Returns
 - A `NamedTuple` with the fields:
    + `∂ₜn_frz`: Specific number freezing rate [kg⁻¹(air) s⁻¹].
    + `∂ₜq_frz`: Specific mass freezing rate [kg(water) kg⁻¹(air) s⁻¹].
"""
function liquid_freezing_rate(opt::CMP.RainFreezing, pdf, tps, q, ρ, N, T)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    # Bigg (1953) volumetric freezing rate [m⁻³(water) s⁻¹]
    return _liquid_freezing_rate_from_J(pdf, opt(T, T_freeze), q, ρ, N, T, T_freeze)
end

"""
    _liquid_freezing_is_active(::Type{FT}, q, n, T, T_freeze)

Whether heterogeneous freezing of a liquid population can proceed: the population is present
in both moments and the state is supercooled. Shared by every liquid-freezing entry point so
that they agree on where freezing occurs.
"""
@inline function _liquid_freezing_is_active(::Type{FT}, q, n, T, T_freeze) where {FT}
    ϵₘ, ϵₙ = UT.ϵ_numerics_2M_M(FT), UT.ϵ_numerics_2M_N(FT)
    return (n > ϵₙ) & (q > ϵₘ) & (T < T_freeze)
end

"""
    _liquid_freezing_rate_from_J(pdf, J, q, ρ, N, T, T_freeze)

The bulk raindrop freezing rate driven by an arbitrary VOLUMETRIC freezing rate coefficient `J`
[m⁻³(water) s⁻¹], integrated over the SB2006 exponential raindrop PSD.

Factored out of [`liquid_freezing_rate`](@ref) so that more than one nucleation pathway can share
one PSD treatment. Nucleation pathways act in PARALLEL on the same drop, so their volumetric
coefficients ADD, and because the bulk rate is exactly linear in `J` the sum can be formed before
the PSD integral rather than after:

```
∂ₜ(n, q)|total = (J_bigg + J_koop) ⋅ f(PSD)
```

That linearity is the reason this is a factoring and not a rewrite. It also rules out the
tempting shortcut of scaling one pathway's answer by `J_koop / J_bigg`: that divides by a
quantity which overflows `Float32` near `ΔT = 132` K, and it would silently reintroduce the
Bigg exponential into a term that must not depend on it.

The arithmetic and its association are preserved exactly from the original, so the Bigg-only
path through [`liquid_freezing_rate`](@ref) is bit-identical to what it was before the factoring.

# Arguments
 - `pdf`: the [`CMP.RainParticlePDF_SB2006`](@ref) raindrop size distribution.
 - `J`: total volumetric freezing rate coefficient [m⁻³(water) s⁻¹].
 - `q`: rain specific content [kg(water) kg⁻¹(air)].
 - `ρ`: air density [kg(air) m⁻³(air)].
 - `N`: raindrop number concentration [m⁻³(air)].
 - `T`: air temperature [K].
 - `T_freeze`: freezing temperature [K].

# Returns
 - `(; ∂ₜn_frz, ∂ₜq_frz)`, specific number and mass freezing rates.
"""
@inline function _liquid_freezing_rate_from_J(pdf, J, q, ρ, N, T, T_freeze)
    FT = float(UT.promote_typeof(q, ρ, N, T))
    (; ρw) = pdf  # [kg(water) m⁻³(water)]
    n = N / ρ     # specific number concentration [kg⁻¹(air)]

    # Solve for the pdf parameters
    (; Dr_mean) = CM2.pdf_rain_parameters(pdf, q, ρ, N)

    # Diameter PSD moments via exponential_Mⁿ:
    #   M_D^k = ∫ D^k n(D) dD
    # The freezing probability per unit time for a single drop of diameter D
    # is J_drop(D) = J * V(D) = J * (π/6) * D³  [s⁻¹].

    # Number freezing rate per kg air:  ∂n/∂t = ∫ J_drop(D) n(D) dD
    #                                         = J * (π/6) * M_D³
    M_D³ = DT.exponential_Mⁿ(Dr_mean, n, 3)  # [m³ · kg⁻¹(air)]

    # Mass freezing rate per kg air:    ∂q/∂t = ∫ x(D) * J_drop(D) n(D) dD
    #                                         = J * ρw * (π/6)² * M_D⁶
    M_D⁶ = DT.exponential_Mⁿ(Dr_mean, n, 6)  # [m⁶ · kg⁻¹(air)]

    # Association preserved from the Bigg-only original, so that path stays bit-identical.
    V_1 = FT(π) / 6
    ∂ₜn_frz = J * V_1 * M_D³         # [kg⁻¹(air) s⁻¹]   — specific
    ∂ₜq_frz = J * ρw * V_1^2 * M_D⁶  # [kg(water) kg⁻¹(air) s⁻¹]  — specific

    # Return the computed rate only if N and L are (essentially) non-zero, and T is colder than -4°C.
    # Otherwise, return zero.
    cond = _liquid_freezing_is_active(FT, q, n, T, T_freeze)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(FT))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(FT))

    return (; ∂ₜn_frz, ∂ₜq_frz)
end

"""
    liquid_freezing_rate(
        opt::CMP.RainFreezing, pdf::CMP.CloudParticlePDF_SB2006,
        tps, q, ρ, N, T,
    )

Compute the rate of cloud-droplet immersion freezing into ice using the same
Bigg (1953) kinetics as the rain version, but integrated over the
generalized-gamma cloud-droplet PSD (SB2006).

The cloud PSD in diameter is

```
n(D) = N₀c · D^νcD · exp(-λc · D^μcD)   ,   νcD = 3νc + 2,  μcD = 3μc.
```

Bigg's per-drop freezing probability is `J_bigg(T) · (π/6) · D³`. Integrating
against the PSD gives closed-form number- and mass-freezing rates:

```
∂ₜn_frz = J_bigg ·      (π/6)  · M_D³(N₀c, λc, νcD, μcD)
∂ₜq_frz = J_bigg · ρw · (π/6)² · M_D⁶(N₀c, λc, νcD, μcD)
```

with `M_Dᵏ` the kth diameter moment computed by `generalized_gamma_Mⁿ`.
Volume-weighting → bigger drops freeze first

# Arguments
 - `opt`: The [`CMP.RainFreezing`](@ref) parameterization
   (Bigg / Barklie-Gokhale parameters; despite the type name, the kinetics
   apply to any liquid-drop PSD).
 - `pdf`: The [`CMP.CloudParticlePDF_SB2006`](@ref) cloud-droplet PSD.
 - `tps`: Thermodynamics parameters.
 - `q`: Cloud-liquid specific content [kg(water) kg⁻¹(air)].
 - `ρ`: Air density [kg(air) m⁻³(air)].
 - `N`: Cloud-droplet number concentration [m⁻³(air)].
 - `T`: Air temperature [K].

# Returns
 - A `NamedTuple` with the fields:
    + `∂ₜn_frz`: Specific number freezing rate [kg⁻¹(air) s⁻¹].
    + `∂ₜq_frz`: Specific mass freezing rate [kg(water) kg⁻¹(air) s⁻¹].
"""
function liquid_freezing_rate(
    opt::CMP.RainFreezing, pdf::CMP.CloudParticlePDF_SB2006,
    tps, q, ρ, N, T,
)
    FT = float(UT.promote_typeof(q, ρ, N, T))
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    (; ρw) = pdf
    n = N / ρ     # specific number concentration [kg⁻¹(air)]

    # Solve for the diameter-space PDF parameters.
    (; λc, νcD, μcD) = CM2.pdf_cloud_parameters(pdf, q, ρ, N)

    # Bigg (1953) volumetric freezing rate per unit drop water volume.
    J_bigg = opt(T, T_freeze)  # [m⁻³(water) s⁻¹]

    # Diameter moments via the closed-form generalized-gamma formula.
    # When N or L is essentially zero, `pdf_cloud_parameters` returns
    # `λc = Inf`, which makes `B^(-k/μ) → 0`, so the moments vanish
    M_D³ = DT.generalized_gamma_Mⁿ(νcD, μcD, λc, n, 3)  # [m³ · kg⁻¹(air)]
    M_D⁶ = DT.generalized_gamma_Mⁿ(νcD, μcD, λc, n, 6)  # [m⁶ · kg⁻¹(air)]

    V_1 = FT(π) / 6
    ∂ₜn_frz = J_bigg * V_1 * M_D³          # [kg⁻¹(air) s⁻¹]
    ∂ₜq_frz = J_bigg * ρw * V_1^2 * M_D⁶   # [kg(water) kg⁻¹(air) s⁻¹]

    # Check non-trivial number/mass and that T < -4 °C.
    cond = _liquid_freezing_is_active(FT, q, n, T, T_freeze)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(FT))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(FT))

    return (; ∂ₜn_frz, ∂ₜq_frz)
end

"""
    PRANDTL_NUMBER_AIR(FT)

The Prandtl number of air, `ν_air / α_therm` with `α_therm = K_therm / (ρₐ cp)` the thermal
diffusivity. It is the dimensionless group belonging in the CONDUCTION half of a ventilated
transfer balance, where [`CMP.AirProperties`](@ref) otherwise offers only the Schmidt number
`ν_air / D_vapor`, which belongs in the VAPOUR half.

**Why it is a constant and not computed from `aps`.** Momentum and thermal diffusivity carry
the same `1/ρₐ` scaling, so their ratio is a material property of air and sits at 0.71 to 0.72
over the whole atmospheric range. `AirProperties` stores `ν_air` and `K_therm` as fixed
numbers and carries no air density, so forming `ν_air ρₐ cp_d / K_therm` from them would
attach a spurious `ρₐ` dependence to a quantity that is physically constant: with the shipped
values it returns 0.80 at `ρₐ = 1.2` and 0.47 at `ρₐ = 0.7`, the drift coming entirely from
`ν_air` being held fixed while `α_therm` is not. A constant is the honest reading of the
constants the scheme actually has.

**Measured, because it decides how much retiring the Schmidt-for-Prandtl borrow buys.** The
shipped `ν_air = 1.6e-5` and `D_vapor = 2.26e-5` give `N_sc = 0.708`, so
`cbrt(N_pr)/cbrt(N_sc) = 1.00083`: the borrow the docstrings recorded as an approximation is
inert to 0.08 % in the ventilation factor, not the "under 8 %" previously claimed. Retiring it
is a correctness-of-form change with no measurable magnitude, and it is done because the two
halves of the balance should not share one dimensionless group by accident.

Deliberately not a ClimaParams key while the freezing form is under review, like
[`DENDRITE_GROWTH_VELOCITY`](@ref).
"""
@inline PRANDTL_NUMBER_AIR(::Type{FT}) where {FT} = FT(0.71)

"""
    drop_freezing_heat_timescale(vent, aps, tps, Tₐ, ρₐ, qᵥ, x_drop, ρw, v_drop)

The time a supercooled liquid drop needs to finish freezing once nucleation has occurred,
set by how fast the latent heat of fusion can be disposed of to the air, by conduction and
by evaporation together.

# Arguments
 - `vent`: ventilation coefficients, [`CMP.VentilationFactor`](@ref).
 - `aps`: air properties, [`CMP.AirProperties`](@ref) (uses `K_therm`, `ν_air`, `D_vapor`).
 - `tps`: thermodynamics parameters.
 - `Tₐ`: air temperature [K].
 - `ρₐ`: air density [kg(air) m⁻³(air)].
 - `qᵥ`: ambient water vapour specific humidity [kg(water) kg⁻¹(air)].
 - `x_drop`: drop mass [kg(water)].
 - `ρw`: liquid water density [kg(water) m⁻³(water)].
 - `v_drop`: drop fall speed [m s⁻¹], for the ventilation Reynolds number.

# Returns
 - `τ_heat`: the heat-dissipation-limited freezing time of one drop [s].

Nucleating a drop and freezing it through are different things. Bigg (1953) kinetics give the
first; this gives the second, in two stages.

**Recalescence.** The drop sits at `Tₐ = T_freeze - ΔT` and its own heat capacity absorbs
latent heat until the drop reaches `T_freeze`. Freezing a mass fraction `x_ad` releases
`x_ad Lf` per unit drop mass against a warming cost of `c_w ΔT`, so

```
x_ad = c_w ΔT / Lf(Tₐ).
```

This stage is dendritic, takes milliseconds, and is treated as instantaneous.

**Heat-transfer-limited freezing of the remainder.** The drop surface is now at `T_freeze`,
mixed-phase equilibrium, while the air is `ΔT` colder AND drier. Latent heat leaves across the
surface by two ventilated channels in parallel: conduction down `ΔT`, and evaporation down the
vapour-density deficit `Δρᵥ` between the melting-point surface and the environment. With
`h = Nu K_therm / D`, `Nu = 2 F_v` and the matching mass-transfer coefficient, the export is
`Q = 2π D (F_v,heat K_therm ΔT + F_v,vap Lᵥ D_vapor Δρᵥ)`. The heat still to be shed is
`x_drop (1 - x_ad) Lf`, giving

```
τ_heat = x_drop max(Lf(Tₐ) - c_w ΔT, 0) / (2π D (F_v,h K_therm ΔT + F_v,v Lᵥ D_vapor Δρᵥ)),
Δρᵥ    = max(ρᵥ,sat(T_freeze, over liquid) - ρₐ qᵥ, 0),
```

with `D` the volume-equivalent drop diameter. `Nu = 2 F_v` is this codebase's convention:
`Common.ventilation_factor` tends to `aᵥ ≈ 0.78` for a still drop and `2π D K_therm ΔT`
is the conduction solution for a sphere. This is the full Musil (1970, Eq. A7) balance that
`P3Scheme.compute_max_freeze_rate` already applies to riming, here in its self-freezing
counterpart. The two halves take their own dimensionless group -
[`PRANDTL_NUMBER_AIR`](@ref) for conduction, the Schmidt number for vapour - retiring the
Schmidt-for-Prandtl borrow the conduction-only form carried; that retirement is measured to be
worth 0.08 %, so the magnitude here comes from the vapour term and not from the ventilation.

**The evaporated mass is NOT routed to vapour, deliberately.** The term enters the TIMESCALE
only: the conversion still sends 100 % of the freezing mass from rain to ice, with no
rain-to-vapour side channel. The neglected fraction is `s (1-x_ad) Lf / Lᵥ`, about 3 to 4 % of
the frozen mass with `s ≈ 0.3` the vapour share of the export, so the cell's sensible heating
is overstated by `s (1-x_ad) Lf` per kilogram frozen - about 0.007 K at `Δq = 1e-4` - and the
vapour budget misses ~4e-6 kg/kg against an ambient 1e-3 to 1e-2. That misattribution is the
same magnitude class as the single-temperature recalescence transient the scheme already
accepts, it is documented rather than hidden, and the complete alternative (a 96/4 split flux)
is recorded as a later refinement.

The numerator carries the whole physics content: `τ_heat` DECREASES with supercooling and
reaches exactly zero at

```
ΔT* = Lf(T_freeze) / (2 c_w - c_i) ≈ 53 K   (Tₐ ≈ 220 K),
```

where the drop's own cold absorbs all of the latent heat and freezing runs at dendrite speed
with nothing to export. The homogeneous-freezing limit falls out of the balance instead of
being a threshold branch. `ΔT*` is 53 K and not `Lf/c_w ≈ 80` K because `Lf` is itself
temperature dependent, `Lf(Tₐ) = Lf(T_freeze) - (c_w - c_i) ΔT`;
`P3Scheme.compute_max_freeze_rate` carries the same 53 K in its own denominator. The vapour
term sits in the DENOMINATOR, so it moves no threshold: `ΔT*` and the deep-supercooled
behaviour are untouched by it, bit for bit.

Read as a bound this is a physics change inside the ordinary mixed-phase range, not a
pathological-tail guard: `τ_heat` overtakes the Bigg nucleation time near `ΔT ≈ 20` K for
millimetre drops. It is also, and for the same reason, no bound at all colder than ~220 K.
`diag/rainfrz_tauheat_magnitudes.jl` measures both statements.

Total on degenerate states, with no new division hazard. `x_drop / D` is evaluated as
`cbrt(x_drop² π ρw / 6)` so that it is zero rather than `0/0` on an empty rain state, and
`ΔT ≤ 0` returns zero, which is the no-limiting value and where the freezing rate is zero
regardless. Whenever `ΔT > 0` the conduction term alone is already strictly positive, and
`Δρᵥ` is clamped at zero and only added, so the denominator cannot be driven to zero or made
to flip the sign of the timescale. The clamp itself is for pathological inputs only: the
surface holds `ρᵥ,sat(0 °C) = 4.85e-3` kg m⁻³ against a liquid-saturated environment that
carries less at every `Tₐ < T_freeze`, so an unclamped `Δρᵥ < 0` needs an ambient liquid
supersaturation ratio above 2.1 at `ΔT = 10` K and above 6.9 at `ΔT = 25` K - unreachable in
cloud, but reachable through a corrupt `q_tot`, where a negative `Δρᵥ` would turn evaporation
into a heat SOURCE and shorten the freezing time without limit.
"""
@inline function drop_freezing_heat_timescale(vent, aps, tps, Tₐ, ρₐ, qᵥ, x_drop, ρw, v_drop)
    FT = UT.promote_typeof(Tₐ, ρₐ, qᵥ, x_drop, ρw, v_drop)
    (; K_therm, D_vapor) = aps
    T_freeze = TDI.T_freeze(tps)
    c_w = TDI.cp_l(tps)
    L_f = TDI.Lf(tps, Tₐ)
    L_v = TDI.Lᵥ(tps, Tₐ)
    ΔT = T_freeze - Tₐ

    D = cbrt(6 * x_drop / (FT(π) * ρw))  # volume-equivalent drop diameter [m]
    v_term = _ -> v_drop
    F_v_heat = CO.ventilation_factor(vent, aps, v_term, PRANDTL_NUMBER_AIR(FT))(D)
    F_v_vap = CO.ventilation_factor(vent, aps, v_term)(D)  # Schmidt number

    # Vapour density deficit between the melting-point drop surface and the environment,
    # built the same way `P3Scheme.compute_max_freeze_rate` builds its own Δρᵥ_sat: a
    # saturation pressure converted through `p2q` at the air density and multiplied back by
    # it. Over LIQUID at the surface, because the wet-growth surface is an ice-water mixture;
    # against the AMBIENT vapour rather than a saturation value, which is what gives the
    # balance its humidity dependence - drier air freezes drops faster.
    ρᵥ_sfc = ρₐ * TDI.p2q(tps, T_freeze, ρₐ, TDI.saturation_vapor_pressure_over_liquid(tps, T_freeze))
    Δρᵥ = max(ρᵥ_sfc - ρₐ * qᵥ, zero(FT))

    # Latent heat left to shed after recalescence, per unit drop mass [J kg⁻¹]. Zero at and
    # below ΔT*, where the drop's own cold covers all of it.
    ΔL = max(L_f - c_w * ΔT, zero(FT))
    # x_drop / D, written to vanish rather than divide zero by zero on an empty rain state
    x_over_D = cbrt(x_drop^2 * FT(π) * ρw / 6)
    ∂ₜQ = F_v_heat * K_therm * ΔT + F_v_vap * L_v * D_vapor * Δρᵥ
    τ_heat = x_over_D * ΔL / (2 * FT(π) * ∂ₜQ)
    return ifelse(ΔT > 0, τ_heat, zero(FT))
end

"""
    DENDRITE_GROWTH_VELOCITY(FT)

Dendritic ice growth velocity in supercooled water [m s⁻¹], the speed at which the freezing
front crosses a drop once nucleation has occurred.

**This is the one tunable constant of the heat-limited freezing prototype**, and it is
deliberately not a ClimaParams key while the form is under review.

Measured dendrite tip velocities in supercooled water rise steeply with supercooling and then
plateau at the kinetic limit, of order 0.1 to 1 m s⁻¹. `0.5` m s⁻¹ sits in the middle of that
plateau. The value matters far less than it looks:
[`drop_freezing_dendrite_timescale`](@ref) only becomes the largest of the three composed
timescales at deep supercooling, past `ΔT*`, where
[`drop_freezing_heat_timescale`](@ref) has already vanished. Anywhere in the plausible range
it gives `τ_dend` of order 1 to 10 ms for a raindrop, so `h/τ_dend` is 200 to 2000 at a 2 s
step and the linearized-implicit update converts 99.5 to 99.95 % of the drop population
either way. Its job is to remove roughly thirteen orders of absurdity from `f`, not to time
the conversion, and no choice inside the plateau changes the converted fraction meaningfully.
"""
@inline DENDRITE_GROWTH_VELOCITY(::Type{FT}) where {FT} = FT(0.5)

"""
    drop_freezing_dendrite_timescale(x_drop, ρw)

The time the freezing front needs to cross a drop of mass `x_drop` [kg] at the dendritic
growth velocity, `τ_dend = D / v_dend`, with `D` the volume-equivalent diameter and `v_dend`
from [`DENDRITE_GROWTH_VELOCITY`](@ref).

Even a drop cold enough to absorb all of its own latent heat of fusion does not become ice
instantaneously: the dendrites still have to propagate through it. This is the third and
irreducible sequential stage of freezing, and past `ΔT*` - where
[`drop_freezing_heat_timescale`](@ref) is exactly zero because there is no heat left to
export - it is the ONLY stage bounding the conversion. That is precisely the deep-supercooling
regime in which the unbounded Bigg rate reaches 2.9e13 kg/kg/s.

Exactly zero on an empty rain state, where `x_drop` and hence `D` are zero.

Computes its own `D` rather than sharing one with `drop_freezing_heat_timescale`, which costs
a second `cbrt` (libm `cbrt` is not reliably common-subexpression-eliminated across uses - see
the note in `Microphysics2M.rain_evaporation`). Kept separate deliberately: this is a
prototype whose two stages have to be readable and testable in isolation, and the rain
freezing branch already pays for a PSD solve and an `exp`. Threading one `D` through both is
the obvious optimization if the form is adopted.
"""
@inline function drop_freezing_dendrite_timescale(x_drop, ρw)
    FT = UT.promote_typeof(x_drop, ρw)
    D = cbrt(6 * x_drop / (FT(π) * ρw))
    return D / DENDRITE_GROWTH_VELOCITY(FT)
end

"""
    homogeneous_freezing_rate_coefficient(hom, tps, T)

The Koop (2000) homogeneous freezing rate coefficient `J_koop` [m⁻³(water) s⁻¹] for DILUTE
water, with the fitted window handled at both ends.

# Arguments
 - `hom`: the [`CMP.Koop2000`](@ref) parameterization.
 - `tps`: thermodynamics parameters.
 - `T`: air temperature [K].

# Returns
 - `J_koop`: volumetric homogeneous freezing rate coefficient [m⁻³(water) s⁻¹], zero warmer
   than the fitted window.

Koop's rate is a function of the water-activity difference `Δa_w`. For a dilute drop the
solution water activity is 1, so

```
Δa_w(T) = 1 - e_sat,ice(T) / e_sat,liq(T),
```

both from the thermodynamics interface. Raindrops are dilute, which is what makes the pure-water
limit the right one here.

**Units are already SI and must not be converted again.** `homogeneous_J_cubic` returns
`10^(logJ + 6)`, and the `+6` is exactly the cm⁻³ to m⁻³ conversion of Koop's original fit; its
docstring and `docs/src/plots/linear_HOM_J.jl` both say so. Multiplying by another `1e6` here
would be invisible until someone compared against the fitted window.

**`homogeneous_J_cubic` THROWS a `DomainError` outside `[Δa_w_min, Δa_w_max]`.** That is not a
silent clamp to guard against, it is a GPU-kernel abort of exactly the family this campaign has
been chasing, and the error message interpolates values. So `Δa_w` is clamped into the window
BEFORE the call and the call can never throw, rather than being guarded after the fact.

The two ends of the window are treated differently, and neither choice is arbitrary:

  - **Warm side (`Δa_w < Δa_w_min`): hard ZERO, not the edge value.** Holding a fit flat outside
    its fitted range is the precise defect this thread exists to remove - it is what Bigg does at
    82 K of supercooling. Physically, homogeneous freezing genuinely vanishes at warm
    temperatures rather than plateauing. Measured, the choice is also numerically free: the
    lower-edge `J` is `4.2e2` m⁻³ s⁻¹, which for a 1.15 mm drop is `τ_nuc ≈ 3e6` s, about a
    month, against `τ_nuc ≈ 2` s from Bigg at the same temperature. Held flat it would overtake
    Bigg only warmer than `ΔT ≈ 1.2` K, inside the `T < T_freeze - 4` gate where all freezing is
    already zero. So the hard zero is made on principle, and it is inert either way.

  - **Cold side (`Δa_w > Δa_w_max`): CLAMP at the upper edge,** and this is provably inert rather
    than merely small. The clamped `J` is `2.9e24` m⁻³ s⁻¹, giving `τ_nuc ≈ 4e-16` s for a 1.15 mm
    drop against `τ_dend ≈ 2.3e-3` s - thirteen orders apart. Any larger `J` only shortens a
    `τ_nuc` that is already negligible in `τ_eff = τ_nuc + τ_heat + τ_dend`, so the dendrite stage
    owns the answer no matter what the extrapolation would have said. Same self-protecting
    structure as `τ_heat` vanishing at `ΔT*`: the composition is insensitive exactly where the
    parameterization is least trustworthy.

`e_sat,liq` underflows to zero before `e_sat,ice` does at `Float32`, so the guarded quotient
gives `0/floatmin = 0` and `Δa_w = 1`, which clamps to the cold edge - the correct branch.
"""
@inline function homogeneous_freezing_rate_coefficient(hom::CMP.Koop2000, tps, T)
    (; Δa_w_min, Δa_w_max) = hom
    e_i = TDI.saturation_vapor_pressure_over_ice(tps, T)
    e_l = TDI.saturation_vapor_pressure_over_liquid(tps, T)
    # `e_l` underflows to zero before `e_i` does, so the select gives Δa_w = 1 there, which clamps
    # to the cold edge - the correct branch. Written as a select rather than `max(e_l, floatmin)`
    # because `T` is a `Dual` under the temperature-coupled substep and `floatmin` of a dual type
    # is not something to rely on.
    Δa_w = 1 - ifelse(e_l > 0, e_i / e_l, zero(e_l))
    # Clamped BEFORE the call: `homogeneous_J_cubic` throws outside the window.
    J = CM_HomIce.homogeneous_J_cubic(hom, clamp(Δa_w, Δa_w_min, Δa_w_max))
    return ifelse(Δa_w < Δa_w_min, zero(J), J)
end

"""
    _composed_liquid_freezing_moments(vent, aps, tps, T, ρ, qᵥ, ρw, J, n, D_at, v_at, quad)

The two freezing moments of ANY liquid drop population, composed per size and integrated against
that population's own size distribution.

# Arguments
 - `vent`, `aps`, `tps`: ventilation coefficients, air properties, thermodynamics parameters.
 - `T`, `ρ`, `qᵥ`: air temperature [K], density [kg m⁻³] and vapour specific humidity.
 - `ρw`: liquid water density [kg m⁻³].
 - `J`: the SUMMED volumetric nucleation coefficient [m⁻³(water) s⁻¹], heterogeneous plus
   homogeneous, with any budget already applied to the heterogeneous part.
 - `n`: specific number concentration of the population [kg⁻¹(air)].
 - `D_at`: maps a quadrature node to a drop diameter [m]. This is where the size distribution
   enters, and it is the ONLY thing that differs between categories.
 - `v_at`: `(D, x_drop) -> v` [m s⁻¹], the fall speed law for this size regime, for the
   ventilation Reynolds number.
 - `quad`: the rule whose weight function IS the distribution; see the note below.

# Returns
 - `(; ∂ₜn_frz, ∂ₜq_frz)`, specific number and mass freezing rates, ungated.

**This is the categorization-blindness principle in code.** Cloud drops and rain drops are not
different substances, they are different sizes, so the same composition governs both and size
does the differentiating work:

```
r(D)     = 1 / (τ_nuc(D) + τ_heat(D) + τ_dend(D)),    1/τ_nuc(D) = J (π/6) D³
∂ₜn_frz  = ∫ r(D) n(D) dD
∂ₜq_frz  = ∫ x(D) r(D) n(D) dD,                       x(D) = (π/6) ρw D³
```

Only `D_at`, `v_at` and the rule change between categories. Nothing here is category-aware, and
in particular the kinetic series is never omitted by hand for small drops: at cloud sizes it
evaluates to `1` to within `1e-6` on its own wherever the heterogeneous pathway drives freezing,
and it correctly stops evaluating to `1` where the homogeneous pathway takes over (see
[`cloud_freezing_rate`](@ref)).

**Both distributions normalize their own quadrature weight to one.** For the SB2006 exponential
rain PSD, `u = D/Dr_mean` gives `n(D) dD = n exp(-u) du`; for the SB2006 generalized-gamma cloud
PSD, `t = λc D^μcD` gives `n(D) dD = n t^α exp(-t) dt` with `α = (νc+1)/μc - 1`. In both cases
the weight integrates to one, so the moments are exactly `n Σᵢ wᵢ g(uᵢ)` with no normalization
factor - which is also why a mistyped node table would show up as a broken zero-limiting limit
rather than as a silent rescaling. The tests assert `Σ wᵢ = 1` for both rules.

`r(D)` is written as one reciprocal of the summed timescales rather than as a nucleation rate
times a retained fraction, because that form takes the right limit where `J` has overflowed:
`τ_nuc` is zero and `r` is the finite `1/τ_freeze`, while `∂ₜn|nuc/(1 + τ_f/τ_nuc)` would be
`Inf/Inf`.

**It needs exactly one select, and the select is a PRESENCE gate on the nucleation pathway.**
`1/(J (π/6) D³)` has the right VALUE at `J (π/6) D³ = 0` - `Inf`, so `r` is zero - and NaN
PARTIALS, because the derivative of `1/x` at zero is `-∂x/x²` and a zero over a zero is not a
number. The value is guarded and the partials are not, which is the failure this select exists
to remove. Two state families reach it, and they are both physical rather than pathological:

  - **no nucleation pathway is open.** `J_het` is exactly zero once the INP budget is exhausted
    (`n_active ≥ INPC/ρ`, which at trace ice loadings is the ordinary case, not the edge case -
    at the `ad_compat_tests` cloud-edge state `n_ice = 30` kg⁻¹ against a budget of 1.6 kg⁻¹),
    and `J_koop` is a hard zero warmer than its fitted window. Warm of that window with the
    budget spent, the summed coefficient is exactly zero: no drop nucleates, so the freezing
    rate is exactly zero, in value AND in partials.
  - **no population is present.** `D` is zero on an empty state, so `J (π/6) D³` is zero for the
    same reason, and the fall-speed power law `α x^β` additionally has an infinite derivative
    at `x = 0`. Both are discarded by the same select.

Gated on presence rather than floored in magnitude: a floor under `J (π/6) D³` would make the
rate merely small instead of absent, and it would hand the state the partials of a CONSTANT,
which is a wrong derivative rather than a missing one. `zero(inv_τ_nuc)` carries a zero value
and zero partials, which is the correct derivative of a rate that is identically zero over a
neighbourhood of the state.

Each node rebuilds `D` from `x` inside the two timescale functions, and each rebuilds the
node-invariant parts of the wet-growth balance (the surface vapour density, the latent heats, the
two dimensionless groups' cube roots). Deliberate: this is a prototype whose stages have to stay
readable and testable in isolation, and the timescale functions are the tested objects. Hoisting
the node-invariant work and threading one `D` per node is the obvious optimization if the form is
adopted.
"""
@inline function _composed_liquid_freezing_moments(
    vent, aps, tps, T, ρ, qᵥ, ρw, J, n, D_at::DF, v_at::VF, quad,
) where {DF, VF}
    FT = UT.promote_typeof(T, ρ, qᵥ, ρw, n)
    V_1 = FT(π) / 6
    ∂ₜn_frz = zero(FT)
    ∂ₜq_frz = zero(FT)
    for i in eachindex(quad.nodes)
        u = FT(quad.nodes[i])
        w = FT(quad.weights[i])
        D = D_at(u)
        # ONE presence predicate for the whole node, read on the VALUE LANE. `inv_τ_nuc > 0`
        # would NOT do: `ForwardDiff` orders `Dual`s lexicographically, so a value of zero with
        # live partials compares GREATER than zero and the guard selects the present branch,
        # running the very reciprocal it exists to avoid. See `UT.guarded_quotient`.
        inv_τ_nuc = J * V_1 * D^3
        present = FD.value(inv_τ_nuc) > zero(FD.value(inv_τ_nuc))
        # Benign substitutes BEFORE the arithmetic, so nothing non-finite is ever CONSTRUCTED on
        # the absent branch - selecting away from a NaN is only as safe as the predicate. Two are
        # needed and they are not the same hazard: `1/inv_τ_nuc` is `Inf` with NaN partials where
        # no pathway is open, and the rain fall-speed law `α x^β` has an INFINITE derivative at
        # `x = 0` where no population is present, which would poison `τ_freeze` on its own.
        # On the present branch both substitutes are the identity, so the rate is unchanged.
        Ds = ifelse(present, D, one(D))
        inv_s = ifelse(present, inv_τ_nuc, one(inv_τ_nuc))
        x_drop = V_1 * ρw * Ds^3
        v_drop = v_at(Ds, x_drop)
        τ_freeze =
            drop_freezing_heat_timescale(vent, aps, tps, T, ρ, qᵥ, x_drop, ρw, v_drop) +
            drop_freezing_dendrite_timescale(x_drop, ρw)
        r = ifelse(present, 1 / (1 / inv_s + τ_freeze), zero(inv_τ_nuc))
        ∂ₜn_frz += w * r
        ∂ₜq_frz += w * r * x_drop
    end
    return (; ∂ₜn_frz = n * ∂ₜn_frz, ∂ₜq_frz = n * ∂ₜq_frz)
end

"""
    RAIN_FREEZING_QUADRATURE_ORDER

Number of Gauss-Laguerre nodes [`rain_freezing_rate`](@ref) uses for the per-size composition
inside the raindrop PSD integrals.

**Eight, and the choice is measured rather than picked.** Four already makes the zero-limiting
limit exact - with `τ_heat, τ_dend → 0` the two integrands are exactly `u³` and `u⁶`, and an
`n`-node Gauss-Laguerre rule is exact to polynomial degree `2n-1` - so the exactness sibling of
the analytic identity alone does not discriminate. What discriminates is the composed
integrand. Against an adaptive reference over the band where the composition changes the answer
(`ΔT ≤ 30` K), eight nodes hold the mass moment to 1.5 % and the number moment to 9 %, the best
of 4, 5, 6 and 8; the errors are not monotone in the node count, because Gauss-Laguerre error
alternates in sign, so "more is better" is not available as an argument and the count has to be
measured.

**Where no count in this range converges, and why that is acceptable.** Deeper than
`ΔT ≈ 35` K the number integrand's peak moves below the smallest node (`u₁ = 0.170` here) and
keeps going: past `ΔT*` it sits near `u ≈ 9e-4`. Fixed-node rules under-resolve it, by about a
factor three at the trigger cell, and adding nodes inside the affordable range does not fix it.
The error is one-signed - always an under-estimate of the number rate - and it is inert through
the substep, because `h/τ ≫ 1` at both the resolved and the under-resolved value and the
implicit update converts essentially the whole population either way. The MASS moment, which is
what the conservation bounds are stated on, converges everywhere and is exact to 1e-10 past
`ΔT*`, because its `D³` weight suppresses exactly the small-drop end the number moment cannot
resolve.

See `diag/rainfrz_tauheat_magnitudes.jl`, which regenerates the convergence table, and the
quadrature testsets in `test/rosenbrock_mode_tests.jl`.
"""
const RAIN_FREEZING_QUADRATURE_ORDER = 8

"""
    rain_freezing_quadrature()

The nodes and weights of the [`RAIN_FREEZING_QUADRATURE_ORDER`](@ref)-point Gauss-Laguerre
rule, which approximates `∫₀^∞ exp(-u) g(u) du ≈ Σᵢ wᵢ g(uᵢ)`.

The `exp(-u)` weight is the SB2006 exponential raindrop PSD written in `u = D/Dr_mean`, so this
is the PSD's own quadrature rather than a general-purpose rule applied to it, and the
zero-limiting limit reduces to the analytic moments exactly.

Returned as `Float64` tuples and converted at use, so a `Float32` kernel gets correctly rounded
nodes and a `ForwardDiff.Dual` one promotes without a separate table. Written as literals rather
than built from `FastGaussQuadrature` so that the values are compile-time constants inside a GPU
kernel with no global to capture; `test/rosenbrock_mode_tests.jl` asserts them against the
recurrence they satisfy, so the table cannot silently drift.

The order is fixed for production and passed as the `quad` keyword of
[`rain_freezing_rate`](@ref) only by the convergence tests, which build other orders host-side.
"""
@inline function rain_freezing_quadrature()
    nodes = (
        0.17027963230510093, 0.90370177679938, 2.251086629866131, 4.266700170287659,
        7.0459054023934655, 10.758516010180996, 15.740678641278004, 22.863131736889265,
    )
    weights = (
        0.36918858934163495, 0.4187867808143447, 0.17579498663717255, 0.033343492261215794,
        0.0027945362352256834, 9.076508773358139e-5, 8.48574671627257e-7, 1.0480011748715153e-9,
    )
    return (; nodes, weights)
end

"""
    rain_freezing_rate(opt, hom, vent, aps, tps, evap, pdf_r, q, ρ, N, T, qᵥ; quad)

Compose a Bigg nucleation freezing rate with the heat-dissipation limit on finishing the
freeze, and return the bulk rate the composition supports.

# Arguments
 - `opt`: the [`CMP.RainFreezing`](@ref) parameterization, for the Bigg volumetric rate.
 - `hom`: the [`CMP.Koop2000`](@ref) parameterization, for the homogeneous volumetric rate.
 - `vent`: ventilation coefficients, [`CMP.VentilationFactor`](@ref).
 - `aps`: air properties, [`CMP.AirProperties`](@ref).
 - `tps`: thermodynamics parameters.
 - `evap`: [`CMP.EvaporationSB2006`](@ref), for the mean-drop fall speed `α x̄^β √(ρ0/ρ)`.
 - `pdf_r`: the [`CMP.RainParticlePDF_SB2006`](@ref) raindrop size distribution.
 - `q`: rain specific content [kg(water) kg⁻¹(air)].
 - `ρ`: air density [kg(air) m⁻³(air)].
 - `N`: raindrop number concentration [m⁻³(air)].
 - `T`: air temperature [K].
 - `qᵥ`: ambient water vapour specific humidity [kg(water) kg⁻¹(air)], for the evaporative
   half of the wet-growth balance in [`drop_freezing_heat_timescale`](@ref).

# Keyword arguments
 - `quad`: the Gauss-Laguerre rule for the PSD integrals, [`rain_freezing_quadrature`](@ref).
   Production uses the default; the convergence tests pass other orders.

# Returns
 - A `NamedTuple` with the fields:
    + `∂ₜn_frz`: limited specific number freezing rate [kg⁻¹(air) s⁻¹].
    + `∂ₜq_frz`: limited specific mass freezing rate [kg(water) kg⁻¹(air) s⁻¹].
    + `τ_nuc`: the combined nucleation timescale of the mean-mass drop [s].
    + `τ_heat`: its heat-dissipation freezing timescale [s],
      [`drop_freezing_heat_timescale`](@ref).
    + `τ_dend`: its dendrite-propagation timescale [s],
      [`drop_freezing_dendrite_timescale`](@ref).
    + `J_bigg`, `J_koop`: the two volumetric nucleation coefficients [m⁻³(water) s⁻¹].

**Nucleation pathways are PARALLEL; conversion stages are SEQUENTIAL.** Two independent ways to
nucleate the same drop give additive RATES, so their volumetric coefficients add; nucleating and
then finishing the freeze are one after the other, so those TIMESCALES add:

```
1/τ_nuc = (J_bigg(T) + J_koop(T)) ⋅ V_drop = (J_bigg + J_koop) ⋅ x̄ / ρw
τ_eff   = τ_nuc + τ_heat + τ_dend
```

both stated per drop, so the sum is formed BEFORE the PSD integral: both pathways go through one
PSD treatment and nothing is ever divided by `J_bigg`.

Adding the homogeneous pathway is what stops the scheme relying on Bigg's out-of-range
extrapolation being accidentally right. Bigg is a fit over roughly 0 to 40 K of supercooling;
past that it is the only thing the scheme had, and it is the term that reached `2.9e13` kg/kg/s.
Measured handoff, drop-size INDEPENDENT because the same `V_drop` multiplies both:

  - `J_koop` is identically zero warmer than `ΔT ≈ 30.9` K (`Δa_w < 0.26`);
  - `J_bigg = J_koop` at `ΔT ≈ 35.2` K (`T ≈ 238` K, `Δa_w ≈ 0.29`), and Koop dominates below;
  - by `ΔT = 40` K Koop leads Bigg by five orders, which is the point: the physically supported
    pathway takes over before the extrapolated one runs away.

At the trigger cell (`ΔT = 82` K) the clamped `J_koop = 2.9e24` is one order BELOW the
extrapolated `J_bigg = 2.8e25`, so the sum is only 10 % above Bigg alone and the vertex behaviour
is set by `τ_dend`, unchanged.

Nucleate, then shed the latent heat, then let the front cross the drop. Each stage dominates in
its own regime and the sum needs no branch:

  - weak supercooling: `τ_nuc` dominates, the Bigg exponential being small and Koop zero;
  - the ordinary mixed-phase range from roughly 20 K: `τ_heat` overtakes it for large drops;
  - past `ΔT* ≈ 53` K: `τ_heat` is exactly zero, and `τ_dend` alone bounds the conversion.

The third stage is what reaches the state this whole thread exists for. At the box's trigger
cell (`T = 191` K, `ΔT = 82` K, `x̄ ≈ 5e-7` kg) `τ_nuc ≈ 7e-17` s and `τ_heat = 0`, so the
two-stage form passed the raw Bigg rate through untouched; with `τ_dend ≈ 2e-3` s the rate
falls by about thirteen orders of magnitude, and per size the mass moment falls by another
factor ~30 on top of that, because `τ_dend ∝ D` and the mass integral lives on the large drops
the mean-drop factor under-limited.

`τ_dend` is included UNCONDITIONALLY, not gated on `τ_heat` having vanished. Gating would make
`τ_eff` jump by `τ_dend` as `ΔT` crosses `ΔT*`, and since `τ_nuc(ΔT*) ≈ 1e-8` s against
`τ_dend ≈ 2e-3` s that is a five-order DISCONTINUITY in a tendency at one temperature. This
scheme has just finished removing a `C⁰` kink from the shape solver for making the bracketing
solve worse; introducing a jump discontinuity to buy a bit-identity property would be moving
backwards. Dendrite propagation is also a real stage at every supercooling, merely subdominant
where the other two are large.

**The composition is done PER SIZE, inside the PSD integrals.** The limiting ratio scales as
`τ_heat/τ_nuc ∝ D⁵/F_v(D)`, so a factor two in diameter is a factor ~30 in the ratio and the
population STRADDLES the crossover at band temperatures. A single mean-drop factor multiplying
the whole bulk rate therefore over-limits the sub-mean drops, which nucleate slowly but freeze
fast, and under-limits the super-mean ones. Each drop instead carries its own conversion rate,

```
r(D)     = 1 / (τ_nuc(D) + τ_heat(D) + τ_dend(D)),   1/τ_nuc(D) = (J_bigg + J_koop) (π/6) D³
∂ₜn_frz  = ∫ r(D) n(D) dD
∂ₜq_frz  = ∫ x(D) r(D) n(D) dD,                      x(D) = (π/6) ρw D³
```

evaluated by [`rain_freezing_quadrature`](@ref) against the SB2006 exponential PSD
`n(D) dD = n exp(-u) du` with `u = D/Dr_mean`, which is the Gauss-Laguerre weight exactly.

**This resolves the two known caveats of the mean-drop form at once, with no branch.** The
first is the size-selectivity above. The second is the EVENT MASS: `∂ₜq_frz/∂ₜn_frz` is the
mean mass converted per freezing event, and the mean-drop form preserves it at `20 x̄` at every
temperature, because both moments took the same divisor. `20 x̄` is the right answer for RARE,
volume-selective freezing - only the biggest drops go, and they carry twenty times the mean
mass - and the wrong one in the every-drop-freezing limit, where the converted mass must relax
toward `x̄`. Per size it slides down on its own: `r(D)` saturates at `1/(τ_heat + τ_dend)`
where `J` is large, the `D³` volume selectivity drops out of it, and the measured event mass
falls monotonically from `18.6 x̄` at `ΔT = 5` K through `x̄` near `ΔT = 25` K. Deeper than
that it keeps falling, to `~0.05 x̄`, and that undershoot is real rather than an artifact: with
`τ_heat = 0` past `ΔT*` the surviving `τ_dend = D/v_dend` makes the smallest drops the fastest
converters, so the number integral is dominated by drops far below `x̄`.

**The quadrature's accuracy is measured, and it is not uniform.** The rule is exact in the
zero-limiting limit: with `τ_heat, τ_dend → 0` the integrands are exactly `u³` and `u⁶`, which
`RAIN_FREEZING_QUADRATURE_ORDER = 8` integrates to round-off (Gauss-Laguerre is exact to degree
`2n-1`, so four nodes would already do it), recovering `J (π/6) M_D³` and `J ρw (π/6)² M_D⁶`.
Against an adaptive reference the mass moment holds to 1.5 % or better at every supercooling
and to 1e-10 past `ΔT*`. The NUMBER moment holds to a few percent only while the integrand's
peak stays above the smallest node `u₁ = 0.170`, which is true out to `ΔT ≈ 35` K and false
beyond it: past `ΔT*` the peak sits at `u ≈ 9e-4` and no fixed-node rule in this range resolves
it. The resulting error there is a systematic UNDER-estimate of the number rate, roughly a
factor three at the trigger cell, and it is inert through the substep - the implicit update
converts essentially the whole population at either value, since `h/τ ≫ 1` both ways.

Evaluating `τ_nuc` at the mean-mass drop was exact for the number moment of the UNLIMITED rate,
and that identity is what `_liquid_freezing_rate_from_J` still carries: for the SB2006
exponential rain PSD `x̄ = π ρw Dr_mean³`, so

```
∂ₜn_frz|nuc = J (π/6) M_D³ = J (π/6) Γ(4) n Dr_mean³ = J π n Dr_mean³ = n / τ_nuc
```

exactly. The returned `τ_nuc`, `τ_heat` and `τ_dend` are the mean-mass drop's and are
DIAGNOSTICS only: the composition no longer uses them, and the reduction is no longer
`1/(1 + (τ_heat + τ_dend)/τ_nuc)` for either moment.

**Bit-identity is NOT claimed, and the earlier two-stage form's claim is withdrawn.** With
`τ_dend` always positive on a populated state, `r(D)` is strictly below the nucleation rate
everywhere, so the returned rate is never bit-identical to the unlimited one except where that
rate is already exactly zero (an empty rain state, or above `T_freeze - 4`). At weak
supercooling the NUMBER moment is within about 2 % of the unlimited one, but the MASS moment is
not: it is weighted by `D⁶`, so it is carried by drops several times the PSD scale, and those
are heat-limited far warmer than the mean drop is. Measured, at `ΔT = 5` K the mass rate is
about 9 % below the mean-drop form's. That is a real effect the mean-drop factor was hiding, not
a regression, and the tests bracket it rather than asserting an equality that does not hold.

**Where the boundedness lives, and it is not in the rate.** Past `ΔT*` the returned rate is
still far larger than the donor can supply over a step - at the trigger cell it is about
1 kg/kg/s against `q_rai ≈ 1e-4`, so `rate ⋅ dt` exceeds `q_rai` by four orders. That is
expected and correct for a relaxation whose timescale is short compared with `dt`, and it is
why `rate ⋅ dt ≤ q_rai` is the WRONG thing to require. The substep donor-linearizes this
transfer onto the `q_rai` diagonal at `-1/τ_eff`, and the linearized-implicit update of a decay
`q/τ_eff` is

```
Δq = -h (q/τ_eff) / (1 + h/τ_eff),
```

bounded by `q` for every `h` and tending to complete conversion as `h/τ_eff → ∞`. At the
trigger cell `h/τ_eff ≈ 1e3`, so the step converts about 99.9 % of the rain and no more. The
conversion is therefore donor-bounded THROUGH the implicit update, and timestep-independent:
no `dt` appears anywhere in the rate. The tests assert boundedness by driving the substep entry,
not by bounding the rate.

An `∂ₜ(n, q)|nuc` that has already overflowed no longer produces a `NaN`, and that is a
consequence of the per-size form rather than a guard added to it. `J_bigg` grows like
`exp(0.65 ΔT)` and reaches `floatmax(Float32)` near `ΔT = 132` K; the mean-drop form divided
one overflowed quantity by another and returned `Inf/Inf`, while `r(D) = 1/(τ_nuc + τ_f)` sends
`τ_nuc → 0` and returns the finite `1/τ_f`, which is the correct limit of an
infinitely-fast nucleation. Bounding `J_bigg` itself, or the INP budget, is still a separate
matter.
"""
@inline function rain_freezing_rate(opt, hom, vent, aps, tps, evap, pdf_r, q, ρ, N, T, qᵥ;
    quad = rain_freezing_quadrature(),
)
    FT = UT.promote_typeof(q, ρ, N, T)
    T_freeze = TDI.T_freeze(tps)
    (; ρw) = pdf_r
    (; α, β, ρ0) = evap
    (; Dr_mean, xr_mean) = CM2.pdf_rain_parameters(pdf_r, q, ρ, N)
    n = N / ρ                # specific number concentration [kg⁻¹(air)]

    # Parallel nucleation pathways: the coefficients add, and the sum is formed BEFORE the PSD
    # integral so both share one treatment and `J_bigg` is never a denominator.
    J_bigg = opt(T, T_freeze)
    J_koop = homogeneous_freezing_rate_coefficient(hom, tps, T)
    J = J_bigg + J_koop

    # Per-size composition, integrated against the exponential PSD by Gauss-Laguerre in
    # `u = D/Dr_mean`, whose `exp(-u)` weight IS the PSD. Shared with the cloud category through
    # `_composed_liquid_freezing_moments`; only the node-to-diameter map and the fall-speed law
    # differ, which is the categorization-blindness principle in code.
    v_at = (D, x_drop) -> α * x_drop^β * sqrt(ρ0 / ρ)
    (; ∂ₜn_frz, ∂ₜq_frz) = _composed_liquid_freezing_moments(
        vent, aps, tps, T, ρ, qᵥ, ρw, J, n, u -> Dr_mean * u, v_at, quad,
    )

    # The same presence and temperature gate `_liquid_freezing_rate_from_J` applies, so the two
    # entry points agree on where freezing happens at all.
    cond = _liquid_freezing_is_active(FT, q, n, T, T_freeze)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(FT))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(FT))

    # Mean-mass-drop timescales, reported as DIAGNOSTICS. The composition above uses the
    # per-node values, not these; they are what `diag/rainfrz_tauheat_magnitudes.jl` tabulates
    # and what the crossover and ΔT* claims are stated in terms of.
    v_bar = α * xr_mean^β * sqrt(ρ0 / ρ)
    τ_heat = drop_freezing_heat_timescale(vent, aps, tps, T, ρ, qᵥ, xr_mean, ρw, v_bar)
    τ_dend = drop_freezing_dendrite_timescale(xr_mean, ρw)
    inv_τ_nuc = J * xr_mean / ρw  # = (J_bigg + J_koop) V_drop
    return (;
        ∂ₜn_frz,
        ∂ₜq_frz,
        # `Inf` where no pathway is open or no drop is present, which is the physical answer for a
        # nucleation waiting time; guarded so it is `Inf` with ZERO partials rather than NaN ones.
        τ_nuc = UT.guarded_quotient(one(inv_τ_nuc), inv_τ_nuc, oftype(inv_τ_nuc, Inf)),
        τ_heat,
        τ_dend,
        J_bigg,
        J_koop,
    )
end

"""
    immersion_limit_rate(opt::CMP.Frostenberg2023, T, ρ; τ, inpc_log_shift, n_active)

Compute the F23-INPC-imposed upper limit on the cloud-droplet immersion freezing number rate.

!!! note
    This bound is not part of the default immersion path, which is Bigg alone
    (see [`cloud_freezing_rate`](@ref)). It is retained as the ingredient of a
    selectable INP-budgeted immersion option.

The Frostenberg 2023 climatology specifies a target INP concentration `INPC(T)`
in air [m⁻³(air)]. Treating that target as a budget that should be activated
on a relaxation timescale `τ`, the maximum number of crystals nucleated per
kg of air per second is

```
∂ₜn_lim = max(0, INPC(T)/ρ - n_active) / τ.
```

`n_active` is the depletion proxy; see [`n_active`](@ref).

# Arguments
 - `opt`: The [`CMP.Frostenberg2023`](@ref) parameters.
 - `T`: Air temperature [K].
 - `ρ`: Air density [kg(air) m⁻³(air)].

# Keyword arguments
 - `τ`: Relaxation timescale `[s]` (default `300`).
 - `inpc_log_shift`: Additive shift to `log(INPC)` (default `0`).
 - `n_active`: Depletion proxy [kg⁻¹(air)].

# Returns
 - A `NamedTuple` `(; ∂ₜn_frz)` — the specific number freezing-rate cap
   [kg⁻¹(air) s⁻¹]. Zero when `T ≥ T_freeze`.
"""
function immersion_limit_rate(
    opt::CMP.Frostenberg2023, T, ρ;
    τ = oftype(T, 300), inpc_log_shift = zero(T),
    n_active = zero(T),
)
    # The early return must carry the element type the main path produces, or a
    # differentiated (Dual `n_active`) kernel gets a type-unstable NamedTuple and
    # faults at runtime. See UT.promote_typeof.
    FT = UT.promote_typeof(T, ρ, inpc_log_shift, n_active, τ)
    T ≥ opt.T_freeze && return (; ∂ₜn_frz = zero(FT))
    log_inpc = INP_concentration_mean(opt, T) + inpc_log_shift
    INPC_per_kg = exp(log_inpc) / ρ                  # [kg⁻¹(air)]
    ∂ₜn_frz = max(zero(FT), INPC_per_kg - n_active) / τ # [kg⁻¹(air) s⁻¹]
    return (; ∂ₜn_frz)
end

"""
    CLOUD_FREEZING_QUADRATURE_ORDER
    cloud_freezing_quadrature()

The GENERALIZED Gauss-Laguerre rule, weight `t^α exp(-t)` with `α = 1`, that
[`cloud_freezing_rate`](@ref) uses for the per-size composition over the SB2006 cloud-droplet
PSD, and the order it uses. Approximates `∫₀^∞ t^α exp(-t) g(t) dt ≈ Σᵢ wᵢ g(tᵢ)`.

**`α = 1` is derived from the shipped shape parameters, not chosen.** The cloud PSD is
`n(D) = N₀c D^νcD exp(-λc D^μcD)` with `νcD = 3νc + 2` and `μcD = 3μc`; substituting
`t = λc D^μcD` gives `n(D) dD ∝ t^α exp(-t) dt` with

```
α = (νcD + 1)/μcD - 1 = (νc + 1)/μc - 1,
```

which is `1` at the shipped `νc = μc = 1`. The rule is therefore tied to those two parameters,
so `test/rosenbrock_mode_tests.jl` asserts the derived `α` against the rule's, and a ClimaParams
change to either shape parameter fails a test instead of silently mis-integrating the PSD. The
weight integrates to `Γ(α+1) = 1`, matching the rain rule's normalization, so the moments carry
no normalization factor.

**Eight nodes, for the same reason the rain rule uses eight.** The zero-limiting limit needs only
two (the integrands reduce to `t` and `t²`, and an `n`-node rule is exact to degree `2n-1`), so it
does not choose the count. Measured against an adaptive reference the eight-node rule holds the
mass moment to 0.1 % and the number moment to 3 % across the whole range, with the same character
as the rain rule: the mass moment converges everywhere, the number moment is the harder one.

Written as literals rather than built from `FastGaussQuadrature` so the values are compile-time
constants inside a GPU kernel with no global to capture; the tests assert them against the
defining `Σ wᵢ tᵢ^k = Γ(α+1+k)/Γ(α+1)` exactness, so the table cannot silently drift.
"""
const CLOUD_FREEZING_QUADRATURE_ORDER = 8

@inline function cloud_freezing_quadrature()
    nodes = (
        0.4093835732031852, 1.3849631848031398, 2.956254556168862, 5.181943101040071,
        8.161709688145818, 12.070055126837154, 17.24973552614899, 24.58595524365278,
    )
    weights = (
        0.18763254140572333, 0.4389853607311426, 0.2899960707813135, 0.07514138461669735,
        0.007932646648707355, 0.0003086421368133023, 3.348958209797081e-6, 4.721392823193185e-9,
    )
    return (; nodes, weights)
end

"""
    cloud_freezing_rate(opt, hom, vent, aps, tps, pdf_c, q, ρ, N, T, qᵥ; quad)

Freeze cloud droplets with the SAME composition rain gets, evaluated at cloud sizes.

!!! note
    This is the Frostenberg 2023 immersion-mode INP spectrum standing in the
    deposition slot, retained as a non-default option. It is not the default
    deposition parameterization; see the
    [`CMP.AbstractINPTargetSpectrum`](@ref) method above.

# Arguments
 - `opt`: the [`CMP.RainFreezing`](@ref) parameterization, for the Bigg volumetric rate. The
   `Rain` in the name is historical; the kinetics apply to any liquid drop.
 - `hom`: the [`CMP.Koop2000`](@ref) parameterization, for the homogeneous volumetric rate.
 - `vent`, `aps`, `tps`: ventilation coefficients, air properties, thermodynamics parameters.
 - `pdf_c`: the [`CMP.CloudParticlePDF_SB2006`](@ref) cloud-droplet size distribution.
 - `q`, `ρ`, `N`, `T`, `qᵥ`: cloud liquid specific content, air density, droplet number
   concentration [m⁻³], air temperature and vapour specific humidity.

# Keyword arguments
 - `quad`: the rule, [`cloud_freezing_quadrature`](@ref) by default.

# Returns
 - `(; ∂ₜn_frz, ∂ₜq_frz, J_het, J_koop, τ_nuc, τ_heat, τ_dend)`. The three timescales are the
   mean-mass DROPLET's and are diagnostics; the composition uses the per-node values.

**Liquid is liquid.** Cloud drops and rain drops are not different substances, they are different
sizes, so this is [`rain_freezing_rate`](@ref)'s composition with the cloud PSD substituted, not a
cloud-specific parameterization. Tendencies differ because sizes differ, never because of the
chosen categorization - which matters because a single-category warm scheme should get the same
physics, and because the alternative asymmetry (cloud freezing bounded and two-stage-free, rain
freezing raw and kinetically limited) is indefensible in either direction.

**The immersion coefficient is Bigg alone.** Cloud-droplet immersion freezing runs the
Barklie-Gokhale form of Bigg (1953) with no ice-nucleating-particle budget above it, which is what
the reference P3 implementation does: its `qcheti` is `exp(a_imm (T₀ - T))` with no budget, no
bound and no aerosol dependence. An INP-spectrum bound on the heterogeneous coefficient is
available as a non-default option through [`immersion_limit_rate`](@ref), which converts an INP
budget from a number RATE to the coefficient that would produce it through this population's own
PSD integral; the bound belongs on the heterogeneous COEFFICIENT and never on the summed rate, so
that homogeneous freezing is never throttled by a budget it does not draw on.

**Uncapped Bigg is bounded by the composition, not left unbounded.** The Bigg coefficient grows
without limit at deep supercooling, and the composition is what holds the realized rate finite:
`τ_eff = τ_nuc + τ_heat + τ_dend` with `τ_nuc = ρw / (J x̄)`, so a larger `J` shrinks `τ_nuc` while
`τ_dend` is a floor INDEPENDENT of `J`. The conversion is therefore bounded by `1/τ_dend` however
large the coefficient becomes.

**The kinetic series is not omitted for small drops, and it is not always inert.** The expectation
was that at cloud sizes it evaluates to `1` automatically, so uniformity costs nothing. Measured,
that is true wherever the HETEROGENEOUS pathway drives freezing - the
retained fraction is above `0.995` at every cloud size warmer than the Koop window, and above
`0.9999` for a 20 μm droplet - and it is FALSE once the homogeneous pathway takes over. At 233 K
`J_koop` collapses `τ_nuc` for a 20 μm droplet to about `6e-7` s against `τ_heat ≈ 2.4e-3` s, so
the series retains about `2e-4` of the nucleation rate. That is correct physics rather than a
defect: when every droplet nucleates within a microsecond, the conversion is limited by how fast
the droplets can shed the latent heat of fusion, and even a 20 μm droplet needs milliseconds. It
is also inert THROUGH the substep, where `h/τ ≈ 800` at a 2 s step still converts the population.
The consequence to carry into any evaluation of this change: the structural payoff below 236 K is
large but it is roughly four orders smaller than the coefficient ratio alone would suggest.

Fall speed comes from the Stokes-regime law rather than the SB2006 rain power law, through the
scheme's own [`CO.particle_terminal_velocity`](@ref) with a locally constructed
[`CMP.StokesRegimeVelType`](@ref). This is a size-regime statement, not a category one, and it is
the honest reading of the principle: the rain power law extrapolated to a 20 μm droplet gives
0.17 m/s against a physical 0.015 m/s, which is the same class of out-of-range extrapolation this
thread exists to remove. It enters only through the ventilation Reynolds number, so the price of
getting it wrong would have been about 11 % in `τ_heat`. Constructing the velocity type locally
avoids adding a field to a shipped parameter struct on a prototype branch; threading it from
`CMP.TerminalVelocityParams` is the follow-up if the form is adopted.
"""
@inline function cloud_freezing_rate(
    opt, hom, vent, aps, tps, pdf_c, q, ρ, N, T, qᵥ;
    quad = cloud_freezing_quadrature(),
)
    FT = UT.promote_typeof(q, ρ, N, T)
    T_freeze = TDI.T_freeze(tps)
    (; ρw) = pdf_c
    (; λc, μcD) = CM2.pdf_cloud_parameters(pdf_c, q, ρ, N)
    n = N / ρ                # specific number concentration [kg⁻¹(air)]

    # Bigg (1953) heterogeneous coefficient, and Koop (2000) homogeneous beside it. The two
    # pathways act in parallel on the same drop, so their volumetric coefficients add.
    J_het = opt(T, T_freeze)
    J_koop = homogeneous_freezing_rate_coefficient(hom, tps, T)
    J = J_het + J_koop

    # Stokes-regime fall speed, the size regime cloud droplets are actually in. Built from the
    # scheme's own function so there is no second copy of the expression to drift.
    # Built at the parameters' own float type, not at `FT`: under the exact-AD Jacobian `FT` is a
    # `Dual`, and a velocity law whose constants carried dead partials would cost work in every
    # kernel for nothing. `D` promotes at the multiply.
    stokes = CMP.StokesRegimeVelType(; ρw, ν_air = aps.ν_air, grav = TDI.grav(tps))
    v_stokes = CO.particle_terminal_velocity(stokes, ρ)
    # `t = λc D^μcD`, so `D = (t/λc)^(1/μcD)`. Written as a general power rather than a `cbrt` so
    # that only `α` is assumed about the shape, and `α` is asserted in the tests.
    D_at = t -> (t / λc)^(1 / μcD)
    (; ∂ₜn_frz, ∂ₜq_frz) = _composed_liquid_freezing_moments(
        vent, aps, tps, T, ρ, qᵥ, ρw, J, n, D_at, (D, x_drop) -> v_stokes(D), quad,
    )

    # The same presence and temperature gate every other liquid freezing entry applies.
    cond = _liquid_freezing_is_active(FT, q, n, T, T_freeze)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(FT))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(FT))

    # Mean-mass-droplet timescales, DIAGNOSTICS only, for the magnitudes probe and the tests.
    x̄ = UT.guarded_quotient(q, n)
    D̄ = cbrt(6 * x̄ / (FT(π) * ρw))
    τ_heat = drop_freezing_heat_timescale(vent, aps, tps, T, ρ, qᵥ, x̄, ρw, v_stokes(D̄))
    τ_dend = drop_freezing_dendrite_timescale(x̄, ρw)
    return (;
        ∂ₜn_frz,
        ∂ₜq_frz,
        J_het,
        J_koop,
        τ_nuc = UT.guarded_quotient(ρw, J * x̄, oftype(J * x̄, Inf)),
        τ_heat,
        τ_dend,
    )
end

# ---------------------------------------------------------------------------
# Deposition nucleation: the INP target-spectrum interface
# ---------------------------------------------------------------------------

"""
    is_active(inp, T, S_i)

Whether the deposition nucleation slot is open at this thermodynamic state, for the INP
target spectrum `inp` (a [`CMP.AbstractINPTargetSpectrum`](@ref)).

# Arguments
 - `inp`: the INP target spectrum.
 - `T`: air temperature [K].
 - `S_i`: ice supersaturation ratio `q_vap/q_sat_ice - 1` [-].

# Returns
 - `Bool`, whether nucleation proceeds.

Each spectrum states its own activation window, because the window belongs to the
measurement the spectrum is fitted to and not to the relaxation form the slot shares.
"""
function is_active end

"""
    delivery_rate(inp, mp, tps, T, S_i)

The inverse of the activation time of the INP target spectrum `inp` [s⁻¹]: how fast the
deficit between the target concentration and the already-activated concentration is
converted into ice crystals.

# Arguments
 - `inp`: the INP target spectrum, an [`CMP.AbstractINPTargetSpectrum`](@ref).
 - `mp`: the [`CMP.Microphysics2MParams`](@ref) parameter container.
 - `tps`: thermodynamics parameters.
 - `T`: air temperature [K].
 - `S_i`: ice supersaturation ratio [-].

# Returns
 - The inverse activation time [s⁻¹].

Returned as an inverse time rather than a time because that is how the slot consumes it:
the deficit is multiplied by it. A spectrum whose activation time is a fixed constant
returns the reciprocal of that constant.
"""
# One method, on the abstract spectrum: the activation WINDOW belongs to the measurement a
# spectrum is fitted to, the relaxation FORM belongs to the slot, and this is the form.
function delivery_rate end

"""
    is_active(inp::CMP.ExponentialSupercoolingINP, T, S_i)

The two activation conditions of the exponential-in-supercooling deposition spectrum:
cold enough, and supersaturated with respect to ice by at least the threshold. Both come
from the reference P3 implementation, quoted in
[`CMP.ExponentialSupercoolingINP`](@ref).
"""
@inline is_active(inp::CMP.ExponentialSupercoolingINP, T, S_i) =
    (T < inp.T_thr) & (S_i >= inp.S_thr)

"""
    delivery_rate(inp::CMP.AbstractINPTargetSpectrum, mp, tps, T, S_i)

The seed delivery rate: the reciprocal of the time diffusional growth needs to build a
crystal of the nascent radius `r_nuc` out of vapor,

```
1/τ_act = 2 G_ice(T) max(S_i, 0) / (ρ_i r_nuc²).
```

It follows from the diffusional growth law `r ∂ₜr = G_ice S_i / ρ_i` integrated from zero
to `r_nuc` at fixed supersaturation, so the deposition slot activates ice no faster than
vapor can make an ice crystal of the size the slot claims to make. `G_ice` is the standard
combined conductivity and diffusivity coefficient the package already evaluates
([`CO.G_func_ice`](@ref)); the same coefficient builds the ice term of the activation
supersaturation balance in `AerosolActivation.max_supersaturation`, so the two places that
need a diffusional relaxation rate read it from one function.

The delivery time is state dependent and falls steeply with supersaturation. At the
activation threshold `S_i = S_thr` in cold cirrus conditions it is of order 10 s, and it is
shorter everywhere the slot fires harder, so it is the reciprocal, not the time, that is
smooth: only `max(S_i, 0)` appears, no division by `S_i`, and there is no singular branch
to mask. Every target spectrum takes this form. A spectrum whose window carries a
supersaturation threshold makes the `max` inert wherever its slot is active; one whose
window admits `S_i = 0` relies on the `max`, and the rate falls continuously to zero
there rather than being cut off by the gate.

The seed geometry enters through `r_nuc` and `ρ_i` from [`CMP.ice_seed`](@ref), which is the
single source for the nascent crystal shared with the starter mass, the ice
number-adjustment bound, the orphan-ice drain, the shape-solve bracket, and the melt limit.
The equivalent solid-sphere form `8π G_ice S_i r_nuc / (3 m_nuc)` is the same quantity
rewritten through `m_nuc = (π/6) ρ_i D_nuc³`; the form above is the fundamental one and is
what is evaluated.
"""
@inline function delivery_rate(inp::CMP.AbstractINPTargetSpectrum, mp, tps, T, S_i)
    (; r_nuc, ρ_i) = CMP.ice_seed(mp.ice.scheme)
    return 2 * CO.G_func_ice(mp.warm_rain.air_properties, tps, T) * max(S_i, 0) /
           (ρ_i * r_nuc^2)
end

"""
    (inp::CMP.Frostenberg2023)(T)

The target ice nucleating particle concentration [m⁻³] of the Frostenberg 2023
climatology, the mean of its lognormal INPC(T) distribution
([`INP_concentration_mean`](@ref)) exponentiated back out of log space.

The future stochastic reading of the spectrum adds an additive shift to this
log-mean before exponentiating (`inpc_log_shift`, deferred to that commit); the
value here is that shift at its present default of zero, so no extra term
appears.
"""
@inline (inp::CMP.Frostenberg2023)(T) = exp(INP_concentration_mean(inp, T))

"""
    is_active(inp::CMP.Frostenberg2023, T, S_i)

The Frostenberg 2023 activation window: below freezing, and not subsaturated with
respect to ice.

The two literals this carried, colder than 15 K below freezing and supersaturated by
more than 5 percent, are GONE. They were the default closure's values and not this
spectrum's, which has no threshold at either quantity, and the docstring that called
them the closure's own defaults was wrong. Their removal is what makes the rate
continuous in both `T` and `S_i` across the two former thresholds.

The subsaturation floor states the window rather than doing the limiting: the shared
[`delivery_rate`](@ref) carries `max(S_i, 0)`, so the rate already falls continuously to
zero as saturation is approached from above, which is the role the retired vapor-excess
cap played and which a gate cannot play continuously. The boundary
point `S_i == 0` is inside the window by the letter of the ruling; nothing physical
turns on it, and a strict `S_i > 0` would differ only there.
"""
@inline is_active(inp::CMP.Frostenberg2023, T, S_i) =
    (T < inp.T_freeze) & (S_i >= zero(S_i))

"""
    deposition_rate(inp::CMP.AbstractINPTargetSpectrum, mp, tps, micro, thermo)
    deposition_rate(::Nothing, mp, tps, micro, thermo)

Nucleate pristine ice from the vapor phase onto ice nucleating particles.

# Arguments
 - `inp`: the INP target spectrum, an [`CMP.AbstractINPTargetSpectrum`](@ref); `nothing`
   disables the process.
 - `mp`: the [`CMP.Microphysics2MParams`](@ref) parameter container. Reads
   `mp.ice.scheme` (the nascent crystal), `mp.ice.inp_depletion_model` (the depletion
   proxy) and `mp.warm_rain.air_properties` (the growth coefficient).
 - `tps`: thermodynamics parameters.
 - `micro`: microphysics state `NamedTuple`; reads the specific contents `q_tot`, `q_lcl`,
   `q_rai` and `q_ice` [kg kg⁻¹] and the specific ice number `n_ice` [kg⁻¹].
 - `thermo`: thermodynamic state `NamedTuple`; reads `ρ` [kg m⁻³] and `T` [K].

# Returns
 - A `NamedTuple` `(; ∂ₜn_frz, ∂ₜq_frz)`, the specific number rate [kg⁻¹ s⁻¹] and the
   specific mass rate [kg kg⁻¹ s⁻¹]. Both are zero outside the activation window.

The slot is a target-limited relaxation. The spectrum `inp` supplies a target number
concentration `N_t(T)` [m⁻³], the depletion model supplies how much of it is already
activated ([`n_active`](@ref)), and the deficit between them is delivered at the spectrum's
own rate ([`delivery_rate`](@ref)):

```
∂ₜn_frz = max(0, N_t(T)/ρ - n_active) · delivery_rate(inp, ...)   [kg⁻¹ s⁻¹]
∂ₜq_frz = m_nuc · ∂ₜn_frz                                          [kg kg⁻¹ s⁻¹]
```

Every crystal is created at the nascent mass `m_nuc` of [`CMP.ice_seed`](@ref), so the pair
is internally consistent at one starter mass per crystal by construction, and there is no
bound on the mass moment that could break that consistency. The vapor the injection draws
is small: at the nascent size the demand is two orders below the vapor excess that the
activation threshold already guarantees is present, and the implicit substep update the
slot is evaluated inside absorbs what remains.

The liquid grouping is owned here rather than at the caller: `q_vap` is formed from
`q_tot` less `q_lcl + q_rai` less `q_ice`, so one site decides which condensate enters the
saturation ratio.

The activation window is the spectrum's own and is applied to both moments together, so a
closed slot returns an exactly zero pair rather than a small inconsistent one.
"""
@inline function deposition_rate(inp::CMP.AbstractINPTargetSpectrum, mp, tps, micro, thermo)
    (; ρ, T) = thermo
    q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    q_vap = TDI.q_vap(micro.q_tot, micro.q_lcl + micro.q_rai, micro.q_ice)
    S_i = q_vap / q_sat_ice - 1
    n_act = n_active(mp.ice.inp_depletion_model, micro.n_ice)
    (; m_nuc) = CMP.ice_seed(mp.ice.scheme)
    ∂ₜn_frz = max(0, inp(T) / ρ - n_act) * delivery_rate(inp, mp, tps, T, S_i)
    ∂ₜn_frz = ifelse(is_active(inp, T, S_i), ∂ₜn_frz, zero(∂ₜn_frz))
    return (; ∂ₜn_frz, ∂ₜq_frz = m_nuc * ∂ₜn_frz)
end

@inline deposition_rate(::Nothing, mp, tps, micro, thermo) =
    (; ∂ₜn_frz = zero(thermo.ρ), ∂ₜq_frz = zero(thermo.ρ))

# ---------------------------------------------------------------------------
# INP-activation memory dispatch
# ---------------------------------------------------------------------------

"""
    n_active(model::CMP.NIceProxyDepletion, n_ice)

Return the depletion proxy `n_active` to subtract from the INP target
concentration in [`deposition_rate`](@ref) and any analogous budgeted rate.
For `NIceProxyDepletion` (the only model currently provided) this is the
in-cell ice number `n_ice`. (A prognostic activation-memory model returning a
host-supplied tracer is deferred to a follow-up PR.)
"""
@inline n_active(::CMP.NIceProxyDepletion, n_ice) = n_ice

end # end module
