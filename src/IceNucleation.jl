# This file contains the `HetIceNucleation` and `HomIceNucleation` modules.

"""
Parameterization for heterogenous cloud ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP
import CloudMicrophysics.Common as CO
import CloudMicrophysics.ThermodynamicsInterface as TDI
import CloudMicrophysics.Microphysics2M as CM2
import CloudMicrophysics.DistributionTools as DT
import CloudMicrophysics.Utilities as UT

export dust_activated_number_fraction
export MohlerDepositionRate
export deposition_J
export ABIFM_J
export P3_deposition_N_i
export P3_het_N_i
export INP_concentration_frequency
export INP_concentration_mean
export liquid_freezing_rate
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
    FT = float(UT.promote_typeof(q, ρ, N, T))
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    (; ρw) = pdf  # [kg(water) m⁻³(water)]
    n = N / ρ     # specific number concentration [kg⁻¹(air)]

    # Solve for the pdf parameters
    (; Dr_mean) = CM2.pdf_rain_parameters(pdf, q, ρ, N)

    # Bigg (1953) volumetric freezing rate:
    J_bigg = opt(T, T_freeze)  # [m⁻³(water) s⁻¹]

    # Diameter PSD moments via exponential_Mⁿ:
    #   M_D^k = ∫ D^k n(D) dD
    # The freezing probability per unit time for a single drop of diameter D
    # is J_drop(D) = J_bigg * V(D) = J_bigg * (π/6) * D³  [s⁻¹].

    # Number freezing rate per kg air:  ∂n/∂t = ∫ J_drop(D) n(D) dD
    #                                         = J_bigg * (π/6) * M_D³
    M_D³ = DT.exponential_Mⁿ(Dr_mean, n, 3)  # [m³ · kg⁻¹(air)]

    # Mass freezing rate per kg air:    ∂q/∂t = ∫ x(D) * J_drop(D) n(D) dD
    #                                         = J_bigg * ρw * (π/6)² * M_D⁶
    M_D⁶ = DT.exponential_Mⁿ(Dr_mean, n, 6)  # [m⁶ · kg⁻¹(air)]

    V_1 = FT(π) / 6
    ∂ₜn_frz = J_bigg * V_1 * M_D³         # [kg⁻¹(air) s⁻¹]   — specific
    ∂ₜq_frz = J_bigg * ρw * V_1^2 * M_D⁶  # [kg(water) kg⁻¹(air) s⁻¹]  — specific

    # Return the computed rate only if N and L are (essentially) non-zero, and T is colder than -4°C.
    # Otherwise, return zero.
    ϵ = UT.ϵ_numerics(FT)
    cond = (n > ϵ) & (q > ϵ) & (T < T_freeze - 4)
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
    ϵ = UT.ϵ_numerics(FT)
    cond = (n > ϵ) & (q > ϵ) & (T < T_freeze - 4)  # TODO: Make a parameter.
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
