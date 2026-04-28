# This file contains the `HetIceNucleation` and `HomIceNucleation` modules.

"""
Parameterization for heterogenous cloud ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP
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
export f23_immersion_limit_rate
export f23_deposition_rate
export n_active, INP_relaxation_tendency

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
    MohlerDepositionRate(dust, ip, S_i, T, dSi_dt, N_aer)

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

Calculate the deposition nucleation rate coefficient, `J` [m⁻² s⁻¹]
for different minerals in liquid droplets.

# Arguments
  - `dust`: a struct with dust parameters
  - `Δa_w`: change in water activity [unitless].

See [China2017](@cite) for details on the parameterization.
Returns zero for unsupported aerosol types.
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
 - `Nᵢ`: number of ice crystals nucleated via deposition nucleation with units of m^-3.

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
    + `het_a`: empirical parameter [C⁻¹]
    + `het_B`: water-type dependent parameter [cm⁻³ s⁻¹]
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
    Vₗ_cm³ = Vₗ * 1_000_000  # converted from m^3 to cm^3
    Tₛ = T₀ - T
    return Nₗ * (1 - exp(-het_B * Vₗ_cm³ * Δt * exp(het_a * Tₛ)))
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
  - `T`: air temperature [K]

The mean `log(INPC)` is given by:
```
μ(T) = log(-(T_celsius / 10)^9)
```
with the corresponding INPC obtained by exponentiating: `exp(μ(T))`.

For details see: [Frostenberg2023](@cite), doi.org/10.5194/acp-23-10883-2023
"""
function INP_concentration_mean((; T_freeze)::CMP.Frostenberg2023, T)
    T_celsius = min(T - T_freeze, 0)
    return 9log(-T_celsius / 10)  # = log((-T_celsius / 10)^9)
end

"""
    liquid_freezing_rate(parameterization, pdf, tps, q_rai, ρ, N_rai, T)

Compute the rate of liquid water freezing into ice.

# Arguments
 - `parameterization`: The [`CMP.RainFreezing`](@ref) parameterization.
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
function liquid_freezing_rate(parameterization::CMP.RainFreezing, pdf, tps, q, ρ, N, T)
    FT = eltype(q)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    (; ρw) = pdf  # [kg(water) m⁻³(water)]
    n = N / ρ     # specific number concentration [kg⁻¹(air)]

    # Solve for the pdf parameters
    (; Dr_mean) = CM2.pdf_rain_parameters(pdf, q, ρ, N)

    # Bigg (1953) volumetric freezing rate:
    J_bigg = parameterization(T, T_freeze)  # [m⁻³(water) s⁻¹]

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
    ϵₘ, ϵₙ = UT.ϵ_numerics_2M_M(FT), UT.ϵ_numerics_2M_N(FT)
    cond = (n > ϵₙ) & (q > ϵₘ) & (T < T_freeze - 4)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(n))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(q))

    return (; ∂ₜn_frz, ∂ₜq_frz)
end

"""
    liquid_freezing_rate(parameterization, pdf::CMP.CloudParticlePDF_SB2006, tps, q, ρ, N, T)

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
∂ₜn_frz = J_bigg · (π/6)   · M_D³(N₀c, λc, νcD, μcD)
∂ₜq_frz = J_bigg · ρw · (π/6)² · M_D⁶(N₀c, λc, νcD, μcD)
```

with `M_Dᵏ` the kth diameter moment computed by
[`DT.generalized_gamma_Mⁿ`](@ref). Volume-weighting → bigger drops freeze
first, captured analytically; no numerical PSD integration required.

# Arguments
 - `parameterization`: The [`CMP.RainFreezing`](@ref) parameterization
   (Bigg / Barklie–Gokhale parameters; the `Rain` in the name is historical —
   the kinetics apply to any liquid-drop PSD).
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
    parameterization::CMP.RainFreezing, pdf::CMP.CloudParticlePDF_SB2006,
    tps, q, ρ, N, T,
)
    FT = eltype(q)
    T_freeze = TDI.TD.Parameters.T_freeze(tps)
    (; ρw) = pdf
    n = N / ρ     # specific number concentration [kg⁻¹(air)]

    # Solve for the diameter-space PDF parameters.
    (; λc, νcD, μcD) = CM2.pdf_cloud_parameters(pdf, q, ρ, N)

    # Bigg (1953) volumetric freezing rate per unit drop water volume.
    J_bigg = parameterization(T, T_freeze)  # [m⁻³(water) s⁻¹]

    # Diameter moments via the closed-form generalized-gamma formula.
    # When N or L is essentially zero, `pdf_cloud_parameters` returns
    # `λc = Inf`, which makes `B^(-k/μ) → 0`, so the moments vanish — no
    # extra guard needed.
    M_D³ = DT.generalized_gamma_Mⁿ(νcD, μcD, λc, n, 3)  # [m³ · kg⁻¹(air)]
    M_D⁶ = DT.generalized_gamma_Mⁿ(νcD, μcD, λc, n, 6)  # [m⁶ · kg⁻¹(air)]

    V_1 = FT(π) / 6
    ∂ₜn_frz = J_bigg * V_1 * M_D³          # [kg⁻¹(air) s⁻¹]
    ∂ₜq_frz = J_bigg * ρw * V_1^2 * M_D⁶   # [kg(water) kg⁻¹(air) s⁻¹]

    # Same gates as the rain branch: non-trivial number/mass and T < -4 °C.
    ϵₘ, ϵₙ = UT.ϵ_numerics_2M_M(FT), UT.ϵ_numerics_2M_N(FT)
    cond = (n > ϵₙ) & (q > ϵₘ) & (T < T_freeze - 4)
    ∂ₜn_frz = ifelse(cond, ∂ₜn_frz, zero(n))
    ∂ₜq_frz = ifelse(cond, ∂ₜq_frz, zero(q))

    return (; ∂ₜn_frz, ∂ₜq_frz)
end

"""
    f23_immersion_limit_rate(f23_params, T, ρ; τ, inpc_log_shift, n_active)

Compute the **F23-INPC-imposed upper limit** on the cloud-droplet immersion
freezing number rate.

The Frostenberg 2023 climatology specifies a target INP concentration `INPC(T)`
in air [m⁻³(air)]. Treating that target as a budget that should be activated
on a relaxation timescale `τ`, the maximum number of crystals nucleated per
kg of air per second is

```
∂ₜn_lim = max(0, INPC(T)/ρ - n_active) / τ.
```

`n_active` is the depletion proxy supplied by the host (see
[`AbstractINPDepletion`](@ref) and [`n_active`](@ref)): with the default
`NIceProxyDepletion` model the host passes `n_ice` (or anything; pass
zero to recover the no-depletion form), and with `PrognosticINPDecay`
the host passes `n_INP_used`.

This is the "INPC-only" leg of Pathway 2 — the cap that forces the realized
immersion-freezing rate to fall below pure Bigg kinetics in clean-air or
low-q_lcl regimes. Combined with [`liquid_freezing_rate`](@ref) on the cloud
PSD, the realized Pathway-2 rate is `min(bigg, ∂ₜn_lim)`.

# Arguments
 - `f23_params`: The [`CMP.Frostenberg2023`](@ref) parameters.
 - `T`: Air temperature [K].
 - `ρ`: Air density [kg(air) m⁻³(air)].

# Keyword arguments
 - `τ`: Relaxation timescale [s] (default `300`).
 - `inpc_log_shift`: Additive shift to `log(INPC)` (e.g. an OU-SIF stochastic
   excursion). Default `0`.
 - `n_active`: Depletion proxy [kg⁻¹(air)] (default `0`, i.e. no depletion;
   recovers the no-memory cap).

# Returns
 - A `NamedTuple` `(; ∂ₜn_frz)` — the specific number freezing-rate cap
   [kg⁻¹(air) s⁻¹]. Zero when `T ≥ T_freeze`.
"""
function f23_immersion_limit_rate(
    f23_params::CMP.Frostenberg2023, T, ρ;
    τ = oftype(T, 300), inpc_log_shift = zero(T),
    n_active = zero(T),
)
    T ≥ f23_params.T_freeze && return (; ∂ₜn_frz = zero(T))
    log_inpc = INP_concentration_mean(f23_params, T) + inpc_log_shift
    INPC_per_kg = exp(log_inpc) / ρ                  # [kg⁻¹(air)]
    ∂ₜn_frz = max(zero(T), INPC_per_kg - n_active) / τ # [kg⁻¹(air) s⁻¹]
    return (; ∂ₜn_frz)
end

"""
    f23_deposition_rate(f23_params, tps, T, ρ, q_vap, n_ice;
                        τ, inpc_log_shift, T_thresh, S_i_thresh, ρ_i, D_nuc)

Compute the **F23 deposition nucleation rate** (Pathway 1).

The rate is the Frostenberg 2023 INP concentration (treated as a budget)
relaxed toward depletion at `n_ice` over timescale `τ`:

```
∂ₜn_frz = max(0, INPC(T)/ρ - n_ice) / τ
```

Each newly nucleated crystal is assigned a starter mass

```
m_starter = (π/6) · ρ_i · D_nuc³
```

(default ≈ 4.8e-13 kg for a 10 μm solid-ice crystal). The mass tendency is
the implied mass injection, capped by half the local vapor excess over
ice saturation per relaxation window — this prevents the channel from
sourcing more vapor than physically available when OU-SIF excursions push
INPC above the local supply:

```
q_excess = max(0, q_vap - q_sat_ice)
∂ₜq_frz  = min(m_starter · ∂ₜn_frz,  ½ q_excess / τ)
```

`q_sat_ice` is computed internally from `tps` via
`TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)`.

Two physical gates close both rates outside the activation window:

  - `T < T_thresh` (default `T_freeze − 15 K`, i.e. −15 °C): below this
    temperature deposition nucleation is active.
  - `S_i ≡ q_vap/q_sat_ice − 1 > S_i_thresh` (default 5 %): ice
    supersaturation is required for vapor to nucleate onto INPs.

The `n_ice` argument is a proxy for already-activated INPs; in a Phillips-
style framework it would be replaced by an explicit prognostic
`N_INP_active`.

# Arguments
 - `f23_params`: The [`CMP.Frostenberg2023`](@ref) parameters.
 - `tps`: Thermodynamics parameters (used for the ice-saturation curve).
 - `T`: Air temperature [K].
 - `ρ`: Air density [kg(air) m⁻³(air)].
 - `q_vap`: Water-vapor specific content [kg(vap) kg⁻¹(air)].
 - `n_ice`: Specific ice-crystal number concentration [kg⁻¹(air)] (proxy
   for already-activated INPs).

# Keyword arguments
 - `τ`: Relaxation timescale [s] (default `300`).
 - `inpc_log_shift`: Additive shift to `log(INPC)` (default `0`).
 - `T_thresh`: Activation temperature threshold [K] (default
   `f23_params.T_freeze − 15`).
 - `S_i_thresh`: Activation ice-supersaturation threshold (default `0.05`).
 - `ρ_i`: Solid-ice density used for the starter mass [kg/m³] (default
   `916.7`).
 - `D_nuc`: Nascent crystal diameter [m] (default `10e-6`, the small-D
   tail of the P3 distribution).

# Returns
 - A `NamedTuple` `(; ∂ₜn_frz, ∂ₜq_frz)` with the specific number rate
   [kg⁻¹(air) s⁻¹] and specific mass rate [kg(ice) kg⁻¹(air) s⁻¹]. Zero
   outside the activation window.
"""
function f23_deposition_rate(
    f23_params::CMP.Frostenberg2023, tps, T, ρ, q_tot, q_liq, q_ice, n_ice; m_nuc,
    T_thresh, S_i_thresh, τ_act = 300, inpc_log_shift = 0,
)
    q_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T, ρ)
    q_vap = TDI.q_vap(q_tot, q_liq, q_ice)
    S_i = q_vap / q_sat_ice - 1
    # Activation gates.
    gated = (T < T_thresh) & (S_i > S_i_thresh)
    # Nucleation rate, limited by ambient INP availability and relaxation time τ_act.
    log_inpc = INP_concentration_mean(f23_params, T) + inpc_log_shift
    INPC_per_kg = exp(log_inpc) / ρ
    ∂ₜn_frz = max(0, INPC_per_kg - n_ice) / τ_act
    ∂ₜn_frz = ifelse(gated, ∂ₜn_frz, zero(∂ₜn_frz))
    # Vapor-excess cap on the implied mass injection.
    q_excess = max(0, q_vap - q_sat_ice)
    # Implied mass injection, capped by half the local vapor excess per relaxation window.
    ∂ₜq_frz = min(m_nuc * ∂ₜn_frz, q_excess / (2τ_act))
    return (; ∂ₜn_frz, ∂ₜq_frz)
end

# ---------------------------------------------------------------------------
# F23 INP-activation memory dispatch
# ---------------------------------------------------------------------------

"""
    n_active(model::CMP.AbstractINPDepletion, n_ice, n_INP_used)

Return the depletion proxy `n_active` to subtract from the F23 INPC
target in [`f23_deposition_rate`](@ref) and any analogous
INPC-budgeted rate. Dispatches on the host's chosen depletion model:

- `NIceProxyDepletion`: returns `n_ice` (legacy "n_ice as proxy" form).
- `PrognosticINPDecay`: returns `n_INP_used` (the prognostic activation-
  memory tracer).

The host always passes both `n_ice` and `n_INP_used`; only one is used
per call, but keeping both at the call site lets us swap models without
restructuring caller code.
"""
@inline n_active(::CMP.NIceProxyDepletion,    n_ice, n_INP_used) = n_ice
@inline n_active(::CMP.PrognosticINPDecay,    n_ice, n_INP_used) = n_INP_used

"""
    INP_relaxation_tendency(model, n_INP_used)

Return the relaxation contribution to `∂ₜ n_INP_used` (the host's
prognostic activation-memory tracer):

```
∂ₜ n_INP_used (relax) = -n_INP_used / τ_INP_decay
```

This is the **decay** half of the n_INP_used budget; the **source**
half (just-activated INPs from F23 deposition + cloud-immersion) is
added at the BMT call site where those rates are available.

For `NIceProxyDepletion`, returns 0 — the host does not maintain
`n_INP_used` in that mode.
"""
@inline INP_relaxation_tendency(::CMP.NIceProxyDepletion, n_INP_used) = zero(n_INP_used)
@inline INP_relaxation_tendency(model::CMP.PrognosticINPDecay, n_INP_used) =
    -n_INP_used / model.τ_INP_decay

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
