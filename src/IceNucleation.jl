# This file contains the `HetIceNucleation` and `HomIceNucleation` modules.

"""
Parameterization for heterogenous cloud ice nucleation.
"""
module HetIceNucleation

import ..Parameters as CMP

export dust_activated_number_fraction
export MohlerDepositionRate
export deposition_J
export ABIFM_J
export P3_deposition_N_i
export P3_het_N_i
export INP_concentration_frequency
export INP_concentration_mean

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
