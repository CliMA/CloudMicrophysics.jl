export IceNucleationParameters
export Frostenberg2023

"""
    Mohler2006{FT}

Parameters for ice nucleation from Mohler et al 2006
DOI: 10.5194/acp-6-3007-2006

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Mohler2006{FT} <: ParametersType
    "max allowed supersaturation [-]"
    Sᵢ_max::FT
    "threshold temperature [K]"
    T_thr::FT
end

function Mohler2006(td::CP.ParamDict)
    name_map = (;
        :Mohler2006_maximum_allowed_Si => :Sᵢ_max,
        :Mohler2006_threshold_T => :T_thr,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Mohler2006(; parameters...)
end

"""
    Koop2000{FT}

Parameters for ice nucleation from Koop et al 2000
DOI: 10.1038/35020537

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Koop2000{FT} <: ParametersType
    "min Δaw [-]"
    Δa_w_min::FT
    "max Δaw [-]"
    Δa_w_max::FT
    "coefficient [-]"
    c₁::FT
    "coefficient [-]"
    c₂::FT
    "coefficient [-]"
    c₃::FT
    "coefficient [-]"
    c₄::FT
    "coefficient [-]"
    linear_c₁::FT
    "coefficient [-]"
    linear_c₂::FT
end

function Koop2000(td::CP.ParamDict)
    name_map = (;
        :Koop2000_min_delta_aw => :Δa_w_min,
        :Koop2000_max_delta_aw => :Δa_w_max,
        :Koop2000_J_hom_coeff1 => :c₁,
        :Koop2000_J_hom_coeff2 => :c₂,
        :Koop2000_J_hom_coeff3 => :c₃,
        :Koop2000_J_hom_coeff4 => :c₄,
        :Linear_J_hom_coeff1 => :linear_c₁,
        :Linear_J_hom_coeff2 => :linear_c₂,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Koop2000(; parameters...)
end

"""
    MorrisonMilbrandt2014{FT}

Parameters for ice nucleation from Morrison & Milbrandt (2015)
DOI: 10.1175/JAS-D-14-0065.1

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct MorrisonMilbrandt2014{FT} <: ParametersType
    "Cutoff temperature for deposition nucleation `[K]`; maps to the ClimaParams homogeneous nucleation temperature, 233 K"
    T_dep_thres::FT
    "coefficient [-]"
    c₁::FT
    "coefficient [-]"
    c₂::FT
    "freezing temperature of water [K]"
    T₀::FT
    "heterogeneous freezing parameter a [K⁻¹]"
    het_a::FT
    "heterogeneous freezing parameter B [m⁻³ s⁻¹]"
    het_B::FT
end

function MorrisonMilbrandt2014(td::CP.ParamDict)
    name_map = (;
        :temperature_homogenous_nucleation => :T_dep_thres,
        :Thompson2004_c1_Cooper => :c₁,
        :Thompson2004_c2_Cooper => :c₂,
        :temperature_water_freeze => :T₀,
        :BarklieGokhale1959_a_parameter => :het_a,
        :BarklieGokhale1959_B_parameter => :het_B,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return MorrisonMilbrandt2014(; parameters...)
end

export RainFreezing

"""
    RainFreezing{FT}

Parameters for heterogeneous (Bigg-type) immersion freezing of rain drops.

Stores the empirical Barklie-Gokhale (1959) / Bigg (1953) parameters
used by Morrison & Milbrandt (2015).

# Fields
$(DocStringExtensions.FIELDS)

# Callable interface

    (rf::RainFreezing)(T, T₀) → het_B * exp(het_a * (T₀ - T))

Compute the volumetric freezing rate [m⁻³ s⁻¹]
"""
@kwdef struct RainFreezing{FT} <: ParametersType
    "empirical parameter [K⁻¹]"
    het_a::FT
    "water-type dependent parameter [m⁻³ s⁻¹]"
    het_B::FT
end

function RainFreezing(td::CP.ParamDict)
    name_map = (;
        :BarklieGokhale1959_a_parameter => :het_a,
        :BarklieGokhale1959_B_parameter => :het_B,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return RainFreezing(; parameters...)
end

# Callable: returns the Bigg (1953) volumetric freezing rate [m⁻³(water) s⁻¹]
(rf::RainFreezing)(T, T₀) = rf.het_B * exp(rf.het_a * (T₀ - T))

"""
    IceNucleationParameters{DEP, HOM, P3_type}

Parameters for ice nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct IceNucleationParameters{DEP, HOM, P3_type} <: ParametersType
    "deposition nucleation parameters, e.g. [`Mohler2006`](@ref)"
    deposition::DEP
    "homogeneous nucleation parameters, e.g. [`Koop2000`](@ref)"
    homogeneous::HOM
    "P3 ice nucleation parameters, e.g. [`MorrisonMilbrandt2014`](@ref)"
    p3::P3_type
end

IceNucleationParameters(toml_dict::CP.ParamDict) =
    IceNucleationParameters(;
        deposition = Mohler2006(toml_dict),
        homogeneous = Koop2000(toml_dict),
        p3 = MorrisonMilbrandt2014(toml_dict),
    )


"""
    Frostenberg2023{FT}

Parameters for frequency distribution of INP concentration
DOI: 10.5194/acp-23-10883-2023

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Frostenberg2023{FT} <: ParametersType
    "standard deviation"
    σ::FT
    "coefficient"
    a::FT
    "coefficient"
    b::FT
    "freezing temperature [K]"
    T_freeze::FT
    "log of the coefficient `a`"
    log_a::FT = log(a)
end

function Frostenberg2023(td::CP.ParamDict)
    name_map = (;
        :Frostenberg2023_standard_deviation => :σ,
        :Frostenberg2023_a_coefficient => :a,
        :Frostenberg2023_b_coefficient => :b,
        :temperature_water_freeze => :T_freeze,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Frostenberg2023(; parameters...)
end

# ---------------------------------------------------------------------------
# INP target spectra for the deposition nucleation slot
# ---------------------------------------------------------------------------

export AbstractINPTargetSpectrum
export ExponentialSupercoolingINP

"""
    AbstractINPTargetSpectrum

The super-type for the parameterizations that supply a TARGET ice nucleating particle
concentration to the deposition nucleation slot.

A concrete subtype answers three questions, and nothing else:

  - what the target concentration is, through the callable `(inp)(T)` → `N_t` [m⁻³];
  - where the closure is active, through
    `HetIceNucleation.is_active(inp, T, S_i)`;
  - how fast the deficit is delivered, through
    `HetIceNucleation.delivery_rate(inp, mp, tps, T, S_i)` [s⁻¹].

The deficit-relaxation rate form itself belongs to the slot and is shared, so a new
spectrum is added by writing those three methods rather than a fourth rate body.
The slot is `HetIceNucleation.deposition_rate`.
"""
abstract type AbstractINPTargetSpectrum <: ParametersType end

"""
    ExponentialSupercoolingINP{FT}

An ice nucleating particle spectrum that is exponential in supercooling, with a ceiling:

```math
N_t(T) = \\min\\big(a \\exp(b (T_0 - T)),\\ N_{max}\\big)
```

The shipped values are those of Cooper (1986) in the form given by Thompson et al. (2004)
and used by Morrison and Milbrandt (2015) appendix C(a). Fletcher (1962) is the same law
with different coefficients; Meyers et al. (1992) is a different law, exponential in
supersaturation rather than in supercooling.

`T_thr` and `S_thr` are the activation window: the closure is active where the air is
colder than `T_thr` and supersaturated with respect to ice by at least `S_thr`. Both
belong to the spectrum rather than to the slot, because they state where the fit is
claimed to hold.

The window and the ceiling are the reference P3 implementation's, which applies the same
law as

```fortran
! module_mp_p3.f90 (WRF 4.6.0)
if (t(i,k).lt.258.15 .and. supi_cld.ge.0.05) then
   dum = 0.005*exp(0.304*(273.15-t(i,k)))*1000.*inv_rho(i,k)
   dum = min(dum,100.e3*inv_rho(i,k)*SCF(k))
```

so `a = 0.005 · 1000 = 5` m⁻³, `b = 0.304` K⁻¹, `T_thr = 258.15` K, `S_thr = 0.05` and
`N_max = 1.0e5` m⁻³, which is the Fortran's 100 per liter. `SCF` is the subgrid cloud
fraction and is one where subgrid cloud fraction is not used. The ceiling distinguishes
this closure from
[`MorrisonMilbrandt2014`](@ref)'s `P3_deposition_N_i`, which clamps the temperature instead
and so tops out an order of magnitude higher.

# Fields
$(DocStringExtensions.FIELDS)

# Callable interface

    (inp::ExponentialSupercoolingINP)(T) → min(a exp(b (T₀ - T)), N_max)

The target ice nucleating particle concentration [m⁻³].
"""
@kwdef struct ExponentialSupercoolingINP{FT} <: AbstractINPTargetSpectrum
    "prefactor [m⁻³]"
    a::FT
    "exponent coefficient [K⁻¹]"
    b::FT
    "reference (freezing) temperature [K]"
    T₀::FT
    "ceiling on the target concentration [m⁻³]"
    N_max::FT
    "activation temperature threshold [K]"
    T_thr::FT
    "activation ice-supersaturation threshold [-]"
    S_thr::FT
end

function ExponentialSupercoolingINP(td::CP.ParamDict)
    name_map = (;
        :P3_cooper_deposition_prefactor => :a,
        :P3_cooper_deposition_exponent_coefficient => :b,
        :temperature_water_freeze => :T₀,
        :P3_cooper_deposition_max_concentration => :N_max,
        :P3_cooper_deposition_temperature_threshold => :T_thr,
        :P3_cooper_deposition_ice_supersaturation_threshold => :S_thr,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return ExponentialSupercoolingINP(; parameters...)
end

# Callable: the target ice nucleating particle concentration [m⁻³]
((; a, b, T₀, N_max)::ExponentialSupercoolingINP)(T) = min(a * exp(b * (T₀ - T)), N_max)

ShowMethods.field_units(::ExponentialSupercoolingINP) =
    (; a = "m⁻³", b = "K⁻¹", T₀ = "K", N_max = "m⁻³", T_thr = "K")

# ---------------------------------------------------------------------------
# INP-activation memory models
# ---------------------------------------------------------------------------

export NIceProxyDepletion

"""
    NIceProxyDepletion

Use the in-cell ice number `n_ice` as the depletion proxy for INP
activation. This is the legacy / always-on form: a column with no ice
sees the full INP target; activation events do not by themselves
deplete the budget on a memory timescale, but the ice they create
proxies "INPs already used" downstream until that ice sublimates,
sediments out, or melts.

Conflates two physically distinct counts: "ice in column" and
"INPs already activated in this air parcel". Drop a fresh anvil into
clean air below it and the deposition channel artificially shuts off.

The delivery rate is the closure's, not this model's: see
`HetIceNucleation.delivery_rate`.
"""
struct NIceProxyDepletion end
