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


# ---------------------------------------------------------------------------
# INP target spectra for deposition nucleation
# ---------------------------------------------------------------------------

export INPTargetSpectrum

"""
    INPTargetSpectrum

The supertype of the target ice nucleating particle (INP) concentrations of deposition nucleation.

They enter `HetIceNucleation.deposition_rate`. A subtype defines

  - the target concentration, `(inp)(T)` [m⁻³];
  - the activation conditions, `HetIceNucleation.is_active(inp, T, S_i)`;
  - the relaxation rate, `HetIceNucleation.delivery_rate(inp, mp, tps, T, S_i)` [s⁻¹].

A subtype with an INP budget for immersion freezing defines a method of `HetIceNucleation.immersion_limit_rate`.
"""
abstract type INPTargetSpectrum <: ParametersType end

"""
    Frostenberg2023{FT}

Parameters for frequency distribution of INP concentration
DOI: 10.5194/acp-23-10883-2023

As a target spectrum for deposition nucleation, its target is the median of the climatology, and it has an INP
budget for immersion freezing, `HetIceNucleation.immersion_limit_rate`.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct Frostenberg2023{FT} <: INPTargetSpectrum
    "standard deviation"
    σ::FT
    "coefficient"
    b::FT
    "freezing temperature [K]"
    T_freeze::FT
    "relaxation timescale of deposition nucleation and of the immersion budget [s]"
    τ::FT
    "temperature below which deposition nucleation is active [K]"
    T_thr::FT
    "ice supersaturation at and above which deposition nucleation is active [-]"
    S_thr::FT
end

function Frostenberg2023(td::CP.ParamDict)
    name_map = (;
        :Frostenberg2023_standard_deviation => :σ,
        :Frostenberg2023_b_coefficient => :b,
        :temperature_water_freeze => :T_freeze,
        :Frostenberg2023_deposition_relaxation_timescale => :τ,
        :P3_cooper_deposition_temperature_threshold => :T_thr,
        :P3_cooper_deposition_ice_supersaturation_threshold => :S_thr,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    return Frostenberg2023(; parameters...)
end

export ExponentialSupercoolingINP

"""
    ExponentialSupercoolingINP{FT}

An ice nucleating particle spectrum that is exponential in supercooling, with a ceiling:

```math
N_t(T) = \\min\\big(a \\exp(b (T_0 - T)),\\ N_{max}\\big)
```

The default values are those of Cooper (1986) in the form of Thompson et al. (2004), used by Morrison and
Milbrandt (2015), appendix C(a). Deposition nucleation is active where the air is colder than `T_thr` and its ice
supersaturation is at least `S_thr`. The ceiling `N_max` and the activation conditions follow the P3 reference
implementation.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ExponentialSupercoolingINP{FT} <: INPTargetSpectrum
    "prefactor [m⁻³]"
    a::FT
    "exponent coefficient [K⁻¹]"
    b::FT
    "reference (freezing) temperature [K]"
    T₀::FT
    "ceiling on the target concentration [m⁻³]"
    N_max::FT
    "temperature below which deposition nucleation is active [K]"
    T_thr::FT
    "ice supersaturation at and above which deposition nucleation is active [-]"
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

ShowMethods.field_units(::ExponentialSupercoolingINP) =
    (; a = "m⁻³", b = "K⁻¹", T₀ = "K", N_max = "m⁻³", T_thr = "K", S_thr = "-")

export StochasticINPSpectrum

"""
    StochasticINPSpectrum{FT, S}

A target spectrum that shifts the log concentration of the base spectrum `base`, with a ceiling and a lowest
temperature:

    (inp::StochasticINPSpectrum)(T) → min(exp(log(base(T)) + log_shift), N_max)

Deposition nucleation is active where the base spectrum is active and `T ≥ T_min`. With `log_shift = 0` and an
inactive ceiling and lowest temperature, the spectrum is the base spectrum.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct StochasticINPSpectrum{FT, S <: INPTargetSpectrum} <: INPTargetSpectrum
    "the base spectrum, whose log concentration is shifted"
    base::S
    "the shift of the log concentration [-]"
    log_shift::FT
    "the lowest temperature at which the spectrum is active [K]"
    T_min::FT
    "the ceiling on the target concentration [m⁻³]"
    N_max::FT
end

"""
    StochasticINPSpectrum(toml_dict::CP.ParamDict; log_shift, T_min, N_max)

A [`StochasticINPSpectrum`](@ref) on the [`Frostenberg2023`](@ref) base spectrum.
"""
function StochasticINPSpectrum(td::CP.ParamDict; log_shift, T_min, N_max)
    FT = CP.float_type(td)
    return StochasticINPSpectrum(;
        base = Frostenberg2023(td), log_shift = FT(log_shift), T_min = FT(T_min), N_max = FT(N_max),
    )
end

ShowMethods.field_units(::StochasticINPSpectrum) = (; log_shift = "-", T_min = "K", N_max = "m⁻³")

# ---------------------------------------------------------------------------
# Depletion models of the INP target
# ---------------------------------------------------------------------------

export NIceProxyDepletion

"""
    NIceProxyDepletion

Use the specific ice number `n_ice` of the cell state as the activated INP concentration.

In this form, a column with no ice sees the full target, and the ice that nucleation creates counts as activated
INPs until it sublimates, sediments out, or melts. It also counts ice that sedimented in from above, so ice falling
into clean air below suppresses nucleation there.
"""
struct NIceProxyDepletion end
