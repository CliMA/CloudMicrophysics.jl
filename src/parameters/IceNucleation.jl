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

# ---------------------------------------------------------------------------
# F23 INP-activation memory models
# ---------------------------------------------------------------------------

export NIceProxyDepletion

"""
    NIceProxyDepletion{FT}

Use the in-cell ice number `n_ice` as the depletion proxy for F23
activation. In this form, a column with no ice
sees the full INPC target; activation events do not by themselves
deplete the budget on a memory timescale, but the ice they create
proxies "INPs already used" downstream until that ice sublimates,
sediments out, or melts.

Conflates two physically distinct counts: "ice in column" and
"INPs already activated in this air parcel". Drop a fresh anvil into
clean air below it and the F23 channel artificially shuts off.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct NIceProxyDepletion{FT}
    "F23 activation relaxation timescale `[s]` (default `300`)"
    τ_act::FT
end
NIceProxyDepletion(; τ_act = 300) = NIceProxyDepletion(τ_act)
