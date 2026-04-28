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

Parameters for ice nucleation from  Morrison & Milbrandt 2014
DOI: 10.1175/JAS-D-14-0065.1

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct MorrisonMilbrandt2014{FT} <: ParametersType
    "Cutoff temperature for deposition nucleation [K]"
    T_dep_thres::FT
    "coefficient [-]"
    c₁::FT
    "coefficient [-]"
    c₂::FT
    "T₀"
    T₀::FT
    "heterogeneous freezing parameter a [°C^-1]"
    het_a::FT
    "heterogeneous freezing parameter B [cm^-3 s^-1]"
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

Stores the empirical Barklie–Gokhale (1959) / Bigg (1953) parameters
used by Morrison & Milbrandt (2015).

# Fields
$(DocStringExtensions.FIELDS)

# Callable interface

    (rf::RainFreezing)(T, T₀) → het_B * 10⁶ * exp(het_a * (T₀ − T))

Returns the volumetric freezing rate in SI units [m⁻³ s⁻¹].
The stored `het_B` is in [cm⁻³ s⁻¹]; the factor 10⁶ converts cm³ → m³.
"""
@kwdef struct RainFreezing{FT} <: ParametersType
    "empirical parameter [°C⁻¹]"
    het_a::FT
    "water-type dependent parameter [cm⁻³ s⁻¹]"
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
# het_B is stored in [cm⁻³ s⁻¹]; multiply by 10⁶ to convert to [m⁻³(water) s⁻¹].
(rf::RainFreezing)(T, T₀) = rf.het_B * 1_000_000 * exp(rf.het_a * (T₀ - T))

"""
    IceNucleationParameters{FT, DEP, HOM, P3_type}

Parameters for ice nucleation

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct IceNucleationParameters{DEP, HOM, P3_type} <: ParametersType
    deposition::DEP
    homogeneous::HOM
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
# F23 INP-activation memory models
# ---------------------------------------------------------------------------

export AbstractINPDepletion, NIceProxyDepletion, PrognosticINPDecay

"""
    AbstractINPDepletion

Abstract type for the model of how F23 INP-activation budgets are
"depleted" within an air parcel. Choose a concrete subtype to control
the value subtracted from the F23 INPC target in the F23 deposition and
immersion-cap rates:

```
∂ₜn_frz = max(0, INPC(T)/ρ - n_active) / τ_act
```

where `n_active` is supplied by the host. The depletion model also
carries `τ_act`, the F23 activation relaxation timescale, so the host
doesn't need to wire that knob separately. The model selects how the
host sources `n_active`:

- [`NIceProxyDepletion`](@ref): use existing `n_ice` (zero memory).
- [`PrognosticINPDecay`](@ref): use a prognostic `n_INP_used` tracer
  with a finite memory timescale.
"""
abstract type AbstractINPDepletion end

"""
    NIceProxyDepletion{FT}

Use the in-cell ice number `n_ice` as the depletion proxy for F23
activation. This is the legacy / always-on form: a column with no ice
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
struct NIceProxyDepletion{FT} <: AbstractINPDepletion
    "F23 activation relaxation timescale [s] (default `300`)"
    τ_act::FT
end
NIceProxyDepletion(; τ_act = 300) = NIceProxyDepletion(τ_act)

"""
    PrognosticINPDecay{FT}

Phillips-style activation memory: the host maintains a passive scalar
`n_INP_used` whose budget is

```
∂t (ρ n_INP_used) = -∇·(ρ u n_INP_used)                  # advection
                    + ρ (∂ₜn_frz_F23_dep + ∂ₜn_imm_F23)  # source: just-activated INPs
                    - ρ n_INP_used / τ_INP_decay          # relaxation
```

so a parcel that activated INPs at low altitude carries that memory as
it lofts to colder levels, but the memory decays with timescale
`τ_INP_decay` representing aerosol re-supply from the unmodelled
mixing/scavenging cycle. The activation/relaxation knob is `τ_INP_decay`;
the F23 activation timescale `τ_act` is also carried for consistency:

| `τ_INP_decay`        | Behaviour                                                |
|----------------------|----------------------------------------------------------|
| → 0                  | collapses to current always-on `(INPC − 0)/τ_act`        |
| ~ τ_act (~minutes)   | rapid replenishment, F23 keeps firing in steady state    |
| ~ hours              | activation memory persists across a cumulus lifetime     |
| → ∞                  | strict Phillips depletion (INPs gone forever per parcel) |

# Fields
$(DocStringExtensions.FIELDS)

# Boundary conditions

A reasonable first cut is `n_INP_used = 0` at the BL floor (fresh
aerosols enter with the full INPC budget) and zero-gradient at the top
of domain. Rigorous treatment couples the floor source to the surface
aerosol flux, which we don't model.

# Compared to full Phillips-Demott

A full multi-species Phillips depletion (4–10 aerosol tracers) would
distinguish dust / soot / biological / sulphate INPs and track each
species' activation/scavenging budget. This single-scalar form captures
the dominant physical effect (activation memory + finite re-supply
timescale) without aerosol speciation. It can be replaced by the full
framework later when the model gains an aerosol prognostic.
"""
struct PrognosticINPDecay{FT} <: AbstractINPDepletion
    "F23 activation relaxation timescale [s] (default `300`)"
    τ_act::FT
    "Decay timescale of activation memory [s] (τ_INP_decay → 0 recovers no-memory; τ → ∞ recovers strict depletion)"
    τ_INP_decay::FT
end
function PrognosticINPDecay(; τ_INP_decay, τ_act = 300)
    FT = typeof(τ_INP_decay)
    return PrognosticINPDecay{FT}(FT(τ_act), τ_INP_decay)
end
