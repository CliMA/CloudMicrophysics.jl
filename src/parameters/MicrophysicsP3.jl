export ParametersP3
export MassPowerLaw, AreaPowerLaw, SlopePowerLaw, SmoothSlopePowerLaw, SlopeConstant, VentilationFactor
export DEFAULT_SLOPE_LAW, DEFAULT_ASPECT_RATIO
export ice_seed

### ----------------------------- ###
### --- SUB-PARAMETERIZATIONS --- ###
### ----------------------------- ###

"""
    MassPowerLaw{FT}

Parameters for mass(size) relation.

From measurements of mass grown by vapor diffusion and aggregation in midlatitude cirrus
by Brown and Francis (1995) [BrownFrancis1995](@cite)

A part of the [`ParametersP3`](@ref) parameter set.

!!! note
    The `BF1995_mass_coeff_alpha` parameter is provided in units of [`g μm^(-β_va)`]
    but the `α_va` field is stored in SI-like units of [`kg m^(-β_va)`] for
    consistency with the rest of the code.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct MassPowerLaw{FT} <: ParametersType
    "Coefficient in mass(size) relation [`kg m^(-β_va)`]"
    α_va::FT
    "Coefficient in mass(size) relation [`-`]"
    β_va::FT
end
function MassPowerLaw(toml_dict::CP.ParamDict)
    name_map = (;
        :BF1995_mass_coeff_alpha => :α_va,
        :BF1995_mass_exponent_beta => :β_va,
    )
    (; β_va) = p = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    α_va = p.α_va * 10^(6 * β_va - 3)
    FT = CP.float_type(toml_dict)
    return MassPowerLaw{FT}(; α_va, β_va)
end

"""
    AreaPowerLaw{FT}

Parameters for area(size) relation.

```math
A(D) = γ D^σ
```

where `γ` and `σ` are coefficients in area(size) for ice side plane, column, bullet,
and planar polycrystal aggregates. Values are from Mitchell (1996) [Mitchell1996](@cite)

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AreaPowerLaw{FT} <: ParametersType
    "Scale [`m^(2-σ)`]"
    γ::FT
    "Power [`-`]"
    σ::FT
end
function AreaPowerLaw(toml_dict::CP.ParamDict)
    name_map = (; :M1996_area_coeff_gamma => :γ, :M1996_area_exponent_sigma => :σ)
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return AreaPowerLaw{FT}(; params...)
end

"""
    SlopeLaw

The top-level super-type for slope parameterizations.

See [`SlopePowerLaw`](@ref) and [`SlopeConstant`](@ref) for concrete implementations.
"""
abstract type SlopeLaw <: ParametersType end

"""
    SlopePowerLaw{FT}

Slope parameter μ as a power law in shape parameter λ:

```math
μ(λ) = a λ^b - c
```

and is limited to:

```math
0 ≤ μ ≤ μ_{max}
```

See also Eq. 3 in Morrison and Milbrandt (2015) [MorrisonMilbrandt2015](@cite)

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SlopePowerLaw{FT} <: SlopeLaw
    "Scale [`m^b`]"
    a::FT
    "Power [`-`]"
    b::FT
    "Offset [`-`]"
    c::FT
    "Upper limiter [`-`]"
    μ_max::FT
end
function SlopePowerLaw(toml_dict::CP.ParamDict)
    name_map = (;
        :Heymsfield_mu_coeff1 => :a,
        :Heymsfield_mu_coeff2 => :b,
        :Heymsfield_mu_coeff3 => :c,
        :Heymsfield_mu_cutoff => :μ_max,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return SlopePowerLaw(; params...)
end

"""
    SmoothSlopePowerLaw{FT}

Slope parameter μ as a power law in slope λ, identical to [`SlopePowerLaw`](@ref)
but with the `0 ≤ μ ≤ μ_max` limiters applied as smooth transitions of sharpness `κ`:

```math
μ(λ) = M_{μ_{max}}\\big(M_0(a λ^b - c)\\big)
```

where `M_0(x) = κ^{-1} \\log(1 + e^{κ x})` is a smooth `max(x, 0)` and
`M_{μ_{max}}(x) = μ_{max} - κ^{-1} \\log(1 + e^{κ (μ_{max} - x)})` is a smooth
`min(x, μ_{max})`. The corner width in μ is of order `1/κ`; as `κ → ∞`, `μ(λ)`
converges to the hard-clamped [`SlopePowerLaw`](@ref).

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SmoothSlopePowerLaw{FT} <: SlopeLaw
    "Scale [`m^b`]"
    a::FT
    "Power [`-`]"
    b::FT
    "Offset [`-`]"
    c::FT
    "Upper limiter [`-`]"
    μ_max::FT
    "Corner sharpness [`-`]"
    κ::FT
end
function SmoothSlopePowerLaw(toml_dict::CP.ParamDict)
    name_map = (;
        :Heymsfield_mu_coeff1 => :a,
        :Heymsfield_mu_coeff2 => :b,
        :Heymsfield_mu_coeff3 => :c,
        :Heymsfield_mu_cutoff => :μ_max,
        :P3_mu_smoothing_sharpness => :κ,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return SmoothSlopePowerLaw(; params...)
end

"""
    SlopeConstant{FT}

Slope parameter μ as a constant:

```math
μ(λ) = μ_{const}
```

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SlopeConstant{FT} <: SlopeLaw
    "Slope parameter μ [`-`]"
    μ::FT
end
function SlopeConstant(toml_dict::CP.ParamDict)
    name_map = (; :P3_constant_slope_parameterization_value => :μ)
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return SlopeConstant(; params...)
end

"""
    VentilationFactor{FT}

Parameters for ventilation factor:

```math
F(D) = a_{v} + b_{v}  N_{Sc}^{1/3} N_{Re}(D)^{1/2}
```
where `N_{Sc}` is the Schmidt number and `N_{Re}(D)` is the Reynolds number for a particle with diameter `D`.

From Seifert and Beheng (2006) [SeifertBeheng2006](@cite),
see also Eq. (13-61) in Pruppacher and Klett (2010) [PruppacherKlett2010](@cite)

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct VentilationFactor{FT} <: ParametersType
    "Constant coefficient in ventilation factor [`-`]"
    aᵥ::FT
    "Linear coefficient in ventilation factor [`-`]"
    bᵥ::FT
end
function VentilationFactor(toml_dict::CP.ParamDict)
    name_map = (;
        :SB2006_ventilation_factor_coeff_av => :aᵥ,
        :SB2006_ventilation_factor_coeff_bv => :bᵥ,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return VentilationFactor(; params...)
end

"""
    LocalRimeDensity{FT}
    (ρ′_rim::LocalRimeDensity)(Rᵢ)

Local rime density parameterization based on Cober and List (1993) [CoberList1993](@cite),
Eq. 16 and 17.

Given an instance `ρ′_rim::LocalRimeDensity`, obtain the local rime density
for a given Rᵢ [μm m s⁻¹ °C⁻¹] by calling `ρ′_rim(Rᵢ)`.

The parameterization is given by:

```math
ρ'_{rim} = a + b R_i + c R_i^2, \\quad 1 ≤ R_i ≤ 8,
```
The range is extended to `R_i ≤ 12`, by linearly interpolating between
`ρ′_rim(8)` and `ρ_ice = 916.7 kg/m³`. The latter is the solid bulk ice density.

For calculating Rᵢ, see [`compute_local_rime_density`](@ref CloudMicrophysics.P3Scheme.compute_local_rime_density).

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct LocalRimeDensity{FT} <: ParametersType
    "Constant coefficient"
    a::FT
    "Linear coefficient"
    b::FT
    "Quadratic coefficient"
    c::FT
    "Density of solid bulk ice [`kg m⁻³`]"
    ρ_ice::FT
end
function LocalRimeDensity(toml_dict::CP.ParamDict)
    name_map = (;
        :CL1993_local_rime_density_constant_coeff => :a,
        :CL1993_local_rime_density_linear_coeff => :b,
        :CL1993_local_rime_density_quadratic_coeff => :c,
        :density_ice_water => :ρ_ice,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    return LocalRimeDensity(; params...)
end
# Bounds of the Cober and List (1993) riming index Rᵢ, following the P3 fortran code,
# `microphy_p3.f90`, Line 3315. At RIME_DENSITY_Rᵢ_MAX the rime density equals the solid bulk
# ice density, so the Rᵢ → ∞ limit is solid ice.
const RIME_DENSITY_Rᵢ_MIN = 1
const RIME_DENSITY_Rᵢ_MAX = 12

function ((; a, b, c, ρ_ice)::LocalRimeDensity)(Rᵢ)
    Rᵢ = clamp(Rᵢ, RIME_DENSITY_Rᵢ_MIN, RIME_DENSITY_Rᵢ_MAX)

    # Eq. 17 in Cober and List (1993), in [kg / m³], valid for 1 ≤ Rᵢ ≤ 8
    ρ′_rim_CL93(Rᵢ) = a + b * Rᵢ + c * Rᵢ^2

    ρ′_rim = if Rᵢ ≤ 8
        ρ′_rim_CL93(Rᵢ)
    else
        # following P3 fortran code, microphy_p3.f90, Line 3323
        #   https://github.com/P3-microphysics/P3-microphysics/blob/main/src/microphy_p3.f90#L3323
        # for 8 < Rᵢ ≤ 12, linearly interpolate between ρ′_rim(8) ≡ 611 kg/m³ and ρ_ice = 916.7 kg/m³
        ρ′_rim8 = ρ′_rim_CL93(8)
        f_ρ_ice = (Rᵢ - 8) / (RIME_DENSITY_Rᵢ_MAX - 8)
        (1 - f_ρ_ice) * ρ′_rim8 + f_ρ_ice * ρ_ice  # Linear interpolation beyond 8.
    end
    return ρ′_rim
end

"""
    AspectRatio

Aspect-ratio treatment for the ice terminal-velocity correction. Each subtype
is a functor `(state, D)` returning the multiplicative velocity factor:
`Oblate` returns `cbrt(ϕᵢ(state, D))`, `NoAspectRatio` returns `1`.
The functor methods are defined in `P3Scheme`, where `ϕᵢ` is available.
"""
abstract type AspectRatio end
struct Oblate <: AspectRatio end
struct NoAspectRatio <: AspectRatio end

### ----------------------------- ###
### --- TOP-LEVEL CONSTRUCTOR --- ###
### ----------------------------- ###

"""
    DEFAULT_SLOPE_LAW

The single source of truth for the `ParametersP3` slope law default. Referenced by the
constructors that forward to it ([`P3IceParams`](@ref), [`Microphysics2MParams`](@ref))
so that changing the default is one edit rather than one per forwarding layer. See also
[`DEFAULT_ASPECT_RATIO`](@ref).
"""
const DEFAULT_SLOPE_LAW = :smooth_powerlaw

"""
    DEFAULT_ASPECT_RATIO

The single source of truth for the `ParametersP3` aspect-ratio default, forwarded the
same way as [`DEFAULT_SLOPE_LAW`](@ref).
"""
const DEFAULT_ASPECT_RATIO = Oblate()

"""
    IceStickingEfficiency{FT}

The efficiency with which colliding ice particles stick, as two dimensionless factors on the
ice self-collection (aggregation) rate.

Without them the scheme collects every geometric encounter - an efficiency of exactly one at every
temperature and rime fraction - which is the one value that is certainly wrong. Aggregation
efficiency is far below unity and falls steeply as the air cools, the usual account being the
quasi-liquid surface layer that lets crystals adhere near the melting point and is absent when cold.

`eii(T)` ramps linearly from `e_cold` at or below `T_cold` to `e_warm` at and above the freezing
point. `Eii_fact(F_rim)` is one at or below `F_rim_lo`, falls linearly to zero at `F_rim_hi`, and is
exactly zero above it: heavily rimed ice is smooth and dense and does not aggregate.

PROVENANCE, stated exactly because a search did not support the obvious claim. These values are the
P3 REFERENCE IMPLEMENTATION's own, adopted for consistency with the scheme this one implements.
They are NOT documented in Morrison and Milbrandt (2015), whose appendix has no ice self-collection
section and which lists "aggregation and riming efficiencies" among the uncertainties it does not
address, nor in the 2025 three-moment P3 paper; and the reference source carries no citation for
them, alongside a commented-out predecessor with values differing by up to a factor of a hundred.
Treat them as an uncertain tuned parameter on the reference's own account.

Note the ramp slope is DERIVED from the endpoints here, where the reference hard-codes its
reciprocal span - so changing `T_cold` keeps the ramp consistent rather than silently breaking it.
"""
@kwdef struct IceStickingEfficiency{FT} <: ParametersType
    "sticking efficiency at and below `T_cold` [-]"
    e_cold::FT
    "sticking efficiency at and above the freezing point [-]"
    e_warm::FT
    "temperature at and below which the efficiency is `e_cold` [K]"
    T_cold::FT
    "rime mass fraction at and below which collection is unreduced [-]"
    F_rim_lo::FT
    "rime mass fraction at and above which collection is shut off [-]"
    F_rim_hi::FT
end

IceStickingEfficiency(toml_dict::CP.ParamDict) = IceStickingEfficiency(;
    CP.get_parameter_values(toml_dict,
        (;
            :P3_ice_sticking_efficiency_cold => :e_cold,
            :P3_ice_sticking_efficiency_warm => :e_warm,
            :P3_ice_sticking_efficiency_T_cold => :T_cold,
            :P3_ice_collection_rime_shutoff_start => :F_rim_lo,
            :P3_ice_collection_rime_shutoff_end => :F_rim_hi,
        ), "CloudMicrophysics")...)

"""
    ParametersP3

Parameters for P3 bulk microphysics scheme.

From Morrison and Milbrandt (2015) [MorrisonMilbrandt2015](@cite)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParametersP3{FT, SLOPELAW <: SlopeLaw, AR <: AspectRatio} <: ParametersType
    "Mass-size relation, e.g. [`MassPowerLaw`](@ref)"
    mass::MassPowerLaw{FT}
    "Area-size relation, e.g. [`AreaPowerLaw`](@ref)"
    area::AreaPowerLaw{FT}
    "Slope relation, e.g. [`SlopePowerLaw`](@ref) or [`SlopeConstant`](@ref)"
    slope::SLOPELAW
    "Ventilation relation, e.g. [`VentilationFactor`](@ref)"
    vent::VentilationFactor{FT}
    "Local rime density, e.g. [`LocalRimeDensity`](@ref)"
    ρ_rim_local::LocalRimeDensity{FT}
    "Ice aggregation sticking efficiency, an [`IceStickingEfficiency`](@ref)"
    sticking::IceStickingEfficiency{FT}
    "Wet growth time scale [`s`]"
    τ_wet::FT
    "Ice number-concentration adjustment timescale [`s`]"
    τ_numadj::FT
    "Diameter of a nascent ice crystal [`m`]; see [`ice_seed`](@ref)"
    D_nuc::FT
    "Cloud ice density [`kg m⁻³`]"
    ρ_i::FT
    "Cloud liquid water density [`kg m⁻³`]"
    ρ_l::FT
    "Water freeze temperature [`K`]"
    T_freeze::FT
    "Terminal-velocity aspect-ratio treatment, an [`AspectRatio`](@ref)"
    aspect_ratio::AR = DEFAULT_ASPECT_RATIO
end

"""
    ParametersP3(toml_dict::CP.ParamDict; [slope_law = DEFAULT_SLOPE_LAW], [aspect_ratio = DEFAULT_ASPECT_RATIO])

Create a `ParametersP3` object from a `ClimaParams` TOML dictionary.

# Arguments
- `toml_dict::CP.ParamDict`: A `ClimaParams` TOML dictionary
- `slope_law`: Slope law to use (`:constant`, `:powerlaw`, or, by default, `:smooth_powerlaw`).
  The default is the SMOOTH law: the hard-clamped [`SlopePowerLaw`](@ref) is only `C^0` at its two
  kinks, and the lower (`μ = 0`) kink makes `log(L/N)` locally increase with `λ`, so the shape map
  `λ ↦ log(L/N)` loses monotonicity and `get_distribution_logλ` admits multiple solutions there.
  [`SmoothSlopePowerLaw`](@ref) removes that by construction. `:powerlaw` is retained for
  reproducing earlier results, not because it is preferred.
- `aspect_ratio`: an [`AspectRatio`](@ref); by default, [`DEFAULT_ASPECT_RATIO`](@ref)

"""
function ParametersP3(toml_dict::CP.ParamDict;
    slope_law = DEFAULT_SLOPE_LAW,
    aspect_ratio = DEFAULT_ASPECT_RATIO,
)
    @assert slope_law in (:constant, :powerlaw, :smooth_powerlaw)
    params = CP.get_parameter_values(toml_dict,
        (;
            :density_ice_water => :ρ_i,  # TODO: Use `WaterProperties` struct for ice and liquid water density
            :density_liquid_water => :ρ_l,
            :temperature_water_freeze => :T_freeze,
            :P3_wet_growth_timescale => :τ_wet,
            :P3_ice_number_adjustment_timescale => :τ_numadj,
            :P3_ice_nucleation_diameter => :D_nuc,
        ), "CloudMicrophysics")
    slope = if slope_law == :powerlaw
        SlopePowerLaw(toml_dict)
    elseif slope_law == :smooth_powerlaw
        SmoothSlopePowerLaw(toml_dict)
    else
        SlopeConstant(toml_dict)
    end
    return ParametersP3(;
        mass = MassPowerLaw(toml_dict),
        area = AreaPowerLaw(toml_dict),
        slope,
        vent = VentilationFactor(toml_dict),
        ρ_rim_local = LocalRimeDensity(toml_dict),
        sticking = IceStickingEfficiency(toml_dict),
        aspect_ratio,
        params...,
    )
end

### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

"""
    ice_seed(p3)

The nascent ice crystal: the smallest particle the scheme creates, as the three quantities
its consumers need.

# Arguments
 - `p3`: the [`ParametersP3`](@ref) scheme parameters.

# Returns
 - A `NamedTuple` `(; m_nuc, r_nuc, ρ_i)`:
    + `m_nuc`: the mass of one nascent crystal [kg], `ρ_i (π/6) D_nuc³`, a solid ice sphere
      of the nascent diameter.
    + `r_nuc`: the nascent radius [m], `D_nuc/2`.
    + `ρ_i`: the solid ice density [kg m⁻³] the mass is built from.

This is the SINGLE SOURCE for the nascent crystal, and it is deliberately shared rather
than restated per process. Its consumers are the deposition nucleation starter mass and
seed delivery time, the lower bound of the ice number adjustment and the presence predicate
that goes with it, the orphan-ice drain timescale, the upper edge of the shape-solve
bracket, and the conduction-limited melt-fraction bound. Each of those is an argument about
the smallest particle in the population, so they must move together: a number adjustment
whose lower mass bound differed from the mass at which crystals are created would adjust a
freshly nucleated population on the step it was created.

Only the diameter is a free parameter (`P3_ice_nucleation_diameter`); the mass and the
radius are derived here so that no consumer can hold its own copy of the relation.
"""
@inline function ice_seed((; D_nuc, ρ_i)::ParametersP3)
    r_nuc = D_nuc / 2
    m_nuc = ρ_i * π * D_nuc^3 / 6
    return (; m_nuc, r_nuc, ρ_i)
end

# Unit annotations for verbose show (used by ShowMethods.verbose_show_type_and_fields)
ShowMethods.field_units(::MassPowerLaw) = (; α_va = "kg m^(-β_va)")
ShowMethods.field_units(::AreaPowerLaw) = (; γ = "m^(2-σ)")
ShowMethods.field_units(::SlopePowerLaw) = (; a = "m^b")
ShowMethods.field_units(::SmoothSlopePowerLaw) = (; a = "m^b")
ShowMethods.field_units(::LocalRimeDensity) = (; ρ_ice = "kg m⁻³")
ShowMethods.field_units(::ParametersP3) =
    (; τ_wet = "s", τ_numadj = "s", D_nuc = "m", ρ_i = "kg m⁻³", ρ_l = "kg m⁻³", T_freeze = "K")
