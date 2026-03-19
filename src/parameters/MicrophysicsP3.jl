export ParametersP3
export MassPowerLaw, AreaPowerLaw, SlopePowerLaw, SlopeConstant, VentilationFactor

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
    "Scale [`μm^(2-σ)`]"
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
for a given Rᵢ [m² s⁻¹ °C⁻¹] by calling `ρ′_rim(Rᵢ)`.

The parameterization is given by:

```math
ρ'_{rim} = a + b R_i + c R_i^2, \\quad 1 ≤ R_i ≤ 8,
```
The range is extended to `R_i ≤ 12`, by linearly interpolating between 
`ρ′_rim(8)` and `ρ_ice = 900 kg/m³`. The latter is the solid bulk ice density.

For calculating Rᵢ, see [`compute_local_rime_density`](@ref CloudMicrophysics.P3Scheme.compute_local_rime_density).
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
function ((; a, b, c, ρ_ice)::LocalRimeDensity)(Rᵢ)
    Rᵢ = clamp(Rᵢ, 1, 12)  # P3 fortran code, microphy_p3.f90, Line 3315 clamps to 1 ≤ Rᵢ ≤ 12

    # Eq. 17 in Cober and List (1993), in [kg / m³], valid for 1 ≤ Rᵢ ≤ 8
    ρ′_rim_CL93(Rᵢ) = a + b * Rᵢ + c * Rᵢ^2

    ρ′_rim = if Rᵢ ≤ 8
        ρ′_rim_CL93(Rᵢ)
    else
        # following P3 fortran code, microphy_p3.f90, Line 3323
        #   https://github.com/P3-microphysics/P3-microphysics/blob/main/src/microphy_p3.f90#L3323
        # for 8 < Rᵢ ≤ 12, linearly interpolate between ρ′_rim(8) ≡ 611 kg/m³ and ρ_ice = 916.7 kg/m³
        ρ′_rim8 = ρ′_rim_CL93(8)
        f_ρ_ice = (Rᵢ - 8) / (12 - 8)
        (1 - f_ρ_ice) * ρ′_rim8 + f_ρ_ice * ρ_ice  # Linear interpolation beyond 8.
    end
    return ρ′_rim
end

### ----------------------------- ###
### --- TOP-LEVEL CONSTRUCTOR --- ###
### ----------------------------- ###

"""
    ParametersP3

Parameters for P3 bulk microphysics scheme.

From Morrison and Milbrandt (2015) [MorrisonMilbrandt2015](@cite)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParametersP3{FT, SLOPELAW <: SlopeLaw} <: ParametersType
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
    "Wet growth time scale [`s`]"
    τ_wet::FT
    "Cloud ice density [`kg m⁻³`]"
    ρ_i::FT
    "Cloud liquid water density [`kg m⁻³`]"
    ρ_l::FT
    "Water freeze temperature [`K`]"
    T_freeze::FT
end

"""
    ParametersP3(toml_dict::CP.ParamDict; [slope_law = :powerlaw])

Create a `ParametersP3` object from a `ClimaParams` TOML dictionary.

# Arguments
- `toml_dict::CP.ParamDict`: A `ClimaParams` TOML dictionary
- `slope_law`: Slope law to use (`:constant` or, by default, `:powerlaw`)

"""
function ParametersP3(toml_dict::CP.ParamDict; slope_law = :powerlaw)
    @assert slope_law in (:constant, :powerlaw)
    params = CP.get_parameter_values(toml_dict,
        (;
            :density_ice_water => :ρ_i,
            :density_liquid_water => :ρ_l,
            :temperature_water_freeze => :T_freeze,
            :P3_wet_growth_timescale => :τ_wet,
        ), "CloudMicrophysics")
    return ParametersP3(;
        mass = MassPowerLaw(toml_dict),
        area = AreaPowerLaw(toml_dict),
        slope = slope_law == :powerlaw ? SlopePowerLaw(toml_dict) : SlopeConstant(toml_dict),
        vent = VentilationFactor(toml_dict),
        ρ_rim_local = LocalRimeDensity(toml_dict),
        params...,
    )
end

### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

# Unit annotations for verbose show (used by ShowMethods.verbose_show_type_and_fields)
ShowMethods.field_units(::MassPowerLaw) = (; α_va = "kg m^(-β_va)")
ShowMethods.field_units(::AreaPowerLaw) = (; γ = "μm^(2-σ)")
ShowMethods.field_units(::SlopePowerLaw) = (; a = "m^b")
ShowMethods.field_units(::LocalRimeDensity) = (; ρ_ice = "kg m⁻³")
ShowMethods.field_units(::ParametersP3) =
    (; τ_wet = "s", ρ_i = "kg m⁻³", ρ_l = "kg m⁻³", T_freeze = "K")
