export ParametersP3
export MassPowerLaw, AreaPowerLaw, SlopePowerLaw, SlopeConstant, VentilationSB2005

### ----------------------------- ###
### --- SUB-PARAMETERIZATIONS --- ###
### ----------------------------- ###

"""
    MassPowerLaw{FT}

Parameters for mass(size) relation.

From measurements of mass grown by vapor diffusion and aggregation in midlatitude cirrus by Brown and Francis (1995)
doi: 10.1175/1520-0426(1995)012<0410:IMOTIW>2.0.CO;2

A part of the [`ParametersP3`](@ref) parameter set.

!!! note
    The `BF1995_mass_coeff_alpha` parameter is provided in units of [`g μm^(-β_va)`]
    but the `α_va` field is stored in SI-like units of [`kg m^(-β_va)`] for 
    consistency with the rest of the code.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct MassPowerLaw{FT} <: ParametersType{FT}
    "Coefficient in mass(size) relation [`kg m^(-β_va)`]"
    α_va::FT
    "Coefficient in mass(size) relation [`-`]"
    β_va::FT
end
function MassPowerLaw(toml_dict::CP.AbstractTOMLDict)
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

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct AreaPowerLaw{FT} <: ParametersType{FT}
    "Coefficient in area(size) for ice side plane, column, bullet, and planar polycrystal aggregates from Mitchell 1996 [`μm^(2-σ)`]"
    γ::FT
    "Coefficient in area(size) for ice side plane, column, bullet, and planar polycrystal aggregates from Mitchell 1996 [`-`]"
    σ::FT
end
function AreaPowerLaw(toml_dict::CP.AbstractTOMLDict)
    name_map =
        (; :M1996_area_coeff_gamma => :γ, :M1996_area_exponent_sigma => :σ)
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return AreaPowerLaw{FT}(; params...)
end

"""
    SlopeLaw{FT}

The top-level super-type for slope parameterizations.

See [`SlopePowerLaw`](@ref) and [`SlopeConstant`](@ref) for concrete implementations.
"""
abstract type SlopeLaw{FT} <: ParametersType{FT} end

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

See also Eq. 3 in Morrison and Milbrandt 2015.

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct SlopePowerLaw{FT} <: SlopeLaw{FT}
    "Scale [`m^b`]"
    a::FT
    "Power [`-`]"
    b::FT
    "Offset [`-`]"
    c::FT
    "Upper limiter [`-`]"
    μ_max::FT
end
function SlopePowerLaw(toml_dict::CP.AbstractTOMLDict)
    name_map = (;
        :Heymsfield_mu_coeff1 => :a,
        :Heymsfield_mu_coeff2 => :b,
        :Heymsfield_mu_coeff3 => :c,
        :Heymsfield_mu_cutoff => :μ_max,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return SlopePowerLaw{FT}(; params...)
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
@kwdef struct SlopeConstant{FT} <: SlopeLaw{FT}
    "Slope parameter μ [`-`]"
    μ::FT
end
function SlopeConstant(toml_dict::CP.AbstractTOMLDict)
    name_map = (; :P3_constant_slope_parameterization_value => :μ)
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return SlopeConstant{FT}(; params...)
end

"""
    VentilationSB2005{FT}

Parameters for ventilation factor:

```math
F(r) = a_{vent} + b_{vent}  N_{Sc}^{1/3} N_{Re}(r)^{1/2}
```

From Seifert and Beheng 2005, doi: 10.1007/s00703-005-0112-4

A part of the [`ParametersP3`](@ref) parameter set.

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct VentilationSB2005{FT} <: ParametersType{FT}
    "Ventilation factor a [`-`]"
    vent_a::FT
    "Ventilation factor b [`-`]"
    vent_b::FT
end
function VentilationSB2005(toml_dict::CP.AbstractTOMLDict)
    name_map = (;
        :p3_ventillation_a => :vent_a, # TODO: fix typo in TOML
        :p3_ventiallation_b => :vent_b,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return VentilationSB2005{FT}(; params...)
end

### ----------------------------- ###
### --- TOP-LEVEL CONSTRUCTOR --- ###
### ----------------------------- ###

"""
    ParametersP3{FT}

Parameters for P3 bulk microphysics scheme.
    
From Morrison and Milbrandt 2015 [MorrisonMilbrandt2015](@cite)

# Fields
$(DocStringExtensions.FIELDS)
"""
@kwdef struct ParametersP3{FT, SLOPELAW <: SlopeLaw{FT}} <: ParametersType{FT}
    "Mass-size relation"
    mass::MassPowerLaw{FT}
    "Area-size relation"
    area::AreaPowerLaw{FT}
    "Slope relation"
    slope::SLOPELAW
    "Ventilation relation"
    vent::VentilationSB2005{FT}
    "Cloud ice density [`kg m⁻³`]"
    ρ_i::FT
    "Cloud liquid water density [`kg m⁻³`]"
    ρ_l::FT
    "Water freeze temperature [`K`]"
    T_freeze::FT
end

"""
    ParametersP3(FT)

Create a `ParametersP3` object from a floating point type `FT`.

# Examples

```jldoctest
julia> import CloudMicrophysics.Parameters as CMP

julia> CMP.ParametersP3(Float64)
ParametersP3{Float64}
├── mass: MassPowerLaw
│   ├── α_va = 0.018537721864540644 [kg μm^(-β_va)]
│   └── β_va = 1.9 [-]
├── area: AreaPowerLaw
│   ├── γ = 0.2285 [μm^(2-σ)]
│   └── σ = 1.88 [-]
├── slope: SlopePowerLaw
│   ├── a = 0.00191 [m^b]
│   ├── b = 0.8 [-]
│   ├── c = 2.0 [-]
│   └── μ_max = 6.0 [-]
├── vent: VentilationSB2005
│   ├── vent_a = 0.78 [-]
│   └── vent_b = 0.308 [-]
├── ρ_i = 916.7 [kg m⁻³]
├── ρ_l = 1000.0 [kg m⁻³]
└── T_freeze = 273.15 [K]
```
"""
ParametersP3(::Type{FT}; kw...) where {FT} = ParametersP3(CP.create_toml_dict(FT); kw...)

"""
    ParametersP3(toml_dict::CP.AbstractTOMLDict; [slope_law = :powerlaw])

Create a `ParametersP3` object from a `ClimaParams` TOML dictionary.

# Arguments
- `toml_dict::CP.AbstractTOMLDict`: A `ClimaParams` TOML dictionary
- `slope_law`: Slope law to use (`:constant` or, by default, `:powerlaw`)

"""
function ParametersP3(toml_dict::CP.AbstractTOMLDict; slope_law = :powerlaw)
    @assert slope_law in (:constant, :powerlaw)
    mass = MassPowerLaw(toml_dict)
    area = AreaPowerLaw(toml_dict)
    slope = if slope_law == :powerlaw
        SlopePowerLaw(toml_dict)
    else
        SlopeConstant(toml_dict)
    end
    vent = VentilationSB2005(toml_dict)
    name_map = (;
        :density_ice_water => :ρ_i,
        :density_liquid_water => :ρ_l,
        :temperature_water_freeze => :T_freeze,
    )
    params = CP.get_parameter_values(toml_dict, name_map, "CloudMicrophysics")
    FT = CP.float_type(toml_dict)
    return ParametersP3{FT, typeof(slope)}(; mass, area, slope, vent, params...)
end

### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

function Base.show(
    io::IO,
    p::Union{
        ParametersP3,
        MassPowerLaw,
        AreaPowerLaw,
        SlopePowerLaw,
        SlopeConstant,
        VentilationSB2005,
    },
)
    indent = get(io, :indent, "")
    prefix = get(io, :prefix, "")

    # Get type information
    T = typeof(p)
    FT = eltype(p)
    type = get(io, :typeinfo, true) ? "{$FT}" : ""

    # Print type name, handling nested types
    type_name = T.name.name
    println(io, "$(prefix)$(type_name)$type")

    # Print each field
    fields = fieldnames(T)
    for (i, field) in enumerate(fields)
        value = getfield(p, field)
        is_last = i == length(fields)
        prefix_char = is_last ? "└" : "├"

        if typeof(value) <: Number
            # Simple value - print directly with unit
            unit = _get_parameter_unit(field)
            println(
                io,
                "$(indent)$(prefix_char)── $(field) = $(value) [$(unit)]",
            )
        else
            # Nested struct - recursively show with proper context
            nested_io = IOContext(
                io,
                :prefix => "$(prefix_char)── $(field): ",
                :indent => "$(indent)│   ",
                :typeinfo => false,
            )
            show(nested_io, value)
        end
    end
end

function _get_parameter_unit(field::Symbol)
    units = Dict(
        :α_va => "kg μm^(-β_va)",
        :γ => "μm^(2-σ)",
        :a => "m^b",
        :ρ_i => "kg m⁻³",
        :ρ_l => "kg m⁻³",
        :T_freeze => "K",
    )
    # unitless parameters: β_va, σ, b, c, μ_max, μ, vent_a, vent_b
    return get(units, field, "-")
end
