export ParametersP3
export MassPowerLaw, AreaPowerLaw, SlopePowerLaw, SlopeConstant, VentilationSB2005

"""
    MassPowerLaw{FT}

Parameters for mass(size) relation.

From measurements of mass grown by vapor diffusion and aggregation in midlatitude cirrus by Brown and Francis (1995)
doi: 10.1175/1520-0426(1995)012<0410:IMOTIW>2.0.CO;2

!!! note
    The `BF1995_mass_coeff_alpha` parameter is provided in units of [g μm^(-β_va)], but
    the `α_va` field is stored in SI-like units of [kg m^(-β_va)] for consistency with the rest of the code.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct MassPowerLaw{FT} <: ParametersType{FT}
    "Coefficient in mass(size) relation [kg m^(-β_va)]"
    α_va::FT
    "Coefficient in mass(size) relation [-]"
    β_va::FT
end
function MassPowerLaw(td::CP.AbstractTOMLDict)
    name_map = (;
        :BF1995_mass_coeff_alpha => :α_va,
        :BF1995_mass_exponent_beta => :β_va,
    )
    (; β_va) = p = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    α_va = p.α_va * 10^(6 * β_va - 3)
    FT = CP.float_type(td)
    return MassPowerLaw{FT}(; α_va, β_va)
end

"""
    AreaPowerLaw{FT}

Parameters for area(size) relation.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct AreaPowerLaw{FT} <: ParametersType{FT}
    "Coefficient in area(size) for ice side plane, column, bullet, and planar polycrystal aggregates from Mitchell 1996 [μm^(2-σ)]"
    γ::FT
    "Coefficient in area(size) for ice side plane, column, bullet, and planar polycrystal aggregates from Mitchell 1996 [-]"
    σ::FT
end
function AreaPowerLaw(td::CP.AbstractTOMLDict)
    name_map = (;
        :M1996_area_coeff_gamma => :γ,
        :M1996_area_exponent_sigma => :σ,
    )
    params = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return AreaPowerLaw{FT}(; params...)
end

"""
    SlopeLaw{FT}

The top-level super-type for slope parameterizations.
"""
abstract type SlopeLaw{FT} <: ParametersType{FT} end

"""
    SlopePowerLaw{FT}

Slope parameter μ as a power law in shape parameter λ:
    μ(λ) = a λ^b - c
and is limited to:
    0 ≤ μ ≤ μ_max

See also Eq. 3 in Morrison and Milbrandt 2015.

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SlopePowerLaw{FT} <: SlopeLaw{FT}
    "Scale, a [m^0.8]"
    a::FT
    "Power, b [-]"
    b::FT
    "Offset, c [-]"
    c::FT
    "Upper limiter, μ_max [-]"
    μ_max::FT
end
function SlopePowerLaw(td::CP.AbstractTOMLDict)
    name_map = (;
        :Heymsfield_mu_coeff1 => :a,
        :Heymsfield_mu_coeff2 => :b,
        :Heymsfield_mu_coeff3 => :c,
        :Heymsfield_mu_cutoff => :μ_max,
    )
    params = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return SlopePowerLaw{FT}(; params...)
end

"""
    SlopeConstant{FT}

Slope parameter μ as a constant:
    μ(λ) = μ_const

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct SlopeConstant{FT} <: SlopeLaw{FT}
    "Slope parameter μ"
    μ::FT
end
function SlopeConstant(td::CP.AbstractTOMLDict)
    name_map = (;
        :Heymsfield_mu_coeff1 => :μ,
    )
    params = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return SlopeConstant{FT}(; params...)
end

"""
    VentilationSB2005{FT}

Parameters for ventilation factor:
    F(r) = a_{vent} + b_{vent}  N_{Sc}^{1/3} N_{Re}(r)^{1/2}

From Seifert and Beheng 2005, doi: 10.1007/s00703-005-0112-4

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct VentilationSB2005{FT} <: ParametersType{FT}
    "ventilation factor a"
    vent_a::FT
    "ventilation factor b"
    vent_b::FT
end
function VentilationSB2005(td::CP.AbstractTOMLDict)
    name_map = (;
        :p3_ventillation_a => :vent_a, # TODO: fix typo in TOML
        :p3_ventiallation_b => :vent_b,
    )
    params = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return VentilationSB2005{FT}(; params...)
end

"""
    ParametersP3{FT}

Parameters for P3 bulk microphysics scheme.
    
From Morrison and Milbrandt 2015, doi: 10.1175/JAS-D-14-0065.1

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ParametersP3{FT} <: ParametersType{FT}
    "Mass-size relation"
    mass::MassPowerLaw{FT}
    "Area-size relation"
    area::AreaPowerLaw{FT}
    "Slope relation"
    slope::SlopeLaw{FT}
    "Ventilation relation"
    vent::VentilationSB2005{FT}
    "Cloud ice density [kg/m3]"
    ρ_i::FT
    "Cloud liquid water density [kg/m3]"
    ρ_l::FT
    "Water freeze temperature [K]"
    T_freeze::FT
end

ParametersP3(::Type{FT}) where {FT <: AbstractFloat} =
    ParametersP3(CP.create_toml_dict(FT))

function ParametersP3(td::CP.AbstractTOMLDict; slope_law = :powerlaw)
    @assert slope_law in (:constant, :powerlaw)
    mass = MassPowerLaw(td)
    area = AreaPowerLaw(td)
    slope = if slope_law == :powerlaw
        SlopePowerLaw(td)
    else
        SlopeConstant(td)
    end
    vent = VentilationSB2005(td)
    name_map = (;
        :density_ice_water => :ρ_i,
        :density_liquid_water => :ρ_l,
        :temperature_water_freeze => :T_freeze,
    )
    params = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return ParametersP3{FT}(; mass, area, slope, vent, params...)
end


### ----------------- ###
### ----- UTILS ----- ###
### ----------------- ###

function Base.show(io::IO, p::ParametersP3{FT}) where {FT}
    println(io, "ParametersP3{$FT}")
    io_opts = ( :indent => "│   ", :typeinfo => false)
    show(IOContext(io, :prefix => "├── mass: ", io_opts...), p.mass)
    show(IOContext(io, :prefix => "├── area: ", io_opts...), p.area)
    show(IOContext(io, :prefix => "├── slope: ", io_opts...), p.slope)
    show(IOContext(io, :prefix => "├── vent: ", io_opts...), p.vent)
    println(io, "├── ρ_i = $(p.ρ_i) [kg m⁻³]")
    println(io, "├── ρ_l = $(p.ρ_l) [kg m⁻³]")
    println(io, "└── T_freeze = $(p.T_freeze) [K]")
end

function Base.show(io::IO, p::MassPowerLaw{FT}) where {FT}
    indent = get(io, :indent, "")
    prefix = get(io, :prefix, "")
    type = get(io, :typeinfo, true) ? "{$FT}" : ""
    println(io, "$(prefix)MassPowerLaw$type")
    println(io, "$(indent)├── α_va = $(p.α_va) [kg μm^(-β_va)]")
    println(io, "$(indent)└── β_va = $(p.β_va) [-]")
end

function Base.show(io::IO, p::AreaPowerLaw{FT}) where {FT}
    indent = get(io, :indent, "")
    prefix = get(io, :prefix, "")
    type = get(io, :typeinfo, true) ? "{$FT}" : ""
    println(io, "$(prefix)AreaPowerLaw$type")
    println(io, "$(indent)├── γ = $(p.γ) [μm^(2-σ)]")
    println(io, "$(indent)└── σ = $(p.σ) [-]")
end

function Base.show(io::IO, p::SlopePowerLaw{FT}) where {FT}
    indent = get(io, :indent, "")
    prefix = get(io, :prefix, "")
    type = get(io, :typeinfo, true) ? "{$FT}" : ""
    println(io, "$(prefix)SlopePowerLaw$type")
    println(io, "$(indent)├── a = $(p.a) [m^b]")
    println(io, "$(indent)├── b = $(p.b) [-]")
    println(io, "$(indent)├── c = $(p.c) [-]")
    println(io, "$(indent)└── μ_max = $(p.μ_max) [-]")
end

function Base.show(io::IO, p::SlopeConstant{FT}) where {FT}
    indent = get(io, :indent, "")
    prefix = get(io, :prefix, "")
    type = get(io, :typeinfo, true) ? "{$FT}" : ""
    println(io, "$(prefix)SlopeConstant$type")
    println(io, "$(indent)└── μ = $(p.μ) [-]")
end

function Base.show(io::IO, p::VentilationSB2005{FT}) where {FT}
    indent = get(io, :indent, "")
    prefix = get(io, :prefix, "")
    type = get(io, :typeinfo, true) ? "{$FT}" : ""
    println(io, "$(prefix)VentilationSB2005$type")
    println(io, "$(indent)├── vent_a = $(p.vent_a) [-]")
    println(io, "$(indent)└── vent_b = $(p.vent_b) [-]")
end
