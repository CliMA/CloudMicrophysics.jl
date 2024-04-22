export ParametersP3

"""
    ParametersP3{FT}

Parameters for P3 bulk microphysics scheme from
Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

# Fields
$(DocStringExtensions.FIELDS)
"""
Base.@kwdef struct ParametersP3{FT} <: ParametersType{FT}
    "Coefficient in mass(size) relation [g/μm^β_va]"
    α_va::FT
    "Coefficient in mass(size) relation [-]"
    β_va::FT
    "Cloud ice density [kg/m3]"
    ρ_i::FT
    "Cloud liquid water density [kg/m3]"
    ρ_l::FT
    "Coefficient in area(size) for ice side plane, column, bullet, and planar polycrystal aggregates from Mitchell 1996 [μm^(2-σ)]"
    γ::FT
    "Coefficient in area(size) for ice side plane, column, bullet, and planar polycrystal aggregates from Mitchell 1996 [-]"
    σ::FT
    "Coefficient for shape parameter mu for ice. See eq 3 in Morrison and Milbrandt 2015. Units: [m^0.8]"
    a::FT
    "Coefficient for shape parameter mu for ice. See eq 3 in Morrison and Milbrandt 2015. Units: [-]"
    b::FT
    "Coefficient for shape parameter mu for ice. See eq 3 in Morrison and Milbrandt 2015. Units: [-]"
    c::FT
    "Limiter for shape parameter mu for ice. See eq 3 in Morrison and Milbrandt 2015. Units: [-]"
    μ_max::FT
    "Water freeze temperature [K]"
    T_freeze::FT
    "ventillation factor a"
    vent_a::FT
    "ventillation factor b"
    vent_b::FT
end

ParametersP3(::Type{FT}) where {FT <: AbstractFloat} =
    ParametersP3(CP.create_toml_dict(FT))

function ParametersP3(td::CP.AbstractTOMLDict)
    name_map = (;
        :BF1995_mass_coeff_alpha => :α_va,
        :BF1995_mass_exponent_beta => :β_va,
        :density_ice_water => :ρ_i,
        :density_liquid_water => :ρ_l,
        :M1996_area_coeff_gamma => :γ,
        :M1996_area_exponent_sigma => :σ,
        :Heymsfield_mu_coeff1 => :a,
        :Heymsfield_mu_coeff2 => :b,
        :Heymsfield_mu_coeff3 => :c,
        :Heymsfield_mu_cutoff => :μ_max,
        :temperature_water_freeze => :T_freeze,
        :p3_ventillation_a => :vent_a,
        :p3_ventiallation_b => :vent_b,
    )
    parameters = CP.get_parameter_values(td, name_map, "CloudMicrophysics")
    FT = CP.float_type(td)
    return ParametersP3{FT}(; parameters...)
end
