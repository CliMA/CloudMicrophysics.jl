export ParametersP3

"""
    ParametersP3{FT}

Parameters for P3 bulk microphysics scheme from
Morrison and Milbrandt 2015 doi: 10.1175/JAS-D-14-0065.1

# Fields
$(DocStringExtensions.FIELDS)
"""
struct ParametersP3{FT} <: ParametersType{FT}
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
end

function ParametersP3(
    ::Type{FT},
    toml_dict::CP.AbstractTOMLDict = CP.create_toml_dict(FT),
) where {FT}
    (; data) = toml_dict
    return ParametersP3(
        FT(data["BF1995_mass_coeff_alpha"]["value"]),
        FT(data["BF1995_mass_exponent_beta"]["value"]),
        FT(data["density_ice_water"]["value"]),
        FT(data["density_liquid_water"]["value"]),
        FT(data["M1996_area_coeff_gamma"]["value"]),
        FT(data["M1996_area_exponent_sigma"]["value"]),
        FT(data["Heymsfield_mu_coeff1"]["value"]),
        FT(data["Heymsfield_mu_coeff2"]["value"]),
        FT(data["Heymsfield_mu_coeff3"]["value"]),
        FT(data["Heymsfield_mu_cutoff"]["value"]),
    )
end
