"""
    CloudMicrophysicsParametersP3{FT}

Parameters for P3 bulk microphysics scheme

# Fields
$(DocStringExtensions.FIELDS)
"""
struct CloudMicrophysicsParametersP3{FT} <: AbstractCloudMicrophysicsParameters
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
end

CloudMicrophysicsParametersP3(param_set) = CloudMicrophysicsParametersP3(
    param_set.α_va_BF1995,
    param_set.β_va_BF1995,
    param_set.ρ_cloud_ice,
    param_set.ρ_cloud_liq,
    param_set.γ_M1996,
    param_set.σ_M1996,
)
