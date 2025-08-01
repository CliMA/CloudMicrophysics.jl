module PrecipitationSusceptibility

import ..Microphysics2M as CM2
import ..Parameters as CMP

using ForwardDiff

export precipitation_susceptibility_autoconversion
export precipitation_susceptibility_accretion

"""
A structure containing the logarithmic derivatives of the production of
precipitation with respect to the specific contents and number
densities of cloud liquid water and rain water.
"""
Base.@kwdef struct precip_susceptibility_rates{FT}
    d_ln_pp_d_ln_q_lcl::FT = FT(0)
    d_ln_pp_d_ln_q_rai::FT = FT(0)
    d_ln_pp_d_ln_N_lcl::FT = FT(0)
    d_ln_pp_d_ln_N_rai::FT = FT(0)
end

"""
    precipitation_susceptibility_autoconversion(param_set, scheme, q_lcl, q_rai, ρ, N_lcl)

- `scheme` - type for 2-moment rain autoconversion parameterization
- `q_lcl` - cloud liquid water specific content
- `q_rai` - rain water specific content
- `ρ` - air density
- `N_lcl` - cloud droplet number density

Returns the precipitation susceptibility rates due to autoconversion as a `precip_susceptibility_rates`
object, using automatic differentiation.
Works for any 2-moment scheme, as long as autoconversion is defined for it.
"""
function precipitation_susceptibility_autoconversion(
    scheme::CMP.SB2006{FT},
    q_lcl::FT,
    q_rai::FT,
    ρ::FT,
    N_lcl::FT,
) where {FT}
    grad = ForwardDiff.gradient(
        x -> log(
            CM2.autoconversion(scheme.acnv, scheme.pdf_c, exp.(x)...).dq_rai_dt,
        ),
        log.(abs.([q_lcl, q_rai, ρ, N_lcl])),
    )
    return precip_susceptibility_rates(
        d_ln_pp_d_ln_q_lcl = grad[1],
        d_ln_pp_d_ln_q_rai = grad[2],
        d_ln_pp_d_ln_N_lcl = grad[4],
        d_ln_pp_d_ln_N_rai = FT(0),
    )
end

"""
    precipitation_susceptibility_accretion(param_set, scheme, q_lcl, q_rai, ρ, N_lcl)

- `scheme` - type for 2-moment rain autoconversion parameterization
- `q_lcl` - cloud liquid water specific content
- `q_rai` - rain water specific content
- `ρ` - air density
- `N_lcl` - cloud droplet number density

Returns the precipitation susceptibility rates due to accretion as a `precip_susceptibility_rates`
object, using automatic differentiation.
Works for any 2-moment scheme, as long as accretion is defined for it.
"""
function precipitation_susceptibility_accretion(
    scheme::CMP.SB2006{FT},
    q_lcl::FT,
    q_rai::FT,
    ρ::FT,
    N_lcl::FT,
) where {FT}

    grad = ForwardDiff.gradient(
        x -> log(CM2.accretion(scheme, exp.(x)...).dq_rai_dt),
        log.(abs.([q_lcl, q_rai, ρ, N_lcl])),
    )
    return precip_susceptibility_rates(
        d_ln_pp_d_ln_q_lcl = grad[1],
        d_ln_pp_d_ln_q_rai = grad[2],
        d_ln_pp_d_ln_N_lcl = grad[4],
        d_ln_pp_d_ln_N_rai = FT(0),
    )
end

end # module
