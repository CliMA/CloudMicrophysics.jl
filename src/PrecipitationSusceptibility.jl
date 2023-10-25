module PrecipitationSusceptibility

import ..Microphysics2M as CM2
import ..Parameters as CMP

using ForwardDiff

export precipitation_susceptibility_autoconversion
export precipitation_susceptibility_accretion

"""
A structure containing the logarithmic derivatives of the production of
precipitation with respect to the specific humidities and number
densities of liquid and rain water.
"""
@kwdef struct precip_susceptibility_rates{FT}
    d_ln_pp_d_ln_q_liq::FT = FT(0)
    d_ln_pp_d_ln_q_rai::FT = FT(0)
    d_ln_pp_d_ln_N_liq::FT = FT(0)
    d_ln_pp_d_ln_N_rai::FT = FT(0)
end

"""
    precipitation_susceptibility_autoconversion(param_set, scheme, q_liq, q_rai, ρ, N_liq)

- `scheme` - type for 2-moment rain autoconversion parameterization
- `q_liq` - cloud water specific humidity
- `q_rai` - rain water specific humidity
- `ρ` - air density
- `N_liq` - cloud droplet number density

Returns the precipitation susceptibility rates due to autoconversion as a `precip_susceptibility_rates`
object, using automatic differentiation.
Works for any 2-moment scheme, as long as autoconversion is defined for it.
"""
function precipitation_susceptibility_autoconversion(
    scheme::CMP.SB2006{FT},
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT}
    grad = ForwardDiff.gradient(
        x -> log(CM2.autoconversion(scheme.acnv, exp.(x)...).dq_rai_dt),
        log.(abs.([q_liq, q_rai, ρ, N_liq])),
    )
    return precip_susceptibility_rates(
        d_ln_pp_d_ln_q_liq = grad[1],
        d_ln_pp_d_ln_q_rai = grad[2],
        d_ln_pp_d_ln_N_liq = grad[4],
        d_ln_pp_d_ln_N_rai = FT(0),
    )
end

"""
    precipitation_susceptibility_accretion(param_set, scheme, q_liq, q_rai, ρ, N_liq)

- `scheme` - type for 2-moment rain autoconversion parameterization
- `q_liq` - cloud water specific humidity
- `q_rai` - rain water specific humidity
- `ρ` - air density
- `N_liq` - cloud droplet number density

Returns the precipitation susceptibility rates due to accretion as a `precip_susceptibility_rates`
object, using automatic differentiation.
Works for any 2-moment scheme, as long as accretion is defined for it.
"""
function precipitation_susceptibility_accretion(
    scheme::CMP.SB2006{FT},
    q_liq::FT,
    q_rai::FT,
    ρ::FT,
    N_liq::FT,
) where {FT}

    grad = ForwardDiff.gradient(
        x -> log(CM2.accretion(scheme, exp.(x)...).dq_rai_dt),
        log.(abs.([q_liq, q_rai, ρ, N_liq])),
    )
    return precip_susceptibility_rates(
        d_ln_pp_d_ln_q_liq = grad[1],
        d_ln_pp_d_ln_q_rai = grad[2],
        d_ln_pp_d_ln_N_liq = grad[4],
        d_ln_pp_d_ln_N_rai = FT(0),
    )
end

end # module
