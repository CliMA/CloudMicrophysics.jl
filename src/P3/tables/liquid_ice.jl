"""
    collision_point(params, vel, pdf, generator, (logλ, w, ρ_rim, logρₐ, logx̄))

Return the collected rates of one liquid species at a point of its collision table, without the freeze and shed partition.

The value is an `SVector` of the first three elements of [`collection_integrals`](@ref): the
collected number rate per unit ice number and per unit liquid number [m³/s], the collected mass
rate per unit ice number and per unit liquid mass [m³/s], and the impact scale `log κ̄`. The table
stores the logarithms of the two rates. Per unit liquid number, the liquid size distribution
depends only on its mean mass, see [`liquid_table_scales`](@ref).
"""
function collision_point(params, vel, pdf, generator, (logλ, w, ρ_rim, logρₐ, logx̄))
    r = collection_integrals(params, vel, pdf, generator, logλ, w, ρ_rim, logρₐ, logx̄, ())
    return SA.SVector(r[1], r[2], r[3])
end

"""
    rime_density_point(params, vel, pdf, generator, logR, (logλ, w, ρ_rim, logρₐ, logx̄))

Return the mean local rime density [kg/m³] of the mass of one liquid species that the ice collects, at each point of `logR`.

The mean rime density is the ratio of the collected mass rate to the rime volume rate of
[`collection_integrals`](@ref). The local rime density depends on the temperature `T` [°C]
through `Rᵢ = κ / |T|`, see [`rime_density_at`](@ref), so a point of `logR` stands for the
temperature `T = -exp(log κ̄ - logR)` at the impact scale `log κ̄` of the state.
"""
function rime_density_point(params, vel, pdf, generator, logR, (logλ, w, ρ_rim, logρₐ, logx̄))
    K = Val(length(logR))
    args = (params, vel, pdf, generator, logλ, w, ρ_rim, logρₐ, logx̄)
    log_impact = collection_integrals(args..., ())[3]
    r = collection_integrals(args..., ntuple(k -> -exp(log_impact - logR[k]), K))
    return SA.SVector(ntuple(k -> r[2] / r[k + 3], K))
end

"""
    collection_integrals(params, vel, pdf, generator, logλ, w, ρ_rim, logρₐ, logx̄, T°Cs)

Return the collected rates of one liquid species by the ice of [`table_state`](@ref), without the freeze and shed partition.

The result is an `SVector` of
1. the number rate per unit ice number and per unit liquid number [m³/s],
2. the mass rate per unit ice number and per unit liquid mass [m³/s],
3. the impact scale `log κ̄`, the mean of `log κ` over the collected mass, see [`impact_scale`](@ref),
4. the rime volume rate per unit ice number and per unit liquid mass at each temperature of
   `T°Cs` [m⁶/(kg s)].

The integrals are those of [`∫liquid_ice_collisions`](@ref) with the rule `generator`, with the
rates of one ice particle of [`particle_collection`](@ref).

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `vel`: [`CMP.Chen2022VelType`](@ref)
- `pdf`: the size distribution of the liquid species
- `generator`: `QuadratureRule`
- `logλ`, `w`, `ρ_rim`: the ice state, see [`table_state`](@ref)
- `logρₐ`: log air density [log(kg/m³)]
- `logx̄`: log mean mass of the liquid species [log(kg)], see [`liquid_table_scales`](@ref)
- `T°Cs`: `Tuple` of temperatures [°C]
"""
function collection_integrals(params, vel, pdf, generator, logλ, w, ρ_rim, logρₐ, logx̄, T°Cs)
    (; particles, dist) = table_state(params, generator, logλ, w, ρ_rim)
    liquid = table_liquid(vel, particles, pdf, exp(logρₐ), logx̄)
    rates_at = particle_collection(params, particles, liquid, T°Cs, generator)
    n_i = DT.size_distribution(dist)
    p = liquid_ice_collision_quantile(eltype(particles))
    rates = integrate(Dᵢ -> n_i(Dᵢ) * rates_at(Dᵢ), velocity_integral_bounds(particles, dist, liquid.v_i; p), generator)
    return SA.SVector(
        rates[1] / liquid.number,
        rates[2] / liquid.mass,
        rates[3] / rates[2],
        ntuple(k -> rates[k + 3] / liquid.mass, Val(length(T°Cs)))...,
    )
end

"""
    particle_collection(params, particles, liquid, T°Cs, generator)

Return the function of the ice diameter `Dᵢ` [m] that gives the rates at which one ice particle collects one liquid species.

The rates are an `SVector` of the collected number rate [1/s], the collected mass rate [kg/s],
the collected mass rate times `log κ` (see [`impact_scale`](@ref)) [kg/s], and the rime volume
rate at each temperature of `T°Cs` [m³/s], each integrated over the liquid size distribution
with the rule `generator`, split at the fall-speed crossing (see
[`crossing_integral_bounds`](@ref)).

# Arguments
- `params`: [`CMP.ParametersP3`](@ref)
- `particles`: [`P3Particles`](@ref)
- `liquid`: the parts of the integrals of [`table_liquid`](@ref)
- `T°Cs`: `Tuple` of temperatures [°C]
- `generator`: `QuadratureRule`
"""
function particle_collection(params, particles, liquid, T°Cs, generator)
    (; v_i, v_l, n_l, bounds_l, m_liq) = liquid
    K = Val(length(T°Cs))
    ρ′s = ntuple(k -> CollisionRimeDensity(params.ρ_rim_local, T°Cs[k], nothing, nothing), K)
    return function (Dᵢ)
        v_i_at_Dᵢ = v_i(Dᵢ)
        coeffs = collision_cross_section_ice_liquid_coeffs(particles, Dᵢ)
        inner = D -> begin
            Δv = abs(v_i_at_Dᵢ - v_l(D))
            ∂ₜN = evalpoly(D, coeffs) * Δv * n_l(D)
            ∂ₜM = ∂ₜN * m_liq(D)
            log_κ = log(max(impact_scale(D, Δv), floatmin(Δv)))
            B = ntuple(k -> ∂ₜM / rime_density_at(ρ′s[k], Δv, D), K)
            return SA.SVector(∂ₜN, ∂ₜM, ∂ₜM * log_κ, B...)
        end
        return integrate(inner, crossing_integral_bounds(bounds_l, v_l, v_i_at_Dᵢ), generator)
    end
end

"""
    impact_scale(D, Δv)

Return the product `κ = r Δv` of the radius `r = D / 2` [μm] of a drop of diameter `D` [m] and its impact speed `Δv` [m/s], see [`rime_density_at`](@ref).
"""
@inline impact_scale(D, Δv) = D * 1_000_000 * Δv / 2

"""
    table_liquid(vel, particles, pdf, ρₐ, logx̄)

Return the parts of the liquid-ice collision integrals with one liquid species of mean mass `exp(logx̄)`.

The result is a `NamedTuple` with the ice and liquid terminal velocities `v_i` and `v_l`, the
liquid size distribution `n_l`, its integration bounds `bounds_l`, the liquid particle mass
`m_liq`, and the `mass` and `number` of the distribution, see [`liquid_table_scales`](@ref).

# Arguments
- `vel`: [`CMP.Chen2022VelType`](@ref)
- `particles`: [`P3Particles`](@ref)
- `pdf`: the size distribution of the liquid species
- `ρₐ`: air density [kg/m³]
- `logx̄`: log mean mass of the liquid species [log(kg)]
"""
function table_liquid(vel, particles, pdf, ρₐ, logx̄)
    (L, N) = reference_liquid(pdf, exp(logx̄))
    q = L / ρₐ
    scales = liquid_table_scales(pdf, L, N, ρₐ)
    p = liquid_ice_collision_quantile(eltype(particles))
    return (;
        v_i = ice_particle_terminal_velocity(vel, particles, ρₐ),
        v_l = CO.particle_terminal_velocity(vel.rain, ρₐ),
        n_l = DT.size_distribution(pdf, q, ρₐ, N),
        bounds_l = CM2.get_size_distribution_bounds(pdf, q, ρₐ, N, p),
        m_liq = D -> pdf.ρw * CO.volume_sphere_D(D),
        scales.mass,
        scales.number,
    )
end

function reference_liquid(::CMP.CloudParticlePDF_SB2006, x̄)
    N = oftype(x̄, 1e8)
    return (N * x̄, N)
end
function reference_liquid((; ρw, N0_min, N0_max)::CMP.RainParticlePDF_SB2006_limited, x̄)
    Dr_mean = cbrt(x̄ / (π * ρw))
    N₀ = sqrt(N0_min * N0_max)
    return (π * ρw * N₀ * Dr_mean^4, N₀ * Dr_mean)
end

"""
    liquid_table_scales(pdf, L, N, ρₐ)

Return the scales of the liquid-ice collision tables for a liquid species of content `L` and number `N`.

A liquid size distribution is its number times a shape that depends only on its mean mass. The
result is a `NamedTuple` with
- `present`: whether the distribution is nonzero, from primal values,
- `logx̄`: the log of its mean mass [log(kg)], the table coordinate,
- `mass`: its mass [kg/m³], which scales the mass rates,
- `number`: its number [1/m³], which scales the number rates.

For the limited rain distribution, these are the values of the distribution with its clamped
intercept and slope.

# Arguments
- `pdf`: [`CMP.CloudParticlePDF_SB2006`](@ref) or [`CMP.RainParticlePDF_SB2006_limited`](@ref)
- `L`: liquid content [kg/m³]
- `N`: number concentration [1/m³]
- `ρₐ`: air density [kg/m³]
"""
@inline function liquid_table_scales(::CMP.CloudParticlePDF_SB2006, L, N, ρₐ)
    L_p, N_p = primal_value(L), primal_value(N)
    ϵ = UT.ϵ_numerics(typeof(L_p))
    present = (N_p ≥ ϵ) & (L_p / primal_value(ρₐ) ≥ ϵ)
    logx̄ = log(max(L, floatmin(L_p))) - log(max(N, floatmin(N_p)))
    return (; present, logx̄, mass = ifelse(present, L, zero(L)), number = ifelse(present, N, zero(N)))
end
@inline function liquid_table_scales(pdf::CMP.RainParticlePDF_SB2006_limited, L, N, ρₐ)
    (; N₀r, Dr_mean) = CM2.pdf_rain_parameters(pdf, L / ρₐ, ρₐ, N)
    present = primal_value(N₀r) > 0
    # The number and mass of the exponential distribution `N₀r exp(-D / Dr_mean)`
    number = N₀r * Dr_mean
    mass = π * pdf.ρw * number * Dr_mean^3
    logx̄ = log(π * pdf.ρw) + 3 * log(max(Dr_mean, floatmin(primal_value(Dr_mean))))
    return (; present, logx̄, mass = ifelse(present, mass, zero(mass)), number = ifelse(present, number, zero(number)))
end

"""
    table_log_temperature(T, T_freeze)

Return `log|T - T_freeze|` [log(°C)] of the temperature `T` [K], at least the log of `floatmin`.

The rime density tables are read at `logR = log κ̄ - log|T - T_freeze|`, see
[`tabulated_collection`](@ref).
"""
@inline function table_log_temperature(T, T_freeze)
    T°C = T - T_freeze
    return log(max(-T°C, floatmin(primal_value(T°C))))
end

"""
    tabulated_collection(collision, rime_density, key, pdf, L, N, ρₐ, logT)

Return the collected rates of one liquid species per unit ice number from its liquid-ice tables, without the partition.

The result is a `NamedTuple` with the mass rate `M` [kg/s], the number rate `N` [1/s], the mean
local rime density `ρ′` [kg/m³], and the `mass` and `logx̄` of [`liquid_table_scales`](@ref).

# Arguments
- `collision`: the collision table of the species, see [`collision_point`](@ref)
- `rime_density`: the rime density table of the species, see [`rime_density_point`](@ref)
- `key`: the ice coordinates of [`ice_table_key`](@ref)
- `pdf`: the size distribution of the species
- `L`: liquid content [kg/m³]
- `N`: number concentration [1/m³]
- `ρₐ`: air density [kg/m³]
- `logT`: the temperature coordinate of [`table_log_temperature`](@ref)
"""
@inline function tabulated_collection(collision, rime_density, key, pdf, L, N, ρₐ, logT)
    (; present, logx̄, mass, number) = liquid_table_scales(pdf, L, N, ρₐ)
    (log_number, log_mass, log_impact) = lookup_or(present, 0, collision, key..., logx̄)
    return (;
        M = exp(log_mass) * mass,
        N = exp(log_number) * number,
        ρ′ = lookup_or(present, 1, rime_density, log_impact - logT, key..., logx̄)[1],
        mass,
        logx̄,
        present,
    )
end

"""
    ∫liquid_ice_collisions(
        psd_c, psd_r, aps, tps, vel, m_liq, state::TabulatedP3State,
        L_c, N_c, L_r, N_r, ρₐ, T,
    )

Return the liquid-ice collision rates of [`∫liquid_ice_collisions`](@ref) from the tables, as an `SVector` in the same order.

The collected rates come from [`tabulated_collection`](@ref), and the freeze and shed partition
from [`tabulated_wet_set`](@ref). Above the freezing temperature every collected particle sheds,
and below `T_freeze - L_f / cp_l` every one freezes.
"""
@inline function ∫liquid_ice_collisions(
    psd_c, psd_r, aps, tps, vel, m_liq, state::TabulatedP3State,
    L_c, N_c, L_r, N_r, ρₐ, T,
)
    (; particles, dist) = state
    tables = state.integration
    (; groups) = tables
    key = ice_table_key(state, ρₐ)
    (; heat, denom, T_frz) = max_freeze_rate_scalars(aps, tps, ρₐ, T)
    logT = table_log_temperature(T, particles.params.T_freeze)
    cloud = tabulated_collection(groups.cloud_collision, groups.cloud_rime_density, key, psd_c, L_c, N_c, ρₐ, logT)
    rain = tabulated_collection(groups.rain_collision, groups.rain_rime_density, key, psd_r, L_r, N_r, ρₐ, logT)
    (a, b) = CO.ventilation_coefficients(particles.params.vent, aps)
    A = 2 * (π * heat) / denom
    ρn, M_c, M_r, N_c_col, N_r_col, B_c, B_r, c₁, c₂ = promote(
        dist.ρn, cloud.M, rain.M, cloud.N, rain.N, cloud.M / cloud.ρ′, rain.M / rain.ρ′, A * a, A * b,
    )
    R = typeof(ρn)

    # Above freezing everything sheds, and without a limit to the freezing nothing does
    warm = primal_value(T) ≥ primal_value(T_frz)
    partial = !warm & (primal_value(denom) > 0) & (primal_value(M_c + M_r) > 0)
    wet =
        partial ?
        tabulated_wet_set(tables, state, key, (; cloud..., M = M_c), (; rain..., M = M_r), (c₁, c₂), logT) :
        (;
            cloud_frozen = M_c,
            rain_frozen = M_r,
            cloud_shed = zero(R),
            rain_shed = zero(R),
            cloud_rime = B_c,
            rain_rime = B_r,
        )
    QCFRZ = ifelse(warm, zero(R), ρn * wet.cloud_frozen)
    QRFRZ = ifelse(warm, zero(R), ρn * wet.rain_frozen)
    QCSHD = ifelse(warm, ρn * M_c, ρn * wet.cloud_shed)
    QRSHD = ifelse(warm, ρn * M_r, ρn * wet.rain_shed)
    return SA.SVector(
        QCFRZ, QCSHD, ρn * N_c_col,
        QRFRZ, QRSHD, ρn * N_r_col,
        ρn * M_c + ρn * M_r, ifelse(warm, zero(R), ρn * wet.cloud_rime), ifelse(warm, zero(R), ρn * wet.rain_rime),
    )
end
