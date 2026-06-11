#####
##### 2M+P3 Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState2MP3{FT}

The eight prognostic 2M+P3 species as a `StaticArrays.FieldVector`: vector
algebra for substep updates and linear solves, named fields for readability.
Internal to the [`RosenbrockAverage`](@ref) implementation.
"""
struct MicroState2MP3{FT} <: SA.FieldVector{8, FT}
    q_lcl::FT
    n_lcl::FT
    q_rai::FT
    n_rai::FT
    q_ice::FT
    n_ice::FT
    q_rim::FT
    b_rim::FT
end
SA.similar_type(::Type{<:MicroState2MP3}, ::Type{FT}, ::SA.Size{(8,)}) where {FT} =
    MicroState2MP3{FT}

"""
    _sink_limit(q, dt, n = 3)

Largest admissible sink rate for `q`: deplete at most `1/n` of it over `dt`.
"""
@inline _sink_limit(q, dt, n = 3) = max(zero(q), q) / dt / n

"""
    _coupled_sink_factor(S1, S2, q1, q2, dt)

Common scaling factor in `(0, 1]` for a coupled (mass, number) sink pair
`(S1, S2)`, such that neither `f S1` nor `f S2` exceeds its
[`_sink_limit`](@ref) over `dt`. Scaling the pair by one common factor
preserves the mean particle mass implied by the pair; limiting the members
independently would not.
"""
@inline function _coupled_sink_factor(S1, S2, q1, q2, dt)
    M1 = _sink_limit(q1, dt)
    M2 = _sink_limit(q2, dt)
    f1 = ifelse((S1 < zero(S1)) & (-S1 > M1), M1 / -S1, one(S1))
    f2 = ifelse((S2 < zero(S2)) & (-S2 > M2), M2 / -S2, one(S2))
    return min(f1, f2)
end

"""
    _cap_supersaturation(tend, tps, T, ρ, q_tot, q_lcl, q_rai, q_ice, dt)

Cap the net condensation/deposition tendencies in `tend` so that one substep of
length `dt` cannot overshoot saturation: the liquid channel
(`dq_lcl + dq_rai`) is capped against the analytic saturation-adjustment
condensation increment, then the ice channel (`dq_ice + dq_rim`) against the
deposition increment computed from the vapor remaining after the liquid
step. Each channel is scaled by a common, sign-preserving ratio in `[0, 1]`
(both members of a channel are scaled together, like
[`_coupled_sink_factor`](@ref)). Branch-free: `ifelse` evaluates both arms,
which is safe here and keeps the kernel GPU-friendly.
"""
@inline function _cap_supersaturation(tend, tps, T, ρ, q_tot, q_lcl, q_rai, q_ice, dt)
    (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt,
        dq_ice_dt, dn_ice_dt, dq_rim_dt, db_rim_dt) = tend
    FT = typeof(q_tot)
    T_safe = max(150, T)
    q_vap = max(0, q_tot - q_lcl - q_rai - q_ice)

    Rᵥ = TDI.Rᵥ(tps)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    L_v = TDI.Lᵥ(tps, T_safe)
    L_s = TDI.Lₛ(tps, T_safe)
    qv_sat_liq = TDI.saturation_vapor_specific_content_over_liquid(tps, T_safe, ρ)
    qv_sat_ice = TDI.saturation_vapor_specific_content_over_ice(tps, T_safe, ρ)

    # liquid channel vs the saturation-adjustment condensation increment
    qcon_satadj =
        (q_vap - qv_sat_liq) / (1 + L_v^2 * qv_sat_liq / (cp_d * Rᵥ * T_safe^2)) / dt
    net_liq = dq_lcl_dt + dq_rai_dt
    target_liq = ifelse(
        net_liq > 0, min(net_liq, max(zero(qcon_satadj), qcon_satadj)),
        ifelse(net_liq < 0, max(net_liq, min(zero(qcon_satadj), qcon_satadj)), zero(net_liq)),
    )
    ratio_liq = clamp(
        ifelse(abs(net_liq) > eps(FT), target_liq / net_liq, one(net_liq)),
        zero(net_liq), one(net_liq),
    )
    dq_lcl_dt *= ratio_liq
    dq_rai_dt *= ratio_liq
    dn_lcl_dt *= ratio_liq
    dn_rai_dt *= ratio_liq

    # ice channel vs deposition of the vapor remaining after the liquid step
    qv_after_liq = max(0, q_vap - target_liq * dt)
    qdep_satadj =
        (qv_after_liq - qv_sat_ice) /
        (1 + L_s^2 * qv_sat_ice / (cp_d * Rᵥ * T_safe^2)) / dt
    net_ice = dq_ice_dt + dq_rim_dt
    target_ice = ifelse(
        net_ice > 0, min(net_ice, max(zero(qdep_satadj), qdep_satadj)),
        ifelse(net_ice < 0, max(net_ice, min(zero(qdep_satadj), qdep_satadj)), zero(net_ice)),
    )
    ratio_ice = clamp(
        ifelse(abs(net_ice) > eps(FT), target_ice / net_ice, one(net_ice)),
        zero(net_ice), one(net_ice),
    )
    dq_ice_dt *= ratio_ice
    dq_rim_dt *= ratio_ice
    dn_ice_dt *= ratio_ice
    db_rim_dt *= ratio_ice

    return (; dq_lcl_dt, dn_lcl_dt, dq_rai_dt, dn_rai_dt,
        dq_ice_dt, dn_ice_dt, dq_rim_dt, db_rim_dt)
end

"""
    _limited_2mp3_tendency(mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ, dt)

The instantaneous 2M+P3 tendency projected onto the eight prognostic
species, with [`_cap_supersaturation`](@ref) and the coupled
[`_coupled_sink_factor`](@ref) mass-number limits applied for a substep of
length `dt`.

This is the function the Rosenbrock step differentiates: the limiting sits
inside the Jacobian so the implicit update respects it. Interrogated
empirically: at substep lengths up to tens of seconds the supersaturation
cap improves accuracy substantially (the coupled-sink limits are neutral
there and protect the explicit fallback paths); at very large substeps
(h ≳ 100 s) the limiter switching degrades the single linearization —
increase `nsub` rather than relying on the limiters.
"""
@inline function _limited_2mp3_tendency(mp, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    dt,
)
    full = bulk_microphysics_tendencies(Microphysics2Moment(), mp, tps,
        ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    )
    tend = (; full.dq_lcl_dt, full.dn_lcl_dt, full.dq_rai_dt, full.dn_rai_dt,
        full.dq_ice_dt, full.dn_ice_dt, full.dq_rim_dt, full.db_rim_dt)
    tend = _cap_supersaturation(tend, tps, T, ρ, q_tot, q_lcl, q_rai, q_ice, dt)
    f_liq = _coupled_sink_factor(tend.dq_lcl_dt, tend.dn_lcl_dt, q_lcl, n_lcl, dt)
    f_rai = _coupled_sink_factor(tend.dq_rai_dt, tend.dn_rai_dt, q_rai, n_rai, dt)
    return (;
        dq_lcl_dt = tend.dq_lcl_dt * f_liq, dn_lcl_dt = tend.dn_lcl_dt * f_liq,
        dq_rai_dt = tend.dq_rai_dt * f_rai, dn_rai_dt = tend.dn_rai_dt * f_rai,
        tend.dq_ice_dt, tend.dn_ice_dt, tend.dq_rim_dt, tend.db_rim_dt,
    )
end

"""
    Limited2MP3Tendency(mp, tps, ρ, T, q_tot, logλ, dt)

Callable bundling the frozen per-substep context; applying it to the species
vector evaluates [`_limited_2mp3_tendency`](@ref). A top-level struct rather
than a closure, so `ForwardDiff` differentiates a concretely-typed callable.

The frozen `q_tot` is promoted to the state's element type at the call: a
zero-partial Dual is exact for a constant, and it keeps
`eltype(q_tot)`-keyed sentinel returns in the tendency functions concretely
typed (a no-op in the primal pass). `logλ` stays plain so distribution-shape
slots never carry derivatives (the forward `gamma_inc` has no shape rule);
`T` and `ρ` stay plain because they feed working-type computations that must
remain floats.
"""
struct Limited2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
    logλ::F
    dt::F
end
@inline function (g::Limited2MP3Tendency)(x::SA.StaticVector{8})
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    tend = _limited_2mp3_tendency(g.mp, g.tps,
        g.ρ, g.T, eltype(x)(g.q_tot),
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, g.logλ,
        g.dt,
    )
    return MicroState2MP3(values(tend)...)
end

"""
    _rosenbrock_channel_mask(x)

Diagonal of the channel projection matrix `P` used by
[`_rosenbrock_update`](@ref): 1 for healthy channels, 0 for near-empty ones
(condensed mass below `1e-10`, per channel: liquid, rain, ice+rime).

This is a scheme partition, not a preconditioner: with `P J P` in place of
`J`, a masked channel's row of `(I/h - P J P)` is the identity row, so the
solve returns exactly the forward-Euler update for that channel while
healthy channels stay implicit (an IMEX-style splitting at channel
granularity). The gate is load-bearing: through the coupled-sink limiter,
near-empty channels produce finite but enormous Jacobian entries whose
linearized steady state fabricates phantom number concentrations that
substep refinement cannot remove.
"""
@inline function _rosenbrock_channel_mask(x::MicroState2MP3{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    liq = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    ice = ifelse(x.q_ice < ϵ_empty, zero(FT), one(FT))
    return MicroState2MP3(liq, liq, rai, rai, ice, ice, ice, ice)
end

"""
    _euler_update(x, f, h)

Forward-Euler substep, floored at zero.
"""
@inline _euler_update(x, f, h) = max.(x .+ h .* f, 0)

"""
    _rosenbrock_update(x, f, J, h)

One linearized-implicit (Rosenbrock-Euler) substep: solve

    (I/h - P J P) Δx = f

and return `max.(x + Δx, 0)`, where `P = Diagonal(z)` is the channel
projection of [`_rosenbrock_channel_mask`](@ref). The system is solved in
equilibrated variables: with `S = Diagonal(|x| + h |f| + ϵ)` the similarity
transform `S⁻¹ A S` is O(1)-conditioned — the raw rows span ~9 orders of
magnitude (number vs mass species), and an unscaled Float32 factorization
bleeds roundoff from the large rows into empty species as phantom mass.
Equilibration is exact in exact arithmetic and keeps roundoff relative to
each species' own scale.
"""
@inline function _rosenbrock_update(x::MicroState2MP3{FT}, f, J, h) where {FT}
    I₈ = one(SA.SMatrix{8, 8, FT})
    z = _rosenbrock_channel_mask(x)
    s = abs.(x) .+ h .* abs.(f) .+ eps(FT)
    # dense diagonal matrices: an `SDiagonal` wrapper here defeats the
    # optimizer's static-array stack allocation (heap spills per substep)
    P = I₈ .* z'
    S = I₈ .* s'
    S⁻¹ = I₈ .* inv.(s)'
    A = I₈ / h - S⁻¹ * (P * J * P) * S
    Δx = S * (A \ (S⁻¹ * f))
    return max.(x .+ Δx, 0)
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Compute average 2M+P3 microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Evaluate the limited tendency `f` ([`_limited_2mp3_tendency`](@ref)) and
   its exact 8×8 Jacobian `J` via `ForwardDiff` at the current state.
2. Advance with [`_rosenbrock_update`](@ref): solve `(I/h - P J P) Δx = f`
   in equilibrated variables, where the projection `P` routes near-empty
   channels to forward Euler.
3. Update the local temperature from the latent heating of the realized
   increments.

A non-finite state cannot be differentiated and a non-finite Jacobian (an
exotic state escaping the channel mask) cannot be linearized; both fall back
to a forward-Euler substep of the limited tendency, which is safe to step
explicitly. `logλ` and `q_tot` are held fixed across substeps, matching the
explicit-substepping semantics; the limiting (supersaturation cap and
coupled mass-number sinks) applies the explicit-timestepping semantics
unconditionally and sits inside the differentiated function.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry (without the activation diagnostic).
"""
@inline function bulk_microphysics_tendencies(::RosenbrockAverage, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)

    x = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Limited2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ, h)
        f = g(x)
        x_prev = x
        x = if all(isfinite, x)
            J = FD.jacobian(g, x)
            all(isfinite, J) ? _rosenbrock_update(x, f, J, h) : _euler_update(x, f, h)
        else
            _euler_update(x, f, h)
        end
        Δ = x - x_prev
        T_safe = max(150, Tsub)
        Tsub += (TDI.Lᵥ(tps, T_safe) * (Δ.q_lcl + Δ.q_rai) + TDI.Lₛ(tps, T_safe) * Δ.q_ice) / cp_d
    end

    rates = (x - x₀) / Δt
    return NamedTuple{(
        :dq_lcl_dt, :dn_lcl_dt, :dq_rai_dt, :dn_rai_dt,
        :dq_ice_dt, :dn_ice_dt, :dq_rim_dt, :db_rim_dt,
    )}(
        Tuple(rates),
    )
end
