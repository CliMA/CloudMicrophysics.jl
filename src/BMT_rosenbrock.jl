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
    _instantaneous_2mp3_tendency(mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ)

The raw instantaneous 2M+P3 tendency projected onto the eight prognostic
species: the unlimited process rates of the `Microphysics2Moment` entry,
without timestep-dependent clipping.

This is the function the Rosenbrock step differentiates. It is the model
physics — an analytic time-averaged relaxation (Morrison and Milbrandt
2015, J. Atmos. Sci., 72, 287-311, Appendix C) carries no tendency clip;
the L-stable one-stage Rosenbrock update damps the stiff vapor-exchange
subsystem monotonically, so the overshoot an explicit-Euler supersaturation
cap suppresses cannot arise, while a 1/dt clip would add dt-independent
error and break convergence under refinement (Wan et al. 2020, J. Adv.
Model. Earth Syst., 12, e2019MS001982).
"""
@inline function _instantaneous_2mp3_tendency(mp, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
)
    full = bulk_microphysics_tendencies(Microphysics2Moment(), mp, tps,
        ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    )
    return (; full.dq_lcl_dt, full.dn_lcl_dt, full.dq_rai_dt, full.dn_rai_dt,
        full.dq_ice_dt, full.dn_ice_dt, full.dq_rim_dt, full.db_rim_dt)
end

"""
    Instantaneous2MP3Tendency(mp, tps, ρ, T, q_tot, logλ)

Callable bundling the frozen per-substep context; applying it to the species
vector evaluates [`_instantaneous_2mp3_tendency`](@ref). A top-level struct
rather than a closure, so `ForwardDiff` differentiates a concretely-typed
callable.

The frozen `q_tot` is promoted to the state's element type at the call: a
zero-partial Dual is exact for a constant, and it keeps
`eltype(q_tot)`-keyed fallback returns in the tendency functions concretely
typed (a no-op in the primal pass). `logλ` stays plain so distribution-shape
slots never carry derivatives (the forward `gamma_inc` has no shape rule);
`T` and `ρ` stay plain because they feed working-type computations that must
remain floats.
"""
struct Instantaneous2MP3Tendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
    logλ::F
end
@inline function (g::Instantaneous2MP3Tendency)(x::SA.StaticVector{8})
    (q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim) = x
    tend = _instantaneous_2mp3_tendency(g.mp, g.tps,
        g.ρ, g.T, eltype(x)(g.q_tot),
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, g.logλ,
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
granularity). The gate is load-bearing: near-empty channels produce finite
but enormous Jacobian entries whose linearized steady state fabricates
phantom number concentrations that substep refinement cannot remove.
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
    _rosenbrock_update(x, f, J, z, h)

One linearized-implicit (Rosenbrock-Euler) substep: solve

    (I/h - P J P) Δx = f

and return `max.(x + Δx, 0)`, where `P = Diagonal(z)` is the channel
projection built from the per-scheme channel mask `z` (e.g.
[`_rosenbrock_channel_mask`](@ref) for 2M+P3). The system is solved in
equilibrated variables: with `S = Diagonal(|x| + h |f| + ϵ)` the similarity
transform `S⁻¹ A S` is O(1)-conditioned — the raw rows span ~9 orders of
magnitude (number vs mass species), and an unscaled Float32 factorization
bleeds roundoff from the large rows into empty species as phantom mass.
Equilibration is exact in exact arithmetic and keeps roundoff relative to
each species' own scale.

Dimension-generic over any `StaticArrays.StaticVector{N}` state: the
identity and dense diagonal matrices are sized from the static length `N`,
so any prognostic vector and its matching channel mask plug in unchanged.
"""
@inline function _rosenbrock_update(x::SA.StaticVector{N, FT}, f, J, z, h) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    s = abs.(x) .+ h .* abs.(f) .+ eps(FT)
    # dense diagonal matrices: an `SDiagonal` wrapper here defeats the
    # optimizer's static-array stack allocation (heap spills per substep)
    P = Iₙ .* z'
    S = Iₙ .* s'
    S⁻¹ = Iₙ .* inv.(s)'
    A = Iₙ / h - S⁻¹ * (P * J * P) * S
    Δx = S * (A \ (S⁻¹ * f))
    return max.(x .+ Δx, 0)
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Compute average 2M+P3 microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps of the raw instantaneous
tendency.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Evaluate the raw tendency `f` ([`_instantaneous_2mp3_tendency`](@ref)) and
   its exact 8×8 Jacobian `J` via `ForwardDiff` at the current state.
2. Advance with [`_rosenbrock_update`](@ref): solve `(I/h - P J P) Δx = f`
   in equilibrated variables, where the projection `P` routes near-empty
   channels to forward Euler.
3. Update the local temperature from the latent heating of the realized
   increments.

A non-finite state cannot be differentiated and a non-finite Jacobian (an
exotic state escaping the channel mask) cannot be linearized; both fall back
to a forward-Euler substep of the raw tendency. `logλ` and `q_tot` are held
fixed across substeps, matching the explicit-substepping semantics. The
discrete safeguards are h-free conditioning and projection (channel mask,
equilibration, positivity clamp), not model terms: the differentiated
tendency is the unmodified relaxation physics.

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
        g = Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        f = g(x)
        x_prev = x
        x = if all(isfinite, x)
            J = FD.jacobian(g, x)
            z = _rosenbrock_channel_mask(x)
            all(isfinite, J) ? _rosenbrock_update(x, f, J, z, h) : _euler_update(x, f, h)
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
