#####
##### 2M+P3 Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState2MP3{FT}

The eight prognostic 2M+P3 species as a `StaticArrays.FieldVector`. Internal to
the [`RosenbrockAverage`](@ref) implementation.
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
vector evaluates [`_instantaneous_2mp3_tendency`](@ref). `q_tot` is promoted to
the state's element type at the call; `logλ`, `T`, and `ρ` stay plain.
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
    _rosenbrock_species_mask(x)

Diagonal of the species projection matrix `P` used by
[`_rosenbrock_update`](@ref): 1 for active species, 0 for near-empty ones
(condensed mass below `1e-10`, per species: liquid, rain, ice+rime). A masked
species takes the forward-Euler update while active species stay implicit.
"""
@inline function _rosenbrock_species_mask(x::MicroState2MP3{FT}) where {FT}
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
    _rosenbrock_system(x, f, J, z, h)

Build the equilibrated linear system of one linearized-implicit
(Rosenbrock-Euler) substep at state `x` with raw tendency `f`, Jacobian `J`,
species mask `z`, and substep `h`. Returns `(S, S⁻¹, A)`, the equilibration
matrix `S = Diagonal(|x| + h |f| + ϵ)`, its inverse, and the equilibrated
system matrix `A = I/h - S⁻¹ B S`, where `B = P J P` is the masked Jacobian and
`P = Diagonal(z)` is the species projection built from the per-scheme species
mask `z` (e.g. [`_rosenbrock_species_mask`](@ref) for 2M+P3).
"""
@inline function _rosenbrock_system(
    x::SA.StaticVector{N, FT}, f, J, z, h,
) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    s = abs.(x) .+ h .* abs.(f) .+ eps(FT)
    P = Iₙ .* z'
    S = Iₙ .* s'
    S⁻¹ = Iₙ .* inv.(s)'
    B = P * J * P
    A = Iₙ / h - S⁻¹ * B * S
    return S, S⁻¹, A
end

"""
    _rosenbrock_solve(S, S⁻¹, A, v)

Solve the equilibrated Rosenbrock system from [`_rosenbrock_system`](@ref) for
the unclamped increment of right-hand side `v`: `Δ = S (A \\ (S⁻¹ v))`, the
equilibrated form of `(I/h - P J P)⁻¹ v`.
"""
@inline _rosenbrock_solve(S, S⁻¹, A, v) = S * (A \ (S⁻¹ * v))

"""
    _rosenbrock_update(x, f, J, z, h)

One linearized-implicit (Rosenbrock-Euler) substep: build the equilibrated
system ([`_rosenbrock_system`](@ref)) for `(I/h - B) Δx = f` with the masked
Jacobian `B = P J P`, solve it ([`_rosenbrock_solve`](@ref)), and return
`max.(x + Δx, 0)`.
"""
@inline function _rosenbrock_update(
    x::SA.StaticVector{N, FT}, f, J, z, h,
) where {N, FT}
    S, S⁻¹, A = _rosenbrock_system(x, f, J, z, h)
    Δx = _rosenbrock_solve(S, S⁻¹, A, f)
    return max.(x .+ Δx, 0)
end

# ---- Framework option resolvers shared by the substep driver ----

"""
    _full_species_mask(x)

The all-ones species projection: every species stays in the implicit solve.
"""
@inline _full_species_mask(x::SA.StaticVector{N, FT}) where {N, FT} =
    ones(SA.SVector{N, FT})

"""
    _species_mask(jacobian, growth)

The species projection `x -> z` for a [`Jacobian`](@ref) and
[`GrowthTreatment`](@ref) pair. An explicit growth diagonal removes the
unbounded growth of the exact Jacobian, so its species stay implicit
([`_full_species_mask`](@ref)); the exact Jacobian with implicit growth uses the
near-empty mask ([`_rosenbrock_species_mask`](@ref)).
"""
@inline _species_mask(::ExactJacobian, ::ExplicitGrowthDiagonal) = _full_species_mask
@inline _species_mask(::ExactJacobian, ::ImplicitGrowth) = _rosenbrock_species_mask

"""
    _apply_growth(growth, J)

Apply a [`GrowthTreatment`](@ref) to the Jacobian `J`. [`ImplicitGrowth`](@ref)
returns `J` unchanged; [`ExplicitGrowthDiagonal`](@ref) removes the positive
diagonal entries, leaving the off-diagonals and the negative diagonals.
"""
@inline _apply_growth(::ImplicitGrowth, J) = J
@inline function _apply_growth(::ExplicitGrowthDiagonal, J::SA.SMatrix{N, N, FT}) where {N, FT}
    Iₙ = one(SA.SMatrix{N, N, FT})
    return J - max.(Iₙ .* J, zero(FT))
end

"""
    _saturation_bisection_count(FT)

Number of bisection iterations to resolve a fraction in `[0, 1]` to the
precision of `FT`.
"""
@inline _saturation_bisection_count(::Type{FT}) where {FT} = ceil(Int, -log2(eps(FT)))

"""
    _apply_limiter(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)

Limit the substep increment `d` at state `x`. [`NoLimiter`](@ref) returns `d`.
[`EndStateSaturationAdjustment`](@ref), for a cell at or above ice saturation
whose full-increment end state would drop below it, scales `d` by `s ∈ [0, 1]`
to keep the latent-heated end state at or above ice saturation; otherwise it
returns `d`.
"""
@inline _apply_limiter(::NoLimiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps) = d

@inline function _apply_limiter(::EndStateSaturationAdjustment,
    x::MicroState2MP3{FT}, d::MicroState2MP3{FT},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {FT}
    Sice(xx, TT) = TDI.supersaturation_over_ice(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_ice, ρ, TT)
    latent(dd) = Lv_over_cp * (dd.q_lcl + dd.q_rai) + Ls_over_cp * dd.q_ice
    xf = max.(x .+ d, 0)
    if Sice(x, Tsub) >= 0 && Sice(xf, Tsub + latent(xf .- x)) < 0
        lo = zero(FT)
        hi = one(FT)
        for _ in 1:_saturation_bisection_count(FT)
            s = (lo + hi) / 2
            xs = max.(x .+ s .* d, 0)
            if Sice(xs, Tsub + latent(xs .- x)) >= 0
                lo = s
            else
                hi = s
            end
        end
        return lo .* d
    end
    return d
end

bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment, args...) = throw(
    ArgumentError(
        "RosenbrockAverage on the 2M+P3 model supports only ExactJacobian; use rosenbrock_exact()",
    ),
)

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
   species to forward Euler.
3. Update the local temperature from the latent heating of the realized
   increments.

A non-finite state or Jacobian falls back to a forward-Euler substep of the
raw tendency. `logλ` and `q_tot` are held fixed across substeps.

The 2M+P3 model supports only [`ExactJacobian`](@ref).

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry (without the activation diagnostic).
"""
@inline function bulk_microphysics_tendencies(mode::RosenbrockAverage{ExactJacobian}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, ICE}, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1,
) where {WR, ICE <: CMP.P3IceParams}
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    cp_d = TDI.TD.Parameters.cp_d(tps)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / cp_d
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / cp_d

    x = MicroState2MP3{FT}(q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Instantaneous2MP3Tendency(mp, tps, ρ, Tsub, q_tot, logλ)
        f = g(x)
        x_prev = x
        if all(isfinite, x)
            J = _apply_growth(mode.growth, FD.jacobian(g, x))
            z = _species_mask(mode.jacobian, mode.growth)(x)
            d = if all(isfinite, J)
                _rosenbrock_update(x, f, J, z, h) - x
            else
                _euler_update(x, f, h) - x
            end
            d = _apply_limiter(mode.limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)
            x = max.(x .+ d, 0)
        else
            x = _euler_update(x, f, h)
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
