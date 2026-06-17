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

The system build is separated from the solve so the full-step update and the
per-process attribution can reuse one factorization, solving against the same
`S`, `S⁻¹`, `A`.
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
equilibrated form of `(I/h - P J P)⁻¹ v`. Linear in `v`, so per-process
increments sum to the full-step increment.
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

bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment, args...) = throw(
    ArgumentError(
        "RosenbrockAverage on the 2M+P3 model supports only ExactJacobian; use rosenbrock_exact()",
    ),
)

bulk_microphysics_tendencies(
    ::RosenbrockAverage, ::Microphysics2Moment,
    mp::CMP.Microphysics2MParams{WR, Nothing}, args...,
) where {WR} = throw(
    ArgumentError(
        "RosenbrockAverage on Microphysics2Moment requires P3 ice parameters (with_ice = true)",
    ),
)

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1)

Compute average 2M+P3 microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps of the raw instantaneous
tendency. See the [Rosenbrock-average microphysics substepping](@ref)
documentation page for the substep algorithm.

`logλ` and `q_tot` are held fixed across substeps. The 2M+P3 model supports
only [`ExactJacobian`](@ref); the donor-based matrix is 1M-only.

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

#####
##### 1M Rosenbrock-Euler substepping (`RosenbrockAverage`)
#####

"""
    MicroState1M{FT}

The four prognostic 1M species as a `StaticArrays.FieldVector`. Internal to the
[`RosenbrockAverage`](@ref) implementation, mirroring [`MicroState2MP3`](@ref)
at dimension 4 (1M carries no number species).
"""
struct MicroState1M{FT} <: SA.FieldVector{4, FT}
    q_lcl::FT
    q_icl::FT
    q_rai::FT
    q_sno::FT
end
SA.similar_type(::Type{<:MicroState1M}, ::Type{FT}, ::SA.Size{(4,)}) where {FT} =
    MicroState1M{FT}

"""
    _instantaneous_1m_tendency(mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno)

The raw instantaneous 1M tendency projected onto the four prognostic species:
the unlimited process rates of the `Microphysics1Moment` `Instantaneous` entry
(`_microphysics_source_terms` aggregated by `_aggregate_tendencies`), without
timestep-dependent clipping.
"""
@inline function _instantaneous_1m_tendency(mp, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
)
    src = _microphysics_source_terms(Microphysics1Moment(), mp, tps,
        ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno,
    )
    return _aggregate_tendencies(src)
end

"""
    Raw1MTendency(mp, tps, ρ, T, q_tot)

Callable bundling the frozen per-substep context; applying it to the species
vector evaluates [`_instantaneous_1m_tendency`](@ref), mirroring
[`Instantaneous2MP3Tendency`](@ref). `q_tot` is promoted to the state's element
type at the call; `T` and `ρ` stay plain.
"""
struct Raw1MTendency{P, H, F}
    mp::P
    tps::H
    ρ::F
    T::F
    q_tot::F
end
@inline function (g::Raw1MTendency)(x::SA.StaticVector{4})
    (q_lcl, q_icl, q_rai, q_sno) = x
    tend = _instantaneous_1m_tendency(g.mp, g.tps,
        g.ρ, g.T, eltype(x)(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    return MicroState1M(tend.dq_lcl_dt, tend.dq_icl_dt, tend.dq_rai_dt, tend.dq_sno_dt)
end

"""
    _rosenbrock_species_mask(x::MicroState1M)

Diagonal of the species projection `P` for the 1M state: 1 for active species,
0 for near-empty ones (mass below `1e-10`, per species: liquid, ice, rain,
snow). See the [`MicroState2MP3`](@ref) method for the role of `P` in
[`_rosenbrock_update`](@ref).
"""
@inline function _rosenbrock_species_mask(x::MicroState1M{FT}) where {FT}
    ϵ_empty = FT(1e-10)
    lcl = ifelse(x.q_lcl < ϵ_empty, zero(FT), one(FT))
    icl = ifelse(x.q_icl < ϵ_empty, zero(FT), one(FT))
    rai = ifelse(x.q_rai < ϵ_empty, zero(FT), one(FT))
    sno = ifelse(x.q_sno < ϵ_empty, zero(FT), one(FT))
    return MicroState1M(lcl, icl, rai, sno)
end

"""
    bulk_microphysics_tendencies(::RosenbrockAverage, ::Microphysics1Moment,
        mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1)

Compute average 1M microphysics tendencies over `Δt` using `nsub`
linearized-implicit (Rosenbrock-Euler) substeps of the raw instantaneous
tendency. The [`Jacobian`](@ref), [`GrowthTreatment`](@ref), and
[`TendencyLimiter`](@ref) options of `mode` select the substep matrix, the
growth-diagonal treatment, and the increment limiter.

# Algorithm

For each substep of `h = Δt / nsub`:

1. Build the substep Jacobian `J` from the [`Jacobian`](@ref) provider and apply
   the [`GrowthTreatment`](@ref) to it.
2. Advance with [`_rosenbrock_update`](@ref): solve `(I/h - P J P) Δx = f` in
   equilibrated variables, where the projection `P` is selected by
   [`_species_mask`](@ref).
3. Apply the [`TendencyLimiter`](@ref) to the increment.
4. Update the local temperature from the latent heating of the realized
   increments (constant latent heats, liquid+rain on `L_v`, ice+snow on `L_s`).

A non-finite state or Jacobian falls back to a forward-Euler substep of the raw
tendency; `q_tot` is held fixed across substeps.

Returns the net change in the species over `Δt` divided by `Δt`, in the same
fields as the `Instantaneous` entry.
"""
@inline bulk_microphysics_tendencies(mode::RosenbrockAverage, ::Microphysics1Moment,
    mp::CMP.Microphysics1MParams, tps,
    ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub = 1) =
    _rosenbrock_average_1m(mode, mp, tps, ρ, T, q_tot, q_lcl, q_icl, q_rai, q_sno, Δt, nsub)

# ---- Shared 1M driver + Jacobian/mask/growth/limiter resolvers ----

"""
    _jacobian_provider(jacobian)

The substep Jacobian provider `(g, x) -> SMatrix` for a [`Jacobian`](@ref)
option: [`_jacobian_1m_linearized`](@ref) for [`DonorJacobian`](@ref),
[`_jacobian_1m_relinearized`](@ref) for [`CoupledDonorJacobian`](@ref), and
[`_ad_jacobian_1m`](@ref) for [`ExactJacobian`](@ref).
"""
@inline _jacobian_provider(::DonorJacobian) = _jacobian_1m_linearized
@inline _jacobian_provider(::CoupledDonorJacobian) = _jacobian_1m_relinearized
@inline _jacobian_provider(::ExactJacobian) = _ad_jacobian_1m

"""
    _species_mask(jacobian, growth)

The species projection `x -> z` for a [`Jacobian`](@ref) and
[`GrowthTreatment`](@ref) pair. The donor-based matrices are bounded by their rate
flooring, so every species stays implicit ([`_full_species_mask`](@ref)). An
explicit growth diagonal removes the unbounded growth of the exact Jacobian, so
its species stay implicit too; the exact Jacobian with implicit growth uses the
near-empty mask ([`_rosenbrock_species_mask`](@ref)).
"""
@inline _species_mask(::DonorJacobian, ::GrowthTreatment) = _full_species_mask
@inline _species_mask(::CoupledDonorJacobian, ::GrowthTreatment) = _full_species_mask
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
    _apply_limiter(limiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps)

Limit the substep increment `d` at state `x`. [`NoLimiter`](@ref) returns `d`.
[`EndStateSaturationAdjustment`](@ref), for a cell at or above ice saturation
whose full-increment end state would drop below it, scales `d` by `s ∈ [0, 1]`
to keep the latent-heated end state at or above ice saturation; otherwise it
returns `d`.
"""
@inline _apply_limiter(::NoLimiter, x, d, ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps) = d

"""
    _saturation_bisection_count(FT)

Number of bisection iterations to resolve a fraction in `[0, 1]` to the
precision of `FT`.
"""
@inline _saturation_bisection_count(::Type{FT}) where {FT} = ceil(Int, -log2(eps(FT)))

@inline function _apply_limiter(::EndStateSaturationAdjustment,
    x::MicroState1M{FT}, d::MicroState1M{FT},
    ρ, Tsub, q_tot, Lv_over_cp, Ls_over_cp, tps,
) where {FT}
    Sice(xx, TT) =
        TDI.supersaturation_over_ice(tps, q_tot, xx.q_lcl + xx.q_rai, xx.q_icl + xx.q_sno, ρ, TT)
    latent(dd) = Lv_over_cp * (dd.q_lcl + dd.q_rai) + Ls_over_cp * (dd.q_icl + dd.q_sno)
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

"Exact ForwardDiff Jacobian provider for [`_rosenbrock_average_1m`](@ref)."
@inline _ad_jacobian_1m(g, x) = FD.jacobian(g, x)

"""
    _full_species_mask(x)

The all-ones species projection: every species stays in the implicit
solve.
"""
@inline _full_species_mask(x::SA.StaticVector{N, FT}) where {N, FT} =
    ones(SA.SVector{N, FT})

"""
    _rosenbrock_average_1m(mode, mp, tps, ρ, T, q_tot,
                           q_lcl, q_icl, q_rai, q_sno, Δt, nsub)

Shared 1M Rosenbrock-average substep driver. The substep loop, equilibration,
positivity clamp, and explicit between-substep `T` update are fixed; the
[`RosenbrockAverage`](@ref) `mode` selects the Jacobian, the growth treatment,
and the increment limiter through [`_jacobian_provider`](@ref),
[`_species_mask`](@ref), [`_apply_growth`](@ref), and [`_apply_limiter`](@ref).
"""
@inline function _rosenbrock_average_1m(
    mode::RosenbrockAverage, mp::CMP.Microphysics1MParams, tps, ρ, T, q_tot,
    q_lcl, q_icl, q_rai, q_sno, Δt, nsub,
)
    FT = typeof(q_tot)
    nsub_eff = max(Int(nsub), 1)
    h = Δt / FT(nsub_eff)
    Lv_over_cp = TDI.TD.Parameters.LH_v0(tps) / TDI.TD.Parameters.cp_d(tps)
    Ls_over_cp = TDI.TD.Parameters.LH_s0(tps) / TDI.TD.Parameters.cp_d(tps)
    jacobian = _jacobian_provider(mode.jacobian)
    mask = _species_mask(mode.jacobian, mode.growth)

    x = MicroState1M{FT}(q_lcl, q_icl, q_rai, q_sno)
    x₀ = x
    Tsub = T
    for _ in 1:nsub_eff
        g = Raw1MTendency(mp, tps, ρ, Tsub, q_tot)
        f = g(x)
        x_prev = x
        if all(isfinite, x)
            J = _apply_growth(mode.growth, jacobian(g, x))
            z = mask(x)
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
        Tsub += Lv_over_cp * (Δ.q_lcl + Δ.q_rai) + Ls_over_cp * (Δ.q_icl + Δ.q_sno)
    end

    rates = (x - x₀) / Δt
    return (;
        dq_lcl_dt = rates.q_lcl, dq_icl_dt = rates.q_icl,
        dq_rai_dt = rates.q_rai, dq_sno_dt = rates.q_sno,
    )
end

"""
    _jacobian_1m(FT; lcl_lcl, lcl_icl, ..., sno_sno)

Assemble the 4×4 Jacobian of the 1M tendency over the state
`(q_lcl, q_icl, q_rai, q_sno)` from named species-pair entries: keyword
`<receiver>_<donor>` is `∂(dq_receiver/dt)/∂q_donor`. Entries not supplied are
zero. Centralizes the index layout so a hand-built Jacobian is filled by physical
coupling name rather than numeric `[i, j]` position.
"""
@inline _jacobian_1m(::Type{FT};
    lcl_lcl = zero(FT), lcl_icl = zero(FT), lcl_rai = zero(FT), lcl_sno = zero(FT),
    icl_lcl = zero(FT), icl_icl = zero(FT), icl_rai = zero(FT), icl_sno = zero(FT),
    rai_lcl = zero(FT), rai_icl = zero(FT), rai_rai = zero(FT), rai_sno = zero(FT),
    sno_lcl = zero(FT), sno_icl = zero(FT), sno_rai = zero(FT), sno_sno = zero(FT),
) where {FT} = SA.SMatrix{4, 4, FT}(
    # column-major: each column is a donor q_*, each row a receiver dq_*/dt
    lcl_lcl, icl_lcl, rai_lcl, sno_lcl,  # ∂/∂q_lcl
    lcl_icl, icl_icl, rai_icl, sno_icl,  # ∂/∂q_icl
    lcl_rai, icl_rai, rai_rai, sno_rai,  # ∂/∂q_rai
    lcl_sno, icl_sno, rai_sno, sno_sno,  # ∂/∂q_sno
)

"""
    _jacobian_1m_linearized(g::Raw1MTendency, x::MicroState1M)

Donor-based Jacobian provider for [`_rosenbrock_average_1m`](@ref): the donor-based
linearized system matrix `M` that [`LinearizedAverage`](@ref) uses (built by
`_linearize`), assembled by name via [`_jacobian_1m`](@ref).

`M` is the donor-based linearization, not the exact derivative of the raw
tendency: each donor to receiver transfer is linearized only in its donor species,
so collector-species couplings are absent. It is evaluated at the same
`(ρ, Tsub, q_tot, x)` as the raw tendency `f = g(x)`.
"""
@inline function _jacobian_1m_linearized(g::Raw1MTendency, x::MicroState1M{FT}) where {FT}
    (; q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(
        Microphysics1Moment(), g.mp, g.tps, g.ρ, g.T, FT(g.q_tot),
        q_lcl, q_icl, q_rai, q_sno,
    )
    q_min = TDI.TD.Parameters.q_min(g.tps)
    M = _linearize(src, q_lcl, q_icl, q_rai, q_sno, q_min)
    return _jacobian_1m(FT;
        lcl_lcl = M.M11, lcl_icl = M.M12,
        icl_icl = M.M22,
        rai_lcl = M.M31, rai_rai = M.M33, rai_sno = M.M34,
        sno_lcl = M.M41, sno_icl = M.M42, sno_rai = M.M43, sno_sno = M.M44,
    )
end

"""
    _vapor_exchange_rates(g::Raw1MTendency, x::MicroState1M, q_tot)

The four vapor-to-species phase-change rates `(S_phase_change_vap_lcl,
S_phase_change_vap_icl, S_phase_change_vap_rai, S_phase_change_vap_sno)` of the
1M source terms at state `x` and total water `q_tot`, as a `MicroState1M` in the
`(lcl, icl, rai, sno)` receiver order. Used to differentiate the vapor-exchange
coupling with respect to `q_tot`.
"""
@inline function _vapor_exchange_rates(g::Raw1MTendency, x::MicroState1M, q_tot)
    (; q_lcl, q_icl, q_rai, q_sno) = x
    src = _microphysics_source_terms(
        Microphysics1Moment(), g.mp, g.tps, g.ρ, g.T, q_tot,
        q_lcl, q_icl, q_rai, q_sno,
    )
    return MicroState1M(
        src.S_phase_change_vap_lcl, src.S_phase_change_vap_icl,
        src.S_phase_change_vap_rai, src.S_phase_change_vap_sno,
    )
end

"""
    _jacobian_1m_relinearized(g::Raw1MTendency, x::MicroState1M)

Coupled donor-based Jacobian provider for [`_rosenbrock_average_1m`](@ref): the
donor-based matrix [`_jacobian_1m_linearized`](@ref) with the vapor-competition
(Wegener-Bergeron-Findeisen) coupling restored.

The coupling is approximated by its vapor part: each rate depends on the
vapor specific content `q_vap = q_tot - q_lcl - q_icl - q_rai - q_sno`, giving a
condensate derivative `-∂(rate)/∂q_vap`, obtained from the derivative of the
vapor-to-species rates ([`_vapor_exchange_rates`](@ref)) with respect to `q_tot`
and added to each receiver row. The direct condensate dependence of the rates
(for example rain ventilation and the condensate availability terms) is not
recovered; use [`ExactJacobian`](@ref) for the full derivative.
"""
@inline function _jacobian_1m_relinearized(g::Raw1MTendency, x::MicroState1M{FT}) where {FT}
    Jdonor = _jacobian_1m_linearized(g, x)
    dS_dq_tot = FD.derivative(qt -> _vapor_exchange_rates(g, x, qt), FT(g.q_tot))
    wbf = -dS_dq_tot
    coupling = SA.SMatrix{4, 4, FT}(
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
        wbf[1], wbf[2], wbf[3], wbf[4],
    )
    return Jdonor + coupling
end
