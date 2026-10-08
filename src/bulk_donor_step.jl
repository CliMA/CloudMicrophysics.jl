#####
##### Donor-based linearized tendency and its backward Euler step
#####

"""
    LinearizedTendency{Q, MT, ET}
    LinearizedTendency(Q)

Linearized tendency `dq/dt ≈ M q + e` of the species of the state type `Q`.

The `M` matrix and `e` vector can be accessed using the species names from `Q` as symbols,
e.g., if `Q = Condensates1M` (which has fields: `q_lcl`, `q_icl`, `q_rai`, `q_sno`), then
`lin[:q_rai, :q_lcl]` is the entry of `M` that transfers mass from `q_lcl` to `q_rai`, and
`lin[:q_rai]` is the entry of `e` that adds mass to `q_rai`.

Both `M` and `e` can be assigned, e.g., `lin[:q_rai, :q_lcl] += D`.

`LinearizedTendency(Q)` returns a linearized tendency with `M = 0` and `e = 0`.

# Fields
- `M::MMatrix`: matrix of the model [1/s].
- `e::MVector`: source vector of the model [kg/kg/s].
"""
struct LinearizedTendency{Q, MT, ET}
    M::MT
    e::ET
end
@inline function LinearizedTendency(::Type{Q}) where {Q}
    N = length(Q)
    M = zero(SA.MMatrix{N, N, eltype(Q), N * N})
    e = zero(SA.MVector{N, eltype(Q)})
    return LinearizedTendency{Q, typeof(M), typeof(e)}(M, e)
end

# name-based access and assignment of the `M` matrix and `e` vector.
@inline Base.getindex(lin::LinearizedTendency{Q}, row::Symbol, col::Symbol) where {Q} =
    lin.M[Base.fieldindex(Q, row), Base.fieldindex(Q, col)]
@inline Base.getindex(lin::LinearizedTendency{Q}, species::Symbol) where {Q} =
    lin.e[Base.fieldindex(Q, species)]
@inline function Base.setindex!(lin::LinearizedTendency{Q}, v, row::Symbol, col::Symbol) where {Q}
    lin.M[Base.fieldindex(Q, row), Base.fieldindex(Q, col)] = v
    return lin
end
@inline function Base.setindex!(lin::LinearizedTendency{Q}, v, species::Symbol) where {Q}
    lin.e[Base.fieldindex(Q, species)] = v
    return lin
end

"""
    _relaxation_transfer(S, τ, Δt)

Return the transfer of a relaxation over a substep of width `Δt` [kg/kg].

The relaxation `dq/dt = (q⋆ - q) / τ`, with the timescale `τ` [s], starts with the rate
`S = (q⋆ - q) / τ` [kg/kg/s] and transfers
`S τ (1 - exp(-Δt/τ))` over the substep (Morrison and Milbrandt, 2015, Appendix C):
- For `Δt ≪ τ`, the transfer tends to `S Δt`, the transfer at the instantaneous rate.
- For `Δt ≫ τ`, the transfer tends to `S τ = q⋆ - q`,
  so the substep does not pass the equilibrium `q⋆`.

A disabled process has `S = 0` and `τ = Inf`, and its transfer is zero.
"""
@inline function _relaxation_transfer(S, τ, Δt)
    τ_c = clamp(τ, eps(typeof(τ)), floatmax(typeof(τ)))
    return S * τ_c * -expm1(-Δt / τ_c)
end

"""
    _is_ice(species)

Whether the 1-moment species `species` is an ice species (`:q_icl`, `:q_sno`); the other
species (`:q_lcl`, `:q_rai`) are liquid.
"""
@inline _is_ice(species::Symbol) = species === :q_icl || species === :q_sno

"""
    _changes_phase(term)

Whether `term` releases or absorbs latent heat: every exchange with the vapor does, and a
[`Transfer`](@ref) does when its donor and receiver are of different phases (freezing, melting,
riming); a transfer within one phase (autoconversion, accretion, shedding) does not.
"""
@inline _changes_phase(::Transfer{Donor, Receiver}) where {Donor, Receiver} = _is_ice(Donor) != _is_ice(Receiver)
@inline _changes_phase(::Union{VaporExchange, VaporRelaxation, JointVaporTransfer}) = true

# decay scale of a term: the limiter factor `f_k` of its donor for the terms that change phase, 1 otherwise
@inline _decay_scale(t::Transfer{Donor}, f_k) where {Donor} = ifelse(_changes_phase(t), f_k[Donor], one(f_k[Donor]))
@inline _decay_scale(
    ::Union{VaporExchange{Condensate}, VaporRelaxation{Condensate}, JointVaporTransfer{Condensate}},
    f_k,
) where {Condensate} =
    f_k[Condensate]

"""
    donor_coefficients(term, q, q_min, Δt, f_k = one.(q))

Return the donor-based linearization `(s, D)` of `term` at the species `q`,
for a substep of width `Δt`.

The term changes its donor, the species that loses mass, by `s - D q_donor`,
with the source `s` [kg/kg/s] and the decay `D` [1/s]:
- `Transfer(:Donor => :Receiver, S)`: `s = 0` and `D = S / max(q_min, q[Donor])`.
- `VaporExchange(:Condensate, S)`: `s = S` if `S ≥ 0`, and
  `D = -S / max(q_min, q[Condensate])` otherwise.
- `VaporRelaxation(:Condensate, S, τ)`: the transfer over the substep is
  `Δq = S τ (1 - exp(-Δt/τ))`; `s = Δq / Δt` if `Δq ≥ 0`, and
  `D = -Δq / (max(q[Condensate] + Δq, q_min) Δt)` otherwise.
- `JointVaporTransfer(:Condensate, Δq)`: the same with the transfer `Δq` solved by the joint relaxation.

The optional `f_k`, a state like `q`, scales the decay `D` of a term that changes phase
(`_changes_phase`) by the entry of its donor; it holds the per-donor factors of the
latent-heating limiter of `linearized_step_1m`. The default leaves the decays unscaled.
"""
@inline donor_coefficients(t::Transfer{Donor}, q, q_min, Δt, f_k = one.(q)) where {Donor} =
    (zero(t.S), _decay_scale(t, f_k) * t.S / max(q_min, q[Donor]))
@inline function donor_coefficients(
    t::VaporExchange{Condensate, FT}, q, q_min, Δt, f_k = one.(q),
) where {Condensate, FT}
    is_source = t.S >= zero(FT)
    D = _decay_scale(t, f_k) * (-t.S) / max(q_min, q[Condensate])
    return (ifelse(is_source, t.S, zero(FT)), ifelse(is_source, zero(FT), D))
end
@inline function donor_coefficients(t::VaporRelaxation{Condensate}, q, q_min, Δt, f_k = one.(q)) where {Condensate}
    Δq = _relaxation_transfer(t.S, t.τ, Δt)
    return _transfer_coefficients(Δq, q[Condensate], q_min, Δt, _decay_scale(t, f_k))
end
@inline donor_coefficients(t::JointVaporTransfer{Condensate}, q, q_min, Δt, f_k = one.(q)) where {Condensate} =
    _transfer_coefficients(t.Δq, q[Condensate], q_min, Δt, _decay_scale(t, f_k))
# source and matched decay of a transfer `Δq` over the substep, the decay scaled by `f`: the unscaled
# decay removes exactly `|Δq|` when acting alone (the `q_min` floor keeps it finite when the whole pool is removed)
@inline function _transfer_coefficients(Δq, q_cond, q_min, Δt, f = one(Δq))
    D = ifelse(Δq < zero(Δq), -Δq / (max(q_cond + Δq, q_min) * Δt), zero(Δq))
    return (max(zero(Δq), Δq) / Δt, f * D)
end

@inline function _add_donor_term!(lin, t::Transfer{Donor, Receiver}, q, q_min, Δt, f_k) where {Donor, Receiver}
    (_, D) = donor_coefficients(t, q, q_min, Δt, f_k)
    lin[Donor, Donor] -= D
    lin[Receiver, Donor] += D
    return nothing
end
@inline function _add_donor_term!(
    lin, t::Union{VaporExchange{Condensate}, VaporRelaxation{Condensate}, JointVaporTransfer{Condensate}}, q, q_min, Δt,
    f_k,
) where {Condensate}
    (s, D) = donor_coefficients(t, q, q_min, Δt, f_k)
    lin[Condensate] += s
    lin[Condensate, Condensate] -= D
    return nothing
end

"""
    donor_linearization(terms, q, q_min, Δt, f_k = one.(q))

Return the donor-based [`LinearizedTendency`](@ref) of `terms` at the species `q`.

Each term adds the coefficients `(s, D)` of [`donor_coefficients`](@ref): `-D` to the
diagonal entry of its donor, `D` to the entry in the row of its receiver if the receiver is
one of the species, and `s` to the entry of `e` of its receiver. The optional `f_k` is passed
to [`donor_coefficients`](@ref): it scales the decays of the terms that change phase by the
entry of their donor (the per-donor factors of the latent-heating limiter).

# Arguments
- `terms`: `NamedTuple` of process terms, such as the output of `_microphysics_source_terms`.
- `q`: species, such as a [`Condensates1M`](@ref) [kg/kg].
- `q_min`: lower bound of the donor in the decay coefficients [kg/kg].
- `Δt`: width of the substep [s].
- `f_k`: decay scale per donor, a state like `q` [-].
"""
@inline function donor_linearization(terms, q::Q, q_min, Δt, f_k = one.(q)) where {Q}
    lin = LinearizedTendency(Q)
    UU.unrolled_foreach(t -> (@inline; _add_donor_term!(lin, t, q, q_min, Δt, f_k)), values(terms))
    return lin
end

"""
    limiter_decay_totals(terms, q, q_min, Δt)

Return per donor the total decay of the terms that change phase and the total decay of the
other terms, `(Dpc, Dcol)`, as two states like `q` [1/s]; the inputs of
[`_donor_limiter_scale`](@ref).
"""
@inline function limiter_decay_totals(terms, q::Q, q_min, Δt) where {Q}
    return UU.unrolled_reduce(
        (acc, t) -> (@inline; _add_decay_total(acc, t, q, q_min, Δt)), values(terms), (zero(q), zero(q)),
    )
end
@inline function _add_decay_total((Dpc, Dcol), t::Transfer{Donor}, q, q_min, Δt) where {Donor}
    (_, D) = donor_coefficients(t, q, q_min, Δt)
    return ifelse(
        _changes_phase(t),
        (Base.setindex(Dpc, Dpc[Donor] + D, Donor), Dcol),
        (Dpc, Base.setindex(Dcol, Dcol[Donor] + D, Donor)),
    )
end
@inline function _add_decay_total(
    (Dpc, Dcol), t::Union{VaporExchange{Condensate}, VaporRelaxation{Condensate}, JointVaporTransfer{Condensate}}, q,
    q_min,
    Δt,
) where {Condensate}
    (_, D) = donor_coefficients(t, q, q_min, Δt)
    return (Base.setindex(Dpc, Dpc[Condensate] + D, Condensate), Dcol)
end

"""
    _donor_limiter_scale(f, Dpc, Dcol, Δt)

Factor `f_k` for the phase-change decays of one donor such that their realized transfer
over the substep, `q f_k Dpc Δt / (1 + (Dcol + f_k Dpc) Δt)`, is `f` times the unscaled one
(`f ∈ (0, 1]`); `f_k ≤ f`, with `f_k = 1` for `f = 1`. `Dpc` is the total decay of the donor's
terms that change phase and `Dcol` that of its other terms ([`limiter_decay_totals`](@ref)).
"""
@inline function _donor_limiter_scale(f, Dpc, Dcol, Δt)
    FT = typeof(f)
    return f * (one(FT) + Dcol * Δt) / (one(FT) + Dcol * Δt + (one(FT) - f) * Dpc * Δt)
end

"""
    donor_rates(terms, q, q_new, α, q_min, Δt, f_k = one.(q))

Return the rate of each term over a donor-based substep of width `Δt` from `q` to `q_new`.

The rates are positive in the direction of the term [kg/kg/s].
With the coefficients of [`donor_coefficients`](@ref) (scaled by the optional `f_k` as in
[`donor_linearization`](@ref)) and the factor `α` that scales the source `e` [-], the rate is:
- `Transfer(:Donor => :Receiver, S)`: `D q_new[Donor]`.
- `VaporExchange(:Condensate, S)` and `VaporRelaxation(:Condensate, S, τ)`:
  `α s - D q_new[Condensate]`.

[`species_tendency`](@ref) of the terms with these rates is `(q_new - q) / Δt`.
"""
@inline donor_rates(terms, q, q_new, α, q_min, Δt, f_k = one.(q)) =
    UU.unrolled_map(t -> _donor_rate(t, q, q_new, α, q_min, Δt, f_k), terms)

@inline function _donor_rate(t::Transfer{Donor}, q, q_new, _, q_min, Δt, f_k) where {Donor}
    (_, D) = donor_coefficients(t, q, q_min, Δt, f_k)
    return D * q_new[Donor]
end
@inline function _donor_rate(
    t::Union{VaporExchange{Condensate}, VaporRelaxation{Condensate}, JointVaporTransfer{Condensate}}, q, q_new, α,
    q_min, Δt,
    f_k,
) where {Condensate}
    (s, D) = donor_coefficients(t, q, q_min, Δt, f_k)
    return α * s - D * q_new[Condensate]
end

"""
    _solve_2x2(A, b)

Return the solution `x` of `A x = b` for a `2 × 2` diagonal block `A` of the backward Euler
step in [`backward_euler_solve`](@ref), by Cramer's rule.

The block is `I/Δt - M` restricted to two species, with the decays `D ≥ 0` of
[`donor_coefficients`](@ref):

    A = [ 1/Δt + D₁    -D₂₁      ]
        [ -D₁₂         1/Δt + D₂ ]

Here `Dⱼ` is the total decay of species `j`, and `Dⱼᵢ ≤ Dⱼ` is the part of it transferred to
species `i`. Then `det A = (1/Δt + D₁)(1/Δt + D₂) - D₁₂ D₂₁ > 0`, so the closed form needs
no pivoting or check for a singular block.
"""
@inline function _solve_2x2(A, b)
    det = muladd(-A[1, 2], A[2, 1], A[1, 1] * A[2, 2])  # one rounding for the difference
    return SA.SVector(
        (b[1] * A[2, 2] - A[1, 2] * b[2]) / det,
        (A[1, 1] * b[2] - A[2, 1] * b[1]) / det,
    )
end

"""
    backward_euler_solve(lin, q, α, Δt)

Return the species `q_new` at the end of a substep of width `Δt` that starts from the species
`q`, from the backward Euler step `(I/Δt - M) q_new = q/Δt + α e` of the
[`LinearizedTendency`](@ref) `lin`, with the source scaled by `α` [-].

For the 1-moment species, with the cloud species `c = (q_lcl, q_icl)` and the precipitation
species `p = (q_rai, q_sno)`, the step is block lower triangular:

    [ A_cc    0    ] [ q_new_c ]   [ b_c ]
    [ -M_pc   A_pp ] [ q_new_p ] = [ b_p ]

with `A = I/Δt - M` and `b = q/Δt + α e`. `M_pc` holds the transfers from the cloud to the
precipitation species, such as autoconversion and accretion. The step is then solved by
[`_solve_2x2`](@ref), first for `q_new_c` and then for `q_new_p`:

    A_cc q_new_c = b_c
    A_pp q_new_p = b_p + M_pc q_new_c

The solve assumes that the block `M_cp` from the precipitation to the cloud species is zero.
"""
@inline function backward_euler_solve((; M, e)::LinearizedTendency{Q}, q::Q, α, Δt) where {Q <: Condensates1M}
    to_index(s) = Base.fieldindex(Q, s)
    c = map(to_index, SA.SVector((:q_lcl, :q_icl)))  # cloud  species indices
    p = map(to_index, SA.SVector((:q_rai, :q_sno)))  # precip species indices
    invΔt = one(eltype(q)) / Δt
    A = invΔt * SA.I - M
    b = α * e + invΔt * q

    q_cloud = _solve_2x2(A[c, c], b[c])
    M_pc = M[p, c]
    r = muladd.(M_pc[:, 1], q_cloud[1], muladd.(M_pc[:, 2], q_cloud[2], b[p]))  # b_p + M_pc q_new_c
    q_precip = _solve_2x2(A[p, p], r)

    q_new = similar(q)
    q_new[c], q_new[p] = q_cloud, q_precip
    return Q(q_new)
end
