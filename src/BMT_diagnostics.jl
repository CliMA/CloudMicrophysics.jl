#####
##### Substep record sinks and the diagnostic entries (`Verbose`, `Trace`)
#####

"""
    RecordSink

Supertype of the substep record sinks.

A sink is the whole of the difference between a diagnostic run and a production one.
The substep body is one function: it evaluates the per-process primal, sums it to the
one tendency, builds one system matrix, takes one solve, rescales the increment, floors
it, and hands what it did to its sink as one context object
([`_record_context`](@ref)). The sink observes and returns nothing, so no diagnostic
mode can move the state the model marches: `Verbose`, `Trace` and production produce
bit-identical states by construction, not by a test comparing two implementations.

Because the sink is an ordinary argument and the physics never branches on it, the
diagnostic modes carry no knowledge of the [`Jacobian`](@ref) option, of the state
size, or of the microphysics scheme. They serve [`ExactJacobian`](@ref) and
[`ManualJacobian`](@ref) alike, the 8- and 9-component states alike, and 1M and 2M+P3
alike, with no mode-specific code anywhere in this file.

Production's sink is [`NullSink`](@ref), which compiles away. [`TraceSink`](@ref) and
[`VerboseSink`](@ref) collect one record per substep and are reached only through
[`Trace`](@ref) and [`Verbose`](@ref).
"""
abstract type RecordSink end

"""
    NullSink()

The sink production runs with: [`record!`](@ref) does nothing.

Nothing about a production substep changes for being recordable. `NullSink` has no
fields, its `record!` method has an empty body, and the per-process rates a
[`VerboseSink`](@ref) would attribute are never evaluated, because
[`_recorded_processes`](@ref) dispatches the request away without calling its argument.
The context object is then built from values the substep already holds and immediately
discarded, so the compiler removes the whole recording layer: no arithmetic, no
allocation, and no branch on the mode.
"""
struct NullSink <: RecordSink end

"""
    TraceSink()

The sink of [`Trace`](@ref): one record per substep of what that substep realized.

Each record is `(; δ_acc, α_w, α_l)`:

- `δ_acc`: the increment the substep accepted, `x_new - x`, so the state after substep
  `k` is the entry state plus the first `k` increments.
- `α_w`: the water bound's rescale of the increment, one wherever the bound is inactive.
- `α_l`: the [`TendencyLimiter`](@ref)'s rescale, one wherever the limiter is inactive.

Records are collected in a `Vector` and narrowed to a concrete element type by the
entry that owns the sink, so an element is boxed while the trace is being built. That
cost belongs to the diagnostic path alone, which is also the only path that constructs
a sink with storage.
"""
struct TraceSink <: RecordSink
    records::Vector{Any}
end
TraceSink() = TraceSink(Any[])

"""
    VerboseSink()

The sink of [`Verbose`](@ref): one record per substep attributing the increment that
substep accepted to the processes that produced it.

Each record is `(; pieces, correction, total, α_w, α_l)`:

- `pieces`: a `NamedTuple` over the primal's process slots, each the increment that
  process's own rate produces through the substep's own system matrix and its own
  scalars.
- `correction`: `total - sum(values(pieces))`, the part of the accepted increment that
  no process is responsible for. The linear solve and the two scalar rescales are
  attributed exactly, so what remains is the positivity floor and the rime-pair
  projection, plus the roundoff of re-summing the solve.
- `total`: the increment the substep accepted, the same quantity a
  [`TraceSink`](@ref) records as `δ_acc`.
- `α_w`, `α_l`: the water bound's and the limiter's rescale, the same two scalars a
  [`TraceSink`](@ref) records, carried alongside the attribution so a caller can check
  the `correction` claim rather than take it on faith: at `α_w = α_l = 1` neither
  rescale touched this substep, so a nonzero `correction` there is attributable to the
  positivity floor and the rime-pair projection alone, without needing a floor-engaged
  flag this sink does not carry.

The reconstruction `total == sum(values(pieces)) + correction` holds by construction,
because `correction` is defined as that difference. It is a statement about
attribution, not a claim that the floors did nothing: `correction` vanishes only where
no floor engaged, and asserting that it vanishes everywhere would assert the absence of
the effect the slot exists to expose.

Records are stored as in [`TraceSink`](@ref).
"""
struct VerboseSink <: RecordSink
    records::Vector{Any}
end
VerboseSink() = VerboseSink(Any[])

"""
    SubstepOperator

Supertype of the substep's system matrix, carried in the record context as the substep
actually applied it: `W \\ f` is the increment production takes for the right-hand side
`f`.

Two operators cover the substep's two branches. [`EquilibratedSolve`](@ref) is the
linearized-implicit system `(I/h - P J P)` in the equilibrated variables of
[`_rosenbrock_system`](@ref); [`ExplicitStep`](@ref) is that system reduced to `I/h`,
whose solve is the forward-Euler increment the substep falls back to. A sink therefore
attributes whichever increment production actually took, and keeps no branch of its own
that could fall out of step with the substep's.

`\\` is linear in `f` under both operators, which is what makes the attribution add up:
solving against each process's rate and summing returns the solve against their sum.
"""
abstract type SubstepOperator end

"""
    EquilibratedSolve(S, S⁻¹, A)

The substep's linearized-implicit system, as [`_rosenbrock_system`](@ref) built it and
[`_rosenbrock_solve`](@ref) inverts it.

A sink solves against the factorization the substep already holds, so a per-process
increment costs a back substitution and not a second system build. This is what makes
the per-process attribution the decomposition of production's own increment rather than
a set of independent solves that would have to be argued back into agreement with it.
"""
struct EquilibratedSolve{M} <: SubstepOperator
    S::M
    S⁻¹::M
    A::M
end

"""
    ExplicitStep(h)

The substep's system reduced to `I/h`, the operator of the forward-Euler branch, whose
solve is `h f`.

Written as that multiplication rather than as a solve of `I/h`, so a recorded increment
is bit-identical to the fallback production takes; `v / (1/h)` and `h * v` are not the
same floating-point number.
"""
struct ExplicitStep{FT} <: SubstepOperator
    h::FT
end

@inline Base.:\(W::EquilibratedSolve, f) = _rosenbrock_solve(W.S, W.S⁻¹, W.A, f)
@inline Base.:\(W::ExplicitStep, f) = W.h .* f

"""
    _record_context(pp, W, h, α_w, α_l, δ_acc)

The one object a substep hands its sink: the per-process rates `pp` it summed into its
tendency, the system matrix `W` it solved ([`SubstepOperator`](@ref)), the substep
width `h`, the water bound's rescale `α_w`, the limiter's rescale `α_l`, and the
increment `δ_acc` it accepted after the positivity floor and the pair projection.

Every field is a value the substep computed for its own use, so a record cannot
describe a different substep from the one that produced it, and a sink cannot
desynchronize from the march.

`pp` is `nothing` for a sink that attributes no processes, and lives in the substep's
own state space: the 9-component march hands per-process rates whose temperature row is
that process's own latent heating, which is exact because the row is linear in the mass
rates.
"""
@inline _record_context(pp, W, h, α_w, α_l, δ_acc) = (; pp, W, h, α_w, α_l, δ_acc)

"""
    _recorded_processes(sink, per_process)

The per-process rates for `sink`'s context: `per_process()` where the sink attributes
them, and `nothing` otherwise, without calling `per_process`.

The choice is by dispatch on the sink type, so a production substep neither evaluates
the decomposition nor tests a flag to decide not to.
"""
@inline _recorded_processes(::RecordSink, per_process) = nothing
@inline _recorded_processes(::VerboseSink, per_process) = per_process()

"""
    record!(sink, ctx)

Observe one substep through the context [`_record_context`](@ref) built for it.

A `record!` method reads the context and writes to its sink. It returns `nothing` and
touches no state, which is what makes the diagnostic modes free of any influence on the
march.

[`VerboseSink`](@ref)'s method attributes the accepted increment: it solves the
substep's own system `W` against each process's rate and applies the substep's own
`α_w` and `α_l` in the substep's own order, so the pieces are the terms of production's
increment and not a second opinion about it. Whatever the floor and the pair projection
then did to that increment lands in `correction`, defined as the difference, which is
why the reconstruction is exact and why `correction` is not asserted to be zero.
"""
@inline record!(::NullSink, ctx) = nothing

function record!(sink::TraceSink, ctx)
    push!(sink.records, (; ctx.δ_acc, ctx.α_w, ctx.α_l))
    return nothing
end

function record!(sink::VerboseSink, ctx)
    (; pp, W, α_w, α_l, δ_acc) = ctx
    pieces = map(f_p -> oftype(f_p, α_l .* (α_w .* (W \ f_p))), pp)
    correction = δ_acc - sum(values(pieces))
    push!(sink.records, (; pieces, correction, total = δ_acc, α_w, α_l))
    return nothing
end

"""
    _narrow_records(records)

The collected records with a concrete element type, so a caller iterating them is not
reading a `Vector{Any}`. The records themselves are unchanged.
"""
_narrow_records(records::Vector) = map(identity, records)

"""
    _trace_extras(sink::TraceSink)

The [`Trace`](@ref) fields merged onto the wrapped mode's own return: `substeps`, one
record per substep in the order the march took them.
"""
_trace_extras(sink::TraceSink) = (; substeps = _narrow_records(sink.records))

"""
    _verbose_extras(sink::VerboseSink, Δt)

The [`Verbose`](@ref) fields merged onto the wrapped mode's own return.

- `processes`: the per-process averaged tendencies, each process's increments summed
  over the substeps and divided by `Δt`, in the same slots the primal returns.
- `correction`: the non-attributable part of the averaged tendency, the substeps'
  corrections summed and divided by `Δt`.
- `substeps`: the per-substep records, for a caller that needs the attribution substep
  by substep rather than averaged.

`sum(values(processes)) + correction` reproduces the net averaged tendency to the
roundoff of re-summing across substeps; within a single record the identity is exact
(see [`VerboseSink`](@ref)).
"""
function _verbose_extras(sink::VerboseSink, Δt)
    substeps = _narrow_records(sink.records)
    pieces_sum = reduce((a, b) -> map(+, a, b), (r.pieces for r in substeps))
    processes = map(δ -> δ / Δt, pieces_sum)
    correction = sum(r.correction for r in substeps) / Δt
    return (; processes, correction, substeps)
end

"""
    bulk_microphysics_tendencies(v::Verbose{<:RosenbrockAverage}, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1, w = zero(ρ), p = zero(ρ))

The wrapped mode's own averaged 2M+P3 tendency, with the per-process attribution of
every substep increment that produced it.

This entry runs the production march, with a [`VerboseSink`](@ref) in place of
production's [`NullSink`](@ref) and no other difference, so the net tendency it returns
is bit-identical to the wrapped mode's own: the limiter, the water bound, the
positivity floor, the rime-pair projection and the per-substep shape refresh are all
the ones production runs, because they are the same code. The attribution works for
either [`Jacobian`](@ref) option and for either state size with no mode-specific code.

Returns everything the wrapped mode returns, plus `processes`, `correction` and
`substeps`; see [`_verbose_extras`](@ref) for what each holds and
[`VerboseSink`](@ref) for the reconstruction identity they satisfy.

This is a diagnostic path: it allocates its records on the host and is not meant for
the model time step.
"""
function bulk_microphysics_tendencies(
    v::Verbose{<:RosenbrockAverage}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1, w = zero(ρ), p = zero(ρ),
)
    sink = VerboseSink()
    net = bulk_microphysics_tendencies(
        v.mode, cm, mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub, w, p; sink,
    )
    return merge(net, _verbose_extras(sink, Δt))
end

"""
    bulk_microphysics_tendencies(t::Trace{<:RosenbrockAverage}, ::Microphysics2Moment,
        mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub = 1, w = zero(ρ), p = zero(ρ))

The wrapped mode's own averaged 2M+P3 tendency, with the substep sequence that produced
it: one [`TraceSink`](@ref) record per substep, holding the increment that substep
accepted and the two scalars it was rescaled by.

As for [`Verbose`](@ref), the march is production's, so a traced run cannot drift from
what the model does, and the net tendency is bit-identical to the wrapped mode's. The
number of substeps is the mode's ordinary runtime `nsub`, and tracing costs nothing on
any entry that constructs no `Trace`.

Returns everything the wrapped mode returns, plus `substeps`.
"""
function bulk_microphysics_tendencies(
    t::Trace{<:RosenbrockAverage}, cm::Microphysics2Moment,
    mp::CMP.Microphysics2MParams, tps,
    ρ, T, q_tot,
    q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
    Δt, nsub = 1, w = zero(ρ), p = zero(ρ),
)
    sink = TraceSink()
    net = bulk_microphysics_tendencies(
        t.mode, cm, mp, tps, ρ, T, q_tot,
        q_lcl, n_lcl, q_rai, n_rai, q_ice, n_ice, q_rim, b_rim, logλ,
        Δt, nsub, w, p; sink,
    )
    return merge(net, _trace_extras(sink))
end
