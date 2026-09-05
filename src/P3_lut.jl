"""
    IceKernelTableGrid(; nx, nf, nr, na, logx̄_lo, logx̄_hi, ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi)

Sizes and axis bounds of an [`IceKernelTable`](@ref) grid, carrying no data.

A caller filling a table on the device needs the axis decode before any table exists, so the
geometry is a standalone object rather than a field only [`IceKernelTable`](@ref) carries. Axes
are uniform in their own coordinate. `F_rim` spans its full physical range, `[0, 1]`, by
construction.

# Fields
$(DocStringExtensions.FIELDS)
"""
struct IceKernelTableGrid{FT}
    "grid points along the mean-particle-mass axis"
    nx::Int
    "grid points along the rime-mass-fraction axis"
    nf::Int
    "grid points along the rime-density axis"
    nr::Int
    "grid points along the air-density axis"
    na::Int
    "log10 of the lower mean-particle-mass bound [kg]"
    logx̄_lo::FT
    "log10 of the upper mean-particle-mass bound [kg]"
    logx̄_hi::FT
    "lower rime-mass-fraction bound, always `0`"
    F_rim_lo::FT
    "upper rime-mass-fraction bound, always `1`"
    F_rim_hi::FT
    "lower rime density bound [kg/m³]"
    ρ_rim_lo::FT
    "upper rime density bound [kg/m³]"
    ρ_rim_hi::FT
    "log10 of the lower air-density bound [kg/m³]"
    logρₐ_lo::FT
    "log10 of the upper air-density bound [kg/m³]"
    logρₐ_hi::FT
end

function IceKernelTableGrid(;
    nx::Int = 75, nf::Int = 13, nr::Int = 20, na::Int = 18,
    logx̄_lo, logx̄_hi, ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi,
)
    FT = typeof(logx̄_lo)
    IceKernelTableGrid{FT}(
        nx, nf, nr, na, logx̄_lo, FT(logx̄_hi), zero(FT), one(FT),
        FT(ρ_rim_lo), FT(ρ_rim_hi), FT(logρₐ_lo), FT(logρₐ_hi),
    )
end

# Scalar into a device broadcast, matching the quadrature rules' own methods; without it a caller
# filling a table by broadcasting an entry function has the grid collected as a container.
Base.broadcastable(g::IceKernelTableGrid) = (g,)

"""
    grid_size(grid)

The array shape an [`IceKernelTableGrid`](@ref) describes, for a caller allocating a table.
"""
@inline grid_size(g::IceKernelTableGrid) = (g.nx, g.nf, g.nr, g.na)
@inline grid_length(g::IceKernelTableGrid) = g.nx * g.nf * g.nr * g.na

# Coordinate of index `k` (zero-based) on a uniform axis of `n` points spanning `[lo, hi]`.
@inline _axis_value(lo::FT, hi::FT, n::Int, k::Int) where {FT} =
    n == 1 ? lo : lo + (hi - lo) * FT(k) / FT(n - 1)

"""
    IceKernelTable

Tabulated `log` values of a P3 ice kernel, on a uniform grid in mean particle mass, rime mass
fraction, rime density, and air density.

The ice self-collection rate, the ventilation integral, and the two mean ice fall speeds each
depend on the ice size distribution only through these four state parameters, so one table
type and one [`lookup`](@ref) method serve all of them; the constructors below differ only in
which entry function filled the array. Air density is a table axis because it enters the
Chen (2022) velocity exponents rather than as a size-distribution prefactor, so it changes the
shape of the terminal-velocity curve and cannot be factored out of the kernel.

`FT` is `isbits` only when `A` is a concrete array type, so a table parameterized on a device
array type ships to a kernel like any other quadrature rule.
"""
struct IceKernelTable{FT, A <: AbstractArray{FT, 4}}
    logI::A
    logx̄_lo::FT
    logx̄_hi::FT
    F_rim_lo::FT
    F_rim_hi::FT
    ρ_rim_lo::FT
    ρ_rim_hi::FT
    logρₐ_lo::FT
    logρₐ_hi::FT
end

# The table reaches a kernel inside a `P3TabulatedQuadrature`, captured by closures at many
# host call sites, rather than as a broadcast argument of its own. A carrier holding a bare
# device array is refused at compile time (not `isbits`); one holding a manually converted
# reference whose array has since been collected reads freed memory, which returns zeros
# rather than faulting - `exp(0) = 1` at every log-valued entry, a finite and plausible but
# wrong answer. Converting through this rule, inside the launch machinery that owns the
# array's lifetime, removes that failure class.
Adapt.adapt_structure(to, t::IceKernelTable) = IceKernelTable(
    Adapt.adapt(to, t.logI),
    t.logx̄_lo,
    t.logx̄_hi,
    t.F_rim_lo,
    t.F_rim_hi,
    t.ρ_rim_lo,
    t.ρ_rim_hi,
    t.logρₐ_lo,
    t.logρₐ_hi,
)

@inline function _table_state(scheme, logx̄, F_rim, ρ_rim, logρₐ)
    FT = typeof(logx̄)
    ρₐ = exp10(logρₐ)
    x̄ = exp10(logx̄)
    # Every stored kernel is independent of the number concentration once its own call-site
    # scaling is divided back out, so any positive value serves as the generating number here.
    n = FT(1e4)
    q = x̄ * n
    q_rim = F_rim * q
    b_rim = ρ_rim > 0 ? q_rim / ρ_rim : zero(FT)
    state = state_from_prognostic(scheme, ρₐ * q, ρₐ * n, ρₐ * q_rim, ρₐ * b_rim)
    (state, get_distribution_logλ(state), ρₐ, ρₐ * n)
end

"""
    ice_self_collection_table_entry(scheme, vel, quad, logx̄, F_rim, ρ_rim, logρₐ)

One grid entry: the number-factored self-collection kernel `I = dNdt / ρn_ice²`, evaluated
through [`ice_self_collection`](@ref) so the table cannot drift from the physics it replaces.
"""
@inline function ice_self_collection_table_entry(scheme, vel, quad, logx̄, F_rim, ρ_rim, logρₐ)
    (state, logλ, ρₐ, ρn) = _table_state(scheme, logx̄, F_rim, ρ_rim, logρₐ)
    ice_self_collection(state, logλ, vel, ρₐ; quad).dNdt / ρn^2
end

"""
    ice_ventilation_table_entry(scheme, vel, aps, quad, logx̄, F_rim, ρ_rim, logρₐ)

One grid entry of the ventilation integral `∫ D · F_v(D) · N′(D) dD`, divided by `ρn_ice` so
the stored value is independent of the number concentration.

The table serves both [`ice_deposition_timescale`](@ref) and [`ice_melt`](@ref): they build
this integrand identically over the same bounds, and temperature enters neither of them
through the integral itself, only through the scalar prefactors each applies afterward.
"""
@inline function ice_ventilation_table_entry(scheme, vel, aps, quad, logx̄, F_rim, ρ_rim, logρₐ)
    (state, logλ, ρₐ, ρn) = _table_state(scheme, logx̄, F_rim, ρ_rim, logρₐ)
    ice_ventilation_integral(vel, aps, ρₐ, state, logλ; quad) / ρn
end

"""
    ice_velocity_n_table_entry(scheme, vel, quad, logx̄, F_rim, ρ_rim, logρₐ)
    ice_velocity_m_table_entry(scheme, vel, quad, logx̄, F_rim, ρ_rim, logρₐ)

One grid entry of the number- and mass-weighted mean ice fall speeds. These are already means,
so the stored value carries no call-site scaling at all.
"""
@inline function ice_velocity_n_table_entry(scheme, vel, quad, logx̄, F_rim, ρ_rim, logρₐ)
    (state, logλ, ρₐ, _) = _table_state(scheme, logx̄, F_rim, ρ_rim, logρₐ)
    ice_terminal_velocity_number_weighted(vel, ρₐ, state, logλ; quad)
end

@inline function ice_velocity_m_table_entry(scheme, vel, quad, logx̄, F_rim, ρ_rim, logρₐ)
    (state, logλ, ρₐ, _) = _table_state(scheme, logx̄, F_rim, ρ_rim, logρₐ)
    ice_terminal_velocity_mass_weighted(vel, ρₐ, state, logλ; quad)
end

"""
    _axes_at(grid, lin)

The four axis coordinates of linear entry `lin` (one-based).
"""
@inline function _axes_at(grid::IceKernelTableGrid, lin::Integer)
    i = Int(lin) - 1
    ix = i % grid.nx
    i ÷= grid.nx
    jf = i % grid.nf
    i ÷= grid.nf
    kr = i % grid.nr
    i ÷= grid.nr
    ma = i
    (_axis_value(grid.logx̄_lo, grid.logx̄_hi, grid.nx, ix),
        _axis_value(grid.F_rim_lo, grid.F_rim_hi, grid.nf, jf),
        _axis_value(grid.ρ_rim_lo, grid.ρ_rim_hi, grid.nr, kr),
        _axis_value(grid.logρₐ_lo, grid.logρₐ_hi, grid.na, ma))
end

"""
    ice_self_collection_table_logI(scheme, vel, quad, grid, lin)
    ice_ventilation_table_logK(scheme, vel, aps, quad, grid, lin)
    ice_velocity_n_table_logK(scheme, vel, quad, grid, lin)
    ice_velocity_m_table_logK(scheme, vel, quad, grid, lin)

The stored value of grid entry `lin` (linear, one-based): decode the axes, evaluate the entry,
and take the log, floored to stay finite where the kernel underflows. Scalar and
allocation-free, so a caller may broadcast one of these over a device array to fill a table
where it will be read; [`fill_ice_kernel_table`](@ref) is the host-threaded equivalent.
"""
@inline function ice_self_collection_table_logI(
    scheme, vel, quad, grid::IceKernelTableGrid{FT}, lin::Integer,
) where {FT}
    (lx, fr, rr, la) = _axes_at(grid, lin)
    log(max(ice_self_collection_table_entry(scheme, vel, quad, lx, fr, rr, la), floatmin(FT)))
end

@inline function ice_ventilation_table_logK(
    scheme, vel, aps, quad, grid::IceKernelTableGrid{FT}, lin::Integer,
) where {FT}
    (lx, fr, rr, la) = _axes_at(grid, lin)
    log(max(ice_ventilation_table_entry(scheme, vel, aps, quad, lx, fr, rr, la), floatmin(FT)))
end

@inline function ice_velocity_n_table_logK(
    scheme, vel, quad, grid::IceKernelTableGrid{FT}, lin::Integer,
) where {FT}
    (lx, fr, rr, la) = _axes_at(grid, lin)
    log(max(ice_velocity_n_table_entry(scheme, vel, quad, lx, fr, rr, la), floatmin(FT)))
end

@inline function ice_velocity_m_table_logK(
    scheme, vel, quad, grid::IceKernelTableGrid{FT}, lin::Integer,
) where {FT}
    (lx, fr, rr, la) = _axes_at(grid, lin)
    log(max(ice_velocity_m_table_entry(scheme, vel, quad, lx, fr, rr, la), floatmin(FT)))
end

"""
    IceKernelTable(logI, grid)

Wrap a filled array in its table. This is how a caller-side fill, on a device or otherwise,
returns its result; the shape assertion keeps a wrongly shaped array from being read with its
axes silently reinterpreted.
"""
function IceKernelTable(
    logI::AbstractArray{FT, 4}, grid::IceKernelTableGrid{FT},
) where {FT}
    size(logI) == grid_size(grid) || error(
        "table array is $(size(logI)) but its grid describes $(grid_size(grid))")
    IceKernelTable(
        logI, grid.logx̄_lo, grid.logx̄_hi, grid.F_rim_lo, grid.F_rim_hi,
        grid.ρ_rim_lo, grid.ρ_rim_hi, grid.logρₐ_lo, grid.logρₐ_hi,
    )
end

"""
    fill_ice_kernel_table(logK, grid)

Threaded host fill from a per-entry function `logK(grid, lin)`. Every host constructor below
calls this function with its own entry function, so the four tables share one fill and differ
only in what they store.
"""
function fill_ice_kernel_table(logK, grid::IceKernelTableGrid{FT}) where {FT}
    a = Array{FT}(undef, grid_size(grid))
    Threads.@threads for lin in 1:grid_length(grid)
        @inbounds a[lin] = logK(grid, lin)
    end
    IceKernelTable(a, grid)
end

"""
    IceKernelTable(scheme, vel, grid; quad)
    ice_ventilation_table(scheme, vel, aps, grid; quad)
    ice_velocity_n_table(scheme, vel, grid; quad)
    ice_velocity_m_table(scheme, vel, grid; quad)

Generate a table on the host by evaluating the live integrand at every grid entry with
`Threads.@threads`, on the same grid geometry for all four outputs.

`quad` is the generating quadrature rule and should be of higher order than production's, so
the generator's own quadrature error stays well below the interpolation budget the table is
measured against.
"""
function IceKernelTable(
    scheme, vel, grid::IceKernelTableGrid{FT}; quad,
) where {FT}
    fill_ice_kernel_table((g, lin) -> ice_self_collection_table_logI(scheme, vel, quad, g, lin), grid)
end

ice_ventilation_table(scheme, vel, aps, grid::IceKernelTableGrid; quad) =
    fill_ice_kernel_table((g, lin) -> ice_ventilation_table_logK(scheme, vel, aps, quad, g, lin), grid)
ice_velocity_n_table(scheme, vel, grid::IceKernelTableGrid; quad) =
    fill_ice_kernel_table((g, lin) -> ice_velocity_n_table_logK(scheme, vel, quad, g, lin), grid)
ice_velocity_m_table(scheme, vel, grid::IceKernelTableGrid; quad) =
    fill_ice_kernel_table((g, lin) -> ice_velocity_m_table_logK(scheme, vel, quad, g, lin), grid)

function IceKernelTable(
    scheme, vel;
    quad,
    nx::Int = 75, nf::Int = 13, nr::Int = 20, na::Int = 18,
    logx̄_lo, logx̄_hi, ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi,
)
    grid = IceKernelTableGrid(;
        nx, nf, nr, na, logx̄_lo, logx̄_hi, ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi)
    IceKernelTable(scheme, vel, grid; quad)
end

# Fractional index on a uniform axis of `n` points spanning `[lo, hi]`: arithmetic, not a
# search. `t` clamps into the axis before truncation, because `unsafe_trunc` on a non-finite or
# out-of-range value is undefined behavior and a degenerate state can reach here with a
# non-finite coordinate. Values outside the axis clamp to the edge.
@inline function _uniform_index(v, lo, hi, n::Int)
    t = (v - lo) / (hi - lo) * (n - 1)
    t = isnan(t) ? zero(t) : clamp(t, zero(t), oftype(t, n - 1))
    i = clamp(unsafe_trunc(Int, t), 0, n - 2)
    (i + 1, clamp(t - i, zero(t), one(t)))
end

"""
    lookup(tab::IceKernelTable, x̄, F_rim, ρ_rim, ρₐ)

Quadrilinear interpolation of `log I`, returned exponentiated. Sixteen fetches share one set of
index offsets and weights, which is what makes a second output co-gridded with the first
nearly free.
"""
@inline function lookup(tab::IceKernelTable{FT}, x̄, F_rim, ρ_rim, ρₐ) where {FT}
    nx, nf, nr, na = size(tab.logI)
    (ix, ax) = _uniform_index(log10(x̄), tab.logx̄_lo, tab.logx̄_hi, nx)
    (jf, af) = _uniform_index(F_rim, tab.F_rim_lo, tab.F_rim_hi, nf)
    (kr, ar) = _uniform_index(ρ_rim, tab.ρ_rim_lo, tab.ρ_rim_hi, nr)
    (ma, aa) = _uniform_index(log10(ρₐ), tab.logρₐ_lo, tab.logρₐ_hi, na)
    v = zero(FT)
    @inbounds for di in 0:1, dj in 0:1, dk in 0:1, dm in 0:1
        w = (di == 0 ? 1 - ax : ax) * (dj == 0 ? 1 - af : af) *
            (dk == 0 ? 1 - ar : ar) * (dm == 0 ? 1 - aa : aa)
        v += w * tab.logI[ix + di, jf + dj, kr + dk, ma + dm]
    end
    exp(v)
end

"""
    _ice_self_collection(state, logλ, vel, ρₐ, tab::IceKernelTable)

Tabulated self-collection: the rate through a lookup rather than the nested quadrature that
[`ice_self_collection`](@ref) otherwise evaluates.

A quadrature rule integrates an integrand; a table replaces the whole integral and is keyed on
the state, which by the call site in `P3_processes.jl` has already been captured into closures.
The dispatch therefore sits at this entry, the only level where a table can be substituted.
"""
@inline function _ice_self_collection(state, logλ, vel, ρₐ, tab::IceKernelTable)
    (; ρq_ice, ρn_ice) = state
    FT = typeof(ρₐ)
    (ρq_ice > 0) & (ρn_ice > 0) || return (; dNdt = zero(FT))
    I = lookup(tab, ρq_ice / ρn_ice, state.F_rim, state.ρ_rim, ρₐ)
    (; dNdt = I * ρn_ice^2)
end

"""
    P3TabulatedQuadrature(rule; selfcol = nothing, vent = nothing, vel_n = nothing, vel_m = nothing)

The table mode's carrier: a real quadrature rule, together with whichever
[`IceKernelTable`](@ref)s have been adopted.

The parameter struct holds one quadrature object and every P3 integral consumes it, so
tabulating a single term cannot be expressed by storing a table there directly: the
untabulated terms still need a rule, and read `quad.n` and call [`node`](@ref) /
[`weight`](@ref), which a table alone does not have. A carrier is therefore needed from the
first adopted table onward, and it delegates every `QuadratureRule` method to the rule it
carries, so an untabulated term integrates exactly as before.

Each keyword is a table or `nothing`; `nothing` means the carried rule serves that term, so
which terms are tabulated is exactly which fields are populated. Adding a term to the mode adds
a field here and a dispatch method on this type; nothing else changes.
"""
struct P3TabulatedQuadrature{Q, S, V, N, M} <: QuadratureRule
    "the carried rule, used unchanged by every term that is not tabulated"
    rule::Q
    "the ice self-collection table, or `nothing`"
    selfcol::S
    "the ventilation integral, read by deposition and melt, or `nothing`"
    vent::V
    "the number-weighted mean fall speed, or `nothing`"
    vel_n::N
    "the mass-weighted mean fall speed, or `nothing`"
    vel_m::M
    "mirrors `rule.n`, read directly by the quadrature sums in `Quadrature.integrate`"
    n::Int
end

P3TabulatedQuadrature(rule; selfcol = nothing, vent = nothing, vel_n = nothing,
    vel_m = nothing) = P3TabulatedQuadrature(rule, selfcol, vent, vel_n, vel_m, rule.n)

# Every field is adapted; an absent table is `nothing`, which `Adapt.adapt` returns unchanged,
# so this one method covers every combination of adopted tables, including the all-`nothing`
# carrier that behaves exactly like the rule it carries.
Adapt.adapt_structure(to, q::P3TabulatedQuadrature) = P3TabulatedQuadrature(
    Adapt.adapt(to, q.rule),
    Adapt.adapt(to, q.selfcol),
    Adapt.adapt(to, q.vent),
    Adapt.adapt(to, q.vel_n),
    Adapt.adapt(to, q.vel_m),
    q.n,
)

"""
    P3_LUT_OUTPUTS

The terms [`p3_lut_carrier`](@ref) can tabulate, and the only names it accepts.

A name is a PROCESS rather than a stored array. It happens that each of these four is one table
today, and that is not a property to rely on: a process whose integrand splits by species or by
regime needs several table objects and still takes one name, because a caller selects the work it
wants done differently, not the arrays that get it done.
"""
const P3_LUT_OUTPUTS = (:selfcol, :vent, :vel_n, :vel_m)

"""
    p3_lut_carrier(scheme, vel, aps, grid; quad, generating_quad, outputs)

Build a [`P3TabulatedQuadrature`](@ref) in one call: fill the requested tables on `grid` and
wrap them around `quad`, which remains the rule for every term not requested.

`outputs` states which terms to tabulate, from [`P3_LUT_OUTPUTS`](@ref); a term left out is not
tabulated and its field stays `nothing`.

`quad` is the rule the carrier DELEGATES to for every term it does not tabulate, and
`generating_quad` is the rule the tables are generated FROM. They default to the same rule and
should not be: a generator wants a higher order than production, so that its own quadrature error
stays below the interpolation budget the table is measured against, while an untabulated term must
keep integrating at production's order. Passing one rule for both silently moves every untabulated
term onto the generator's order.

Device residence is the caller's responsibility, as for a single table: pass tables built from
device arrays, or convert them after filling. `CloudMicrophysics` performs no device execution
itself.
"""
function p3_lut_carrier(scheme, vel, aps, grid::IceKernelTableGrid; quad,
    generating_quad = quad, outputs = P3_LUT_OUTPUTS)
    for o in outputs
        o in P3_LUT_OUTPUTS || error(
            "unknown P3 lookup-table output $o; expected one of $(join(P3_LUT_OUTPUTS, ", "))")
    end
    P3TabulatedQuadrature(quad;
        selfcol = :selfcol in outputs ? IceKernelTable(scheme, vel, grid; quad = generating_quad) : nothing,
        vent = :vent in outputs ? ice_ventilation_table(scheme, vel, aps, grid; quad = generating_quad) : nothing,
        vel_n = :vel_n in outputs ? ice_velocity_n_table(scheme, vel, grid; quad = generating_quad) : nothing,
        vel_m = :vel_m in outputs ? ice_velocity_m_table(scheme, vel, grid; quad = generating_quad) : nothing)
end

# The delegation: a carrier is a quadrature rule wherever one is asked for.
@inline node(q::P3TabulatedQuadrature, i, n) = node(q.rule, i, n)
# ... except where a SECOND integral is to be taken with the same rule, which wants the rule
# itself rather than the carrier: handing the carrier on would route that integral's terms through
# the tables as well, which is a different integral from the one asked for.
@inline _plain_rule(q::P3TabulatedQuadrature) = q.rule
@inline weight(q::P3TabulatedQuadrature, i, n) = weight(q.rule, i, n)
@inline inv_weight_fun(q::P3TabulatedQuadrature, y) = inv_weight_fun(q.rule, y)
# Scalar into a device broadcast, matching `GaussLegendre`'s and `ChebyshevGauss`'s own
# methods; without it the carrier is treated as a container and iterated.
Base.broadcastable(q::P3TabulatedQuadrature) = (q,)

# Self-collection reads the table when the carrier holds one and the carried rule when it does
# not, so one carrier expresses both the adopted and the not-adopted state of the term.
@inline _ice_self_collection(state, logλ, vel, ρₐ, q::P3TabulatedQuadrature) =
    _selfcol_from(state, logλ, vel, ρₐ, q.selfcol, q.rule)
@inline _selfcol_from(state, logλ, vel, ρₐ, tab::IceKernelTable, rule) =
    _ice_self_collection(state, logλ, vel, ρₐ, tab)
@inline _selfcol_from(state, logλ, vel, ρₐ, ::Nothing, rule) =
    _ice_self_collection(state, logλ, vel, ρₐ, rule)

# The ventilation integral, read by deposition and melt through `_ventilation_from` in
# `P3_processes.jl`. The scaling is applied here rather than stored: the table holds the
# integral per unit ice number, so the lookup multiplies by `ρn_ice` and the two consumers see
# exactly the integral they used to compute.
@inline function _ice_ventilation(state, logλ, vel, ρₐ, tab::IceKernelTable)
    (; ρq_ice, ρn_ice) = state
    FT = typeof(ρₐ)
    (ρq_ice > 0) & (ρn_ice > 0) || return zero(FT)
    lookup(tab, ρq_ice / ρn_ice, state.F_rim, state.ρ_rim, ρₐ) * ρn_ice
end

@inline _ventilation_from(velocity_params, aps, ρₐ, state, logλ, q::P3TabulatedQuadrature) =
    _vent_from(velocity_params, aps, ρₐ, state, logλ, q.vent, q.rule)
@inline _vent_from(vel, aps, ρₐ, state, logλ, tab::IceKernelTable, rule) =
    _ice_ventilation(state, logλ, vel, ρₐ, tab)
@inline _vent_from(vel, aps, ρₐ, state, logλ, ::Nothing, rule) =
    ice_ventilation_integral(vel, aps, ρₐ, state, logλ; quad = rule)

# The two mean fall speeds are already means, so there is no call-site scaling at all; the
# stored value is the velocity itself. The absent-population return matches what
# `UT.guarded_quotient` returns on the quadrature path.
@inline function _ice_velocity_tab(state, ρₐ, tab::IceKernelTable)
    (; ρq_ice, ρn_ice) = state
    FT = typeof(ρₐ)
    (ρq_ice > 0) & (ρn_ice > 0) || return zero(FT)
    lookup(tab, ρq_ice / ρn_ice, state.F_rim, state.ρ_rim, ρₐ)
end

# `state::P3State` is required here, not decoration: the base method in
# `P3_terminal_velocity.jl` is specific in `state` and this one is specific in the last
# argument, so without it neither method dominates and the call is ambiguous.
@inline _vel_n_from(vel, ρₐ, state::P3State, logλ, p, q::P3TabulatedQuadrature) =
    _veln(vel, ρₐ, state, logλ, p, q.vel_n, q.rule)
@inline _veln(vel, ρₐ, state, logλ, p, tab::IceKernelTable, rule) = _ice_velocity_tab(state, ρₐ, tab)
@inline _veln(vel, ρₐ, state, logλ, p, ::Nothing, rule) = _vel_n_from(vel, ρₐ, state, logλ, p, rule)

@inline _vel_m_from(vel, ρₐ, state::P3State, logλ, p, q::P3TabulatedQuadrature) =
    _velm(vel, ρₐ, state, logλ, p, q.vel_m, q.rule)
@inline _velm(vel, ρₐ, state, logλ, p, tab::IceKernelTable, rule) = _ice_velocity_tab(state, ρₐ, tab)
@inline _velm(vel, ρₐ, state, logλ, p, ::Nothing, rule) = _vel_m_from(vel, ρₐ, state, logλ, p, rule)
