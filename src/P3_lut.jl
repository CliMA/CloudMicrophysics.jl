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

# THE RIME-VOLUME TABLE'S TEMPERATURE AXIS IS GEOMETRIC IN `|T°C|`, not uniform in `1/T°C`.
#
# Temperature enters the rime-volume integrand only through the Cober and List riming index
# `Rᵢ = -(Dₗ·1e6·v_term)/(2·T°C)`, which `LocalRimeDensity` clamps to `[1, 12]`. The integrand's
# active window is therefore exactly a factor of 12 wide in `|T|` for a single drop, but it sits at
# a temperature proportional to `Dₗ·v_term`, so it SLIDES with the ice fall speed, which is another
# axis of this same table. brim-axis-window-slides (40) measured the window's WIDTH to vary by 27
# times across states as well, and the ratio of a state's two thresholds is scale-free, so no
# per-state rescaling of the coordinate can hold every window still. A geometric axis is instead
# scale-INVARIANT: a window of a given width ratio contains the same number of intervals wherever
# it sits, which is the invariance the measurement says is needed.
#
# Measured on this axis alone, at the states' own other five coordinates and over a range wide
# enough that nothing clamps, uniform `1/T°C` reads a p95 of 335.4 percent for cloud at 8 nodes and
# does NOT converge, still 227.9 at 32. Geometric reads 7.53 at 8 and 0.589 at 24, against a 1.29
# percent gate for the whole table.
#
# The endpoints keep both their meaning and their names. `invT_lo` and `invT_hi` are still the two
# endpoint values of `1/T°C`, so a caller's grid construction is unchanged and only the spacing
# between them moves.
@inline _logT_of_invT(v::FT) where {FT} = log(-inv(v))
@inline function _logT_axis_value(invT_lo::FT, invT_hi::FT, n::Int, k::Int) where {FT}
    u = _axis_value(_logT_of_invT(invT_lo), _logT_of_invT(invT_hi), n, k)
    return inv(-exp(u))
end
# A state's coordinate on that axis. `T°C` at or above zero has no logarithm; it is sent to the warm
# end, which is where the clamp already put it and where the rime-volume channels vanish anyway.
@inline _logT_coord(T°C::FT) where {FT} = log(-min(T°C, -floatmin(FT)))

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
    liqice_logx̄_floor(scheme)

The lowest ice mean particle mass a liquid-ice table needs an axis point for, as a base-ten
logarithm: `log10(ice_mean_particle_mass_min(scheme))`.

# Why this is the floor, and why it is derived rather than written down

`_derived_logλ_bracket` brackets the mean-mass target between the scheme's smallest and largest
particle masses, and `get_distribution_logλ` clamps the root into that bracket. Below the smallest
mass the target is unbracketable, `logλ` pins at the endpoint, and every quantity downstream stops
depending on the ice mean mass. Measured on this tree, the entry is BITWISE identical at log10 x̄ of
-13.0 and -12.5 and differs across -12.0, with the release at -12.318774, which is
`log10(ice_mean_particle_mass_min)` to six decimals.

A floor there is therefore optimal in both directions. Below it the lookup's clamp returns the exact
value rather than an interpolated one, because the function is constant. Above it every point is
spent where the function varies. A floor set BELOW the pin spends points on a plateau; a floor set
ABOVE it clamps a band where the function genuinely varies, which is the worse of the two errors and
is silent.

The value is DERIVED and not written as a number on purpose. `ice_mean_particle_mass_min` is
`ice_nucleation_mass`, which is fixed by the nascent-crystal diameter, and that diameter is a
parameter of the scheme rather than a constant of nature: a tree whose nucleation baseline sets it
to 2 micrometres rather than 10 puts the pin at -14.415684 instead, 2.1 decades lower. A floor
written as a literal is correct until that parameter moves and silently wrong afterwards, and the
failure is invisible on the tree it was measured on, because nothing there disagrees with it.
"""
@inline liqice_logx̄_floor(scheme) = log10(ice_mean_particle_mass_min(scheme))

"""
    LiqIceTableGrid(; nx, nf, nr, na, nl, bounds...)

The grid of a liquid-ice collision table: the four axes of [`IceKernelTableGrid`](@ref) plus one in
the mean mass of the LIQUID species, cloud or rain.

# Why a second grid type rather than a rank-generalised one

`IceKernelTable` hard-codes rank four in its type and `_uniform_index`, `lookup` and the bounds
copied onto the table are all written against four axes. Generalising the rank would touch every one
of the four-axis tables, whose throughput is measured and banked, and would put those numbers back
in question to save one struct. Two grid types in one carrier is the cheaper cost by a wide margin.
The axes that both types share keep the same names and the same meaning.
"""
struct LiqIceTableGrid{FT}
    "grid points along the ice mean-particle-mass axis"
    nx::Int
    "grid points along the rime-mass-fraction axis"
    nf::Int
    "grid points along the rime-density axis"
    nr::Int
    "grid points along the air-density axis"
    na::Int
    "grid points along the liquid mean-particle-mass axis"
    nl::Int
    "log10 of the lower ice mean-particle-mass bound [kg]"
    logx̄_lo::FT
    "log10 of the upper ice mean-particle-mass bound [kg]"
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
    "log10 of the lower liquid mean-particle-mass bound [kg]"
    logx̄l_lo::FT
    "log10 of the upper liquid mean-particle-mass bound [kg]"
    logx̄l_hi::FT
end

# Every axis of an interpolated table needs at least two nodes. `_uniform_index` returns the pair
# `(i, i+1)` unconditionally, and the lookup reads both under `@inbounds`, so a one-node axis reads
# outside the array at one end or the other whatever the index is clamped to. The weight on the bad
# read is zero, which makes the result correct whenever the memory happens to hold a finite number
# and silently non-finite when it does not; either way it is an out-of-bounds read. Rejecting it
# here is the only place that can, since the lookup has no size to check against.
@inline function _assert_axis_nodes(ns::NamedTuple)
    for (name, n) in pairs(ns)
        n >= 2 || error(
            "table axis `$name` needs at least 2 nodes, got $n: a linearly " *
            "interpolated axis has no one-node form.",
        )
    end
    return nothing
end

function LiqIceTableGrid(;
    nx::Int = 24, nf::Int = 12, nr::Int = 12, na::Int = 8, nl::Int = 16,
    logx̄_lo, logx̄_hi, ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi, logx̄l_lo, logx̄l_hi,
    F_rim_lo = zero(logx̄_lo), F_rim_hi = one(logx̄_lo),
)
    FT = typeof(logx̄_lo)
    _assert_axis_nodes((; nx, nf, nr, na, nl))
    LiqIceTableGrid{FT}(
        nx, nf, nr, na, nl,
        logx̄_lo, logx̄_hi, F_rim_lo, F_rim_hi,
        ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi, logx̄l_lo, logx̄l_hi,
    )
end

@inline grid_size(g::LiqIceTableGrid) = (g.nx, g.nf, g.nr, g.na, g.nl)

# A grid is a scalar in a broadcast, not a container.
#
# The seam these types exist for is a caller filling the table WHERE IT LIVES, by broadcasting the
# per-entry function over a device array of linear indices with the grid captured as a scalar.
# Without this method the broadcast machinery treats the grid as iterable, falls through to
# `collect`, and throws on the HOST before any kernel is launched, so the failure does not even
# reach the device it was written for. `IceKernelTableGrid` shipped without it and the LUT tier
# found the same hole this week; the only reason it was not seen earlier is that the campaign's own
# probe captured the grid in a closure instead of broadcasting it.
Base.broadcastable(g::LiqIceTableGrid) = (g,)

"""
    _liqice_x̄l_slope(grid, channel)

The exponent `α` such that the stored integral behaves as `x̄_liq^α` at and below the liquid
mean-mass axis floor, DERIVED at fill time rather than written down per species.

# Why the tables carry one

The stored integral is per unit liquid NUMBER. Below the SB2006 mean-mass clamp the two liquid
distributions behave oppositely: the cloud value is constant, so `α = 0` and clamping the axis is
exact, while the windowed rain value is proportional to `x̄_liq`, so `α = 1`, because
`pdf_rain_parameters` builds the amplitude from `Nᵣ_bounded = Lᵣ/xr_mean` rather than from `N` and
that puts a factor of the mean mass into it. With `α = 0` assumed for both, a rain state four
decades below the floor is returned the value AT the floor, which is out by `x̄/xr_min`: measured on
the 2026-08-16 AMIP battery, 272 of 641 states, a median 65 percent and a maximum of 6e+03 times,
and worst where the quantity is largest. See liqice-table-error-anatomy (35).

# Why it is derived and not written down

`α` is a property of the size distribution, not of the species name, so a table built on a different
distribution would need a different constant and nothing would say so. Deriving it makes the table
correct for whatever distribution it was filled from, and the disagreement check below turns a
distribution that does not have a single such exponent into an error rather than a silent bias.

Extending the axis instead is the obvious alternative and is worse: measured, it left the p95
unchanged at 94.8 percent because the clamp then falls BETWEEN nodes and the stored log has a kink
there, which is what axis-kinks-and-coverage (20) requires a node on. Shearing removes the kink
rather than resolving it, and costs no grid points at all.
"""
function _liqice_x̄l_slope(grid, channel)
    FT = typeof(grid.logx̄_lo)
    u0 = grid.logx̄l_lo
    mid(lo, hi) = (lo + hi) / 2
    αs = map((FT(0.05), FT(0.5), FT(0.95))) do t
        lx = grid.logx̄_lo + t * (grid.logx̄_hi - grid.logx̄_lo)
        fr = t
        rr = grid.ρ_rim_lo + t * (grid.ρ_rim_hi - grid.ρ_rim_lo)
        la = mid(grid.logρₐ_lo, grid.logρₐ_hi)
        a = channel(lx, fr, rr, la, u0)
        b = channel(lx, fr, rr, la, u0 - one(FT))
        (a > 0 && b > 0) ? (log(a) - log(b)) / log(FT(10)) : FT(NaN)
    end
    good = filter(isfinite, collect(αs))
    isempty(good) && return zero(FT)
    α = round(sum(good) / length(good))
    maximum(abs.(good .- α)) < FT(0.05) || error(
        "the liquid mean-mass slope below the axis floor is not one exponent: measured $(good). " *
        "A table sheared by a single α cannot represent this size distribution below its clamp.")
    return α
end

"""
    LiqIceTable(logM, logN, x̄l_slope, bounds...)

The full-range part of one liquid species' contribution to the liquid-ice collision entry, on the
five axes of [`LiqIceTableGrid`](@ref).

# What is stored, and why it needs no partition axis

`∫liquid_ice_collisions_split` divides each channel into a full-range integral that carries a
CONSTANT partition and a correction supported on the wet set. Only the first is held here. It
therefore has no dependence on the temperature and none on the liquid magnitudes beyond the two
factors below, because the partition is where both entered.

Two co-gridded outputs, each factored by the two number concentrations it is linear in:

    logM = log( ∫ n_i ∂ₜM_col dD / (ρn_ice · N_liq) )    the collected mass rate
    logN = log( ∫ n_i ∂ₜN_col dD / (ρn_ice · N_liq) )    the collected number rate

The rime volume rate is NOT here. Its integrand divides by the local rime density, which
[`compute_local_rime_density`](@ref) builds from the temperature, so that channel carries `T` and no
table on these five axes can hold it. It takes [`LiqIceBrimTable`](@ref) and its sixth axis instead.

The ice size distribution is linear in `ρn_ice` at fixed mean mass and the liquid one is linear in
`N_liq` at fixed mean mass, so neither number is an axis. That is what keeps the set at five axes,
and it is the property section 9 of `notes/surrogate-options-20260826.md` read off the source.

The two share one set of index offsets and weights, so the second costs thirty-two more fetches and
no more arithmetic.

Storage is logarithmic, which requires every entry to be strictly positive. The generator refuses a
non-positive entry rather than storing its logarithm, because `exp` of a sentinel is a finite wrong
answer and not a detectable one.
"""
struct LiqIceTable{FT, A <: AbstractArray{FT, 5}}
    logM::A
    logN::A
    "the derived `x̄_liq` exponent the stored logs are sheared by; see `_liqice_x̄l_slope`"
    x̄l_slope::FT
    logx̄_lo::FT
    logx̄_hi::FT
    F_rim_lo::FT
    F_rim_hi::FT
    ρ_rim_lo::FT
    ρ_rim_hi::FT
    logρₐ_lo::FT
    logρₐ_hi::FT
    logx̄l_lo::FT
    logx̄l_hi::FT
end

# The same rule, and the same reason, as `IceKernelTable`'s: the conversion belongs inside the launch
# machinery that owns the arrays' lifetime, and the bounds stay `isbits` scalars so that the adapted
# carrier is `isbits` too. See the comment on that rule for the failure it removes.
Adapt.adapt_structure(to, t::LiqIceTable) = LiqIceTable(
    Adapt.adapt(to, t.logM),
    Adapt.adapt(to, t.logN),
    t.x̄l_slope,
    t.logx̄_lo, t.logx̄_hi,
    t.F_rim_lo, t.F_rim_hi,
    t.ρ_rim_lo, t.ρ_rim_hi,
    t.logρₐ_lo, t.logρₐ_hi,
    t.logx̄l_lo, t.logx̄l_hi,
)

"""
    lookup(tab::LiqIceTable, x̄, F_rim, ρ_rim, ρₐ, x̄_liq)

Five-linear interpolation of the two stored logarithms, returned exponentiated as `(M, N)`.
Thirty-two corners, one set of weights, two reads per corner.
"""
@inline function lookup(tab::LiqIceTable{FT}, x̄, F_rim, ρ_rim, ρₐ, x̄_liq) where {FT}
    nx, nf, nr, na, nl = size(tab.logM)
    (ix, ax) = _uniform_index(log10(x̄), tab.logx̄_lo, tab.logx̄_hi, nx)
    (jf, af) = _uniform_index(F_rim, tab.F_rim_lo, tab.F_rim_hi, nf)
    (kr, ar) = _uniform_index(ρ_rim, tab.ρ_rim_lo, tab.ρ_rim_hi, nr)
    (ma, aa) = _uniform_index(log10(ρₐ), tab.logρₐ_lo, tab.logρₐ_hi, na)
    ul = log10(x̄_liq)
    (pl, al) = _uniform_index(ul, tab.logx̄l_lo, tab.logx̄l_hi, nl)
    # The shear is added back with the UNCLAMPED coordinate, which is what makes the region below
    # the axis floor a linear extrapolation in the log rather than a constant. Inside the axis it
    # cancels the shear the fill applied and recovers the stored value exactly at a node.
    shear = tab.x̄l_slope * (ul - tab.logx̄l_lo) * log(FT(10))
    vM = zero(FT)
    vN = zero(FT)
    @inbounds for di in 0:1, dj in 0:1, dk in 0:1, dm in 0:1, dp in 0:1
        w =
            (di == 0 ? 1 - ax : ax) * (dj == 0 ? 1 - af : af) *
            (dk == 0 ? 1 - ar : ar) * (dm == 0 ? 1 - aa : aa) *
            (dp == 0 ? 1 - al : al)
        i, j, k, m, q = ix + di, jf + dj, kr + dk, ma + dm, pl + dp
        vM += w * tab.logM[i, j, k, m, q]
        vN += w * tab.logN[i, j, k, m, q]
    end
    (exp(vM + shear), exp(vN + shear))
end

"""
    LiqIceBrimTableGrid(; nx, nf, nr, na, nl, nt, bounds...)

The grid of the rime-volume table: the five axes of [`LiqIceTableGrid`](@ref) plus one in `1/T°C`.

# Why a sixth axis, and why that coordinate

The rime volume rate divides the collected mass by the LOCAL rime density, which
[`rime_density_at`](@ref) forms as `Rᵢ = -(Dₗ·1e6·v_term)/(2·T°C)` and which uses the temperature
nowhere else. Scaling `v_term` and `T°C` by one factor therefore returns the identical density, an
identity bracket-budget (3) measured to hold to 3e-13 over four temperatures and three scalings. The
temperature enters as a rescaling of one argument rather than as a separable scalar prefactor, so it
cannot be divided out the way the freeze capacity's prefactor can, and one axis in `1/T°C` is what
it costs.
"""
struct LiqIceBrimTableGrid{FT}
    "grid points along the ice mean-particle-mass axis"
    nx::Int
    "grid points along the rime-mass-fraction axis"
    nf::Int
    "grid points along the rime-density axis"
    nr::Int
    "grid points along the air-density axis"
    na::Int
    "grid points along the liquid mean-particle-mass axis"
    nl::Int
    "grid points along the inverse-temperature axis"
    nt::Int
    "log10 of the lower ice mean-particle-mass bound [kg]"
    logx̄_lo::FT
    "log10 of the upper ice mean-particle-mass bound [kg]"
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
    "log10 of the lower liquid mean-particle-mass bound [kg]"
    logx̄l_lo::FT
    "log10 of the upper liquid mean-particle-mass bound [kg]"
    logx̄l_hi::FT
    "lower bound of 1/T°C [1/°C], at the coldest temperature the table covers"
    invT_lo::FT
    "upper bound of 1/T°C [1/°C], at the warmest subfreezing temperature the table covers"
    invT_hi::FT
end

function LiqIceBrimTableGrid(;
    nx::Int = 24, nf::Int = 12, nr::Int = 12, na::Int = 8, nl::Int = 16, nt::Int = 8,
    logx̄_lo, logx̄_hi, ρ_rim_lo, ρ_rim_hi, logρₐ_lo, logρₐ_hi, logx̄l_lo, logx̄l_hi,
    invT_lo, invT_hi, F_rim_lo = zero(logx̄_lo), F_rim_hi = one(logx̄_lo),
)
    FT = typeof(logx̄_lo)
    _assert_axis_nodes((; nx, nf, nr, na, nl, nt))
    LiqIceBrimTableGrid{FT}(
        nx, nf, nr, na, nl, nt,
        logx̄_lo, logx̄_hi, F_rim_lo, F_rim_hi, ρ_rim_lo, ρ_rim_hi,
        logρₐ_lo, logρₐ_hi, logx̄l_lo, logx̄l_hi, invT_lo, invT_hi,
    )
end

@inline grid_size(g::LiqIceBrimTableGrid) = (g.nx, g.nf, g.nr, g.na, g.nl, g.nt)

# Scalar into a broadcast, for the reason given on `LiqIceTableGrid`'s own method.
Base.broadcastable(g::LiqIceBrimTableGrid) = (g,)

"""
    LiqIceBrimTable(logB, x̄l_slope, bounds...)

The full-range rime volume rate of one liquid species, on the six axes of
[`LiqIceBrimTableGrid`](@ref):

    logB = log( ∫ n_i ∂ₜB_col dD / (ρn_ice · N_liq) ) - x̄l_slope · (log10 x̄_liq - logx̄l_lo) · ln10

sheared by the derived liquid mean-mass exponent, for the reason `_liqice_x̄l_slope` gives.

One output rather than two, because nothing else on this grid needs the temperature and a channel
that does not need an axis does not pay for one.
"""
struct LiqIceBrimTable{FT, A <: AbstractArray{FT, 6}}
    logB::A
    "the derived `x̄_liq` exponent the stored log is sheared by; see `_liqice_x̄l_slope`"
    x̄l_slope::FT
    logx̄_lo::FT
    logx̄_hi::FT
    F_rim_lo::FT
    F_rim_hi::FT
    ρ_rim_lo::FT
    ρ_rim_hi::FT
    logρₐ_lo::FT
    logρₐ_hi::FT
    logx̄l_lo::FT
    logx̄l_hi::FT
    invT_lo::FT
    invT_hi::FT
end

Adapt.adapt_structure(to, t::LiqIceBrimTable) = LiqIceBrimTable(
    Adapt.adapt(to, t.logB),
    t.x̄l_slope,
    t.logx̄_lo, t.logx̄_hi,
    t.F_rim_lo, t.F_rim_hi,
    t.ρ_rim_lo, t.ρ_rim_hi,
    t.logρₐ_lo, t.logρₐ_hi,
    t.logx̄l_lo, t.logx̄l_hi,
    t.invT_lo, t.invT_hi,
)

"""
    lookup(tab::LiqIceBrimTable, x̄, F_rim, ρ_rim, ρₐ, x̄_liq, T°C)

Six-linear interpolation of the stored logarithm, returned exponentiated. Sixty-four corners.
The temperature enters as `1/T°C`, which is the coordinate the local rime density is a function of.
"""
@inline function lookup(tab::LiqIceBrimTable{FT}, x̄, F_rim, ρ_rim, ρₐ, x̄_liq, T°C) where {FT}
    nx, nf, nr, na, nl, nt = size(tab.logB)
    (ix, ax) = _uniform_index(log10(x̄), tab.logx̄_lo, tab.logx̄_hi, nx)
    (jf, af) = _uniform_index(F_rim, tab.F_rim_lo, tab.F_rim_hi, nf)
    (kr, ar) = _uniform_index(ρ_rim, tab.ρ_rim_lo, tab.ρ_rim_hi, nr)
    (ma, aa) = _uniform_index(log10(ρₐ), tab.logρₐ_lo, tab.logρₐ_hi, na)
    (pl, al) = _uniform_index(log10(x̄_liq), tab.logx̄l_lo, tab.logx̄l_hi, nl)
    (st, at) = _uniform_index(_logT_coord(T°C),
        _logT_of_invT(tab.invT_lo), _logT_of_invT(tab.invT_hi), nt)
    # added back with the UNCLAMPED coordinate, as in `lookup(::LiqIceTable, ...)`
    shear = tab.x̄l_slope * (log10(x̄_liq) - tab.logx̄l_lo) * log(FT(10))
    v = zero(FT)
    @inbounds for di in 0:1, dj in 0:1, dk in 0:1, dm in 0:1, dp in 0:1, ds in 0:1
        w =
            (di == 0 ? 1 - ax : ax) * (dj == 0 ? 1 - af : af) *
            (dk == 0 ? 1 - ar : ar) * (dm == 0 ? 1 - aa : aa) *
            (dp == 0 ? 1 - al : al) * (ds == 0 ? 1 - at : at)
        v += w * tab.logB[ix + di, jf + dj, kr + dk, ma + dm, pl + dp, st + ds]
    end
    exp(v + shear)
end

@inline grid_length(g::LiqIceTableGrid) = g.nx * g.nf * g.nr * g.na * g.nl
@inline grid_length(g::LiqIceBrimTableGrid) = g.nx * g.nf * g.nr * g.na * g.nl * g.nt

@inline function _axes_at(grid::LiqIceTableGrid, lin::Integer)
    i = Int(lin) - 1
    ix = i % grid.nx
    i ÷= grid.nx
    jf = i % grid.nf
    i ÷= grid.nf
    kr = i % grid.nr
    i ÷= grid.nr
    ma = i % grid.na
    i ÷= grid.na
    pl = i
    (_axis_value(grid.logx̄_lo, grid.logx̄_hi, grid.nx, ix),
        _axis_value(grid.F_rim_lo, grid.F_rim_hi, grid.nf, jf),
        _axis_value(grid.ρ_rim_lo, grid.ρ_rim_hi, grid.nr, kr),
        _axis_value(grid.logρₐ_lo, grid.logρₐ_hi, grid.na, ma),
        _axis_value(grid.logx̄l_lo, grid.logx̄l_hi, grid.nl, pl))
end

@inline function _axes_at(grid::LiqIceBrimTableGrid, lin::Integer)
    i = Int(lin) - 1
    ix = i % grid.nx
    i ÷= grid.nx
    jf = i % grid.nf
    i ÷= grid.nf
    kr = i % grid.nr
    i ÷= grid.nr
    ma = i % grid.na
    i ÷= grid.na
    pl = i % grid.nl
    i ÷= grid.nl
    st = i
    (_axis_value(grid.logx̄_lo, grid.logx̄_hi, grid.nx, ix),
        _axis_value(grid.F_rim_lo, grid.F_rim_hi, grid.nf, jf),
        _axis_value(grid.ρ_rim_lo, grid.ρ_rim_hi, grid.nr, kr),
        _axis_value(grid.logρₐ_lo, grid.logρₐ_hi, grid.na, ma),
        _axis_value(grid.logx̄l_lo, grid.logx̄l_hi, grid.nl, pl),
        _logT_axis_value(grid.invT_lo, grid.invT_hi, grid.nt, st))
end

"""
    liqice_table_entry(scheme, vel, quad, psd_liq, m_liq, logx̄, F_rim, ρ_rim, logρₐ, logx̄l)

One grid entry of a [`LiqIceTable`](@ref): the partition-free collected mass and number rates of one
liquid species, per unit ice number and per unit liquid number.

The rates are evaluated through the same inner and outer integrals the collision entry uses, so the
table cannot drift from the physics it stands in for. The partition is absent by construction rather
than by a choice made here: `∫liquid_ice_collisions_split` applies a CONSTANT partition to the full
range and restores the true one with its correction, and a constant multiplies these integrals
rather than entering them.

# Why the reference numbers are 1

Both integrals are first-degree homogeneous in each number concentration separately: the ice size
distribution is linear in `ρn_ice` at fixed mean mass, and the liquid one is linear in `N_liq` at
fixed mean mass. Evaluating at unit numbers and dividing by their product therefore gives the same
kernel any other pair would, and neither number is an axis.

# Why the rime density is a constant here

The third inner integral divides by the local rime density and is not used. Passing a constant in
its place keeps the temperature out of this entry entirely, which is what makes the result a
function of five arguments. The rime volume rate is [`liqice_brim_table_entry`](@ref)'s.
"""
function liqice_table_entry(scheme, vel, quad, psd_liq, m_liq,
    logx̄, F_rim, ρ_rim, logρₐ, logx̄l)
    FT = typeof(logx̄)
    ρₐ = exp10(logρₐ)
    ρn_ice = one(FT)
    ρq_ice = exp10(logx̄) * ρn_ice
    ρq_rim = F_rim * ρq_ice
    ρb_rim = ρq_rim / ρ_rim
    state = state_from_prognostic(scheme, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    logλ = get_distribution_logλ(state)
    # The stored value is per unit ice number and per unit liquid number, and the integrand is
    # exactly linear in the liquid number, so which number the entry evaluates AT is free in exact
    # arithmetic. It is not free in Float32. At one drop per cubic metre the specific content is
    # four to five decades below anything a run carries, and the rain distribution's own size
    # bounds collapse to an empty interval there, which made the entry exactly zero over 48.4
    # percent of the production grid at Float32 and over 0.0 percent of the same grid at Float64.
    # So the entry evaluates at a representative CONTENT and divides by the number that produces
    # it, which is the same value in exact arithmetic and keeps every evaluation inside the range
    # the working precision handles.
    q_ref = FT(1e-4)
    N_liq = ρₐ * q_ref / exp10(logx̄l)
    L_liq = exp10(logx̄l) * N_liq
    p = FT(0.00001)
    n_liq = DT.size_distribution(psd_liq, L_liq / ρₐ, ρₐ, N_liq)
    bounds_liq = CM2.get_size_distribution_bounds(psd_liq, L_liq / ρₐ, ρₐ, N_liq, p)
    ∂ₜV = volumetric_collision_rate_integrand(vel, ρₐ, state)
    ice_bounds = velocity_integral_bounds(state, logλ, ∂ₜV.v_i; p)
    n_i = DT.size_distribution(state, logλ)
    # The inner integral is written out INSIDE the outer integrand rather than taken from the
    # closure `get_liquid_integrals` returns. That closure captures the size distribution, the
    # velocity functor, the local rime density and the liquid bounds, and calling it from inside
    # this entry's own quadrature materialises it on the heap: this entry allocated 2784 bytes per
    # call and the rime-volume one 5264, while every component of both allocated zero, and the
    # device compiler refused them with `gpu_gc_pool_alloc` so a table could not be filled where it
    # lives. Written out, both allocate NOTHING and both compile for a GPU.
    #
    # A closure is not the problem; a closure crossing a call boundary is. The tuple utilities are
    # kept exactly as they were, because replacing them with loops indexed into the bounds tuple
    # costs far more than the closure did.
    #
    # The rime volume channel is not formed here at all, where before it was computed through a
    # constant local rime density and discarded.
    outer = Dᵢ -> begin
        inner_integrand = D -> begin
            ∂ₜV_D = ∂ₜV(Dᵢ, D)
            ∂ₜN = ∂ₜV_D * n_liq(D)
            SA.SVector(∂ₜN, ∂ₜN * m_liq(D))
        end
        bnds = crossing_integral_bounds(bounds_liq, ∂ₜV, Dᵢ)
        (∂ₜN_col, ∂ₜM_col) = integrate(inner_integrand, bnds, quad)
        n = n_i(Dᵢ)
        SA.SVector(n * ∂ₜM_col, n * ∂ₜN_col)
    end
    (M, N) = integrate(outer, ice_bounds, quad)
    return (M / (ρn_ice * N_liq), N / (ρn_ice * N_liq))
end

"""
    liqice_brim_table_entry(scheme, vel, quad, psd_liq, m_liq, T_freeze,
                            logx̄, F_rim, ρ_rim, logρₐ, logx̄l, invT)

One grid entry of a [`LiqIceBrimTable`](@ref): the partition-free rime volume rate of one liquid
species, per unit ice number and per unit liquid number.

The sixth argument is `1/T°C` rather than a temperature, because that is the coordinate the local
rime density is a function of. The air temperature handed to
[`compute_local_rime_density`](@ref) is reconstructed from it.
"""
function liqice_brim_table_entry(scheme, vel, quad, psd_liq, m_liq, T_freeze,
    logx̄, F_rim, ρ_rim, logρₐ, logx̄l, invT)
    FT = typeof(logx̄)
    ρₐ = exp10(logρₐ)
    ρn_ice = one(FT)
    ρq_ice = exp10(logx̄) * ρn_ice
    ρq_rim = F_rim * ρq_ice
    ρb_rim = ρq_rim / ρ_rim
    state = state_from_prognostic(scheme, ρq_ice, ρn_ice, ρq_rim, ρb_rim)
    logλ = get_distribution_logλ(state)
    # The stored value is per unit ice number and per unit liquid number, and the integrand is
    # exactly linear in the liquid number, so which number the entry evaluates AT is free in exact
    # arithmetic. It is not free in Float32. At one drop per cubic metre the specific content is
    # four to five decades below anything a run carries, and the rain distribution's own size
    # bounds collapse to an empty interval there, which made the entry exactly zero over 48.4
    # percent of the production grid at Float32 and over 0.0 percent of the same grid at Float64.
    # So the entry evaluates at a representative CONTENT and divides by the number that produces
    # it, which is the same value in exact arithmetic and keeps every evaluation inside the range
    # the working precision handles.
    q_ref = FT(1e-4)
    N_liq = ρₐ * q_ref / exp10(logx̄l)
    L_liq = exp10(logx̄l) * N_liq
    p = FT(0.00001)
    T = inv(invT) + T_freeze
    n_liq = DT.size_distribution(psd_liq, L_liq / ρₐ, ρₐ, N_liq)
    bounds_liq = CM2.get_size_distribution_bounds(psd_liq, L_liq / ρₐ, ρₐ, N_liq, p)
    ∂ₜV = volumetric_collision_rate_integrand(vel, ρₐ, state)
    ice_bounds = velocity_integral_bounds(state, logλ, ∂ₜV.v_i; p)
    n_i = DT.size_distribution(state, logλ)
    ρ′_rim = compute_local_rime_density(vel, ρₐ, T, state)
    # Written out inside the outer integrand, for the reason `liqice_table_entry` gives. Only the
    # rime volume channel is formed, since that is all this entry stores.
    outer = Dᵢ -> begin
        v_i_at_Dᵢ = ∂ₜV.v_i(Dᵢ)
        coeffs = collision_cross_section_ice_liquid_coeffs(∂ₜV.state, Dᵢ)
        inner_integrand = D -> begin
            v_term = abs(v_i_at_Dᵢ - ∂ₜV.v_l(D))
            E = one(v_term)  # TODO - Make collision efficiency a function of Dᵢ and Dₗ
            ∂ₜN = E * evalpoly(D, coeffs) * v_term * n_liq(D)
            ∂ₜM = ∂ₜN * m_liq(D)
            return ∂ₜM / rime_density_at(ρ′_rim, v_term, D)
        end
        bnds = crossing_integral_bounds(bounds_liq, ∂ₜV, Dᵢ, v_i_at_Dᵢ)
        return n_i(Dᵢ) * integrate(inner_integrand, bnds, quad)
    end
    B = integrate(outer, ice_bounds, quad)
    return B / (ρn_ice * N_liq)
end

"""
    liqice_table(scheme, vel, grid, psd_liq, m_liq; quad)

Fill a [`LiqIceTable`](@ref) over `grid`, on the HOST with `Threads.@threads`, which is where and
how the four-axis tables are filled. The per-entry function is the seam: a caller that would rather
fill the table where it lives calls [`liqice_table_entry`](@ref) itself and passes the result to the
array constructor, and both routes evaluate the same function, so this fill is the reference the
other would be measured against rather than a second implementation.

Every entry must be strictly positive, because the storage is logarithmic. A non-positive or
non-finite entry is refused with the grid point that produced it rather than replaced by a sentinel:
`exp` of a sentinel is a finite, plausible, wrong answer that passes every `isfinite` check, which
is the same failure the adapt rule exists to prevent at the other end.
"""
function liqice_table(scheme, vel, grid::LiqIceTableGrid{FT}, psd_liq, m_liq; quad) where {FT}
    logM = Array{FT}(undef, grid_size(grid))
    logN = similar(logM)
    # The stored logs are sheared by the derived `x̄_liq` exponent, so that the region below the
    # liquid mean-mass clamp is flat in the stored quantity and the lookup's clamp there returns a
    # linear extrapolation rather than a constant; see `_liqice_x̄l_slope`.
    α = _liqice_x̄l_slope(grid,
        (lx, fr, rr, la, u) ->
            liqice_table_entry(scheme, vel, quad, psd_liq, m_liq, lx, fr, rr, la, u)[1])
    ln10 = log(FT(10))
    Threads.@threads for lin in 1:grid_length(grid)
        pt = _axes_at(grid, lin)
        # Destructured rather than splatted. `f(..., pt...)` splats a tuple into the call, which
        # costs a measured 1152 bytes per entry where passing the components costs nothing, and a
        # fill runs this line once per grid point.
        (Mv, Nv) = liqice_table_entry(scheme, vel, quad, psd_liq, m_liq,
            pt[1], pt[2], pt[3], pt[4], pt[5])
        _refuse_nonpositive(Mv, "M", pt)
        _refuse_nonpositive(Nv, "N", pt)
        sh = α * (pt[5] - grid.logx̄l_lo) * ln10
        @inbounds logM[lin] = log(Mv) - sh
        @inbounds logN[lin] = log(Nv) - sh
    end
    LiqIceTable(logM, logN, α,
        grid.logx̄_lo, grid.logx̄_hi, grid.F_rim_lo, grid.F_rim_hi,
        grid.ρ_rim_lo, grid.ρ_rim_hi, grid.logρₐ_lo, grid.logρₐ_hi,
        grid.logx̄l_lo, grid.logx̄l_hi)
end

"""
    liqice_brim_table(scheme, vel, grid, psd_liq, m_liq, T_freeze; quad)

Fill a [`LiqIceBrimTable`](@ref) over `grid`, on the same host-threaded fill and the same positivity
condition as [`liqice_table`](@ref).
"""
function liqice_brim_table(scheme, vel, grid::LiqIceBrimTableGrid{FT}, psd_liq, m_liq, T_freeze;
    quad) where {FT}
    logB = Array{FT}(undef, grid_size(grid))
    # The same shear as `liqice_table`, derived from THIS table's own channel rather than from the
    # collected-mass one. brim-shear-exponent-transfers (39) measured the two to agree to 3.1e-15
    # at six temperatures, which is what linearity predicts, but the three-point disagreement test
    # inside `_liqice_x̄l_slope` only ever saw the channel it was handed, so deriving it from the
    # mass channel left the agreement asserted rather than checked. Taking it at both ends of the
    # temperature axis and requiring them to match checks the temperature independence too, which
    # the mass channel cannot see at all because it has no temperature axis.
    brim_at(iT) =
        (lx, fr, rr, la, u) ->
            liqice_brim_table_entry(scheme, vel, quad, psd_liq, m_liq, T_freeze,
                lx, fr, rr, la, u, iT)
    α = _liqice_x̄l_slope(grid, brim_at(grid.invT_lo))
    α_hi = _liqice_x̄l_slope(grid, brim_at(grid.invT_hi))
    α == α_hi || error(
        "the rime-volume liquid mean-mass exponent depends on temperature: $α at invT_lo " *
        "$(grid.invT_lo) against $α_hi at invT_hi $(grid.invT_hi). A single scalar shear cannot " *
        "represent this table.")
    ln10 = log(FT(10))
    Threads.@threads for lin in 1:grid_length(grid)
        pt = _axes_at(grid, lin)
        # Destructured rather than splatted, for the reason `liqice_table` gives.
        Bv = liqice_brim_table_entry(scheme, vel, quad, psd_liq, m_liq, T_freeze,
            pt[1], pt[2], pt[3], pt[4], pt[5], pt[6])
        _refuse_nonpositive(Bv, "B", pt)
        @inbounds logB[lin] = log(Bv) - α * (pt[5] - grid.logx̄l_lo) * ln10
    end
    LiqIceBrimTable(logB, α,
        grid.logx̄_lo, grid.logx̄_hi, grid.F_rim_lo, grid.F_rim_hi,
        grid.ρ_rim_lo, grid.ρ_rim_hi, grid.logρₐ_lo, grid.logρₐ_hi,
        grid.logx̄l_lo, grid.logx̄l_hi, grid.invT_lo, grid.invT_hi)
end

@inline function _refuse_nonpositive(v, name, point)
    (isfinite(v) && v > 0) ||
        error(
            "liquid-ice table entry $name is $v at grid point $point; the storage is " *
            "logarithmic and a sentinel there would read back as a finite wrong answer",
        )
    return nothing
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
struct P3TabulatedQuadrature{Q, S, V, N, M, CC, CR, BC, BR} <: QuadratureRule
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
    "the cloud liquid-ice collected mass and number, or `nothing`"
    liqice_col_cloud::CC
    "the rain liquid-ice collected mass and number, or `nothing`"
    liqice_col_rain::CR
    "the cloud liquid-ice rime volume, or `nothing`"
    liqice_brim_cloud::BC
    "the rain liquid-ice rime volume, or `nothing`"
    liqice_brim_rain::BR
    "mirrors `rule.n`, read directly by the quadrature sums in `Quadrature.integrate`"
    n::Int
end

P3TabulatedQuadrature(rule; selfcol = nothing, vent = nothing, vel_n = nothing, vel_m = nothing,
    liqice_col_cloud = nothing, liqice_col_rain = nothing,
    liqice_brim_cloud = nothing, liqice_brim_rain = nothing) =
    P3TabulatedQuadrature(rule, selfcol, vent, vel_n, vel_m,
        liqice_col_cloud, liqice_col_rain, liqice_brim_cloud, liqice_brim_rain, rule.n)

# Every field is adapted; an absent table is `nothing`, which `Adapt.adapt` returns unchanged,
# so this one method covers every combination of adopted tables, including the all-`nothing`
# carrier that behaves exactly like the rule it carries.
Adapt.adapt_structure(to, q::P3TabulatedQuadrature) = P3TabulatedQuadrature(
    Adapt.adapt(to, q.rule),
    Adapt.adapt(to, q.selfcol),
    Adapt.adapt(to, q.vent),
    Adapt.adapt(to, q.vel_n),
    Adapt.adapt(to, q.vel_m),
    Adapt.adapt(to, q.liqice_col_cloud),
    Adapt.adapt(to, q.liqice_col_rain),
    Adapt.adapt(to, q.liqice_brim_cloud),
    Adapt.adapt(to, q.liqice_brim_rain),
    q.n,
)

"""
    P3_LUT_OUTPUTS

The terms [`p3_lut_carrier`](@ref) can tabulate. A name is a PROCESS, not a stored array: `:liqice`
names the liquid-ice collision entry, which needs four table objects holding six arrays, because
the collected mass and number are temperature-free while the rime volume is not, and cloud and rain
have different size distributions over abutting mass ranges.

The four are selected together and cannot be selected apart. Their six integrals come from one
interleaved integrand, so reading some from tables and integrating for the rest would run the same
quadrature and save nothing. The individual constructors [`liqice_table`](@ref) and
[`liqice_brim_table`](@ref) remain available for building and checking one table at a time.
"""
const P3_LUT_OUTPUTS = (:selfcol, :vent, :vel_n, :vel_m, :liqice)

"""
    p3_lut_carrier(scheme, vel, aps, grid; quad, generating_quad, outputs,
                   liqice_grid_cloud, liqice_grid_rain,
                   liqice_brim_grid_cloud, liqice_brim_grid_rain,
                   psd_c, psd_r, m_liq, T_freeze)

Build a [`P3TabulatedQuadrature`](@ref) in one call: fill the requested tables on `grid` and
wrap them around `quad`, which remains the rule for every term not requested.

`outputs` states which terms to tabulate, from [`P3_LUT_OUTPUTS`](@ref); a term left out is not
tabulated and its field stays `nothing`. `:liqice` needs one grid per table object, both liquid size
distributions, `m_liq` and `T_freeze`, none of which this constructor can derive from the ice grid,
and each is refused by name if absent.

`quad` is the rule the carrier DELEGATES to for every term it does not tabulate, and
`generating_quad` is the rule the tables are generated FROM. They default to the same rule and
should not be: a generator wants a higher order than production, so that its own quadrature error
stays below the interpolation budget the table is measured against, while an untabulated term must
keep integrating at production's order. Passing one rule for both silently moves every untabulated
term onto the generator's order. Device residence is the caller's responsibility, as for a single table: pass tables
built from device arrays, or convert them after filling. `CloudMicrophysics` performs no device
execution itself.
"""


function p3_lut_carrier(scheme, vel, aps, grid::IceKernelTableGrid; quad,
    generating_quad = quad, outputs = (:selfcol, :vent, :vel_n, :vel_m),
    liqice_grid_cloud = nothing, liqice_grid_rain = nothing,
    liqice_brim_grid_cloud = nothing, liqice_brim_grid_rain = nothing,
    psd_c = nothing, psd_r = nothing, m_liq = nothing, T_freeze = nothing)
    for o in outputs
        o in P3_LUT_OUTPUTS || error(
            "unknown P3 lookup-table output $o; expected one of $(join(P3_LUT_OUTPUTS, ", "))")
    end
    # Each liquid-ice output needs a grid this constructor cannot derive from the ice grid, so a
    # missing one is refused by name rather than allowed to build a table on a default nobody
    # chose.
    _need(o, g, what) = (o in outputs && g === nothing) &&
                        error("P3 lookup-table output $o requires `$what`")
    # One grid per table object, never one per family. Cloud and rain occupy abutting, DISJOINT
    # liquid mass ranges - `xc_max` is `xr_min` - so a grid shared between them tabulates one
    # species entirely outside its own range, where every lookup clamps to an edge. The rime
    # density also saturates at very different temperatures for the two, so the brim grids differ
    # in their inverse-temperature axis as well.
    _need(:liqice, liqice_grid_cloud, "liqice_grid_cloud")
    _need(:liqice, liqice_grid_rain, "liqice_grid_rain")
    _need(:liqice, liqice_brim_grid_cloud, "liqice_brim_grid_cloud")
    _need(:liqice, liqice_brim_grid_rain, "liqice_brim_grid_rain")
    _need(:liqice, psd_c, "psd_c")
    _need(:liqice, psd_r, "psd_r")
    _need(:liqice, m_liq, "m_liq")
    _need(:liqice, T_freeze, "T_freeze")
    P3TabulatedQuadrature(quad;
        selfcol = :selfcol in outputs ? IceKernelTable(scheme, vel, grid; quad = generating_quad) : nothing,
        vent = :vent in outputs ? ice_ventilation_table(scheme, vel, aps, grid; quad = generating_quad) : nothing,
        vel_n = :vel_n in outputs ? ice_velocity_n_table(scheme, vel, grid; quad = generating_quad) : nothing,
        vel_m = :vel_m in outputs ? ice_velocity_m_table(scheme, vel, grid; quad = generating_quad) : nothing,
        liqice_col_cloud = :liqice in outputs ?
                           liqice_table(scheme, vel, liqice_grid_cloud, psd_c, m_liq; quad = generating_quad) : nothing,
        liqice_col_rain = :liqice in outputs ?
                          liqice_table(scheme, vel, liqice_grid_rain, psd_r, m_liq; quad = generating_quad) : nothing,
        liqice_brim_cloud = :liqice in outputs ?
                            liqice_brim_table(
            scheme,
            vel,
            liqice_brim_grid_cloud,
            psd_c,
            m_liq,
            T_freeze;
            quad = generating_quad,
        ) : nothing,
        liqice_brim_rain = :liqice in outputs ?
                           liqice_brim_table(
            scheme,
            vel,
            liqice_brim_grid_rain,
            psd_r,
            m_liq,
            T_freeze;
            quad = generating_quad,
        ) : nothing)
end


"""
    _liqice_full_range(quad, state, ρₐ, psd_c, L_c, N_c, psd_r, L_r, N_r, T, ∂ₜM_max, ice_bounds)

The nine channels of the split assembly's FULL-RANGE integral, read from the liquid-ice tables, or
`nothing` when they are not all present.

`∫liquid_ice_collisions_split` applies a constant partition to the full range and restores the true
one with a correction on the wet set, so the full range is the only part of the assembly a table can
hold. With that partition constant, every channel is a stored integral times a constant:

    QCFRZ = M_c f      QCSHD = M_c (1-f)     NCCOL = N_c
    QRFRZ = M_r f      QRSHD = M_r (1-f)     NRCOL = N_r
    ∫M_col = M_c + M_r     BCCOL = B_c f     BRCOL = B_r f

with `f` one below freezing and zero at and above it, decided by the same test the split makes.

The stored integrals are per unit ice number and per unit liquid number, both being first-degree
homogeneous, so each is multiplied by `ρn_ice` and its own species' number.

All four tables are required together. The six integrals come from one interleaved integrand, so
reading some from tables and integrating for the rest would run the same quadrature and save
nothing; a partial set therefore leaves the entry on its quadrature path.
"""
@inline _liqice_full_range(quad, state, ρₐ, psd_c, L_c, N_c, psd_r, L_r, N_r, T, ∂ₜM_max,
    ice_bounds) = nothing

@inline function _liqice_full_range(q::P3TabulatedQuadrature, state, ρₐ, psd_c, L_c, N_c,
    psd_r, L_r, N_r, T, ∂ₜM_max, ice_bounds)
    (
        q.liqice_col_cloud === nothing || q.liqice_col_rain === nothing ||
        q.liqice_brim_cloud === nothing || q.liqice_brim_rain === nothing
    ) && return nothing
    (; ρq_ice, ρn_ice) = state
    ((ρq_ice > 0) & (ρn_ice > 0)) || return nothing
    FT = typeof(ρₐ)
    x̄ = ρq_ice / ρn_ice
    x̄_c = L_c / N_c
    x̄_r = L_r / N_r
    T°C = T - state.params.T_freeze
    (Mc, Nc) = lookup(q.liqice_col_cloud, x̄, state.F_rim, state.ρ_rim, ρₐ, x̄_c)
    (Mr, Nr) = lookup(q.liqice_col_rain, x̄, state.F_rim, state.ρ_rim, ρₐ, x̄_r)
    Bc = lookup(q.liqice_brim_cloud, x̄, state.F_rim, state.ρ_rim, ρₐ, x̄_c, T°C)
    Br = lookup(q.liqice_brim_rain, x̄, state.F_rim, state.ρ_rim, ρₐ, x̄_r, T°C)
    sc = ρn_ice * N_c
    sr = ρn_ice * N_r
    M_c = Mc * sc
    N_c_col = Nc * sc
    M_r = Mr * sr
    N_r_col = Nr * sr
    B_c = Bc * sc
    B_r = Br * sr
    # The same scalar test the split makes, so the two paths cannot disagree about which constant
    # partition the full range carries.
    D_lo, D_hi = first(ice_bounds), last(ice_bounds)
    f = iszero(FD.value(∂ₜM_max(sqrt(D_lo * D_hi)))) ? zero(FT) : one(FT)
    return SA.SVector(
        M_c * f, M_c * (1 - f), N_c_col,
        M_r * f, M_r * (1 - f), N_r_col,
        M_c + M_r, B_c * f, B_r * f,
    )
end


"""
    _default_liqice_assembly(quad)

Which assembly [`∫liquid_ice_collisions`](@ref) uses by default.

`PartitionedOuter` evaluates the freeze/shed partition inside the outer integral and so has no
full-range integral for a table to stand in for. A carrier holding the liquid-ice tables is the
statement that the term is tabulated, and `SplitCorrection` is the only assembly that can read
them, so the carrier selects it. A carrier without them keeps the default, and so does a plain
rule.
"""
@inline _default_liqice_assembly(quad) = PartitionedOuter()
@inline _default_liqice_assembly(q::P3TabulatedQuadrature) =
    (
        q.liqice_col_cloud === nothing || q.liqice_col_rain === nothing ||
        q.liqice_brim_cloud === nothing || q.liqice_brim_rain === nothing
    ) ?
    PartitionedOuter() : SplitCorrection()

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
