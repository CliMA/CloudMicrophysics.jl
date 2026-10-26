"""
    TableAxis{FT}

A uniform axis of `n` points from `lo` to `hi`, with coordinates of the float type `FT`.

# Construction

    TableAxis(lo, hi, n)

with `FT` the float type that `lo` and `hi` promote to.

# Fields
$(FIELDS)
"""
struct TableAxis{FT}
    "Coordinate of the first point"
    lo::FT
    "Coordinate of the last point"
    hi::FT
    "Number of points, at least two"
    n::Int
    function TableAxis(lo, hi, n)
        n ≥ 2 || throw(ArgumentError("a table axis needs at least two points, got $n"))
        lo < hi || throw(ArgumentError("a table axis needs lo < hi, got lo = $lo, hi = $hi"))
        (lo, hi) = float.(promote(lo, hi))
        return new{typeof(lo)}(lo, hi, Int(n))
    end
end

@inline axis_point((; lo, hi, n)::TableAxis, k) = lo + (hi - lo) * (k - 1) / (n - 1)

axis_points(axis) = [axis_point(axis, k) for k in 1:(axis.n)]

"""
    grid_point(axes, I)

Return the coordinates of the point with the `CartesianIndex` `I` of the grid spanned by the [`TableAxis`](@ref) objects in `axes`, a tuple in the float types of the axes.
"""
@inline grid_point(axes, I) = ntuple(d -> axis_point(axes[d], I[d]), Val(length(axes)))

struct GridPointFunction{F, A}
    f::F
    axes::A
end
@inline (g::GridPointFunction)(I) = g.f(grid_point(g.axes, I))

"""
    tabulate(f, axes, ArrayType = Array)

Return the array of the values of the point function `f` at the points of the grid spanned by `axes`.

The element at the index `I` is `f` at the coordinates [`grid_point`](@ref)`(axes, I)`, and the
array has one dimension per [`TableAxis`](@ref) of `axes`. For an `Array`, the points are
evaluated by tasks on all threads, in chunks that are scheduled dynamically. For any other array
type, the values are filled by a broadcast over the indices into `ArrayType{T}(undef, dims)`, so
the points are evaluated where the array type places them.

# Arguments
- `f`: point function, which maps an `NTuple` of coordinates to a number or an `SVector`
- `axes`: `NTuple` of [`TableAxis`](@ref)
- `ArrayType`: array type of the result, such as `Array`
"""
function tabulate(f, axes, ArrayType = Array)
    point = GridPointFunction(f, axes)
    return evaluate_points(point, CartesianIndices(map(a -> a.n, axes)), ArrayType)
end

function evaluate_points(point, indices, ::Type{Array})
    values = Array{typeof(point(first(indices)))}(undef, size(indices))
    chunk = cld(length(indices), 16 * Threads.nthreads())
    @sync for linear in Iterators.partition(1:length(indices), chunk)
        Threads.@spawn for i in linear
            @inbounds values[i] = point(indices[i])
        end
    end
    return values
end
function evaluate_points(point, indices, ArrayType)
    values = ArrayType{typeof(point(first(indices)))}(undef, size(indices))
    values .= point.(indices)
    return values
end

"""
    UniformTable{FT, N, K, A}

A function of `N` coordinates with `K` components, stored on a grid that is uniform along each axis.

[`lookup`](@ref) reads all components at once by multilinear interpolation.

# Construction

    UniformTable(values, axes, FT; logs = nothing, inner = false)

stores the array `values` of a point function at the points of the grid spanned by the `N`
[`TableAxis`](@ref) objects in `axes` (see [`tabulate`](@ref)) as `FT`, in an array of the type of
`values`. Each element of `values` is a number or an `SVector` of `K` components. `logs` is an
`NTuple` of `Bool`, whether each component is stored as its logarithm, all `false` by default.
With `inner = true`, `values` spans the axes after the first, and each element holds the values
at the points of the first axis, so the table has one component. The construction throws a
`DomainError` for a value that is not finite or a log component that is not positive.

# Fields
$(FIELDS)
"""
struct UniformTable{FT, N, K, A <: AbstractArray{FT}}
    "Stored values, of size `(K, n₁, …, n_N)` for axes of `n₁` to `n_N` points"
    data::A
    "Coordinate of the first point of each axis"
    lo::NTuple{N, FT}
    "Coordinate of the last point of each axis"
    hi::NTuple{N, FT}
end

UniformTable{FT, N, K}(data, lo, hi) where {FT, N, K} = UniformTable{FT, N, K, typeof(data)}(data, lo, hi)

Adapt.adapt_structure(to, t::UniformTable{FT, N, K}) where {FT, N, K} =
    UniformTable{FT, N, K}(Adapt.adapt(to, t.data), t.lo, t.hi)

function UniformTable(values, axes, ::Type{FT}; logs = nothing, inner = false) where {FT}
    check_table_size(values, inner ? Base.tail(axes) : axes)
    K = component_count(eltype(values))
    inner && K != first(axes).n &&
        throw(DimensionMismatch("$K values per point on an inner axis of $(first(axes).n) points"))
    data = component_array(values, K)
    data = stored_components(inner ? reshape(data, 1, size(data)...) : data, logs)
    return UniformTable{FT, length(axes), inner ? 1 : K}(
        device_array(values, FT.(data)), map(a -> FT(a.lo), axes), map(a -> FT(a.hi), axes),
    )
end

component_count(::Type{<:Real}) = 1
component_count(::Type{<:SA.SVector{K}}) where {K} = K

check_table_size(values, axes) =
    size(values) == map(a -> a.n, axes) ||
    throw(DimensionMismatch("table values of size $(size(values)) on a grid of size $(map(a -> a.n, axes))"))

function component_array(values, K)
    host = host_array(values)
    data = Array{float(eltype(eltype(host)))}(undef, K, size(host)...)
    for I in CartesianIndices(host), k in 1:K
        data[k, I] = host[I][k]
    end
    return data
end

host_array(values::Array) = values
host_array(values) = Array(values)

function stored_components(data, logs)
    for I in CartesianIndices(Base.tail(size(data))), k in axes(data, 1)
        v = data[k, I]
        isfinite(v) || throw(DomainError(v, "a table value must be finite"))
        if !isnothing(logs) && logs[k]
            v > 0 || throw(DomainError(v, "a log component of a table must be positive"))
            data[k, I] = log(v)
        end
    end
    return data
end

device_array(::Array, data) = data
device_array(values, data) = copyto!(similar(values, eltype(data), size(data)), data)

Base.eltype(::UniformTable{FT}) where {FT} = FT
Base.show(io::IO, t::UniformTable{FT, N, K}) where {FT, N, K} =
    print(io, "UniformTable{", FT, "}(", join(Base.tail(size(t.data)), "×"), K == 1 ? "" : ", $K components", ")")

"""
    primal_value(x)

Return the value of `x` with the partials of every nesting level of `ForwardDiff.Dual` dropped.
"""
@inline primal_value(x) = x
@inline primal_value(x::FD.Dual) = primal_value(FD.value(x))

"""
    uniform_index(x, lo, hi, n)

Return the index `i` of the lower point of the cell of a uniform axis that holds `x`, and the fraction of the cell below `x`.

The axis has `n` points from `lo` to `hi`. A coordinate outside the axis, including `±Inf`, is
moved to the nearest end, and `NaN` to the first point. The index follows from the primal value
of `x` (see [`primal_value`](@ref)), and the fraction has the type of `x`.
"""
@inline function uniform_index(x, lo, hi, n)
    t = (x - lo) / (hi - lo) * (n - 1)
    t = ifelse(isnan(t), zero(t), clamp(t, zero(t), oftype(t, n - 1)))
    i = min(unsafe_trunc(Int, primal_value(t)), n - 2)
    return (i + 1, t - i)
end

"""
    axis_cell(t, d, x)

Return the cell of the coordinate `x` on the axis `d` of the [`UniformTable`](@ref) `t`, see [`uniform_index`](@ref).
"""
@inline axis_cell(t, d, x) = uniform_index(x, t.lo[d], t.hi[d], size(t.data, d + 1))

"""
    lookup(t, x...)

Return the `SVector` of the components of the table `t` at the coordinates `x` by multilinear interpolation.

A coordinate outside the grid is moved to the nearest edge, and `NaN` to the first point (see
[`uniform_index`](@ref)). A coordinate with partials, such as a `ForwardDiff.Dual`, gives the
derivative of the interpolation along its axis.

# Arguments
- `t`: [`UniformTable`](@ref)
- `x`: one coordinate per axis of the table
"""
@inline lookup(t::UniformTable{FT, N}, x::Vararg{Real, N}) where {FT, N} =
    interpolate(t, map((d, x_d) -> axis_cell(t, d, x_d), ntuple(identity, Val(N)), x))

"""
    lookup_or(present, default, t, x...)

Return [`lookup`](@ref)`(t, x...)` if `present`, and otherwise, without reading the table, the `SVector` of the same type with every component `default`.
"""
@inline lookup_or(present, default, t::UniformTable{FT, N}, x::Vararg{Real, N}) where {FT, N} =
    interpolate_or(present, default, t, map(@inline((d, x_d) -> axis_cell(t, d, x_d)), ntuple(identity, Val(N)), x))

@inline function interpolate_or(present, default, t::UniformTable{FT, N, K}, cells) where {FT, N, K}
    T = Base.promote_typeof(zero(FT), map(last, cells)...)
    return present ? interpolate(t, cells) : SA.SVector(ntuple(@inline(_ -> T(default)), Val(K)))
end

"""
    interpolate(t, cells)

Return the `SVector` of the components of the [`UniformTable`](@ref) `t` interpolated multilinearly in the cells `cells`, one [`axis_cell`](@ref) per axis.
"""
@inline function interpolate(t::UniformTable{FT, N, K}, cells) where {FT, N, K}
    (; data) = t
    T = Base.promote_typeof(zero(FT), map(last, cells)...)
    # The linear index of the first component at the lower corner, and the stride of each axis
    strides = Base.front(accumulate(*, size(data)))
    base = 1 + sum(ntuple(d -> (first(cells[d]) - 1) * strides[d], Val(N)))
    N == 1 && return interpolate_line(data, base, strides[1], last(cells[1]), Val(K))
    # Interpolate along the first two axes at each corner of the other axes
    a₁, a₂ = last(cells[1]), last(cells[2])
    v = zero(SA.SVector{K, T})
    @inbounds for c in 0:(2 ^ (N - 2) - 1)
        upper = ntuple(d -> isodd(c >> (d - 1)), Val(N - 2))
        offset = base + sum(ntuple(d -> ifelse(upper[d], strides[d + 2], 0), Val(N - 2)); init = 0)
        w = prod(ntuple(d -> ifelse(upper[d], last(cells[d + 2]), 1 - last(cells[d + 2])), Val(N - 2)); init = one(T))
        lower = interpolate_line(data, offset, strides[1], a₁, Val(K))
        higher = interpolate_line(data, offset + strides[2], strides[1], a₁, Val(K))
        v += w * ((1 - a₂) * lower + a₂ * higher)
    end
    return v
end

@inline function interpolate_line(data, offset, stride, a, ::Val{K}) where {K}
    lower = SA.SVector(ntuple(k -> @inbounds(data[offset + k - 1]), Val(K)))
    higher = SA.SVector(ntuple(k -> @inbounds(data[offset + stride + k - 1]), Val(K)))
    return (1 - a) * lower + a * higher
end
