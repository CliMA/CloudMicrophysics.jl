"""
    TableAxis

A uniform axis of `n` points from `lo` to `hi`.

# Construction

    TableAxis(lo, hi, n)

# Fields
$(FIELDS)
"""
struct TableAxis
    "Coordinate of the first point"
    lo::Float64
    "Coordinate of the last point"
    hi::Float64
    "Number of points, at least two"
    n::Int
    function TableAxis(lo, hi, n)
        n ≥ 2 || throw(ArgumentError("a table axis needs at least two points, got $n"))
        lo < hi || throw(ArgumentError("a table axis needs lo < hi, got lo = $lo, hi = $hi"))
        return new(Float64(lo), Float64(hi), Int(n))
    end
end

"""
    axis_point(axis, k)

Return the coordinate of the `k`-th point of the [`TableAxis`](@ref) `axis`.
"""
@inline axis_point((; lo, hi, n)::TableAxis, k) = lo + (hi - lo) * (k - 1) / (n - 1)

"""
    axis_points(axis)

Return the coordinates of the points of the [`TableAxis`](@ref) `axis`.
"""
axis_points(axis) = [axis_point(axis, k) for k in 1:(axis.n)]

"""
    grid_point(axes, I)

Return the coordinates of the point with the `CartesianIndex` `I` of the grid spanned by the [`TableAxis`](@ref) objects in `axes`, an `NTuple` of `Float64`.
"""
@inline grid_point(axes, I) = ntuple(d -> axis_point(axes[d], I[d]), Val(length(axes)))

"""
    GridPointFunction{F, A}

The point function `f` evaluated at the point of the grid spanned by `axes` with a given `CartesianIndex`, see [`grid_point`](@ref).
"""
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
- `f`: point function, which maps an `NTuple` of `Float64` coordinates to a number or an `SVector`
- `axes`: `NTuple` of [`TableAxis`](@ref)
- `ArrayType`: array type of the result, such as `Array`
"""
function tabulate(f, axes, ArrayType = Array)
    point = GridPointFunction(f, axes)
    return evaluate_points(point, CartesianIndices(map(a -> a.n, axes)), ArrayType)
end

"""
    evaluate_points(point, indices, ArrayType)

Return the array of type `ArrayType` of the values of `point` at the `CartesianIndices` `indices`, see [`tabulate`](@ref).
"""
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

    UniformTable(values, axes, FT; logs = nothing)

stores the array `values` of a point function at the points of the grid spanned by the `N`
[`TableAxis`](@ref) objects in `axes` (see [`tabulate`](@ref)) as `FT`, in an array of the type of
`values`. Each element of `values` is a number or an `SVector` of `K` components, and the
components are stored in the first dimension of `data`. `logs` is an `NTuple` of `Bool`, whether
each component is stored as its logarithm, all `false` by default. A lookup interpolates the
logarithm of such a component, which suits a value that scales as a power of the coordinates.
The construction throws a `DomainError` for a value that is not finite or a log component that is
not positive.

    inner_axis_table(values, inner, axes, FT)

stores the values of a point function that evaluates every point of the `inner` axis in one
pass, see [`inner_axis_table`](@ref).

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

function UniformTable(values, axes, ::Type{FT}; logs = nothing) where {FT}
    check_table_size(values, axes)
    K = component_count(eltype(values))
    data = stored_components(component_array(values, K), logs)
    return UniformTable{FT, length(axes), K}(
        device_array(values, FT.(data)), map(a -> FT(a.lo), axes), map(a -> FT(a.hi), axes),
    )
end

"""
    inner_axis_table(values, inner, axes, FT)

Return the [`UniformTable`](@ref) of one component of a point function that evaluates every point of the axis `inner` in one pass.

Each element of the array `values` is an `SVector` of the values at the points of the
[`TableAxis`](@ref) `inner`, at the points of the grid spanned by `axes`. The table is on the axes
`(inner, axes...)`, so the points of `inner` are adjacent in memory.
"""
function inner_axis_table(values, inner, axes, ::Type{FT}) where {FT}
    check_table_size(values, axes)
    K = component_count(eltype(values))
    K == inner.n || throw(DimensionMismatch("$K values per point on an inner axis of $(inner.n) points"))
    data = stored_components(reshape(component_array(values, K), 1, K, size(values)...), nothing)
    all_axes = (inner, axes...)
    return UniformTable{FT, length(all_axes), 1}(
        device_array(values, FT.(data)), map(a -> FT(a.lo), all_axes), map(a -> FT(a.hi), all_axes),
    )
end

"""
    component_count(T)

Return the number of components of a table value of type `T`, a number or an `SVector`.
"""
component_count(::Type{<:Real}) = 1
component_count(::Type{<:SA.SVector{K}}) where {K} = K

"""
    check_table_size(values, axes)

Throw a `DimensionMismatch` unless the array `values` has one element per point of the grid spanned by `axes`.
"""
check_table_size(values, axes) =
    size(values) == map(a -> a.n, axes) ||
    throw(DimensionMismatch("table values of size $(size(values)) on a grid of size $(map(a -> a.n, axes))"))

"""
    component_array(values, K)

Return the `Array` of `Float64` of size `(K, size(values)...)` of the components of the elements of `values`.
"""
function component_array(values, K)
    host = host_array(values)
    data = Array{Float64}(undef, K, size(host)...)
    for I in CartesianIndices(host), k in 1:K
        data[k, I] = host[I][k]
    end
    return data
end

"""
    host_array(values)

Return the array `values` as an `Array`, without a copy for an `Array`.
"""
host_array(values::Array) = values
host_array(values) = Array(values)

"""
    stored_components(data, logs)

Return the components `data` with the components that `logs` marks replaced by their logarithms in place, see [`UniformTable`](@ref).

It throws a `DomainError` for a value that is not finite or a log component that is not positive.
"""
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

"""
    device_array(values, data)

Return the `Array` `data` in an array of the type of `values`.
"""
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
    interpolate(t, ntuple(d -> axis_cell(t, d, x[d]), Val(N)))

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
    v = zero(SA.SVector{K, T})
    @inbounds for c in 0:(2 ^ N - 1)
        upper = ntuple(d -> isodd(c >> (d - 1)), Val(N))
        offset = base + sum(ntuple(d -> ifelse(upper[d], strides[d], 0), Val(N)))
        w = prod(ntuple(d -> ifelse(upper[d], last(cells[d]), 1 - last(cells[d])), Val(N)))
        v += w * SA.SVector(ntuple(k -> data[offset + k - 1], Val(K)))
    end
    return v
end
