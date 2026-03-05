"""
This module contains methods for showing parameter structs in a human-friendly way.

# Main show methods

- [`verbose_show_type_and_fields(io::IO, ::MIME"text/plain", x; kw...)`](@ref)
- [`compact_show_type_and_fields(io::IO, ::MIME"text/plain", x; kw...)`](@ref)
- [`parseable_show_with_fields_no_type_header(io::IO, x; kw...)`](@ref)

The verbose and compact methods display multiline and singleline representations
of the contents of the structs in a human-friendly way. The parseable method
displays a singleline representation of the struct such that it can be copied
and pasted back into the REPL to recreate the struct. See their respective
signatures for additional keyword arguments, which modify how the output is
displayed.

# Use with parameter structs

For convenience all subtypes of [`ParametersType`](@ref CloudMicrophysics.Parameters.ParametersType)
are defined to use:

```julia
Base.show(io::IO, x::ParametersType) =
    ShowMethods.parseable_show_with_fields_no_type_header(io, x; with_module_prefix = false)
Base.show(io::IO, mime::MIME"text/plain", x::ParametersType) =
    ShowMethods.show_type_and_fields(io, mime, x)
```

The second method, [`show_type_and_fields`](@ref), automatically selects between
[`verbose_show_type_and_fields`](@ref) and [`compact_show_type_and_fields`](@ref)
based on the number of fields in the struct (current cutoff is 7 fields).
This minimizes the need to define custom show methods for new parameter structs.

For wrapper-structs (structs that contain other structs, 
e.g. [`Microphysics2MParams`](@ref CloudMicrophysics.Parameters.Microphysics2MParams)),
it may be desireable to define
```julia
Base.show(io::IO, mime::MIME"text/plain", x::WrapperStruct) = 
    ShowMethods.verbose_show_type_and_fields(io, mime, x; with_module_prefix = false)
```
even if the number of fields is small, as the contents of the containing fields
may be large.

# Units

The verbose and compact show methods also display units, if defined.
If the contents of a struct are numbers, it may be useful to define
[`field_units`](@ref) for the struct:
```julia
@kwdef struct MyStruct{FT}
    velocity::FT
    mass::FT
    constant_value::FT
end

Base.show(io::IO, mime::MIME"text/plain", x::MyStruct) = 
    ShowMethods.compact_show_type_and_fields(io, mime, x)

ShowMethods.field_units(x::MyStruct) = (; velocity = "m/s", mass = "kg")
```
If a field is not present in the returned `NamedTuple`, it is assumed to be
dimensionless. In the example above, `constant_value` is dimensionless. When
this struct is displayed in REPL, it will be shown as:
```julia-repl
julia> MyStruct(; velocity = 1.0, mass = 2.0, constant_value = 3.0)
MyStruct(velocity = 1.0 [m/s], mass = 2.0 [kg], constant_value = 3.0 [-])
```
As a reminder, as long as the parameter struct is a subtype of
[`ParametersType`](@ref CloudMicrophysics.Parameters.ParametersType), 
you should not need to define `Base.show` for it.


Further reading: [julia discourse](https://discourse.julialang.org/t/print-vs-two-argument-show-repr-vs-three-argument-show-repr-with-mime-text-plain/117790)
"""
module ShowMethods

public verbose_show_type_and_fields
public compact_show_type_and_fields
public parseable_show_with_fields_no_type_header
public show_type_and_fields


### ----------------------- ###
### Helper methods for show ###
### ----------------------- ###

"""
    verbose_show_type_and_fields(io::IO, ::MIME"text/plain", x; with_module_prefix = false)

Print a verbose representation of the type and fields of `x` to `io`.

# Keyword arguments
- `with_module_prefix = false`: Use module-qualified type name.
    + When true, display e.g. `CloudMicrophysics.Parameters.AirProperties`
    + When false, display e.g. `AirProperties`
"""
function verbose_show_type_and_fields(io::IO, mime::MIME"text/plain", x; with_module_prefix = false)
    compact = get(io, :compact, false)::Bool
    indent = get(io, :indent, "")::String
    typename = _type_name(x; with_module_prefix)
    keys = fieldnames(typeof(x))
    vals = [getfield(x, k) for k in keys]
    units = field_units(x)
    print(io, "$(typename)")
    if compact
        print(io, "(")
        for (k, v) in zip(keys, vals)
            print(io, "$k = ", v, ",")
            k == keys[end] || print(io, " ")
        end
        print(io, ")")
    else
        for (i, (k, v)) in enumerate(zip(keys, vals))
            is_last = i == length(keys)
            connector = is_last ? "└─ " : "├─ "
            continuation = is_last ? "   " : "│  "
            if _is_nested_struct(v)
                print(io, "\n$(indent)$(connector)$(k): ")
                nested_io = IOContext(io, :indent => indent * continuation)
                show(nested_io, mime, v)
            elseif v isa NamedTuple && length(v) > 5
                # Large NamedTuples: print on new lines with tree formatting
                print(io, "\n$(indent)$(connector)$(k):")
                nt_keys = Base.keys(v)
                for (j, nk) in enumerate(nt_keys)
                    nt_last = j == length(nt_keys)
                    nt_connector = nt_last ? "└─ " : "├─ "
                    print(io, "\n$(indent)$(continuation)$(nt_connector)$(nk) = ", v[nk])
                end
            else
                print(io, "\n$(indent)$(connector)$(k) = ", v)
                _print_unit(io, units, k)
            end
        end
    end
end

"""
    compact_show_type_and_fields(io::IO, ::MIME"text/plain", x; kw...)

Print a compact one-line representation of `x`: `TypeName(field = value, ...)`.

# Keyword arguments
- `with_units = true`: Append unit annotations from [`field_units`](@ref).
    + When true, display e.g. `1.0 [m/s]` (if [`field_units`](@ref) is defined for the struct)
    + When false, display only the value, e.g. `1.0`
- `with_kwargs = true`: Print fields as keyword arguments (`field = value`).
    + It is recommended that this is `true` only when the struct has
        a `@kwdef` prefix to its definition, which provides a keyword constructor.
- `with_module_prefix = false`: Use module-qualified type name.
    + When true, display e.g. `CloudMicrophysics.Parameters.AirProperties`
    + When false, display e.g. `AirProperties`
- `skip_fields_by_value = ()`: Skip fields whose values are in this tuple.
    + This is useful for skipping default values when displaying the struct
        in the REPL, for example.

# Example output

```julia-repl
julia> import CloudMicrophysics.Parameters as CMP

julia> ap = CMP.AirProperties(Float32)
AirProperties(K_therm = 0.024 [W/m/K], D_vapor = 2.26e-5 [m²/s], ν_air = 1.6e-5 [m²/s])
```
"""
function compact_show_type_and_fields(io::IO, ::MIME"text/plain", x;
    with_units = true, with_kwargs = true,
    with_module_prefix = false, skip_fields_by_value = (),
)
    typename = _type_name(x; with_module_prefix)
    keys = fieldnames(typeof(x))
    vals = [getfield(x, k) for k in keys]
    units = with_units ? field_units(x) : nothing
    print(io, typename, "(")
    printed_first = false
    for (k, v) in zip(keys, vals)
        v ∈ skip_fields_by_value && continue
        printed_first && print(io, ", ")
        with_kwargs && print(io, k, " = ")
        show(io, v)
        _print_unit(io, units, k)
        printed_first = true
    end
    print(io, ")")
end

"""
	parseable_show_with_fields_no_type_header(io::IO, x; kw...)

Print a parseable (copy-pasteable) one-line representation of `x`.
Thin wrapper around [`compact_show_type_and_fields`](@ref) with
`with_units = false` and `with_module_prefix = true` by default.

Note: This assumes that the type has a `@kwdef` prefix to its definition,
which provides a keyword constructor.

# Examples
```julia-repl
julia> show(stdout, CM.Parameters.NumberAdjustmentHorn2012(τ = 3.0f0))
CloudMicrophysics.Parameters.NumberAdjustmentHorn2012(τ = 3.0f0)
```
"""
parseable_show_with_fields_no_type_header(io::IO, x; with_module_prefix = true, kw...) =
    compact_show_type_and_fields(io, MIME("text/plain"), x;
        with_units = false, with_module_prefix, kw...,
    )


"""
    show_type_and_fields(io::IO, mime::MIME"text/plain", x; threshold = 7, kw...)

Auto-select verbose or compact display based on field count.

Types with more than `threshold` fields are displayed with
[`verbose_show_type_and_fields`](@ref) (multi-line tree).
Types with `threshold` or fewer fields use
[`compact_show_type_and_fields`](@ref) (one-liner with units).
"""
show_type_and_fields(io::IO, mime::MIME"text/plain", x; threshold = 7, kw...) =
    if fieldcount(typeof(x)) > threshold
        verbose_show_type_and_fields(io, mime, x; kw...)
    else
        compact_show_type_and_fields(io, mime, x; kw...)
    end


###
### Internal helper methods
###

"""
    field_units(x)

Return a `NamedTuple` mapping field names to unit strings for type `x`,
or `nothing` if no units are defined.

# Example

To include unit annotations in the verbose or compact `show` output for a type,
define a method that associates field names with units, such as

```julia
@kwdef struct MyStruct{FT}
    velocity::FT
    mass::FT
    constant_value::FT
end

Base.show(io::IO, mime::MIME"text/plain", x::MyStruct) = 
    ShowMethods.compact_show_type_and_fields(io, mime, x)

ShowMethods.field_units(x::MyStruct) = (; velocity = "m/s", mass = "kg")
```

If a field is not present in the returned `NamedTuple`, it is assumed to be
dimensionless. In the example above, `constant_value` is dimensionless. When
this struct is displayed in REPL, it will be shown as:

```julia-repl
julia> MyStruct(; velocity = 1.0, mass = 2.0, constant_value = 3.0)
MyStruct(velocity = 1.0 [m/s], mass = 2.0 [kg], constant_value = 3.0 [-])
```

As a reminder, as long as the parameter struct is a subtype of
[`ParametersType`](@ref CloudMicrophysics.Parameters.ParametersType), 
you should not need to define `Base.show` for it.
"""
field_units(_) = nothing

"Print unit annotation for field `k` based on the units spec."
_print_unit(::IO, _, _) = nothing
function _print_unit(io::IO, units::NamedTuple, k::Symbol)
    unit = haskey(units, k) ? units[k] : "-"
    print(io, " [", unit, "]")
end

function _type_name(x; with_module_prefix = false)
    c = Base.text_colors[:cyan]
    n = Base.text_colors[:normal]
    name = if with_module_prefix
        # name with module prefix, e.g. CloudMicrophysics.Parameters.SB2006
        typeof(x).name.wrapper
    else
        # name without module prefix, e.g. SB2006
        typeof(x).name.name
    end
    return string(c, name, n)
end

"""
    _is_nested_struct(v)

Return `true` if `v` is a composite struct that should be recursively
pretty-printed (i.e. a non-trivial struct with fields, but not a
`Number`, `AbstractString`, `Symbol`, or `Nothing`).
"""
_is_nested_struct(v) =
    isstructtype(typeof(v)) &&
    fieldcount(typeof(v)) > 0 &&
    !(v isa Union{Number, AbstractString, Symbol, Nothing})

end  # module ShowMethods
