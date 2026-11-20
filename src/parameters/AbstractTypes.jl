"""
    ParametersType

The top-level super-type for all cloud microphysics parameters
"""
abstract type ParametersType end


# Temporary fallback until we stop checking eltype from the parameters
Base.eltype(p::ParametersType) =
    isempty(propertynames(p)) ? Any : eltype(getproperty(p, first(propertynames(p))))

"""
    @params struct Name{T...} <: Super ... end

Define a parameter struct whose values are stored in its type.

Each instance has size zero, and the compiler treats every parameter value as a constant.

# When to use
- Use `@params` for a parameter struct that holds numbers.
- A struct whose fields are all `@params` structs needs no macro; it has size zero already.

# Example
```julia
@params struct MyParams{FT} <: ParametersType
    "rate [1/s]"
    k::FT
    "exponent [-]"
    n::FT
end

p = MyParams(; k = 1.0, n = 2.0)  # or MyParams(1.0, 2.0)
p.k                               # 1.0
(; k, n) = p                      # destructuring
rate(p::MyParams, x) = p.k * x^p.n
```

# Type parameters
The values replace the type parameters of the definition: `MyParams(1.0, 2.0)` has type
`MyParams{(k = 1.0, n = 2.0)}`, not `MyParams{Float64}`. Do not rely on the type parameters:
- dispatch on `MyParams`, not on `MyParams{FT}`;
- take the float type from a value, e.g. `typeof(p.k)`.
"""
macro params(ex)
    Meta.isexpr(ex, :struct) || error("@params expects a struct definition")
    head = ex.args[2]
    sup = Meta.isexpr(head, :<:) ? head.args[2] : :Any
    namex = Meta.isexpr(head, :<:) ? head.args[1] : head
    name, tparams = Meta.isexpr(namex, :curly) ? (namex.args[1], namex.args[2:end]) : (namex, Any[])
    lines, fields, kwargs = Any[], Any[], Any[]
    for l in ex.args[3].args
        (l isa LineNumberNode || l isa String) && continue  # skip line numbers and docstrings
        default = nothing
        if Meta.isexpr(l, :(=))  # `field::T = default`
            l, default = l.args
        end
        l isa Symbol && (l = :($l::Any))
        Meta.isexpr(l, :(::)) || error("@params: invalid field declaration `$l`")
        f = l.args[1]
        push!(lines, l)
        push!(fields, f)
        push!(kwargs, isnothing(default) ? f : Expr(:kw, f, default))
    end
    P = gensym(:P)
    return esc(quote
        Core.@__doc__ struct $name{$P} <: $sup end
        @inline $name($(lines...)) where {$(tparams...)} = $name{(; $(fields...))}()
        @inline $name(; $(kwargs...)) = $name($(fields...))
        @inline Base.getproperty(::$name{$P}, s::Symbol) where {$P} = getfield($P, s)
        Base.propertynames(::$name) = $(Tuple(fields))
    end)
end
Base.broadcastable(x::ParametersType) = tuple(x)
Base.show(io::IO, x::ParametersType) =
    ShowMethods.parseable_show_with_fields_no_type_header(io, x; with_module_prefix = false)
Base.show(io::IO, mime::MIME"text/plain", x::ParametersType) =
    ShowMethods.show_type_and_fields(io, mime, x)

"""
    AerosolType

The top-level super-type for all aerosol properties
"""
abstract type AerosolType <: ParametersType end

"""
    AerosolDistributionType

The top-level super-type for all aerosol distribution types
"""
abstract type AerosolDistributionType end

"""
    CloudCondensateType

The top-level super-type for cloud condensate types (liquid and ice)
"""
abstract type CloudCondensateType <: ParametersType end

"""
    PrecipitationType

The top-level super-type for precipitation types (rain and snow)
"""
abstract type PrecipitationType <: ParametersType end

"""
    TerminalVelocityType

The top-level super-type for terminal velocity parameterizations
"""
abstract type TerminalVelocityType <: ParametersType end

"""
    Precipitation2MType

The top-level super-type for 2-moment precipitation parameterizations
"""
abstract type Precipitation2MType <: ParametersType end
