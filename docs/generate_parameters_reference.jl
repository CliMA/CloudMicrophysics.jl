#=
Generate `docs/src/generated/ParametersReference.md`: the mapping between
ClimaParams parameter names and the fields of the CloudMicrophysics parameter
structs.

The page is built from the source code so that it never goes stale:

  1. Every constructor in `src/parameters/*.jl` that takes a `CP.ParamDict` is
     found by parsing the source with `Meta.parseall`. The `:clima_name => :field`
     pairs in its body (the `name_map`s) give the static ClimaParams → field
     mapping, and calls to other constructors give the composition
     (e.g. `Rain` is built from `ParticleMass(Rain)` and `ParticleArea(Rain)`).
  2. Every constructor is then actually called on a fresh `ParamDict`. ClimaParams
     tags each parameter it hands out with a `"used_in"` entry, so the set of
     parameters the constructor *really* reads is known exactly.
  3. The two views are reconciled. A parameter that the static parse attributes
     to a constructor but that is never read is an error (the docs would lie).
     Parameters that are read but not visible statically (dispatch through
     process options, nested constructors) are attributed to the constructors
     that do read them.

Run from `docs/make.jl` before `makedocs`; can also be run standalone with the
docs (or test) environment active.
=#

import ClimaParams as CP
import CloudMicrophysics
import CloudMicrophysics.Parameters as CMP
import InteractiveUtils

const PARAM_DIR = joinpath(pkgdir(CloudMicrophysics), "src", "parameters")
const OUT_DIR = joinpath(@__DIR__, "src", "generated")
const OUT_FILE = joinpath(OUT_DIR, "ParametersReference.md")
const FT = Float64

# Section titles and the source files that belong to them. Files not listed
# here are still documented, under "Other parameter structs".
const GROUPS = [
    "Air and water properties" => ["AirProperties.jl", "WaterProperties.jl"],
    "0-moment microphysics" => ["Microphysics0M.jl", "Microphysics0MParams.jl"],
    "1-moment microphysics" => ["Microphysics1M.jl", "Microphysics1MParams.jl"],
    "1-moment process options" => ["Microphysics1MOptions.jl"],
    "Non-equilibrium cloud formation and 2-moment microphysics" =>
        ["Microphysics2M.jl", "Microphysics2MParams.jl"],
    "P3 scheme" => ["MicrophysicsP3.jl"],
    "Terminal velocity" => ["TerminalVelocity.jl"],
    "Aerosol activation" => ["AerosolActivation.jl"],
    "Ice nucleation" => ["IceNucleation.jl"],
    "Aerosol nucleation" => ["AerosolModalNucleation.jl", "Aerosol_H2SO4_Solution.jl"],
    "Aerosol species" => [
        "AerosolATD.jl", "AerosolAsianDust.jl", "AerosolDesertDust.jl", "AerosolDust.jl",
        "AerosolFeldspar.jl", "AerosolFerrihydrite.jl", "AerosolIllite.jl",
        "AerosolKaolinite.jl", "AerosolMiddleEasternDust.jl", "AerosolSaharanDust.jl",
        "AerosolSeasalt.jl", "AerosolSulfate.jl",
    ],
]

# Process options whose `process_params` slot cannot be inferred from the
# `Microphysics1MOptions` defaults or from a shared abstract supertype.
const OPTION_SLOT_OVERRIDES = Dict(
    :Homogeneous => :cloud_liquid_freezing,
    :Heterogeneous => :cloud_liquid_freezing,
)

# ───────────────────────────── static parse ─────────────────────────────

struct Ctor
    file::String
    name::Symbol
    dispatch::Union{Nothing, Expr, Symbol}  # `X` from a leading `::X` / `::Type{X}` argument
    pairs::Vector{Pair{Symbol, Symbol}}     # ClimaParams name => struct field
    body::Any
end

# A label unique per constructor: `Rain`, `ParticleMass(Rain)`, `Kessler1M`
function label(c::Ctor)
    c.dispatch === nothing && return string(c.name)
    c.name == :process_params_for && return string(c.dispatch)
    d = c.dispatch
    d isa Expr && d.head == :curly && d.args[1] == :Type && (d = d.args[2])
    return string(c.name, "(", d, ")")
end
# The binding the docstring link should point at
docbinding(c::Ctor) = c.name == :process_params_for ? Symbol(string(c.dispatch)) : c.name

isparamdict(ann) = ann == :(CP.ParamDict) || ann == :ParamDict
function argtype(a)
    a isa Expr && a.head == :(::) && return a.args[end]
    a isa Expr && a.head == :kw && return argtype(a.args[1])
    return nothing
end
unwrap_where(sig) = sig isa Expr && sig.head == :where ? unwrap_where(sig.args[1]) : sig

# `:a`, or `Symbol("a-b")` for names that are not valid identifiers
function symof(x)
    x isa QuoteNode && x.value isa Symbol && return x.value
    x isa Expr && x.head == :call && x.args[1] == :Symbol && length(x.args) == 2 &&
        x.args[2] isa String && return Symbol(x.args[2])
    return nothing
end

function collect_pairs!(out, ex)
    ex isa Expr || return
    if ex.head == :call && length(ex.args) == 3 && ex.args[1] == :(=>)
        k, v = symof(ex.args[2]), symof(ex.args[3])
        if k !== nothing && v !== nothing
            push!(out, k => v)
            return
        end
    end
    # CP.get_parameter_values(td, "name" | ["a", "b"], ...) keeps the ClimaParams name
    if ex.head == :call && ex.args[1] == :(CP.get_parameter_values) && length(ex.args) >= 3
        n = ex.args[3]
        n isa String && push!(out, Symbol(n) => Symbol(n))
        n isa Expr && n.head == :vect && all(x -> x isa String, n.args) &&
            foreach(s -> push!(out, Symbol(s) => Symbol(s)), n.args)
    end
    foreach(a -> collect_pairs!(out, a), ex.args)
end

function collect_ctors!(out, ex, file)
    ex isa Expr || return
    if ex.head in (:macrocall, :block, :toplevel)   # docstrings wrap definitions in @doc
        foreach(a -> collect_ctors!(out, a, file), ex.args)
        return
    end
    (ex.head == :function || ex.head == :(=)) && ex.args[1] isa Expr || return
    sig = unwrap_where(ex.args[1])
    sig isa Expr && sig.head == :call && sig.args[1] isa Symbol || return
    args = filter(a -> !(a isa Expr && a.head == :parameters), sig.args[2:end])
    isempty(args) && return
    # the ParamDict must be the last positional argument ...
    isparamdict(argtype(args[end])) || return
    # ... and any other positional argument must be dispatch-only (`::X`)
    lead = args[1:(end - 1)]
    length(lead) <= 1 || return
    dispatch = nothing
    if length(lead) == 1
        a = lead[1]
        a isa Expr && a.head == :(::) && length(a.args) == 1 || return
        dispatch = a.args[1]
    end
    pairs = Pair{Symbol, Symbol}[]
    collect_pairs!(pairs, ex.args[2])
    push!(out, Ctor(file, sig.args[1], dispatch, pairs, ex.args[2]))
end

ctors = Ctor[]
for f in sort(readdir(PARAM_DIR))
    endswith(f, ".jl") || continue
    ex = Meta.parseall(read(joinpath(PARAM_DIR, f), String); filename = f)
    collect_ctors!(ctors, ex, f)
end
const CTOR_NAMES = Set(c.name for c in ctors)

# Does the call `name(arg, td)` select the constructor `c`? `ParticleMass(Rain, td)`
# selects the `::Type{Rain}` method, `process_params_for(Kessler1M(), td)` the
# `::Kessler1M` one. Calls whose first argument is a variable cannot be resolved.
function selects(c::Ctor, call::Expr)
    c.dispatch === nothing && return true
    length(call.args) >= 2 || return false
    arg = call.args[2]
    d = c.dispatch
    d isa Expr && d.head == :curly && d.args[1] == :Type && return arg == d.args[2]
    return arg isa Expr && arg.head == :call && length(arg.args) == 1 && arg.args[1] == d
end

# Labels of the other constructors called in the body (composition)
function called_ctors(c::Ctor)
    out = Set{String}()
    walk(ex) = begin
        ex isa Expr || return
        if ex.head == :call && ex.args[1] isa Symbol && ex.args[1] in CTOR_NAMES
            for other in ctors
                other.name == ex.args[1] && label(other) != label(c) && selects(other, ex) &&
                    push!(out, label(other))
            end
        end
        foreach(walk, ex.args)
    end
    walk(c.body)
    return out
end

# ──────────────────────────── dynamic check ─────────────────────────────

function construct(c::Ctor, td)
    f = getfield(CMP, c.name)
    c.dispatch === nothing && return f(td)
    d = c.dispatch
    d isa Expr && d.head == :curly && d.args[1] == :Type && return f(getfield(CMP, d.args[2]), td)
    d == :Nothing && return f(nothing, td)
    return f(getfield(CMP, d)(), td)
end
used_names(td) = Set(Symbol(k) for (k, v) in td.data if haskey(v, "used_in"))

struct Info
    ctor::Ctor
    used::Set{Symbol}
    calls::Set{String}   # labels of the constructors called in the body
end

infos = Info[]
for c in ctors
    td = CP.create_toml_dict(FT)
    obj = try
        construct(c, td)
    catch e
        # abstract fallbacks such as `process_params_for(::MicrophysicsOption, td)`
        isempty(c.pairs) && continue
        error("Could not construct $(label(c)) ($(c.file)) to verify its parameters: $e")
    end
    used = used_names(td)
    calls = called_ctors(c)
    isempty(c.pairs) && isempty(calls) && continue
    parsed = Set(first.(c.pairs))
    stale = setdiff(parsed, used)
    isempty(stale) || error(
        "$(label(c)) ($(c.file)) maps $(collect(stale)) in its name_map, but the " *
        "constructor never reads them from ClimaParams. Fix the constructor or " *
        "the generator's parser.",
    )
    push!(infos, Info(c, used, calls))
end

const BY_LABEL = Dict(label(i.ctor) => i for i in infos)
used_by_label(l) = haskey(BY_LABEL, l) ? BY_LABEL[l].used : Set{Symbol}()

group_of(file) = (g = findfirst(((_, files),) -> file in files, GROUPS); g === nothing ? "" : GROUPS[g].first)

# Parameters read by `i` that neither its own name_map nor the constructors it
# calls explain (e.g. read through the process options selected by default):
# attribute them to the constructors whose parameters are all among them. When
# one candidate's parameters are a strict subset of another's (`Homogeneous`
# vs `HomogeneousAndHeterogeneous`), keep only the larger one; when two
# candidates read the same parameters (`ConstantTimescale` vs `SubDep2M`),
# prefer the one documented in the same section as `i`.
function unattributed(i::Info)
    rest = setdiff(i.used, Set(first.(i.ctor.pairs)), (used_by_label(l) for l in i.calls)...)
    cands = [j for j in infos if j !== i && !isempty(j.used) && j.used ⊆ rest]
    filter!(j -> !any(k -> k !== j && j.used ⊊ k.used, cands), cands)
    same_group(j) = group_of(j.ctor.file) == group_of(i.ctor.file)
    filter!(j -> same_group(j) || !any(k -> k !== j && k.used == j.used && same_group(k), cands), cands)
    via = sort!([label(j.ctor) for j in cands])
    for j in cands
        setdiff!(rest, j.used)
    end
    return via, rest
end

# ─────────────────────── 1M process option slots ────────────────────────

const OPTION_SLOT = Dict{Symbol, Symbol}()
let opts = CMP.Microphysics1MOptions()
    for slot in fieldnames(typeof(opts))
        T = typeof(getfield(opts, slot))
        OPTION_SLOT[nameof(T)] = slot
        # options sharing a (non-root) abstract supertype share the slot
        S = supertype(T)
        S === CMP.MicrophysicsOption && continue
        for U in InteractiveUtils.subtypes(S)
            OPTION_SLOT[nameof(U)] = slot
        end
    end
    merge!(OPTION_SLOT, OPTION_SLOT_OVERRIDES)
end

function option_slot(c::Ctor)
    opt = Symbol(string(c.dispatch))
    haskey(OPTION_SLOT, opt) && return OPTION_SLOT[opt]
    error(
        "Cannot infer the `process_params` slot of the 1-moment process option `$opt`. " *
        "Add it to `OPTION_SLOT_OVERRIDES` in docs/generate_parameters_reference.jl.",
    )
end

# ───────────────────────────── rendering ─────────────────────────────────

const DEFAULTS = CP.create_toml_dict(FT)
hasdoc(sym) = haskey(Base.Docs.meta(CMP), Base.Docs.Binding(CMP, sym))
reflink(sym) = hasdoc(sym) ? "[`$sym`](@ref CloudMicrophysics.Parameters.$sym)" : "`$sym`"
# `Rain` links to its docstring; `ParticleMass(Rain)` links `ParticleMass` and keeps `(Rain)` as text
function ctor_link(c::Ctor)
    b = docbinding(c)
    (c.dispatch === nothing || c.name == :process_params_for) && return reflink(b)
    return reflink(b) * "`" * label(c)[(length(string(c.name)) + 1):end] * "`"
end
label_link(l::String) = ctor_link(infos[findfirst(i -> label(i.ctor) == l, infos)].ctor)

esc_cell(s) = replace(string(s), "|" => "\\|", "\n" => " ")
function fmt_value(v)
    v isa AbstractVector && length(v) > 4 && return "[" * join(v[1:3], ", ") * ", …] (" * string(length(v)) * " values)"
    v isa AbstractVector && return "[" * join(v, ", ") * "]"
    v isa AbstractString && return "\"" * v * "\""
    return string(v)
end
param_default(n) = fmt_value(DEFAULTS.data[String(n)]["value"])
# ClimaParams descriptions write math as `$...$`; Documenter wants ``...``
param_description(n) = replace(get(DEFAULTS.data[String(n)], "description", ""), r"\$([^\$]+)\$" => s"``\1``")

function render_ctor(io, i::Info)
    c = i.ctor
    println(io, "### ", label(c), "\n")
    if c.name == :process_params_for
        println(io, reflink(docbinding(c)), " option; its parameters are stored in ",
            "`process_params.", option_slot(c), "` of ", reflink(:Microphysics1MParams), ".\n")
    else
        println(io, "Constructor: ", ctor_link(c), "\n")
    end
    if !isempty(i.calls)
        println(io, "Built from: ", join(label_link.(sort(collect(i.calls))), ", "), "\n")
    end
    via, rest = unattributed(i)
    if !isempty(via)
        println(io, "Also reads the parameters of (selected through the default options): ",
            join(label_link.(via), ", "), "\n")
    end
    if !isempty(rest)
        println(io, "Also reads: ", join(("`$n`" for n in sort(collect(rest))), ", "), "\n")
    end
    if !isempty(c.pairs)
        println(io, "| ClimaParams name | Field | Default | Description |")
        println(io, "|:-----------------|:------|:--------|:------------|")
        for (n, f) in c.pairs
            println(io, "| `", n, "` | `", f, "` | ", esc_cell(param_default(n)), " | ",
                esc_cell(param_description(n)), " |")
        end
        println(io)
    end
end

mkpath(OUT_DIR)
open(OUT_FILE, "w") do io
    println(
        io,
        """
```@meta
CurrentModule = CloudMicrophysics
```

# Parameter reference

!!! note
    This page is generated by `docs/generate_parameters_reference.jl` from the
    constructors in `src/parameters/` and the default values in
    ClimaParams v$(pkgversion(CP)). Do not edit it by hand.

The first part lists, for every parameter struct, which ClimaParams parameters
it reads and into which fields they go. The second part is the reverse index:
every ClimaParams parameter used by CloudMicrophysics and the struct fields it
feeds. See the [Parameters interface](../Parameters.md) page for how these
structs are built and overridden.
""",
    )

    grouped = Set{String}()
    for (title, files) in GROUPS
        sel = [i for i in infos if i.ctor.file in files]
        isempty(sel) && continue
        union!(grouped, files)
        println(io, "## ", title, "\n")
        foreach(i -> render_ctor(io, i), sel)
    end
    other = [i for i in infos if !(i.ctor.file in grouped)]
    if !isempty(other)
        @warn "Parameter files not assigned to a section in GROUPS" unique(i.ctor.file for i in other)
        println(io, "## Other parameter structs\n")
        foreach(i -> render_ctor(io, i), other)
    end

    # reverse index
    index = Dict{Symbol, Vector{String}}()
    for i in infos, (n, f) in i.ctor.pairs
        push!(get!(index, n, String[]), "`" * label(i.ctor) * "." * string(f) * "`")
    end
    println(io, "## All ClimaParams parameters used by CloudMicrophysics\n")
    println(io, length(index), " parameters. `Used by` lists `Struct.field` for every ",
        "struct that reads the parameter.\n")
    println(io, "| ClimaParams name | Default | Used by | Description |")
    println(io, "|:-----------------|:--------|:--------|:------------|")
    for n in sort(collect(keys(index)); by = s -> lowercase(string(s)))
        println(io, "| `", n, "` | ", esc_cell(param_default(n)), " | ",
            join(unique(index[n]), ", "), " | ", esc_cell(param_description(n)), " |")
    end
end

@info "Generated $OUT_FILE" constructors = length(infos) parameters = length(
    union(Set{Symbol}(), (Set(first.(i.ctor.pairs)) for i in infos)...),
)
