# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

using CloudMicrophysics, Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

pages = Any[
    "Home" => "index.md",
    "0-moment precipitation microphysics" => "Microphysics0M.md",
    "1-moment precipitation microphysics" => "Microphysics1M.md",
    "Non-equilibrium cloud formation" => "MicrophysicsNonEq.md",
    "Aerosol activation" => "AerosolActivation.md",
    "ARG2000 activation example" => "ARG2000_example.md",
    "API" => "API.md",
    "References" => "References.md",
]

mathengine = MathJax(
    Dict(
        :TeX => Dict(
            :equationNumbers => Dict(:autoNumber => "AMS"),
            :Macros => Dict(),
        ),
    ),
)

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(
    bib,
    sitename = "CloudMicrophysics.jl",
    strict = true,
    format = format,
    checkdocs = :exports,
    clean = true,
    doctest = true,
    modules = [CloudMicrophysics],
    pages = pages,
)

deploydocs(
    repo = "github.com/CliMA/CloudMicrophysics.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
