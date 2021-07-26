using CloudMicrophysics, Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

pages = Any["Home" => "index.md", "References" => "References.md"]

mathengine = MathJax(Dict(
    :TeX => Dict(
        :equationNumbers => Dict(:autoNumber => "AMS"),
        :Macros => Dict(),
    ),
))

format = Documenter.HTML(
    prettyurls = get(ENV, "CI", nothing) == "true",
    mathengine = mathengine,
    collapselevel = 1,
)

makedocs(
    bib,
    sitename = "CloudMicrophysics.jl",
    strict = false, # TODO: make strict = true
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
