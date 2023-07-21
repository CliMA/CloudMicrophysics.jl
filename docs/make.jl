# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

using CloudMicrophysics, Documenter
using DocumenterCitations

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

pages = Any[
    "Home" => "index.md",
    "0-moment precipitation microphysics" => "Microphysics0M.md",
    "1-moment precipitation microphysics" => "Microphysics1M.md",
    "2-moment precipitation microphysics" => "Microphysics2M.md",
    "Non-equilibrium cloud formation" => "MicrophysicsNonEq.md",
    "Smooth transition at thresholds" => "ThresholdsTransition.md",
    "Aerosol activation" => "AerosolActivation.md",
    "ARG2000 activation example" => "ARG2000_example.md",
    "Ice Nucleation" => "IceNucleation.md",
    "Box model ice nucleation example" => "IceNucleationBox.md",
    "Aerosol Nucleation" => "AerosolNucleation.md",
    "Precipitation Susceptibility" => "PrecipitationSusceptibility.md",
    "API" => "API.md",
    "References" => "References.md",
    "Developer's Guide" => "DevelopersGuide.md",
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
