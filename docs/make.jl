# https://github.com/jheinen/GR.jl/issues/278#issuecomment-587090846
ENV["GKSwstype"] = "nul"

using CloudMicrophysics, Documenter
using DocumenterCitations
using Literate

bib = CitationBibliography(joinpath(@__DIR__, "bibliography.bib"))

const GUIDES_DIR = joinpath(@__DIR__, "src/guides/")
const LITERATED_GUIDES_DIR = joinpath(@__DIR__, "src/guides/literated")

if isempty(get(ENV, "CI", ""))
    # only needed when building docs locally; set automatically when built under CI
    # https://fredrikekre.github.io/Literate.jl/v2/outputformats/#Configuration
    extra_literate_config = Dict(
        "repo_root_path" => abspath(joinpath(@__DIR__, "..")),
        "repo_root_url" => "file://" * abspath(joinpath(@__DIR__, "..")),
    )
else
    extra_literate_config = Dict()
end

# To be literated
UnliteratedGuides = ["GettingStarted.jl", "Parameters.jl", "SimpleModel.jl"]
for Guide in UnliteratedGuides
    filepath = joinpath(GUIDES_DIR, Guide)
    Literate.markdown(
        filepath,
        LITERATED_GUIDES_DIR;
        flavor = Literate.DocumenterFlavor(),
        config = extra_literate_config,
    )
end

Parameterizations = [
    "0-moment precipitation microphysics" => "Microphysics0M.md",
    "1-moment precipitation microphysics" => "Microphysics1M.md",
    "2-moment precipitation microphysics" => "Microphysics2M.md",
    "P3 Scheme" => "P3Scheme.md",
    "Terminal Velocity" => "TerminalVelocity.md",
    "Non-equilibrium cloud formation" => "MicrophysicsNonEq.md",
    "Smooth transition at thresholds" => "ThresholdsTransition.md",
    "Aerosol activation" => "AerosolActivation.md",
    "Water Activity" => "WaterActivity.md",
    "Ice Nucleation" => "IceNucleation.md",
    "Aerosol Nucleation" => "AerosolNucleation.md",
    "Cloud diagnostics" => "CloudDiagnostics.md",
    "Precipitation Susceptibility" => "PrecipitationSusceptibility.md",
]

Models = [
    "Adiabatic parcel model" => "IceNucleationParcel0D.md",
    "Box Model" => "IceNucleationBox0D.md",
]

Guides = [
    "Getting started" => "guides/literated/GettingStarted.md",
    "Handling parameters" => "guides/literated/Parameters.md",
    "Building a model" => "guides/literated/SimpleModel.md",
]

pages = Any[
    "Home" => "index.md",
    "Parameterizations" => Parameterizations,
    "How to guides" => Guides,
    "Models" => Models,
    "API" => "API.md",
    "Developer's Guide" => "DevelopersGuide.md",
    "Glossary" => "Glossary.md",
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
    size_threshold_ignore = ["API.md"],
)

makedocs(
    sitename = "CloudMicrophysics.jl",
    format = format,
    checkdocs = :exports,
    modules = [CloudMicrophysics],
    pages = pages,
    plugins = [bib],
)

deploydocs(
    repo = "github.com/CliMA/CloudMicrophysics.jl.git",
    target = "build",
    push_preview = true,
    devbranch = "main",
    forcepush = true,
)
