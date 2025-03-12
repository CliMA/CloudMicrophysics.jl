
println("""
    climaformat.jl has been discontinued in favor of JuliaFormatter.
    To use JuliaFormatter, add it to your base environment with:

        julia -e 'using Pkg; Pkg.add(\"JuliaFormatter\")'

    If you already have it in your local environment, you can use

        julia -e 'import JuliaFormatter; JuliaFormatter.format(\".\")'

    Then, format your current package with

        import JuliaFormatter; JuliaFormatter.format(\".\")

    in a Julia REPL.
    See documentation to read more about this change.
    This file will be removed in future releases.
    """
)
