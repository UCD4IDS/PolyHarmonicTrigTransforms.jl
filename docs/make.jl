using Pkg
Pkg.activate(@__DIR__)  # activate docs environment
Pkg.instantiate()       # install Documenter and dependencies
Pkg.develop(path="..")

using Documenter
using PolyHarmonicTrigTransforms

makedocs(
    modules = [PolyHarmonicTrigTransforms],
    sitename = "PolyHarmonicTrigTransforms.jl",
    pages = [
        "Home" => "index.md",
        "API"  => "api.md",
    ],
    format = Documenter.HTML(),
    clean = true,
    debug = true,
    checkdocs = :none, # disable docstring checks to avoid failing on undocumented symbols
)

