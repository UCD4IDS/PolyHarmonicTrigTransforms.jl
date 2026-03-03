using Pkg
Pkg.activate(@__DIR__)  # activate docs environment
Pkg.instantiate()       # install Documenter and dependencies
Pkg.develop(path="../src")  # ensure the package is available for documentation

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

deploydocs(
    repo = "github.com/UCD4IDS/PolyHarmonicTrigTransforms.git"
)