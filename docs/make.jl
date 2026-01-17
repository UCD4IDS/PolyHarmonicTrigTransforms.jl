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
)

println("Documentation built: $(joinpath(@__DIR__, "build"))")
