using Documenter
using PolyHarmonicTrigTransforms

makedocs(
    modules = [PolyHarmonicTrigTransforms],
    sitename = "PolyHarmonicTrigTransforms.jl",
    pages = [
        "Home" => "index.md",
    ],
    format = Documenter.HTML(),
    clean = true,
)

println("Documentation built: $(joinpath(@__DIR__, "build"))")
