using Documenter
using PolyHarmonicTrigTransforms

makedocs(
    modules = [PolyHarmonicTrigTransforms],
    repo=Documenter.Remotes.GitHub("UCD4IDS", "PolyHarmonicTrigTransforms.jl"),
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
    repo = "github.com/UCD4IDS/PolyHarmonicTrigTransforms.jl.git",
           devbranch="main",
           push_preview=true
)