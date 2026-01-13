using Documenter, PolyharmonicTrigTransforms

makedocs(
 sitename="PolyharmonicTrigTransforms.jl",
    format = Documenter.HTML(collapselevel = 1),
    authors = "",
    clean = true,
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "LLST" => "examples/llst.md",
        ],
        "Functions" => [
            "DST" => "docs/dst.md",
            "Utils" => "functions/Utils.md",
        ],
    ]
)


# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/UCD4IDS/PolyharmonicTrigTransforms.jl.git"
)