using Documenter
using RibosomeAnalysis

makedocs(
    modules = [RibosomeAnalysis],
    authors = "Anthony Egan",
    sitename = "RibosomeAnalysis.jl",
    checkdocs = :none,
    format = Documenter.HTML(
        canonical = "https://eganeganegan.github.io/RibosomeAnalysis.jl/stable/",
        edit_link = "main",
        prettyurls = get(ENV, "CI", "false") == "true",
    ),
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Usage" => "usage.md",
        "Methodology" => "methodology.md",
        "Validation" => "validation.md",
        "Performance" => "performance.md",
        "Architecture" => "architecture.md",
        "API reference" => "reference.md",
        "Releasing" => "releasing.md",
    ],
)

deploydocs(
    repo = "github.com/eganeganegan/RibosomeAnalysis.jl.git",
    devbranch = "main",
    push_preview = false,
)
