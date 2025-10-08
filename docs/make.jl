using Pkg
Pkg.instantiate()
using Documenter
include("../src/DeformationPaths.jl")

if get(ENV, "CI", "false") == "true"
    @info "Running on CI — skipping plotting setup"
end

makedocs(
    sitename = "DeformationPaths.jl",
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "DeformationPath" => "DeformationPath.md", 
            "ConstraintSystem" => "ConstraintSystems.md", 
            "Visualization" => "Visualization.md",
            "Auxiliary Methods" => "Auxiliary.md",
        ],
        "Usage Guide" => "usage.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", sidebar_sitename = true, assets = ["assets/custom.css"]),
    authors = "Matthias Himmelmann",
    modules = [DeformationPaths.DeformationPaths]
)

deploydocs(
    repo = "github.com/matthiashimmelmann/DeformationPaths.jl",
    push_preview = false,
)