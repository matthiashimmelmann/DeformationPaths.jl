using Pkg
Pkg.instantiate()
using Documenter
include("../src/DeformationPaths.jl")

makedocs(
    sitename = "DeformationPaths.jl",
    pages = [
        "Home" => "index.md",
        "API Reference" => [
            "DeformationPath" => "API/DeformationPath.md", 
            "ConstraintSystem" => "API/ConstraintSystems.md", 
            "Visualization" => "API/Visualization.md",
            "Auxiliary Methods" => "API/Auxiliary.md",
        ],
        "Usage Guide" => "usage.md"
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", sidebar_sitename = false, assets = ["assets/custom.css"]),
    authors = "Matthias Himmelmann",
    modules = [DeformationPaths.DeformationPaths],
    remotes = nothing
)