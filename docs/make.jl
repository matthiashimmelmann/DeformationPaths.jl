using Pkg
cd(@__DIR__)
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()
using Documenter
include("../src/DeformationPaths.jl")

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
    modules = [DeformationPaths.DeformationPaths],
    remotes=nothing
)

deploydocs(
    repo = "github.com/matthiashimmelmann/DeformationPaths.jl",
    push_preview = false,
    devbranch="master"
)