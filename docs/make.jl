using Pkg
Pkg.develop(PackageSpec(path=pwd()));
Pkg.instantiate()
using Documenter
include("../src/DeformationPaths.jl")

makedocs(
    sitename = "DeformationPaths.jl",
    pages = [
        "Home" => "index.md",
        "Usage Guide" => "usage.md",
        "API Reference" => [
            "DeformationPath" => "DeformationPath.md", 
            "ConstraintSystem" => "ConstraintSystems.md", 
            "Visualization" => "Visualization.md",
            "Auxiliary Methods" => "Auxiliary.md",
        ]    
    ],
    format = Documenter.HTML(prettyurls = get(ENV, "CI", nothing) == "true", sidebar_sitename = true, assets = ["assets/custom.css"]),
    authors = "Matthias Himmelmann",
    modules = [DeformationPaths.DeformationPaths]
)

deploydocs(
    repo = "github.com/matthiashimmelmann/DeformationPaths.jl",
    push_preview = false,
    devbranch="master",
    versions = nothing
)