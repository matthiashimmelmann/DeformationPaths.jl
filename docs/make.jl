using Pkg
is_ci = get(ENV, "GITHUB_ACTIONS", "false") == "true"
if is_ci
    Pkg.develop(PackageSpec(path=pwd()));
end
Pkg.instantiate()
using Documenter
include("../src/DeformationPaths.jl")

if is_ci
    @info "Running in GitHub Actions CI"
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
else
    @info "Running locally"
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
        modules = [DeformationPaths.DeformationPaths],
        remotes = nothing
    )
end

deploydocs(
    repo = "github.com/matthiashimmelmann/DeformationPaths.jl",
    push_preview = false,
    devbranch="master",
    versions = nothing
)