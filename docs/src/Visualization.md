# Visualization Capabilities

```@meta
CurrentModule = DeformationPaths.DeformationPaths
```

There are various visualization options that are shipped with the package. For instance,
we can plot a geometric constraint system by calling

```@docs
plot(F::AllTypes, filename::Union{String,Nothing}=nothing)
plot!(ax::Union{Axis,Axis3}, F::AllTypes)
```

A specific continuous motion of a geometric constraint system can be animated via
```@docs
animate(D::DeformationPath, F::AllTypes, filename::Union{String,Nothing}=nothing)
```

And we can also visualize a random projection to 2D or 3D of deformation paths:

```@docs
project_deformation_random(D::Union{DeformationPath,Vector{DeformationPath}}, F::AllTypes, projected_dimension::Int, filename::Union{String,Nothing}=nothing)
```

## Plotting Capabilities

```@docs
plot_flexes!(ax, F::AllTypes, flex_Real::Int, flex_color, flex_scale, linewidth, arrowsize)

plot_framework(F::Union{Framework,AngularFramework}, filename::Union{String,Nothing})

plot_frameworkonsurface(F::FrameworkOnSurface, filename::Union{String,Nothing})

plot_spherepacking(F::SpherePacking, filename::Union{String,Nothing})

plot_sphericaldiskpacking(F::SphericalDiskPacking, filename::Union{String,Nothing})

plot_hypergraph(F::VolumeHypergraph, filename::Union{String,Nothing})

plot_polytope(F::Union{Polytope,BodyHinge}, filename::Union{String,Nothing})
```

## Animation Methods

```@docs
animate(F::AllTypes, filename::Union{String,Nothing}=nothing)

animate2D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::Union{String,Nothing})

animate3D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::Union{String,Nothing})

animate3D_frameworkonsurface(D::DeformationPath, F::FrameworkOnSurface, filename::Union{String,Nothing})

animate2D_hypergraph(D::DeformationPath, F::VolumeHypergraph, filename::Union{String,Nothing})

animate3D_polytope(D::DeformationPath, F::Union{Polytope,BodyHinge}, filename::Union{String,Nothing})

animate2D_diskpacking(D::DeformationPath, F::SpherePacking, filename::Union{String,Nothing})

animate3D_spherepacking(D::DeformationPath, F::SpherePacking, filename::Union{String,Nothing})

animate3D_sphericaldiskpacking(D::DeformationPath, F::SphericalDiskPacking, filename::Union{String,Nothing})
```