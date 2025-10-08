# Visualization Capabilities

```@meta
CurrentModule = DeformationPaths.DeformationPaths
```

There are various visualization options that are shipped with the package. For instance,
we can plot a geometric constraint system by calling

```@docs
plot(F::AllTypes, filename::String)
```

A specific continuous motion of a geometric constraint system can be animated via
```@docs
animate(D::DeformationPath, F, filename::String)
```

And we can also visualize a random projection to 2D or 3D of deformation paths:

```@docs
project_deformation_random(D::Union{DeformationPath,Vector{DeformationPath}}, projected_dimension::Int)
```

## Plotting Capabilities

```@docs
plot_flexes!(ax, F::AllTypes, flex_Real::Int, flex_color, flex_scale, linewidth, arrowsize)

plot_framework(F::Union{Framework,AngularFramework}, filename::String)

plot_frameworkonsurface(F::FrameworkOnSurface, filename::String)

plot_spherepacking(F::SpherePacking, filename::String)

plot_sphericaldiskpacking(F::SphericalDiskPacking, filename::String)

plot_hypergraph(F::VolumeHypergraph, filename::String)

plot_polytope(F::Union{Polytope,BodyHinge}, filename::String)
```

## Animation Methods

```@docs
animate(F::AllTypes, filename::String)

animate2D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::String)

animate3D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::String)

animate3D_frameworkonsurface(D::DeformationPath, F::FrameworkOnSurface, filename::String)

animate2D_hypergraph(D::DeformationPath, F::VolumeHypergraph, filename::String)

animate3D_polytope(D::DeformationPath, F::Union{Polytope,BodyHinge}, filename::String)

animate2D_diskpacking(D::DeformationPath, F::SpherePacking, filename::String)

animate3D_spherepacking(D::DeformationPath, F::SpherePacking, filename::String)

animate3D_sphericaldiskpacking(D::DeformationPath, F::SphericalDiskPacking, filename::String)
```