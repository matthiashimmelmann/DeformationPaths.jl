# Constraint System Constructors

```@meta
CurrentModule = DeformationPaths.DeformationPaths
```

```@docs
ConstraintSystem

ConstraintSystem(vertices::Vector{Int}, variables::Vector{Variable}, equations::Vector{Expression}, realization::Matrix{<:Real}, xs)
```

## Types of Geometric Constraint Systems

```@docs
SpherePacking
SphericalDiskPacking
Polytope
FrameworkOnSurface
Framework
VolumeHypergraph
BodyHinge
AngularFramework
```

## Transformations

```@docs
equations!

add_equations!

realization!
```