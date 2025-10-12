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
equations!(G::ConstraintSystem, equations::Vector{Expression})
equations!(F::AllTypes, equations::Vector{Expression})

add_equations!(G::ConstraintSystem, equations::Vector{Expression})
add_equations!(F::AllTypes, equations::Vector{Expression})

realization!(G::ConstraintSystem, realization::Matrix{<:Real})
realization!(F::AllTypes, realization::Matrix{<:Real})
```