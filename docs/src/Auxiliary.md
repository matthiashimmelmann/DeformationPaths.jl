# Auxiliary Functions

## Rigidity Checks

```@meta
CurrentModule = DeformationPaths.DeformationPaths
```

```@docs
is_rigid(F::AllTypes)

is_inf_rigid(F::AllTypes)

is_second_order_rigid(F::AllTypes)
```

## Predictor-Corrector Methods

```@docs
euler_step(G::DeformationPaths.ConstraintSystem, step_size::Real, prev_flex::Vector{<:Real}, point::Vector{<:Real}, K_n::DeformationPaths.ConstraintSystem)

newton_correct(G::DeformationPaths.ConstraintSystem, point::Vector{<:Real})

newton_correct(equations::Vector{Expression}, variables::Vector{Variable}, jac::Matrix{Expression}, point::Vector{<:Real})

symmetric_newton_correct(G::ConstraintSystem, point::Vector{<:Real})

symmetric_newton_correct(equations::Vector{Expression}, variables::Vector{Variable}, jac::Matrix{Expression}, point::Vector{<:Real})
```

## Infinitessimal Flexes

```@docs
compute_nontrivial_inf_flexes(G::ConstraintSystem, point::Vector{<:Real}, K_n::ConstraintSystem)

compute_nonblocked_flex(F::AllTypes)
```

## Transformation Methods

```@docs
to_Array(G::ConstraintSystem, p::Matrix{<:Real})

to_Array(F::Union{SpherePacking,Framework,AngularFramework,FrameworkOnSurface,SphericalDiskPacking,VolumeHypergraph,BodyHinge}, p::Matrix{<:Real})

to_Array(F::Polytope, p::Matrix{<:Real})

to_Matrix(G::ConstraintSystem, q::Vector{<:Real})

to_Matrix(F::AllTypes, q::Vector{<:Real})
```

## Special Polytope Methods

```@docs
fix_antipodals!(F::Polytope)

tetrahedral_symmetry!(F::Polytope)

triangle_shrinking(F::Polytope)
```
