# Auxiliary Functions

## Rigidity Checks

```@meta
CurrentModule = DeformationPaths.DeformationPaths
```

```@docs
is_rigid(F::AllTypes)

is_inf_rigid(F::AllTypes)

is_prestress_stable(F::AllTypes)

is_second_order_rigid(F::AllTypes)

coned_rigidity_phase_space
```

## Predictor-Corrector Methods

```@docs
euler_step

newton_correct

symmetric_newton_correct(G::ConstraintSystem, point::Vector{<:Real})

symmetric_newton_correct(equations::Vector{Expression}, variables::Vector{Variable}, jac::Matrix{Expression}, point::Vector{<:Real})
```

## Infinitessimal Flexes

```@docs
compute_inf_flexes

compute_equilibrium_stresses

compute_trivial_inf_flexes

compute_nontrivial_inf_flexes

compute_nonblocked_flex(F::AllTypes)
```

## Transformation Methods

```@docs
to_Array(G::ConstraintSystem, p::Matrix{<:Real})

to_Array(F::Union{SpherePacking,Framework,AngularFramework,FrameworkOnSurface,SphericalDiskPacking,VolumeHypergraph,BodyHinge,BodyBar}, p::Matrix{<:Real})

to_Array(F::Polytope, p::Matrix{<:Real})

to_Matrix(G::ConstraintSystem, q::Vector{<:Real})

to_Matrix(F::AllTypes, q::Vector{<:Real})

minors
```

## Special Polytope Methods

```@docs
is_in_interior

fix_antipodals!

tetrahedral_symmetry!

triangle_shrinking
```

## I/O methods

```@docs
save_realizations

read_realizations

save_to_Houdini
```