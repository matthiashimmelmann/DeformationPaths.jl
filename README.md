# DeformationPaths.jl

This package allows the approximation of deformation paths corresponding to a variety of geometric constraint systems. To accomplish this, we iteratively apply a combination of Euler predictor steps and Gauss-Newton corrector steps to take a local stpe on the constraint set. This approach is known as homotopy continuation. 

## Installation

```
julia> ]
(@v1.10) pkg> add DeformationPaths
```

## Usage

```julia
using DeformationPaths
```

### Bar-and-Joint Frameworks

Bar-and-Joint frameworks are given by embedded graphs whose edges are equipped with Euclidean distance constraints. We can create a deformation path using the constructor `DeformationPath`. For example, a deformation path starting in a realization of the complete bipartite graph $K_{2,4]$ (giving rise to a `Framework`) can be computed using the following code.

```julia
F = Framework([[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6]], Matrix([0 0; 0 1; 1 -1; 1 0; 1 1; 1 2]'))
D = DeformationPath(F, [1], 500; step_size=0.025)
animate(D,F,"completebipartite_motion")
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/completebipartite_motion.gif)

There is only a single infinitesimal motion, which we select via setting `flex_mult=[1]`. This input selects the initial tangent vector $\sum_{i=1}^r c_i \mathbf{q}_i$ for a basis $(\mathbf{q}_1,\dots,\mathbf{q}_r)$ of nontrivial infinitesimal flexes. We compute `500` predictor corrector steps with a step size of `0.025` and save the resulting animation under the file name `"completebipartite_motion.png"` using the method `animate`.

### Sticky Disk Packings

![](https://github.com/matthiashimmelmann/HomotopyOpt.jl/blob/firstbranch/test/Images/watch1.661497754964e9.gif)

### Polytopes with Coplanarity Constraints

### Volume Hypergraphs

### Spherical Disk Packings
