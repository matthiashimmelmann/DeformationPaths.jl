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

Bar-and-Joint frameworks are given by embedded graphs whose edges are equipped with Euclidean distance constraints. We can create a deformation path using the constructor `DeformationPath`. For example, a deformation path starting in a realization of the complete bipartite graph $K_{2,4}$ (giving rise to a `Framework`) can be computed using the following code. For instance, we can generate the framework on a list of edges and a realization given by a matrix. Here, vertex `i` corresponds to the `i`-th column in the matrix.

```julia
F = Framework([[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6]], Matrix([0 0; 0 1; 1 -1; 1 0; 1 1; 1 2]'))
D = DeformationPath(F, [1], 500; step_size=0.025)
animate(D,F,"completebipartite_motion")
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/completebipartite_motion.gif)

There is only a single infinitesimal motion, which we select via setting `flex_mult=[1]`. This input selects the initial tangent vector $\sum_{i=1}^r c_i \mathbf{q}_i$ for a basis $(\mathbf{q}_1,\dots,\mathbf{q}_r)$ of nontrivial infinitesimal flexes. We compute `500` predictor corrector steps with a step size of `0.025` and save the resulting animation under the file name `"completebipartite_motion.png"` using the method `animate`.

We can also pin vertices in the framework, which are depicted as triangles. This is exemplified through the following animation of the cuspidal double Watt mechanism:

```julia
F = Framework([[1,2],[2,3],[2,4],[3,9],[3,4],[3,5],[4,5],[5,6],[6,7],[7,8],[7,9],[8,9],[8,10],[9,10],[10,11]], Matrix([0 0; 1 0; 2 1; 1 2; 3 2; 4 2; 5 2; 7 2; 6 1; 7 0; 8 0;]'); pinned_vertices=[1,6,11])
D = DeformationPath(F, [-0.5,-0.5], 500; step_size=0.05)
animate(D,F,"double_watt_motion"; padding=0.35, fixed_pair=(1,6), fixed_direction=[4,2])
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/double_watt_motion.gif)

### Sticky Disk Packings

Disk packings are given by a non-overlapping arrangements of disks with fixed radii in $\mathbb{R}^2$. They are called sticky, when existing contacts cannot be broken. In this package, disks that form a contact during the deformation path computation remain in contact. The `DiskPacking` class takes a list of radii and a realization as input. As the optional parameter `pinned_vertices`, we can specify which vertices in the disk packings should be pinned. As an example, we can create an animation using the following code:

```julia
F = DiskPacking([1.,1.,1.,1.], Matrix([0 0; 1.75 -sqrt(2^2-(1.75)^2); 3.5 0; 4.5 sqrt(3)]'); pinned_vertices=[1])
D = DeformationPath(F, [1,1], 250; step_size=0.025)
animate(D,F,"diskpacking_motion")
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/diskpacking_motion.gif)

### Polytopes with Coplanarity Constraints

Polytopes are geometric objects with flat sides. We model them as a collection of vertices, edges and facets. We fix their edge lengths and constrain facets to remain flat by requiring that all edges are normal to the unit outer facet normal. We can compute a deformation of the regular cube using the following code:

```julia
F = Polytope([[1,2,3,4],[5,6,7,8],[1,2,5,6],[2,3,6,7],[3,4,7,8],[1,4,5,8]], Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1]'))
D = DeformationPath(F, [1,1,1], 100; step_size=0.025)
animate(D,F,"cube_motion")
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/cube_motion.gif)

### Volume Hypergraphs

Given a $d+1$-uniform hypergraph $G=(V,E)$ and a realization $p:V\rightarrow \mathbb{R}^d$, we can constrain the volume of its simplicial faces via the determinant of the matrix that can be obtained by adding a row consisting of only ones to the realized vertices belonging to a facet. The `VolumeHypergraph` class allows us to create such structures, for instance a 2-dimensional realization of the octahedral graph with missing faces:

```julia
F = VolumeHypergraph([[1,3,6],[1,2,5],[2,3,4],[1,5,6],[6,4,5]], Matrix([0 0; 3 0; 0 3; 1 1; 1 0.5; 0.5 1]'))
D = DeformationPath(F, [0.333, 1], 350; step_size=0.002)
animate(D, F,"octahedral_decomposition_motion"; fixed_triangle=(6,4,5), skip_stretch=false, target_stretch=0.5, tip_value=0.5)
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/octahedral_decomposition_motion.gif)

In this class, we have additional options. We can fix a triangle using `fixed_triangle=(6,4,5)`. The rescaling can be skipped using `skip_stretch=true`. The first two vertices are fixed to the origin and the $x$-axis, respectively. The third vertex is transformed to the point `(target_stretch, tip_value)`.

### Spherical Disk Packings

Besides planar disk packings, spherical disk packings have also garnered interest in recent times. Instead of the Euclidean distance, the inversive distance is used as a constraint here. The coordinates $(a,b,c)$ of the packing's vertices are represented by affine hyperplanes of the form $ax+by+cy=1$. This choice makes it possible to compute the inversive distance in a particularly elegant way. We can generate a `SphericalDiskPacking` object using its contacts and a realization so that each point has a norm of at least 1. Again, we can pin vertices via the keyword `pinned_vertices`. For example, we can create a deformation of the spherical disk packing corresponding to the flexible Bricard octahedron using the following commands:

```julia
F = SphericalDiskPacking([(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(3,5),(4,5),(2,6),(3,6),(4,6),(5,6)], Matrix([sqrt(2) 0 0; 0 sqrt(2) 0; 0 0 sqrt(2); 0 -sqrt(2) 0; 0 0 -sqrt(2); -sqrt(2) 0 0]'); pinned_vertices=[1])
D = DeformationPath(F, [1], 250; step_size=0.01)
animate(D,F,"sphericaldiskpacking_motion")
```

![](https://github.com/matthiashimmelmann/DeformationPaths.jl/blob/master/data/sphericaldiskpacking_motion.gif)

## Exported functions

- `ConstraintSystem`: Creates a general `ConstraintSystem` object. 
- `Framework`: Creates a bar-and-joint framework object.
- `DiskPacking`: Creates a `DiskPacking` object.
- `Polytope`: Creates a `Polytope` object.
- `VolumeHypergraph`: Creates a volume-constrained hypergraph object.
- `SphericalDiskPacking`: Creates a `SphericalDiskPacking` object. 
- `DeformationPath`: Creates a deformation path corresponding to a given geometric constraint system.
- `animate`: Animates a deformation path. 
- `plot`: Plots the given realization of a geometric constraint system.
- `project_deformation_random`: Visualizes a random projection of the computed deformation path to $\mathbb{R}^2$ or $\mathbb{R}^3$. This can reveal insights into the geometry of the deformation space.
- `to_Matrix`: Transforms an array to a matrix representing the realization. 
- `to_Array`: Transforms a realization to an array in which only non-constant entries appear. 
