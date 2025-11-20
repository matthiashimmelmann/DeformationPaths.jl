module DeformationPaths

import HomotopyContinuation: solve, evaluate, differentiate, newton, Expression, Variable, @var, real_solutions, System, solve, variables
import LinearAlgebra: norm, pinv, nullspace, rank, qr, zeros, inv, cross, det, svd, I, zeros
import GLMakie: NoShading, MultiLightShading, FastShading, GeometryBasics, Vec3f, Vec2f, meshscatter!, surface!, Sphere, mesh!, @lift, poly!, text!, Figure, record, hidespines!, hidedecorations!, lines!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Observable, Point3f, Point2f, connect, faces, Mesh, mesh, save, arrows!
import Combinatorics: powerset
import Colors: distinguishable_colors, red, green, blue, colormap, RGB
import MarchingCubes: MC, march, makemesh
import Polyhedra
import CDDLib
import Base: show
import ProgressMeter: @showprogress

include("GeometricConstraintSystem.jl")
include("InfinitessimalFlexes.jl")
include("PredictorCorrector.jl")
include("Rigidity.jl")
include("GCS_database.jl")

export  ConstraintSystem, 
        Framework,
        AngularFramework,
        DeformationPath,
        VolumeHypergraph,
        animate,
        plot,
        plot!,
        project_deformation_random,
        Polytope,
        to_Matrix,
        to_Array,
        SpherePacking,
        SphericalDiskPacking,
        equations!,
        add_equations!,
        realization!,
        newton_correct,
        FrameworkOnSurface,
        is_rigid,
        is_inf_rigid,
        is_second_order_rigid,
        BodyHinge,
        BodyBar,
        compute_nontrivial_inf_flexes,
        fix_antipodals!,
        tetrahedral_symmetry!,
        compute_nonblocked_flex,
        stich_deformation_paths,
        Dodecahedron,
        add_shadow!

"""
    DeformationPath(G, motion_samples[; tol])
    DeformationPath(G, flex_mult, num_steps, type[; step_size, tol, random_flex, symmetric_newton, start_point])
    DeformationPath(F, flex_mult, num_steps[; random_flex, kwargs...])

Class for constructing approximate deformation paths.

# Attributes
- `G::ConstraintSystem`: The geometric constraint system for which the deformation is computes.
- `step_size::Real`: The step size of the deformation path. 
- `motion_samples::Vector{<:Vector{<:Real}}`: Points in the configuration space representing the approximate motion.
- `motion_matrices::Vector{<:Matrix{<:Real}}`: The points in `motion_samples` as distributed into realizations given by `dxn` matrices for the dimension `d` and the number of vertices `n`.
- `flex_mult::Vector{<:Real}`: The initial infinitesimal flex as a linear combination of the nontrivial infinitesimal flexes at the realization provided by the underlying geometric constraint system.
- `_contacts::Vector`: In some cases (such as sphere packings), the contacts can change during the deformation. This is reflected by this attribute.

# Arguments
- `G::ConstraintSystem`: The underlying geometric constraint system.
- `motion_samples::Vector{<:Vector{<:Real}}`: List of previously computed realizations in array format.
- `flex_mult::Vector`: The initial infinitesimal flex as a linear combination of the nontrivial infinitesimal flexes at the realization provided by the underlying geometric constraint system.
- `num_steps::Int`: Number of steps the algorithm is supposed to take. 
- `type::DataType`: Type of the geometric constraint system for computing the trivial infinitessimal flexes. Possible values: `"framework"`, `"angularframework"`, `"frameworkonsurface"`, `"hypergraph"`, `"polytope"`, `"diskpacking"`, `"bodyhinge"`, `"bodybar"` and `"sphericaldiskpacking"`.
- `step_size::Real`: Step size of the deformation path. 
- `tol::Real` (optional): Numerical tolerance for the approximation that is used for asserting the correctness of the approximation. Default value: `1e-8`.
- `random_flex::Bool` (optional): If `flex_mult` is not provided (e.g. as `[]`), we can instead provide a random linear combination (`true`) or it is uniformly chosen as `[1,...,1]` (`false`). Default value: `false`.
- `symmetric_newton::Bool` (optional): It is possible to use the slightly more efficient symmetric Newton method (`true`) that uses fewer evaluations of the Jacobian, but may take more iterations. Alternatively, the standard Newton's method is used. Default value: `false`.
- `start_point::Union{Nothing, Vector{<:Real}}` (optional): Either, we can explicitly and manually provide a start point for the deformation path or the realization in `G` is automatically selected. Default value: `nothing`.

# Examples
```julia-repl
julia> F = Framework([(1,2)], [0 0; 1 0;])
julia> motion_samples = [[0,0,cos(θ),sin(θ)] for θ in 0:0.025:pi/2]
julia> DeformationPath(F.G, motion_samples)
Deformation Path:
        Constraint System:
                Vertices:       [1, 2]
                Equations:      [-1.0 + (x₁₋₁ - x₁₋₂)^2 + (x₂₋₁ - x₂₋₂)^2]
                Realization:    0.0     1.0
                                0.0     0.0
        Motion:         [
                                [0.0, 0.0, 1.0, ...],
                                [0.0, 0.0, 0.9996875162757026, ...],
                                [0.0, 0.0, 0.9987502603949663, ...],
                                [0.0, 0.0, 0.9971888181122075, ...],
                        ...]

```

```julia-repl
julia> F = Framework([(1,2),(2,3),(3,4),(4,1)], Matrix([0 0; 1 0; 1 1; 0 1;]'))
julia> DeformationPath(F.G, [], 4, typeof(F); step_size=0.1)
Deformation Path:
        Constraint System:
                Vertices:       [1, 2, 3, 4]
                Equations:      [-1.0 + (x₁₋₁ - x₁₋₂)^2 + (x₂₋₁ - x₂₋₂)^2, -1.0 + (x₁₋₂ - x₁₋₃)^2 + (x₂₋₂ - x₂₋₃)^2, ...]
                Realization:    0.0     0.0
                                1.0     0.0
                                1.0     1.0
                                0.0     1.0
        Motion:         [
                                [0.0, 0.0, 1.0, ...],
                                [0.03651261273911354, 0.03651261273911352, 1.0340219488467464, ...],
                                [0.07533400807826349, 0.07533400807826351, 1.0653837593220445, ...],
                                [0.11627080392220492, 0.11627080392220503, 1.0939292079679943, ...],
                        ...]
        Step Size:      0.1
        Flex Selector:  [1.0]

```

```julia-repl
julia> F = BodyHinge([[1,2,3],[1,3,4],[1,4,5],[1,5,2]], Matrix([0 0 1; 1 -1 0; 1 1 0; -1 1 0; -1 -1 0;]'))
julia> DeformationPath(F, [], 5; step_size=0.01)
Deformation Path:
        Constraint System:
                Vertices:       [1, 2, 3, 4, 5]
                Equations:      [-4.0 + (x₁₋₄ - x₁₋₅)^2 + (x₂₋₄ - x₂₋₅)^2 + (x₃₋₄ - x₃₋₅)^2, -3.0 + (x₁₋₁ - x₁₋₂)^2 + (x₂₋₁ - x₂₋₂)^2 + (x₃₋₁ - x₃₋₂)^2, ...]
                Realization:    0.0     0.0     1.0
                                1.0     -1.0    0.0
                                1.0     1.0     0.0
                                -1.0    1.0     0.0
        Motion:         [
                                [0.0, 0.0, 1.0, ...],
                                [4.508862520636648e-17, 9.98697453244557e-18, 1.0001664436067044, ...],
                                [-1.4097240946269324e-17, 4.145853894815855e-17, 1.0006643327504479, ...],
                                [-5.097629777228046e-17, 3.7479573223418395e-17, 1.0014893737348445, ...],
                        ...]
        Step Size:      0.05
        Flex Selector:  [1.0]

```
"""
mutable struct DeformationPath
    G::ConstraintSystem
    step_size::Real    
    motion_samples::Vector{<:Vector{<:Real}}
    motion_matrices::Vector{<:Matrix{<:Real}}
    flex_mult::Vector{<:Real}
    _contacts::Vector
      
    """
        DeformationPath(G, motion_samples[; tol])
        DeformationPath(G, flex_mult, num_steps, type[; step_size, tol, random_flex, symmetric_newton, start_point])

    Constructor for deformation paths when a deformation is already known.

    # Arguments
    - `G::ConstraintSystem`: The underlying geometric constraint system.
    - `motion_samples::Vector{<:Vector{<:Real}}`: List of previously computed realizations in array format.
    - `tol::Real` (optional): Numerical tolerance for the approximation that is used for asserting the correctness of the approximation. Default value: `1e-8`.

    # Returns
    - `DeformationPath` 

    # Examples
    ```julia-repl
    julia> F = Framework([(1,2)], [0 0; 1 0;])
    julia> motion_samples = [[0,0,cos(θ),sin(θ)] for θ in 0:0.025:pi/2]
    julia> DeformationPath(F.G, motion_samples)
    Deformation Path:
            Constraint System:
                    Vertices:       [1, 2]
                    Equations:      [-1.0 + (x₁₋₁ - x₁₋₂)^2 + (x₂₋₁ - x₂₋₂)^2]
                    Realization:    0.0     1.0
                                    0.0     0.0
            Motion:         [
                                    [0.0, 0.0, 1.0, ...],
                                    [0.0, 0.0, 0.9996875162757026, ...],
                                    [0.0, 0.0, 0.9987502603949663, ...],
                                    [0.0, 0.0, 0.9971888181122075, ...],
                            ...]

    ```
    """
    function DeformationPath(G::ConstraintSystem, motion_samples::Vector{<:Vector{<:Real}}; tol::Real=1e-8, skip_check::Bool=false)::DeformationPath
        tol>0 || throw(error("The tolerance `tol` needs to be a positive number, but is $(tol)."))
        skip_check || (all(sample->norm(evaluate(G.equations, G.variables=>sample), Inf) < tol, motion_samples) || throw(error("The `motion_samples` do not satisfy the underlying constraints in the constraint system `G`!")))
        motion_matrices = [to_Matrix(G, Float64.(sample)) for sample in motion_samples]
        new(G, 0., motion_samples, motion_matrices, Vector{Float64}([]), [])
    end

    """
        DeformationPath(G, flex_mult, num_steps, type[; step_size, tol, random_flex, symmetric_newton, start_point])

    Computes an approximate continuous motion of a generic geometric constraint system `G`.

    # Arguments
    - `G::ConstraintSystem`: The underlying geometric constraint system.
    - `flex_mult::Vector`: The initial infinitesimal flex as a linear combination of the nontrivial infinitesimal flexes at the realization provided by the underlying geometric constraint system.
    - `num_steps::Int`: Number of steps the algorithm is supposed to take. 
    - `type::DataType`: Type of the geometric constraint system for computing the trivial infinitessimal flexes. Possible values: `"framework"`, `"angularframework"`, `"frameworkonsurface"`, `"hypergraph"`, `"polytope"`, `"diskpacking"`, `"bodyhinge"`, `"bodybar"` and `"sphericaldiskpacking"`.
    - `step_size::Real`: Step size of the deformation path. 
    - `tol::Real` (optional): Numerical tolerance for the approximation that is used for asserting the correctness of the approximation. Default value: `1e-8`.
    - `random_flex::Bool` (optional): If `flex_mult` is not provided (e.g. as `[]`), we can instead provide a random linear combination (`true`) or it is uniformly chosen as `[1,...,1]` (`false`). Default value: `false`.
    - `symmetric_newton::Bool` (optional): It is possible to use the slightly more efficient symmetric Newton method (`true`) that uses fewer evaluations of the Jacobian, but may take more iterations. Alternatively, the standard Newton's method is used. Default value: `false`.
    - `start_point::Union{Nothing, Vector{<:Real}}` (optional): Either, we can explicitly and manually provide a start point for the deformation path or the realization in `G` is automatically selected. Default value: `nothing`.

    # Returns
    - `DeformationPath` 

        # Examples
    ```julia-repl
    julia> F = Framework([(1,2),(2,3),(3,4),(4,1)], Matrix([0 0; 1 0; 1 1; 0 1;]'))
    julia> DeformationPath(F.G, [], 4, typeof(F); step_size=0.1)
    Deformation Path:
            Constraint System:
                    Vertices:       [1, 2, 3, 4]
                    Equations:      [-1.0 + (x₁₋₁ - x₁₋₂)^2 + (x₂₋₁ - x₂₋₂)^2, -1.0 + (x₁₋₂ - x₁₋₃)^2 + (x₂₋₂ - x₂₋₃)^2, ...]
                    Realization:    0.0     0.0
                                    1.0     0.0
                                    1.0     1.0
                                    0.0     1.0
            Motion:         [
                                    [0.0, 0.0, 1.0, ...],
                                    [0.03651261273911354, 0.03651261273911352, 1.0340219488467464, ...],
                                    [0.07533400807826349, 0.07533400807826351, 1.0653837593220445, ...],
                                    [0.11627080392220492, 0.11627080392220503, 1.0939292079679943, ...],
                            ...]
            Step Size:      0.1
            Flex Selector:  [1.0]

    ```
    """
    function DeformationPath(G::ConstraintSystem, flex_mult::Vector, num_steps::Int, type::DataType; show_progress::Bool=true, step_size::Real=1e-2, tol::Real=1e-13, random_flex::Bool=false, symmetric_newton::Bool=false, start_point::Union{Nothing, Vector{<:Real}}=nothing, time_penalty::Union{Real,Nothing}=2)::DeformationPath
        num_steps>=0 && step_size>=0 && tol>0 || throw(error("The `num_steps`, the `step_size` and `tol` needs to be a nonnegative, but are  $((num_steps, step_size, tol))."))
        if isnothing(start_point)
            start_point = to_Array(G, G.realization)
        else
            norm(evaluate(G.equations, G.variables=>start_point), Inf) < sqrt(tol) || throw(error("The `start_point` does not satisfy the underlying constraints in the constraint system `G`!"))
        end

        if type==Framework
            K_n = Framework([[i,j] for i in eachindex(G.vertices) for j in eachindex(G.vertices) if i<j], G.realization; pinned_vertices=G.pinned_vertices).G    
        elseif type==AngularFramework
            K_n = AngularFramework([[i,j,k] for i in eachindex(G.vertices) for j in eachindex(G.vertices) for k in eachindex(G.vertices) if (i<j && j<k) || (i<k && k<j) || (j<i && i<k)], G.realization; pinned_vertices=G.pinned_vertices).G
        elseif type==FrameworkOnSurface
            K_n = deepcopy(G)
            add_equations!(K_n, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(G.vertices) for j in eachindex(G.vertices) if i<j]])
        elseif type==VolumeHypergraph
            K_n = VolumeHypergraph(collect(powerset(G.vertices, G.dimension+1, G.dimension+1)), G.realization).G
        elseif type==Polytope || type==SpherePacking || type==BodyHinge || type==BodyBar
            K_n = ConstraintSystem(G.vertices, G.variables, vcat(G.equations, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(G.vertices) for j in eachindex(G.vertices) if i<j]]), G.realization, G.xs; pinned_vertices=G.pinned_vertices)
        elseif  type==SphericalDiskPacking
            minkowski_scalar_product(e1,e2) = e1'*e2-1
            inversive_distances = [minkowski_scalar_product(G.realization[:,contact[1]], G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(G.realization[:,contact[1]], G.realization[:,contact[1]]) * minkowski_scalar_product(G.realization[:,contact[2]], G.realization[:,contact[2]])) for contact in powerset(G.vertices, 2, 2)]
            K_n = ConstraintSystem(G.vertices, G.variables, [minkowski_scalar_product(G.xs[:,contact[1]], G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(G.xs[:,contact[1]], G.xs[:,contact[1]]) * minkowski_scalar_product(G.xs[:,contact[2]], G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(G.vertices, 2, 2))], G.realization, G.xs)
        else
            throw("The type must either be 'framework', 'frameworkonsurface', 'diskpacking', 'sphericaldiskpacking', 'hypergraph', 'bodyhinge', 'bodybar' or 'polytope', but is '$(type)'.")
        end
        flex_space = compute_nontrivial_inf_flexes(G, start_point, K_n)
        if flex_mult==[]
            if random_flex
                flex_mult = randn(Float64, size(flex_space)[2])
                flex_mult = flex_mult ./ norm(flex_mult, 1)
            else
                flex_mult = [1/size(flex_space)[2] for _ in axes(flex_space,2)]
            end
        else
            flex_mult = Float64.(flex_mult)
        end
        size(flex_space)[2]==length(flex_mult) || throw("The length of 'flex_mult' match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2]).")
        prev_flex = length(flex_mult)==0 ? [0 for _ in axes(flex_space,1)] : sum(flex_mult[i] .* flex_space[:,i] for i in eachindex(flex_mult))
        prev_flex = prev_flex ./ norm(prev_flex)
        
        motion_samples = [Float64.(start_point)]
        @showprogress enabled=show_progress for i in 1:num_steps
            try
                q, _prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end], K_n; tol=1e-5)
                prev_flex = _prev_flex
                if symmetric_newton
                    q = symmetric_newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                else
                    q = newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                end
                if isapprox(q, motion_samples[end]; atol=1e-10)
                    throw("Slow Progress detected.")
                end
                push!(motion_samples, q)
            catch e
                prev_flex, success = resolve_singularity(G, motion_samples, K_n, prev_flex, step_size; show_progress=show_progress, tol=tol, time_penalty=time_penalty, symmetric_newton=symmetric_newton)
                if !success || e == "The space of nontrivial infinitesimal motions is empty."
                    show_progress && @warn "The approximation of a deformation path ended prematurely."
                    break
                end
            end
        end
        motion_matrices = [to_Matrix(G, Float64.(sample)) for sample in motion_samples]
        new(G, step_size, motion_samples, motion_matrices, Vector{Float64}(flex_mult), [])
    end

    """
        DeformationPath(F::AllTypesWithoutSpherePacking, flex_mult, num_steps[; random_flex=false, kwargs...])

    Create an approximate continuous motion from a geometric constraint system object that is not a sticky sphere packing.

    # Arguments
    - `F::AllTypesWithoutSpherePacking`: The underlying geometric constraint system.
    - `flex_mult::Vector`: The initial infinitesimal flex as a linear combination of the nontrivial infinitesimal flexes at the realization provided by the underlying geometric constraint system.
    - `num_steps::Int`: Number of steps the algorithm is supposed to take. 
    - `random_flex::Bool` (optional): If `flex_mult` is not provided (e.g. as `[]`), we can instead provide a random linear combination (`true`) or it is uniformly chosen as `[1,...,1]` (`false`). Default value: `false`.

    For further arguments, see the base method [`DeformationPath(G::DeformationPaths.ConstraintSystem, motion_samples::Vector{<:Vector{<:Real}})`](@ref).
    """
    function DeformationPath(F::AllTypesWithoutSpherePacking, flex_mult::Vector, num_steps::Int; fast_search::Bool=false, show_progress::Bool=true, random_flex=false, kwargs...)::DeformationPath
        if flex_mult==[] && random_flex
            try
                flex_mult = compute_nonblocked_flex(F; fast_search=fast_search)
            catch e
                flex_mult = []
            end
            flex_mult = isempty(flex_mult) ? [] : flex_mult ./ norm(flex_mult)
        end 
        DeformationPath(F.G, flex_mult, num_steps, typeof(F); show_progress=show_progress, kwargs...)
    end

    """
        DeformationPath(F::SpherePacking, flex_mult, num_steps[; motion_samples, _contacts, step_size, prev_flex, tol, random_flex])

    Constructor method of an approximate motion for sticky sphere packings.

    # Arguments
    - `F::SpherePacking`: The underlying sticky sphere packing.
    - `flex_mult::Vector`: The initial infinitesimal flex as a linear combination of the nontrivial infinitesimal flexes at the realization provided by the underlying geometric constraint system.
    - `num_steps::Int`: Number of steps the algorithm is supposed to take. 
    - `motion_samples::Vector` (optional): Already computed motion samples. Default value: `[]`.
    - `_contacts::Vector` (optional): Current contacts between spheres. Default value: `[]`.
    - `step_size::Real` (optional): Step size of the deformation path. Default value: `1e-2`.
    - `prev_flex::Union{Nothing, Vector}` (optional): Previously computed infinitesimal flex. Default value: `nothing`.
    - `tol::Real` (optional): Numerical tolerance for the approximation that is used for asserting the correctness of the approximation. Default value: `1e-8`.
    - `random_flex::Bool` (optional): If `flex_mult` is not provided (e.g. as `[]`), we can instead provide a random linear combination (`true`) or it is uniformly chosen as `[1,...,1]` (`false`). Default value: `false`.

    # Returns
    - `DeformationPath` 
    """
    function DeformationPath(F::SpherePacking, flex_mult::Vector, num_steps::Int; show_progress::Bool=true, symmetric_newton::Bool=false, motion_samples::Vector=[], _contacts::Vector=[], step_size::Real=1e-2, prev_flex::Union{Nothing, Vector}=nothing, tol::Real=1e-13, random_flex::Bool=false, time_penalty::Union{Real,Nothing}=2)::DeformationPath
        start_point = to_Array(F, F.G.realization)
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(F.G.vertices) for j in eachindex(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
        if isnothing(prev_flex)
            flex_space = compute_nontrivial_inf_flexes(F.G, start_point, K_n)
            if flex_mult==[]
                if random_flex
                    flex_mult = randn(Float64, size(flex_space)[2])
                else
                    flex_mult = [1 for _ in axes(flex_space,2)]
                end
            end
            size(flex_space)[2]==length(flex_mult) || throw("The length of 'flex_mult' must match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2]).")
            prev_flex = sum(flex_mult[i] .* flex_space[:,i] for i in eachindex(flex_mult))
            prev_flex = prev_flex ./ norm(prev_flex)
        end
        if length(motion_samples)==0
            motion_samples = [Float64.(start_point)]
        end
        if length(_contacts)==0
            _contacts = [F.contacts]
        end

        @showprogress enabled=show_progress for i in 1:num_steps
            try
                q, prev_flex = euler_step(F.G, step_size, prev_flex, motion_samples[end], K_n; tol=1e-5)
                if symmetric_newton
                    q = symmetric_newton_correct(F.G, q; tol=tol, time_penalty=time_penalty)
                else
                    q = newton_correct(F.G, q; tol=tol, time_penalty=time_penalty)
                end

                if norm(q-motion_samples[end]) < 1e-6 || norm(q-motion_samples[end]) > step_size*4
                    throw("Slow Progress detected.")
                end

                cur_realization = to_Matrix(F,Float64.(q))
                if any(t->norm(cur_realization[:,t[1]] - cur_realization[:,t[2]]) < F.radii[t[1]] + F.radii[t[2]] - F.tolerance, powerset(F.G.vertices,2,2))
                    _F = SpherePacking(F.G.vertices, F.radii, cur_realization; pinned_vertices=F.G.pinned_vertices, tolerance=step_size)
                    DeformationPath(_F, flex_mult, num_steps-i; motion_samples=motion_samples, _contacts=_contacts, step_size=step_size, prev_flex=prev_flex, tol=tol)
                    break
                end

                push!(motion_samples, q)
                push!(_contacts, F.contacts)
            catch e
                # If Newton's method only diverges once and we are in a singularity,
                # we first try to reverse the previous flex before exiting the routine.
                prev_flex, success = resolve_singularity(F.G, motion_samples, K_n, prev_flex, step_size; show_progress=show_progress, tol=tol, time_penalty=time_penalty, symmetric_newton=symmetric_newton)
                if !success || e == "The space of nontrivial infinitesimal motions is empty."
                    show_progress && @warn "The approximation of a deformation path ended prematurely."
                    break
                end
            end
        end

        motion_matrices = [to_Matrix(F, Float64.(sample)) for sample in motion_samples]
        new(F.G, step_size, motion_samples, motion_matrices, Vector{Float64}(flex_mult), _contacts)
    end
end


"""
    resolve_singularity(G, motion_samples, failure_to_converge, motion_matrices, K_n, prev_flex, step_size[; tol, time_penalty, symmetric_newton])

Attempts to resolve a singularity at `motion_samples[end]`.
"""
function resolve_singularity(G::ConstraintSystem, motion_samples::Vector, K_n::ConstraintSystem, prev_flex::Vector, step_size::Real; show_progress::Bool=true, tol::Real=1e-10, time_penalty::Real=3, symmetric_newton::Bool=false)
    global failure_to_converge = 0
    global success = false
    _prev_flex = copy(prev_flex)
    while failure_to_converge < 3
        global failure_to_converge += 1
        if failure_to_converge==1 && length(motion_samples)>1
            show_progress && @info "Trying smaller step sizes."
            helper_samples = [motion_samples[end]]
            try
                for _ in 1:5
                    q, _prev_flex = euler_step(G, step_size/5, _prev_flex, helper_samples[end], K_n; tol=1e-5)
                    if symmetric_newton
                        q = symmetric_newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                    else
                        q = newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                    end
                    push!(helper_samples, q)
                end
                length(helper_samples)>1 && push!(motion_samples, helper_samples[end])
                if norm(helper_samples[end]-motion_samples[end]) < 1e-6 || norm(helper_samples[end]-motion_samples[end]) > step_size*4
                    _prev_flex = prev_flex
                    continue
                end
                global success = true
                break
            catch
                _prev_flex = prev_flex
                length(helper_samples)>1 && push!(motion_samples, helper_samples[end])
                continue
            end
        elseif failure_to_converge==1
            show_progress && @info "Direction is reversed."
            q = motion_samples[end] - step_size*_prev_flex
            try
                if symmetric_newton
                    q = symmetric_newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                else
                    q = newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                end
                push!(motion_samples, q)
                global success = true
                _prev_flex = -_prev_flex
                break
            catch
                continue
            end
        #=
        elseif failure_to_converge==2 && length(motion_samples)==3
            show_progress && @info "Using quadratic interpolation to predict the next point."
            interpolating_points = motion_samples[end-2:end]
            @var t, a[1:length(interpolating_points[3]), 1:2]
            quadric = interpolating_points[3]+t*a[:,1]+t^2*a[:,2]
            equations = vcat(evaluate(quadric, t=>-step_size) .- interpolating_points[2], evaluate(quadric, t=>-2*step_size) .- interpolating_points[1])
            a_sys = System(equations)
            next_a = real_solutions(solve(a_sys))[1]
            a_variables = variables(a_sys)
            q = evaluate(quadric, vcat(t, a_variables)=>vcat(step_size, next_a))
            try
                if symmetric_newton
                    q = symmetric_newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                else
                    q = newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                end

                if norm(q-motion_samples[end]) < 1e-6 || norm(q-motion_samples[end]) > step_size*5
                    _prev_flex = prev_flex
                    continue
                end

                push!(motion_samples, q)
                push!(motion_matrices, to_Matrix(G, Float64.(q)))
                global success = true
                break
            catch
                continue
            end
        elseif failure_to_converge==2 && length(motion_samples)>=4
            show_progress && @info "Using cubic interpolation to predict the next point."
            interpolating_points = motion_samples[end-3:end]
            @var t, a[1:length(interpolating_points[4]), 1:3]
            cubic = interpolating_points[4]+t*a[:,1]+t^2*a[:,2]+t^3*a[:,3]
            equations = vcat(evaluate(cubic, t=>-step_size) .- interpolating_points[3], evaluate(cubic, t=>-2*step_size) .- interpolating_points[2], evaluate(cubic, t=>-3*step_size) .- interpolating_points[1])
            a_sys = System(equations)
            next_a = real_solutions(solve(a_sys))[1]
            a_variables = variables(a_sys)
            q = evaluate(cubic, vcat(t,a_variables)=>vcat(step_size,next_a))
            try
                if symmetric_newton
                    q = symmetric_newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                else
                    q = newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                end

                if norm(q-motion_samples[end]) < 1e-6 || norm(q-motion_samples[end]) > step_size*5
                    _prev_flex = prev_flex
                    continue
                end

                push!(motion_samples, q)
                global success = true
                break
            catch
                _prev_flex = prev_flex
                continue
            end=#
        else
            show_progress && @info "Acceleration-based cusp method is being used."
            flexes = compute_nontrivial_inf_flexes(G, motion_samples[end], K_n; tol=1e-3)
            if size(flexes)[2]==0
                break
            end
            projection = flexes*pinv(flexes'*flexes)*flexes'
            J = evaluate.(G.jacobian, G.variables=>motion_samples[end])
            acceleration = pinv(J) * (-evaluate.(G.jacobian, G.variables=>_prev_flex)*_prev_flex)
            acceleration, _prev_flex = acceleration ./ norm(acceleration), _prev_flex ./ norm(_prev_flex)
            _prev_flex = - _prev_flex - acceleration
            _prev_flex = projection*_prev_flex
            _prev_flex = _prev_flex ./ norm(_prev_flex)
            q = motion_samples[end] + step_size*_prev_flex
            try
                if symmetric_newton
                    q = symmetric_newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                else
                    q = newton_correct(G, q; tol=tol, time_penalty=time_penalty)
                end

                if norm(q-motion_samples[end]) < 1e-6 || norm(q-motion_samples[end]) > step_size*5
                    _prev_flex = prev_flex
                    continue
                end

                push!(motion_samples, q)
                global success = true
                break
            catch
                _prev_flex = prev_flex
                continue
            end
        end
    end
    return _prev_flex, success
end


"""
    DeformationPath_EdgeContraction(F::Polytope, edge_for_contraction, contraction_target, contraction_start[; kwargs...])

Create an approximate continuous motion from a `Polytope` object induced by continuously perturbing a single edge given by `edge_for_contraction` from the length `contraction_start` to `contraction_target`.

# Arguments
- `F::Polytope`: The underlying polytope.
- `edge_for_contraction::Union{Tuple{Int,Int},Vector{Int}}`: The edge in `F` that is supposed to be contracted.
- `contraction_target::Real`: To what length ratio the edge is suppsed to be contracted.
- `contraction_start::Real`: The edge length where the induced deformation path begins.

For additional keywords, see [`DeformationPath_EdgeContraction(F::Polytope, edge_for_contraction::Union{Tuple{Int,Int},Vector{Int}}, contraction_target::Real)`](@ref).
"""
function DeformationPath_EdgeContraction(F::Polytope, edge_for_contraction::Union{Tuple{Int,Int},Vector{Int}}, contraction_target::Real, contraction_start::Real; kwargs...)::DeformationPath
    D1 = DeformationPath_EdgeContraction(F, edge_for_contraction, contraction_start; kwargs...)
    D2 = DeformationPath_EdgeContraction(F, edge_for_contraction, contraction_target; contraction_start=contraction_start, realization_start=D1.motion_samples[end], kwargs...)
    return D2
end


"""
    stich_deformation_paths(D1, D2)

Combine the `motion_samples` of the two deformation paths `D1` and `D2`.

We assume that both deformation paths start in the same configuration. The `reversed=true`
keyword then reverses the first deformation path and deletes the first `motion_sample` of
the second deformation path. In the case that `reversed=false`, the last realization of 
the first deformation path needs to be the first realization of the second deformation path.

# Arguments
- `D1::DeformationPath`: First `DeformationPath`.
- `D2::DeformationPath`: Second `DeformationPath`.
- `reversed::Bool` (optional): Determines whether the first deformation path is reversed or not. Default value: `true`.
"""
function stich_deformation_paths(D1::DeformationPath, D2::DeformationPath; reversed::Bool=true)::DeformationPath
    D1.G == D2.G || throw(error("The deformation paths need to have the same underlying constraint system."))
    if reversed
        isapprox(norm(D1.motion_samples[1]-D2.motion_samples[1]), 0, atol=1e-4) || throw(error("The motion samples are not compatible."))
        return DeformationPath(D1.G, vcat(reverse(D1.motion_samples), D2.motion_samples[2:end]); skip_check=true)
    else
        isapprox(norm(D1.motion_samples[end]-D2.motion_samples[1]), 0, atol=1e-4) || throw(error("The motion samples are not compatible."))
        return DeformationPath(D1.G, vcat(D1.motion_samples, D2.motion_samples[2:end]); skip_check=true)
    end
end


"""
    DeformationPath_EdgeContraction(F::Polytope, edge_for_contraction, contraction_target[; kwargs...])

Create an approximate continuous motion from a `Polytope` object induced by continuously perturbing a single edge given by `edge_for_contraction` from its original length to the length `contraction_target`.

# Arguments
- `F::Polytope`: The underlying polytope.
- `edge_for_contraction::Union{Tuple{Int,Int},Vector{Int}}`: The edge in `F` that is supposed to be contracted.
- `contraction_target::Real`: To what length ratio the edge is suppsed to be contracted.
- `step_size::Real` (optional): Step size of the deformation path. 
- `tol::Real` (optional): Numerical tolerance for the approximation that is used for asserting the correctness of the approximation. Default value: `1e-8`.
"""
function DeformationPath_EdgeContraction(F::Polytope, edge_for_contraction::Union{Tuple{Int,Int},Vector{Int}}, contraction_target::Real; contraction_start::Union{Real,Nothing}=nothing, realization_start::Union{Nothing,Vector}=nothing, show_progress::Bool=true, step_size::Real=0.002, tol::Real=1e-11, time_penalty::Union{Real,Nothing}=2)::DeformationPath
    edge_for_contraction = [edge_for_contraction[1], edge_for_contraction[2]]
    length(edge_for_contraction)==2 && (edge_for_contraction in [[edge[1],edge[2]] for edge in F.edges] || [edge_for_contraction[2], edge_for_contraction[1]] in [[edge[1],edge[2]] for edge in F.edges]) || throw(error("The `edge_for_contraction` needs to be an edge of the polytope's 1-skeleton!"))
    @var c
    edge_equation = sum( (F.G.xs[:,edge_for_contraction[1]]-F.G.xs[:,edge_for_contraction[2]]) .^2) - sum( (F.G.realization[:,edge_for_contraction[1]]-F.G.realization[:,edge_for_contraction[2]]) .^2)
    edge_variables = variables(edge_equation)
    generic_point = randn(ComplexF64, size(F.G.xs)[1]*2)
    evaluated_edge_equation = evaluate(edge_equation, edge_variables=>generic_point)
    corresponding_equation_index = findfirst(eq->isa(evaluate(eq, edge_variables=>generic_point), ComplexF64) && isapprox(evaluate(eq, edge_variables=>generic_point), evaluated_edge_equation), F.G.equations)
    _G = deepcopy(F.G)
    _G.equations[corresponding_equation_index] = sum( (_G.xs[:,edge_for_contraction[1]]-_G.xs[:,edge_for_contraction[2]]) .^2) - c^2
    if !isnothing(contraction_start) && !isnothing(realization_start)
        start_c_value = contraction_start
        motion_samples = [realization_start]
    else
        start_c_value = sqrt( sum( (_G.realization[:,edge_for_contraction[1]]-_G.realization[:,edge_for_contraction[2]]) .^2) )
        motion_samples = [to_Array(_G, _G.realization)]
    end
    index = 1
    if contraction_target > start_c_value
        local_step_size = step_size
    elseif contraction_target <= start_c_value && contraction_target > 0
        local_step_size = -step_size
    else
        throw(error("The `contraction_target` needs to be bigger than 0, but is $(contraction_target)"))
    end
    while true
        show_progress && println("Trial $index")
        index = index+1
        try
            cur_point = motion_samples[end] + 0.05*(rand(Float64,length(motion_samples[end]))-[0.5 for i in eachindex(motion_samples[end])])
            local_equations = evaluate(_G.equations, c => start_c_value + local_step_size)
            cur_point = newton_correct(local_equations, _G.variables, _G.jacobian, cur_point; tol=tol, time_penalty=time_penalty, armijo_linesearch=false)
            push!(motion_samples, cur_point)
            break
        catch err
            show_progress && println(err)
            continue
        end
    end

    @showprogress enabled=show_progress for step in start_c_value+local_step_size:local_step_size:contraction_target
        local_equations = evaluate(_G.equations, c=>step)
        local_jacobian = _G.jacobian
        try
            cur_point = newton_correct(local_equations, _G.variables, local_jacobian, motion_samples[end]+(motion_samples[end]-motion_samples[end-1]); tol=tol, time_penalty=time_penalty, armijo_linesearch=false)
            push!(motion_samples, cur_point)
        catch e
            if norm(motion_samples[end]-motion_samples[end-1])>step_size*20
                deleteat!(motion_samples, length(motion_samples))
                show_progress && @warn "The approximation of a deformation path ended prematurely. $e"
                break
            else
                continue
            end
            #Does this work?
            interpolating_points = motion_samples[end-2:end]
            @var t, a[1:length(interpolating_points[3]), 1:2]
            quadric = a[:,1]+t*a[:,1]+t^2*a[:,2]
            equations = vcat(evaluate(quadric, t=>step-local_step_size), evaluate(quadric, t=>step-2*local_step_size))
            a_variables = vcat(a...)
            a_matrix = evaluate(differentiate(equations, a_variables), a_variables=>[0 for _ in a_variables])
            a_b_vector = vcat(interpolating_points[2], interpolating_points[1])
            a_solutions = a_matrix \ a_b_vector
            next_point = evaluate(quadric, vcat([t],a_variables)=>vcat([step], a_solutions))
            local_equations = evaluate(_G.equations, c=>step)
            try
                cur_point = newton_correct(local_equations, _G.variables, local_jacobian, next_point; tol=tol, time_penalty=time_penalty)
                push!(motion_samples, cur_point)
            catch err
                show_progress && @warn "The approximation of a deformation path ended prematurely. $err"
                break
            end
        end
    end
    _G.equations = _G.equations[filter(i->i!=corresponding_equation_index, 1:length(_G.equations))]
    return DeformationPath(_G, motion_samples)
end


function Base.show(io::IO, D::DeformationPath)
    """Custom display method for a `DeformationPath`."""
    print(io, "Deformation Path:\n")
    print(io, "\t$(D.G)")
    print(io, "\tMotion:\t\t[\n")
    for sample in D.motion_samples[1:min(4, length(D.motion_samples))]
        print(io,"\t\t\t\t[$(sample[1]), $(length(sample)>=2 ? "$(sample[2]), " : "]")$(length(sample)>=3 ? "$(sample[3]), ...]" : "]"),\n")
    end
    print(io,"\t\t\t...]")
    if D.step_size > 0.0
        print(io, "\n\tStep Size:\t$(D.step_size)")
    end
    if !(isempty(D.flex_mult))
        print(io,"\n\tFlex Selector:\t$(D.flex_mult)")
    end
end

include("Visualization.jl")

end 