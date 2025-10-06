module DeformationPaths

import HomotopyContinuation: evaluate, differentiate, newton, Expression, Variable, @var, real_solutions, System, solve, variables
import LinearAlgebra: norm, pinv, nullspace, rank, qr, zeros, inv, cross, det, svd, I, zeros
import GLMakie: NoShading, GeometryBasics, Vec3f, meshscatter!, surface!, Sphere, mesh!, @lift, poly!, text!, Figure, record, hidespines!, hidedecorations!, lines!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Observable, Point3f, Point2f, connect, faces, Mesh, mesh
import ProgressMeter: @showprogress
import Combinatorics: powerset
import Colors: distinguishable_colors, red, green, blue, colormap, RGB
import MarchingCubes: MC, march, makemesh
import Polyhedra
import CDDLib
import Base:show

include("GeometricConstraintSystem.jl")
using .GeometricConstraintSystem: BodyHinge, ConstraintSystem, Framework, equations!, realization!, to_Array, to_Matrix, VolumeHypergraph, plot, Polytope, SpherePacking, SphericalDiskPacking, FrameworkOnSurface, add_equations!, AngularFramework, compute_nontrivial_inf_flexes, fix_antipodals!

export  ConstraintSystem, 
        Framework,
        AngularFramework,
        DeformationPath,
        VolumeHypergraph,
        animate,
        plot,
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
        compute_nontrivial_inf_flexes,
        fix_antipodals!,
        tetrahedral_symmetry!

"""
Class for constructing approximate deformation paths.

# Attributes
- `G::ConstraintSystem`: The geometric constraint system for which the deformation is computes.
- `step_size::Number`: The step size of the deformation path. 
- `motion_samples::Vector{Vector{Float64}}`: Points in the configuration space representing the approximate motion.
- `motion_matrices::Vector{Matrix{Float64}}`: The points in `motion_samples` as distributed into realizations given by `dxn` matrices for the dimension `d` and the number of vertices `n`.
- `flex_mult::Vector{Float64}`: The initial infinitesimal flex as a linear combination of the nontrivial infinitesimal flexes at the realization provided by the underlying geometric constraint system.
. `_contacts::Vector`: In some cases (such as sphere packings), the contacts can change during the deformation. This is reflected by this attribute.
"""
mutable struct DeformationPath
    G::ConstraintSystem
    step_size::Number
    motion_samples::Vector{Vector{Float64}}
    motion_matrices::Vector{Matrix{Float64}}
    flex_mult::Vector{Float64}
    _contacts::Vector
    
    """
        DeformationPath(G, motion_samples[; tol])

    Constructor for deformation paths when a deformation is already known.

    # Arguments
    - `G::ConstraintSystem`: The underlying geometric constraint system.
    - `motion_samples::Vector{Vector{Float64}}`: List of previously computed realizations in array format.
    - `tol::Float64` (optional): Numerical tolerance for the approximation that is used for asserting the correctness of the approximation. Default value: `1e-8`.

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
                                [0.0, 0.0, 1.0, 0.0, ...],
                                [0.0, 0.0, 0.9996875162757026, 0.024997395914712332, ...],
                                [0.0, 0.0, 0.9987502603949663, 0.04997916927067833, ...],
                                [0.0, 0.0, 0.9971888181122075, 0.07492970727274235, ...],

                        ...]

    ```
    """
    function DeformationPath(G::ConstraintSystem, motion_samples::Vector{Vector{Float64}}; tol::Float64=1e-8)::DeformationPath
        all(sample->norm(evaluate(G.equations, G.variables=>sample), Inf) < tol, motion_samples) || throw(error("The `motion_samples` do not satisfy the underlying constraints in the Constraint System `G`!"))
        motion_matrices = [to_Matrix(G, Float64.(sample)) for sample in motion_samples]
        new(G, 0., motion_samples, motion_matrices, Vector{Float64}([]), [])
    end

    function DeformationPath(G::ConstraintSystem, flex_mult::Vector, num_steps::Int, type::String; step_size::Number=1e-2, newton_tol=1e-14, random_flex=false, symmetric_newton=false, start_point=nothing)::DeformationPath
        println("$step_size, $newton_tol, $flex_mult")
        if start_point == nothing
            start_point = to_Array(G, G.realization)
        end
        if type=="framework"
            K_n = Framework([[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j], G.realization; pinned_vertices=G.pinned_vertices).G
        elseif type=="angularframework"
            K_n = AngularFramework([[i,j,k] for i in 1:length(G.vertices) for j in 1:length(G.vertices) for k in 1:length(G.vertices) if (i<j && j<k) || (i<k && k<j) || (j<i && i<k)], G.realization; pinned_vertices=G.pinned_vertices).G
        elseif type=="frameworkonsurface"
            K_n = deepcopy(G)
            add_equations!(K_n, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j]])
        elseif type=="hypergraph"
            K_n = VolumeHypergraph(collect(powerset(G.vertices, G.dimension+1, G.dimension+1)), G.realization).G
        elseif type=="polytope" || type == "diskpacking" || type == "bodyhinge"
            K_n = ConstraintSystem(G.vertices, G.variables, vcat(G.equations, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j]]), G.realization, G.xs; pinned_vertices=G.pinned_vertices)
        elseif  type=="sphericaldiskpacking"
            minkowski_scalar_product(e1,e2) = e1'*e2-1
            inversive_distances = [minkowski_scalar_product(G.realization[:,contact[1]], G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(G.realization[:,contact[1]], G.realization[:,contact[1]]) * minkowski_scalar_product(G.realization[:,contact[2]], G.realization[:,contact[2]])) for contact in powerset(G.vertices, 2, 2)]
            K_n = ConstraintSystem(G.vertices, G.variables, [minkowski_scalar_product(G.xs[:,contact[1]], G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(G.xs[:,contact[1]], G.xs[:,contact[1]]) * minkowski_scalar_product(G.xs[:,contact[2]], G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(G.vertices, 2, 2))], G.realization, G.xs)
        else
            throw("The type must either be 'framework', 'frameworkonsurface', 'diskpacking', 'sphericaldiskpacking', 'hypergraph', 'bodyhinge' or 'polytope', but is '$(type)'.")
        end

        flex_space = compute_nontrivial_inf_flexes(G, start_point, K_n)
        if flex_mult==[]
            if random_flex
                flex_mult = randn(Float64, size(flex_space)[2])
                flex_mult = flex_mult ./ norm(flex_mult, 1)
            else
                flex_mult = [1/size(flex_space)[2] for _ in 1:size(flex_space)[2]]
            end
        else
            flex_mult = Float64.(flex_mult)
        end
        size(flex_space)[2]==length(flex_mult) || throw("The length of 'flex_mult' match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2]).")
        prev_flex = length(flex_mult)==0 ? [0 for _ in 1:size(flex_space)[1]] : sum(flex_mult[i] .* flex_space[:,i] for i in 1:length(flex_mult))
        prev_flex = prev_flex ./ norm(prev_flex)
        motion_samples, motion_matrices = [Float64.(start_point)], [to_Matrix(G, Float64.(start_point))]
        failure_to_converge = 0
        @showprogress for i in 1:num_steps
            try
                q, prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end], K_n)
                if symmetric_newton
                    q = symmetric_newton_correct(G, q; tol=newton_tol)
                else
                    q = newton_correct(G, q; tol=newton_tol)
                end
                failure_to_converge = 0
                if isapprox(q, motion_samples[end]; atol=1e-12)
                    throw("Slow Progress detected.")
                end
                push!(motion_samples, q)
                push!(motion_matrices, to_Matrix(G, Float64.(q)))                   
            catch e
                i = i - 1
                @warn exception=(e, catch_backtrace())
                if failure_to_converge >= 3 || e == "The space of nontrivial infinitesimal motions is empty."
                    break
                else
                    # If Newton's method only diverges once and we are in a singularity,
                    # we first try to reverse the previous flex before exiting the routine.
                    failure_to_converge += 1
                    if failure_to_converge==1
                        try
                            q, prev_flex = euler_step(G, step_size/4, prev_flex, motion_samples[end], K_n)
                            if symmetric_newton
                                q = symmetric_newton_correct(G, q; tol=newton_tol)
                            else
                                q = newton_correct(G, q; tol=newton_tol)
                            end
                            push!(motion_samples, q)
                            push!(motion_matrices, to_Matrix(G, Float64.(q)))                    
                        catch
                            continue
                        end
                    elseif length(motion_samples)==1
                        @info "Direction was reversed."
                        prev_flex = -prev_flex
                    else 
                        @info "Acceleration-based cusp method is being used."
                        J = evaluate.(G.jacobian, G.variables=>motion_samples[end])
                        acceleration = -pinv(J)*evaluate.(G.jacobian, G.variables=>prev_flex)*prev_flex
                        acceleration = acceleration ./ norm(acceleration)
                        prev_flex = - prev_flex - acceleration
                        prev_flex = prev_flex ./ norm(prev_flex)
                    end
                end
            end
        end
        new(G, step_size, motion_samples, motion_matrices, Vector{Float64}(flex_mult), [])
    end

    function DeformationPath(F::Framework, flex_mult::Vector, num_steps::Int; random_flex=false, kwargs...)::DeformationPath
        if flex_mult==[] && random_flex
            flex_mult = compute_nonblocked_flex(F)
            flex_mult = flex_mult ./ norm(flex_mult)
        end
        DeformationPath(F.G, flex_mult, num_steps, "framework"; kwargs...)::DeformationPath
    end

    function DeformationPath(F::AngularFramework, flex_mult::Vector, num_steps::Int; kwargs...)::DeformationPath
        DeformationPath(F.G, flex_mult, num_steps, "angularframework"; kwargs...)
    end

    function DeformationPath(F::FrameworkOnSurface, flex_mult::Vector, num_steps::Int; kwargs...)::DeformationPath
        DeformationPath(F.G, flex_mult, num_steps, "frameworkonsurface"; kwargs...)
    end

    function DeformationPath(F::VolumeHypergraph, flex_mult::Vector, num_steps::Int; kwargs...)::DeformationPath
        DeformationPath(F.G, flex_mult, num_steps, "hypergraph"; kwargs...)
    end

    function DeformationPath(F::Polytope, flex_mult::Vector, num_steps::Int; random_flex::Bool=false, kwargs...)::DeformationPath
        if flex_mult==[] && random_flex
            try
                flex_mult = compute_nonblocked_flex(F)
            catch e
                @warn e
                flex_mult = []
            end
            flex_mult = isempty(flex_mult) ? [] : flex_mult ./ norm(flex_mult)
        end
        DeformationPath(F.G, flex_mult, num_steps, "polytope"; start_point=to_Array(F, F.G.realization), random_flex=random_flex, kwargs...)
    end

    function DeformationPath(F::Polytope, edge_for_contraction::Union{Tuple{Int,Int},Vector{Int}}, contraction_target::Number; step_size::Number=0.002, tol::Number=1e-10, kwargs...)::DeformationPath
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
        start_c_value = sqrt( sum( (_G.realization[:,edge_for_contraction[1]]-_G.realization[:,edge_for_contraction[2]]) .^2) )
        
        motion_samples = [to_Array(_G, _G.realization)]
        while true
            println("Trial")
            try
                cur_point = motion_samples[end] + 0.01*(rand(Float64,length(motion_samples[end]))-[0.5 for i in 1:length(motion_samples[end])])
                local_equations = evaluate(_G.equations, c => start_c_value - step_size)
                cur_point = newton_correct(local_equations, _G.variables, _G.jacobian, cur_point; tol=tol, time_penalty=1)
                push!(motion_samples, cur_point)
                break
            catch e
                println(e)
                continue
            end
        end

        for step in step_size:step_size:(start_c_value*(1-contraction_target))
            println("Current step: $step")
            local_equations = evaluate(_G.equations, c=>start_c_value-step)
            local_jacobian = _G.jacobian
            try
                cur_point = newton_correct(local_equations, _G.variables, local_jacobian, motion_samples[end]+(step_size/2)*(motion_samples[end]-motion_samples[end-1])/norm(motion_samples[end]-motion_samples[end-1]); tol=tol, time_penalty=1)
                push!(motion_samples, cur_point)
            catch e
                println(e)
                break
            end
        end
        _G.equations = _G.equations[filter(i->i!=corresponding_equation_index, 1:length(_G.equations))]
        DeformationPath(_G, motion_samples)
    end


    function DeformationPath(F::BodyHinge, flex_mult::Vector, num_steps::Int; kwargs...)::DeformationPath
        DeformationPath(F.G, flex_mult, num_steps, "bodyhinge"; kwargs...)
    end

    function DeformationPath(F::SphericalDiskPacking, flex_mult::Vector, num_steps::Int; kwargs...)::DeformationPath
        DeformationPath(F.G, flex_mult, num_steps, "sphericaldiskpacking"; kwargs...)
    end

    """
        DeformationPath(F, flex_mult, num_steps[; motion_samples, _contacts, step_size, prev_flex, newton_tol, random_flex])
    Constructor method of an approximate motion for sticky sphere packings.
    """
    function DeformationPath(F::SpherePacking, flex_mult::Vector, num_steps::Int; motion_samples::Vector=[], _contacts::Vector=[], step_size::Number=1e-2, prev_flex::Union{Nothing, Vector}=nothing, newton_tol=1e-14, random_flex=false)::DeformationPath
        start_point = to_Array(F, F.G.realization)
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
        if prev_flex == nothing
            flex_space = compute_nontrivial_inf_flexes(F.G, start_point, K_n)
            if flex_mult==[]
                if random_flex
                    flex_mult = randn(Float64, size(flex_space)[2])
                else
                    flex_mult = [1 for _ in 1:size(flex_space)[2]]
                end
            end
            size(flex_space)[2]==length(flex_mult) || throw("The length of 'flex_mult' must match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2]).")
            prev_flex = sum(flex_mult[i] .* flex_space[:,i] for i in 1:length(flex_mult))
            prev_flex = prev_flex ./ norm(prev_flex)
        end
        if length(motion_samples)==0
            motion_samples = [Float64.(start_point)]
        end
        if length(_contacts)==0
            _contacts = [F.contacts]
        end

        failure_to_converge = 0
        @showprogress for i in 1:num_steps
            try
                q, prev_flex = euler_step(F.G, step_size, prev_flex, motion_samples[end], K_n)
                q = newton_correct(F.G, q; tol=newton_tol)
                failure_to_converge = 0
                if isapprox(q, motion_samples[end]; atol=1e-12)
                    throw("Slow Progress detected.")
                end

                cur_realization = to_Matrix(F,Float64.(q))
                if any(t->norm(cur_realization[:,t[1]] - cur_realization[:,t[2]]) < F.radii[t[1]] + F.radii[t[2]] - F.tolerance, powerset(F.G.vertices,2,2))
                    _F = SpherePacking(F.G.vertices, F.radii, cur_realization; pinned_vertices=F.G.pinned_vertices, tolerance=step_size)
                    DeformationPath(_F, flex_mult, num_steps-i; motion_samples=motion_samples, _contacts=_contacts, step_size=step_size, prev_flex=prev_flex, newton_tol=newton_tol)
                    break
                end

                push!(motion_samples, q)
                push!(_contacts, F.contacts)    
            catch e
                i = i - 1
                @warn e
                if failure_to_converge == 3 || e == "The space of nontrivial infinitesimal motions is empty."
                    break
                else
                    # If Newton's method only diverges once and we are in a singularity,
                    # we first try to reverse the previous flex before exiting the routine.
                    failure_to_converge += 1
                    if failure_to_converge==1
                        try
                            q, prev_flex = euler_step(F.G, step_size/3, prev_flex, motion_samples[end], K_n)
                            q = newton_correct(F.G, q; tol=newton_tol)
                            push!(motion_samples, q)
                            push!(motion_matrices, to_Matrix(F, Float64.(q)))                    
                        catch
                            continue
                        end
                    elseif length(motion_samples)==1
                        @warn "Direction was reversed."
                        prev_flex = -prev_flex
                    else
                        @info "Acceleration-based cusp method is being used."
                        J = evaluate.(F.G.jacobian, F.G.variables=>motion_samples[end])
                        acceleration = -pinv(J)*evaluate.(F.G.jacobian, F.G.variables=>prev_flex)*prev_flex
                        acceleration = acceleration ./ norm(acceleration)
                        prev_flex = - prev_flex - acceleration
                        prev_flex = prev_flex ./ norm(prev_flex)
                    end
                end
            end
        end

        motion_matrices = [to_Matrix(F, Float64.(sample)) for sample in motion_samples]
        new(F.G, step_size, motion_samples, motion_matrices, Vector{Float64}(flex_mult), _contacts)
    end

    """
        euler_step(G, step_size, prev_flex, point, K_n)

    Euler step predicting the next point along the approximate motion.

    # Returns
    - `predicted_point::Vector{Number}`: The next point predicted by Euler's method.
    - `predicted_inf_flex::Vector{Number}`: The tangent vector predicted by Euler's method.
    """
    function euler_step(G::ConstraintSystem, step_size::Float64, prev_flex::Vector{Number}, point::Vector{Number}, K_n::ConstraintSystem)::Tuple{Vector{Number}, Vector{Number}}
        J = evaluate(G.jacobian, G.variables=>point)
        flex_space = compute_nontrivial_inf_flexes(G, point, K_n; tol=1e-5)
        if size(flex_space)[2]==0
            throw("The space of nontrivial infinitesimal motions is empty.")
        end
        flex_coefficients = pinv(flex_space) * prev_flex
        predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
        predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
        return point+step_size*predicted_inf_flex, predicted_inf_flex
    end
end

function Base.show(io::IO, D::DeformationPath)
    print(io, "Deformation Path:\n")
    print(io, "\t$(D.G)")
    if D.step_size > 0.0
        print(io, "\tStep Size:\t\t$(D.step_size)\n")
    end
    print(io, "\tMotion:\t\t[\n")
    for sample in D.motion_samples[1:4]
        print(io,"\t\t\t\t[$(sample[1]), $(sample[2]), $(sample[3]), $(sample[4]), ...],\n")
    end
    print(io,"\n\t\t\t...]")
    if !(isempty(D.flex_mult))
        print(io,"\n\tFlex Selector:\t\t$(D.flex_mult)")
    end
end


"""
    newton_correct(G, point[; tol, time_penalty])

Apply Newton's method to correct `point` back to the constraints in `G`.

# Arguments
- `G::ConstraintSystem`: The underlying geometric constraint system.
- `point::Vector{Number}`: The initial point that Newton's method is applied to.
- `tol::Number` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Number` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.

# Returns
- `q::Vector{Number}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`

See also [`newton_correct`](@ref)
"""
function newton_correct(G::ConstraintSystem, point::Vector{Number}; tol::Number = 1e-13, time_penalty::Number=2)::Vector{Number}
    return newton_correct(G.equations, G.variables, G.jacobian, point; tol = tol, time_penalty=time_penalty)
end

"""
    newton_correct(equations, variables, jac, point[; tol, time_penalty])

Apply Newton's method to correct `point` back to the constraints in `equations`.

# Arguments
- `equations::Vector{Expression}`: Equations to correct `point` to.
- `variables::Vector{Variable}`: Variables from the affine coordinate ring.
- `jac::Matrix{Expression}`: Jacobian matrix corresponding to `equations` and `variables`.
- `point::Vector{Number}`: The initial point that Newton's method is applied to.
- `tol::Number` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Number` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.

# Returns
- `q::Vector{Number}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`
"""
function newton_correct(equations::Vector{Expression}, variables::Vector{Variable}, jac::Matrix{Expression}, point::Vector{Number}; tol::Number = 1e-13, time_penalty::Number=0.25)::Vector{Number}
    q = Base.copy(point)
    start_time=Base.time()
    global damping = 1
    while(norm(evaluate(equations, variables=>q)) > tol)
        J = evaluate.(jac, variables=>q)
        stress_dimension = size(nullspace(J'; atol=1e-8))[2]
        if stress_dimension > 0
            rand_mat = randn(Float64, length(equations) - stress_dimension, length(equations))
            new_equations = rand_mat*equations
            J = rand_mat*J
        else
            new_equations = equations
        end

        #Armijo Line Search
        r_val = evaluate(new_equations, variables=>q)
        v = -(J \ r_val)
        global damping = 1
        qnew = q + damping*v
        while norm(evaluate(equations, variables=>qnew)) > norm(evaluate(equations, variables=>q)) - 0.1 * damping * v'*v
            global damping = damping*0.6
            qnew = q + damping*v
            if damping < 1e-11 || Base.time()-start_time > length(point)/time_penalty
                throw("Newton's method did not converge in time. damping=$damping and time=$(Base.time()-start_time)")
            end
        end
        q = qnew
    end
    return q
end

"""
    symmetric_newton_correct(G, point[; tol, time_penalty])

Apply symmetric Newton's method to correct `point` back to the constraints in `equations`.

The symmetric Newton corrector evaluates the Jacobian matrix less often.

# Arguments
- `G::ConstraintSystem`: The underlying geometric constraint system.
- `point::Vector{Number}`: The initial point that Newton's method is applied to.
- `tol::Number` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Number` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.

# Returns
- `q::Vector{Number}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`

See also [`symmetric_newton_correct`](@ref)
"""
function symmetric_newton_correct(G::ConstraintSystem, point::Vector{Float64}; tol = 1e-13, time_penalty=2)::Vector{Number}
    return symmetric_newton_correct(G.equations, G.variables, G.jacobian, point; tol = tol, time_penalty=time_penalty)
end


"""
    symmetric_newton_correct(equations, variables, jac, point[; tol, time_penalty])

Apply symmetric Newton's method to correct `point` back to the constraints in `equations`.

The symmetric Newton corrector evaluates the Jacobian matrix less often.

# Arguments
- `equations::Vector{Expression}`: Equations to correct `point` to.
- `variables::Vector{Variable}`: Variables from the affine coordinate ring.
- `jac::Matrix{Expression}`: Jacobian matrix corresponding to `equations` and `variables`.
- `point::Vector{Number}`: The initial point that Newton's method is applied to.
- `tol::Number` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Number` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.

# Returns
- `q::Vector{Number}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`
"""
function symmetric_newton_correct(equations, variables, jacobian, p; tol = 1e-13, time_penalty=2)::Vector{Number}
    global _q = Base.copy(p)
    global qnew = _q
    global damping = 0.15
    J = Matrix{Float64}(evaluate.(jacobian, variables=>_q))
    new_equations, J_new = Base.copy(equations), Base.copy(J)
    index = 0
    # Randomize the linear system of equations
    while length(equations)>0 && norm(evaluate.(equations, variables=>_q), Inf) > tol
        #println("sym: $damping\t $(norm(evaluate(eqns, vars=>q)))")
        if index%6 == 0
            stress_dimension = size(nullspace(J'; atol=1e-8))[2]
            if stress_dimension > 0
                rand_mat = randn(Float64, length(equations) - stress_dimension, length(equations))
                new_equations = rand_mat*equations
                J_new = rand_mat*J
            else
                new_equations = equations
                J_new = J
            end
        end
        #println("$damping\t $(norm(evaluate(equations, variables=>_q)))")
        #qnew = _q - damping * (J_new \ evaluate(new_equations, variables=>_q))
        r_val = evaluate(new_equations, variables=>_q)
        v = -(J_new \ r_val)
        global damping = 1
        qnew = _q + damping*v
        while norm(evaluate(equations, variables=>qnew)) > norm(evaluate(equations, variables=>_q)) - 0.25 * damping * v'*v
            global damping = damping*0.6
            qnew = _q + damping*v
            if damping < 1e-10# || Base.time()-start_time > length(point)/time_penalty
                throw("Newton's method did not converge in time.")
            end
        end

        _q = qnew
        index=index+1
    end
    return _q
end


"""
    is_rigid(F[; tol, newton_tol, tested_random_flexes, symmetric_newton])

Heuristically checks if a geometric constraint system `F` is (continuously) rigid. 
"""
function is_rigid(F::AllTypes; tol::Number=1e-5, newton_tol::Number=1e-13, tested_random_flexes::Int=4, symmetric_newton::Bool=false)::Bool
    if is_inf_rigid(F; tol=tol)
        return true
    end
    for _ in 1:tested_random_flexes
        D = DeformationPath(F, [], 5; step_size=sqrt(tol), newton_tol=newton_tol, random_flex=true, symmetric_newton=symmetric_newton)
        if any(sample->norm(sample-D.motion_samples[1], Inf)>tol, D.motion_samples)
            return false
        end
    end
    return true
end

"""
    is_inf_rigid(F[; tol])

Checks if a geometric constraint system `F` is infinitesimally rigid.
"""
function is_inf_rigid(F::AllTypes; tol::Number=1e-8)::Bool
    if typeof(F)==Framework
        K_n = Framework([[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices).G
    elseif typeof(F)==AngularFramework
        K_n = AngularFramework([[i,j,k] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) for k in 1:length(F.G.vertices) if (i<j && j<k) || (i<k && k<j) || (j<i && i<k)], F.G.realization; pinned_vertices=F.G.pinned_vertices).G
    elseif typeof(F)==FrameworkOnSurface
        K_n = deepcopy(G)
        add_equations!(K_n, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j]])
    elseif typeof(F)==VolumeHypergraph
        K_n = VolumeHypergraph(collect(powerset(F.G.vertices, F.G.dimension+1, F.G.dimension+1)), F.G.realization).G
    elseif typeof(F)==Polytope || typeof(F)==SpherePacking || typeof(F)==BodyHinge
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    elseif typeof(F)==SphericalDiskPacking
        minkowski_scalar_product(e1,e2) = e1'*e2-1
        inversive_distances = [minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[1]]) * minkowski_scalar_product(F.G.realization[:,contact[2]], F.G.realization[:,contact[2]])) for contact in powerset(F.G.vertices, 2, 2)]
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, [minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[1]]) * minkowski_scalar_product(F.G.xs[:,contact[2]], F.G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(F.G.vertices, 2, 2))], F.G.realization, F.G.xs)
    else
        throw("Type of F is not yet supported. It is $(typeof(F)).")
    end
    inf_flexes = nullspace(evaluate(F.G.jacobian, F.G.variables=>to_Array(F, F.G.realization)); atol=tol)
    trivial_inf_flexes = nullspace(evaluate(typeof(K_n)==ConstraintSystem ? K_n.jacobian : K_n.G.jacobian, (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables)=>to_Array(F, F.G.realization)[1:length( (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables))]); atol=tol)
    #println("flexes: $(size(inf_flexes)[2]), nontrivial: $(size(inf_flexes)[2]-size(trivial_inf_flexes)[2])")
    return length(inf_flexes) == length(trivial_inf_flexes)
end

"""
    is_second_order_rigid(F[; tol, newton_tol])

Checks if a geometric constraint system `F` is second-order rigid.

See also [`compute_nonblocked_flex`](@ref)
"""
function is_second_order_rigid(F::AllTypes; kwargs...)::Bool
    flex_mult = compute_nonblocked_flex(F; kwargs...)
    if length(flex_mult)==0
        return true
    else
        return false
    end
end

"""
    compute_nonblocked_flex(F[; tol, newton_tol])

Compute an infinitesimal flex of `F` that is not blocked by an equilibrium stress.
"""
function compute_nonblocked_flex(F::AllTypes; tol::Number=1e-6, newton_tol::Number=1e-13)::Vector
    if typeof(F)==Framework
        K_n = Framework([[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices)
    elseif typeof(F)==Polytope || typeof(F)==SpherePacking || typeof(F)==BodyHinge
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    else
        throw("Type of F is not yet supported. It is $(typeof(F)).")
    end
    flexes = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), K_n; tol=1e-8)
    rigidity_matrix = evaluate.(F.G.jacobian, F.G.variables=>to_Array(F, F.G.realization))
    stresses = nullspace(rigidity_matrix'; atol=1e-8)

    @var λ[1:size(flexes)[2]] ω[1:size(stresses)[2]]
    parametrized_flex = flexes*λ
    parametrized_stress = stresses*ω
    stress_energy = parametrized_stress'*evaluate.(F.G.jacobian, F.G.variables=>Vector{Expression}(parametrized_flex))*parametrized_flex
    stress_poly_system = differentiate(stress_energy, ω)
    projective_stress_system = vcat(stress_poly_system, sum(λ .^ 2) - 1)
    J_stress_energy = Matrix{Expression}(hcat([differentiate(eq, λ) for eq in projective_stress_system]...)')
    for index in 1:size(flexes)[2]*size(stresses)[2]*2
        rand_flex_parameter = randn(Float64, size(flexes)[2])
        rand_flex_parameter = rand_flex_parameter ./ norm(rand_flex_parameter)
        try
            q = newton_correct(projective_stress_system, λ, J_stress_energy, rand_flex_parameter; tol = newton_tol, time_penalty=10)
            return q
        catch e
            continue
        end
    end
    return []
end


"""
    animate(F, filename[; flex_mult, random_flex, num_steps, step_size])

Create an animation of a geometric constraint system `F`.
"""
function animate(F::AllTypes, filename::String; flex_mult::Vector=[], random_flex::Bool=true, num_steps::Int=100, step_size::Float64=1e-2, kwargs...)
    D = DeformationPath(F, flex_mult, num_steps; step_size=step_size, random_flex=random_flex)
    animate(D, F, filename; kwargs...)
end

"""
    animate(D, F, filename)

Create an animation of a geometric constraint system `F` given a previously computed deformation path.

This is a wrapper method for animations for the individual geometric constraint systems. 
It calls one of the following: [`animate2D_framework`](@ref), [`animate3D_framework`](@ref), [`animate3D_frameworkonsurface`](@ref), [`animate2D_hypergraph`](@ref), [`animate3D_polytope`](@ref), [`animate2D_diskpacking`](@ref), [`animate3D_spherepacking`](@ref) and [`animate3D_sphericaldiskpacking`](@ref).
"""
function animate(D::DeformationPath, F, filename::String; kwargs...)
    if typeof(F)==Framework || typeof(F)==AngularFramework
        if F.G.dimension==2
            animate2D_framework(D, F, filename; kwargs...)
        elseif F.G.dimension==3
            animate3D_framework(D, F, filename; kwargs...)
        else
            throw("The dimension of 'F' needs to be either 2 or 3, but is $(F.G.dimension)")
        end
    elseif typeof(F)==FrameworkOnSurface
        animate3D_frameworkonsurface(D, F, filename; kwargs...)
    elseif typeof(F)==VolumeHypergraph
        animate2D_hypergraph(D, F, filename; kwargs...)
    elseif typeof(F)==Polytope || typeof(F)==BodyHinge
        animate3D_polytope(D, F, filename; kwargs...)
    elseif typeof(F)==SpherePacking
        if F.G.dimension==2
            animate2D_diskpacking(D, F, filename; kwargs...)
        elseif F.G.dimension==3
            animate3D_spherepacking(D, F, filename; kwargs...)
        else
            throw("The dimension of 'F' needs to be either 2 or 3, but is $(F.G.dimension)")
        end
    elseif typeof(F)==SphericalDiskPacking
        animate3D_sphericaldiskpacking(D, F, filename; kwargs...)
    else
        throw("The type of 'F' needs to be either Framework, SpherePacking, Polytope, SphericalDiskPacking or VolumeHypergraph, but is $(typeof(F))")
    end
end


"""
    animate2D_framework(D, F, filename)

Compute an animation for a 2-dimensional bar-joint framework.
"""
function animate2D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::String; recompute_deformation_samples::Bool=true, fixed_vertices::Tuple{Int,Int}=(1,2), fixed_direction::Vector{Number}=[1.,0], framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.15, markercolor=:red3, pin_point_offset=0.1, vertex_size::Union{Float64,Int}=55, line_width::Union{Float64,Int}=12, angle_color=:lightgrey, font_color=:lightgrey, angle_size=0.3, edge_color=:steelblue, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    length(fixed_vertices)==2 && fixed_vertices[1] in D.G.vertices && fixed_vertices[2] in D.G.vertices || throw("fixed_vertices is not a vertex of the underlying graph.")
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_vertices[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
    end
    
    if isapprox(norm(fixed_direction),0;atol=1e-6)
        @warn "fixed_direction is $(norm(fixed_direction)) which is too close to 0! We thus set it to [1,0]"
        fixed_direction = [1.,0]
    end
    fixed_direction = fixed_direction ./ norm(fixed_direction)
    for i in 1:length(matrix_coords)
        theta = atan(matrix_coords[i][:,fixed_vertices[2]][2] , matrix_coords[i][:,fixed_vertices[2]][1])
        base_theta = atan(fixed_direction[2], fixed_direction[1])
        theta = theta-base_theta
        rotation_matrix = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        # Rotate the realization to the `fixed_direction`.
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = rotation_matrix*matrix_coords[i][:,j]
        end
    end

    if recompute_deformation_samples
        D.motion_samples = [to_Array(F, matrix_coords[i]) for i in 1:length(matrix_coords)]
    end

    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    limits= [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])]
    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    allVertices=@lift begin
        pointys = matrix_coords[$time]
        [Point2f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end

    if typeof(F)==AngularFramework
        L = angle_size*maximum([minimum([norm(matrix_coords[1][:,bar[1]]-matrix_coords[1][:,bar[2]]) for bar in F.bars]),0])
        angle_points=@lift begin
            output=[]
            for angle in F.angles
                poly_points=[($allVertices)[angle[2]]]
                if norm(($allVertices)[angle[2]]-($allVertices)[angle[1]])==0
                    foreach(_->push!(poly_points, Point2f(($allVertices)[angle[2]]+L*(($allVertices)[angle[2]]-($allVertices)[angle[3]]))/norm(($allVertices)[angle[2]]-($allVertices)[angle[3]])), 1:2)
                elseif norm(($allVertices)[angle[2]]-($allVertices)[angle[3]])==0
                    foreach(_->push!(poly_points, Point2f(($allVertices)[angle[2]]+L*(($allVertices)[angle[2]]-($allVertices)[angle[1]]))/norm(($allVertices)[angle[2]]-($allVertices)[angle[1]])), 1:2)
                else    
                    for t in 0:0.05:1
                        curpt = ($allVertices)[angle[2]] + L*((1-t)*(($allVertices)[angle[1]]-($allVertices)[angle[2]]) + t*(($allVertices)[angle[3]]-($allVertices)[angle[2]])) ./ norm((1-t)*(($allVertices)[angle[1]]-($allVertices)[angle[2]]) + t*(($allVertices)[angle[3]]-($allVertices)[angle[2]]))
                        push!(poly_points, Point2f(curpt))
                    end
                end
                push!(output, poly_points)
            end
            output
        end

        foreach(i->poly!(ax,@lift(($angle_points)[i]), color=(angle_color,0.5), strokewidth=0), 1:length(F.angles))
        foreach(i->lines!(ax,@lift(($angle_points)[i][2:end]), color=angle_color, linewidth=line_width/2), 1:length(F.angles))    
    end

    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=edge_color), F.bars)
    foreach(v->scatter!(ax, @lift([Point2f(($allVertices)[v]-[pin_point_offset,0])]); markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=25, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))

    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
    end
end


"""
    animate3D_framework(D, F, filename)

Compute an animation for a 3-dimensional bar-joint framework.
"""
function animate3D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::String; recompute_deformation_samples::Bool=true, fixed_vertices::Union{Tuple{Int,Int}, Tuple{Int,Int,Int}}=(1,2), fixed_direction=[1.,0,0], framerate::Int=25, animate_rotation=false, azimuth = π / 4, elevation=pi/8, perspectiveness=0., rotation_frames = 240, markercolor=:red3, pin_point_offset=0.05, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_size::Union{Float64,Int}=55, vertex_labels=false, font_color=:lightgrey, line_width::Union{Float64,Int}=12, angle_color=:lightgrey, angle_size=0.3, edge_color=:steelblue, vertex_color=:black, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect = (1, 1, 1), perspectiveness=perspectiveness)
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    length(fixed_vertices)==length(collect(Set(fixed_vertices))) && fixed_vertices[1] in D.G.vertices && fixed_vertices[2] in D.G.vertices && (length(fixed_vertices)==2 || fixed_vertices[3] in D.G.vertices) || throw("The elements of `fixed_vertices`` are not vertices of the underlying graph.")
    
    if isapprox(norm(fixed_direction),0;atol=1e-6)
        @warn "fixed_direction is $(norm(fixed_direction)) which is too close to 0! We thus set it to [1,0,0]"
        fixed_direction = [1.,0,0]
    end
    fixed_direction = fixed_direction ./ norm(fixed_direction)

    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_vertices[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
        edge_vector = Vector(matrix_coords[i][:,fixed_vertices[2]] ./ norm(matrix_coords[i][:,fixed_vertices[2]]))
        rotation_axis = cross(fixed_direction, edge_vector)
        if isapprox(norm(rotation_axis), 0, atol=1e-6)
            if fixed_direction'*edge_vector<0
                rotation_matrix = [-1 0 0; 0 -1 0; 0 0 1;]
            else
                rotation_matrix = [1 0 0; 0 1 0; 0 0 1;]
            end
        else
            rotation_axis = rotation_axis ./ norm(rotation_axis)
            angle = acos(fixed_direction'* edge_vector)
            rotation_matrix = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                                rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                                rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
        end
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = inv(rotation_matrix)*matrix_coords[i][:,j]
        end

        if length(fixed_vertices)==3
            edge_vector_new = Vector(matrix_coords[i][:,fixed_vertices[3]] ./ norm(matrix_coords[i][:,fixed_vertices[3]]))
            target_vector = [0,edge_vector_new[2],edge_vector_new[3]]
            target_vector = target_vector ./ norm(target_vector)
            if isapprox(edge_vector_new[3],0; atol=1e-10)
                angle = 0
            else
                angle = acos(target_vector'* [0,1,0])
            end
            rotation_matrix_new = [ cos(angle)+fixed_direction[1]^2*(1-cos(angle)) fixed_direction[1]*fixed_direction[2]*(1-cos(angle))-fixed_direction[3]*sin(angle) fixed_direction[1]*fixed_direction[3]*(1-cos(angle))+fixed_direction[2]*sin(angle); 
                                fixed_direction[1]*fixed_direction[2]*(1-cos(angle))+fixed_direction[3]*sin(angle) cos(angle)+fixed_direction[2]^2*(1-cos(angle)) fixed_direction[2]*fixed_direction[3]*(1-cos(angle))-fixed_direction[1]*sin(angle); 
                                fixed_direction[1]*fixed_direction[3]*(1-cos(angle))-fixed_direction[2]*sin(angle) fixed_direction[2]*fixed_direction[3]*(1-cos(angle))+fixed_direction[1]*sin(angle) cos(angle)+fixed_direction[3]^2*(1-cos(angle));]
            if edge_vector_new[3]<0
                rotation_matrix_new = inv(rotation_matrix_new)
            end                
            for j in 1:size(matrix_coords[i])[2]
                matrix_coords[i][:,j] = inv(rotation_matrix_new)*matrix_coords[i][:,j]
            end
        end
    end

    if recompute_deformation_samples
        D.motion_samples = [to_Array(F, matrix_coords[i]) for i in 1:length(matrix_coords)]
    end

    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    zlims = [minimum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...))]
    limits = [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]
    xlims!(ax, limits[1]-padding, limits[2]+padding)
    ylims!(ax, limits[1]-padding, limits[2]+padding)
    zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    allVertices=@lift begin
        pointys = matrix_coords[$time]
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end

    if typeof(F)==AngularFramework
        L = angle_size*maximum([minimum([norm(matrix_coords[1][:,bar[1]]-matrix_coords[1][:,bar[2]]) for bar in F.bars]),0])
        angle_points=@lift begin
            output=[]
            for angle in F.angles
                poly_points=[($allVertices)[angle[2]]]
                if norm(($allVertices)[angle[2]]-($allVertices)[angle[1]])==0
                    foreach(_->push!(poly_points, Point3f(($allVertices)[angle[2]]+L*(($allVertices)[angle[2]]-($allVertices)[angle[3]]))/norm(($allVertices)[angle[2]]-($allVertices)[angle[3]])), 1:2)
                elseif norm(($allVertices)[angle[2]]-($allVertices)[angle[3]])==0
                    foreach(_->push!(poly_points, Point3f(($allVertices)[angle[2]]+L*(($allVertices)[angle[2]]-($allVertices)[angle[1]]))/norm(($allVertices)[angle[2]]-($allVertices)[angle[1]])), 1:2)
                else    
                    for t in 0:0.05:1
                        curpt = ($allVertices)[angle[2]] + L*((1-t)*(($allVertices)[angle[1]]-($allVertices)[angle[2]]) + t*(($allVertices)[angle[3]]-($allVertices)[angle[2]])) ./ norm((1-t)*(($allVertices)[angle[1]]-($allVertices)[angle[2]]) + t*(($allVertices)[angle[3]]-($allVertices)[angle[2]]))
                        push!(poly_points, Point3f(curpt))
                    end
                end
                push!(output, poly_points)
            end
            output
        end

        foreach(i->poly!(ax,@lift(($angle_points)[i]), color=(angle_color,0.5), strokewidth=0), 1:length(F.angles))
        foreach(i->lines!(ax,@lift(($angle_points)[i][2:end]), color=angle_color, linewidth=line_width/2), 1:length(F.angles))    
    end

    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=:steelblue), F.bars)
    foreach(v->scatter!(ax, @lift([Point3f(($allVertices)[v]-[pin_point_offset,0,0])]); markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(D.G.vertices))
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=25, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))

    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end

    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.elevation[] = elevation
            ax.azimuth[] = azimuth + 2pi * t / rotation_frames
        end
    end
end


"""
    animate3D_frameworkonsurface(D, F, filename)

Compute an animation for a 3-dimensional bar-joint framework constrained to a surface.
"""
function animate3D_frameworkonsurface(D::DeformationPath, F::FrameworkOnSurface, filename::String; alpha=0.45, framerate::Int=25, animate_rotation=false, azimuth = pi/4, elevation=pi/8, perspectiveness=0., rotation_frames = 480, markercolor=:red3, pin_point_offset=0.05, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_size::Union{Float64,Int}=55, line_width::Union{Float64,Int}=10, edge_color=:steelblue, vertex_labels=true, font_color=:lightgrey, vertex_color=:black, filetype::String="gif", surface_color=:grey80, surface_samples=150)
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect = (1, 1, 1), perspectiveness=perspectiveness)
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]

    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    zlims = [minimum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...))]
    limits = [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]
    xlims!(ax, limits[1]-padding, limits[2]+padding)
    ylims!(ax, limits[1]-padding, limits[2]+padding)
    zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    x = collect(Float64, range(limits[1]-padding, step=(limits[2]-limits[1]+2*padding)/surface_samples, length = surface_samples+1))
    y = collect(Float64, range(limits[1]-padding, step=(limits[2]-limits[1]+2*padding)/surface_samples, length = surface_samples+1))
    z = collect(Float64, range(limits[1]-padding, step=(limits[2]-limits[1]+2*padding)/surface_samples, length = surface_samples+1))
    A = [F.surface([xi,yi,zi]) for xi in x, yi in y, zi in z]
    mc_ranged = MC(A, Int; x, y, z)
    march(mc_ranged, 0.)
    msh = makemesh(GeometryBasics, mc_ranged)
    mesh!(ax, msh, color=(surface_color,alpha), transparency=true)

    time=Observable(1)
    allVertices=@lift begin
        pointys = matrix_coords[$time]
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=:steelblue), F.bars)
    foreach(v->scatter!(ax, @lift([Point3f(($allVertices)[v]-[pin_point_offset,0,0])]); markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(D.G.vertices))
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=25, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end

    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.elevation[] = elevation
            ax.azimuth[] = azimuth + 2pi * t / rotation_frames
        end
    end
end


"""
    animate2D_hypergraph(D, F, filename)

Compute an animation for a 2-dimensional volume hypergraph.
"""
function animate2D_hypergraph(D::DeformationPath, F::VolumeHypergraph, filename::String; alpha=0.2, recompute_deformation_samples::Bool=true, target_stretch::Union{Float64,Int}=1., fixed_triangle::Union{Tuple{Int,Int,Int},Vector{Int},Nothing}=nothing, font_color=:black, skip_stretch::Bool=true, tip_value::Union{Float64,Int}=0.5, framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_size::Union{Float64,Int}=42, line_width::Union{Float64,Int}=6, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    if facet_colors==nothing
        facet_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(F.volumes), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
    end
    skip_scaling = false
    if fixed_triangle==nothing
        fixed_triangle=(1,2,3)
        skip_scaling=true
    else
        all(i->fixed_triangle[i] in D.G.vertices, 1:3) && (Tuple(fixed_triangle) in [Tuple(facet) for facet in F.volumes]) || (Tuple([fixed_triangle[2],fixed_triangle[3],fixed_triangle[1]]) in [Tuple(facet) for facet in F.volumes]) || (Tuple([fixed_triangle[3],fixed_triangle[1],fixed_triangle[2]]) in [Tuple(facet) for facet in F.volumes]) || throw("fixed_triangle is not a vertex of the underlying graph.")
    end
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_triangle[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
    end
    fixed_direction = [1.,0]
    for i in 1:length(matrix_coords)
        theta = atan(matrix_coords[i][:,fixed_triangle[2]][2], matrix_coords[i][:,fixed_triangle[2]][1])
        rotation_matrix = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        # Rotate the realization to the `fixed_direction`.
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = rotation_matrix*matrix_coords[i][:,j]
        end
        # Scale to (1,0)
        if skip_scaling || skip_stretch
            continue
        end
        x_val = matrix_coords[i][1,fixed_triangle[2]]
        scaling_matrix = [target_stretch/x_val 0; 0 x_val/target_stretch]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = scaling_matrix*matrix_coords[i][:,j]
        end
        y_val = matrix_coords[i][:,fixed_triangle[3]]
        shear_matrix = [1 (tip_value-y_val[1])/y_val[2]; 0 1]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = shear_matrix*matrix_coords[i][:,j]
        end
    end

    if recompute_deformation_samples
        D.motion_samples = [to_Array(F, matrix_coords[i]) for i in 1:length(matrix_coords)]
    end

    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    limits= [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])]
    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    allVertices=@lift begin
        pointys = matrix_coords[$time]
        [Point2f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(i->poly!(ax, @lift([($allVertices)[Int64(v)] for v in F.volumes[i]]); color=(facet_colors[i], alpha)), 1:length(F.volumes))
    foreach(i->lines!(ax, @lift([($allVertices)[Int64(v)] for v in vcat(F.volumes[i], F.volumes[i][1])]); linewidth=line_width, linestyle=:dash, color=facet_colors[i]), 1:length(F.volumes))
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(F.G.vertices))
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=26, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))

    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end

    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
    end
end


"""
    animate3D_polytope(D, F, filename)

Compute an animation for a 3-dimensional polytope.
"""
function animate3D_polytope(D::DeformationPath, F::Union{Polytope,BodyHinge}, filename::String; renderEntirePolytope::Bool=true, recompute_deformation_samples::Bool=true, fixed_vertices::Union{Tuple{Int,Int}, Tuple{Int,Int,Int}}=(1,2), alpha=0.6, font_color=:lightgrey, facet_color=:grey98, framerate::Int=25, animate_rotation=false, azimuth = π / 10, elevation=pi/8, perspectiveness=0., rotation_frames = 240, step::Int=1, padding::Union{Float64,Int}=0.1, vertex_size::Union{Float64,Int}=45, line_width::Union{Float64,Int}=8.5, edge_color=:steelblue, special_edge=nothing, special_edge_color=:red3, vertex_color=:black, vertex_labels::Bool=false, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    fixed_vertices[1] in 1:(size(F.G.realization)[2]-length(F.facets)) && fixed_vertices[2] in 1:(size(F.G.realization)[2]-length(F.facets)) && (length(fixed_vertices)==2 || fixed_vertices[3] in 1:(size(F.G.realization)[2]-length(F.facets))) || throw("The elements of `fixed_vertices`` are not vertices of the underlying graph.")
    ax = Axis3(fig[1,1], aspect = (1, 1, 1), perspectiveness=perspectiveness)

    isnothing(special_edge) || (special_edge in [[edge[1],edge[2]] for edge in F.edges] || [special_edge[2], special_edge[1]] in [[edge[1],edge[2]] for edge in F.edges]) || throw(error("The `special_edge` needs to be an edge of the polytope's 1-skeleton!"))

    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_vertices[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
        edge_vector = Vector(matrix_coords[i][:,fixed_vertices[2]] ./ norm(matrix_coords[i][:,fixed_vertices[2]]))
        rotation_axis = cross([1,0,0], edge_vector)
        if isapprox(norm(rotation_axis), 0, atol=1e-4)
            if [1,0,0]'*edge_vector<0
                rotation_matrix = [-1 0 0; 0 -1 0; 0 0 1;]
            else
                rotation_matrix = [1 0 0; 0 1 0; 0 0 1;]
            end
        else
            rotation_axis = rotation_axis ./ norm(rotation_axis)
            angle = acos([1,0,0]' * edge_vector)
            rotation_matrix = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                                rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                                rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
        end

        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = inv(rotation_matrix)*matrix_coords[i][:,j]
        end

        if length(fixed_vertices)==3
            edge_vector_new = Vector(matrix_coords[i][:,fixed_vertices[3]] ./ norm(matrix_coords[i][:,fixed_vertices[3]]))
            target_vector = [0,edge_vector_new[2],edge_vector_new[3]]
            target_vector = target_vector ./ norm(target_vector)
            if isapprox(edge_vector_new[3],0; atol=1e-10)
                angle = 0
            else
                angle = acos(target_vector'* [0,1,0])
            end
            rotation_matrix_new = [ cos(angle)+fixed_direction[1]^2*(1-cos(angle)) fixed_direction[1]*fixed_direction[2]*(1-cos(angle))-fixed_direction[3]*sin(angle) fixed_direction[1]*fixed_direction[3]*(1-cos(angle))+fixed_direction[2]*sin(angle); 
                                fixed_direction[1]*fixed_direction[2]*(1-cos(angle))+fixed_direction[3]*sin(angle) cos(angle)+fixed_direction[2]^2*(1-cos(angle)) fixed_direction[2]*fixed_direction[3]*(1-cos(angle))-fixed_direction[1]*sin(angle); 
                                fixed_direction[1]*fixed_direction[3]*(1-cos(angle))-fixed_direction[2]*sin(angle) fixed_direction[2]*fixed_direction[3]*(1-cos(angle))+fixed_direction[1]*sin(angle) cos(angle)+fixed_direction[3]^2*(1-cos(angle));]
            if edge_vector_new[3]<0
                rotation_matrix_new = inv(rotation_matrix_new)
            end                
            for j in 1:size(matrix_coords[i])[2]
                matrix_coords[i][:,j] = inv(rotation_matrix_new)*matrix_coords[i][:,j]
            end
        end
    end

    if recompute_deformation_samples
        D.motion_samples = [to_Array(F, matrix_coords[i]) for i in 1:length(matrix_coords)]
    end

    centroid = sum([matrix_coords[1][:,i] for i in 1:(size(F.G.realization)[2]-length(F.facets))]) ./ (size(F.G.realization)[2]-length(F.facets))
    for i in 1:length(matrix_coords)
        for j in 1:(size(F.G.realization)[2]-length(F.facets))
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - centroid
        end
    end

    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    zlims = [minimum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...))]
    limits = [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]
    xlims!(ax, limits[1]-padding, limits[2]+padding)
    ylims!(ax, limits[1]-padding, limits[2]+padding)
    zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)
    time=Observable(1)

    allVertices=@lift begin
        pointys = matrix_coords[$time][:,1:(size(F.G.realization)[2]-length(F.facets))]
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end

    allVertices_asLists = @lift begin
        pointys = matrix_coords[$time][:,1:(size(F.G.realization)[2]-length(F.facets))]
        [pointys[:,j] for j in 1:size(pointys)[2]]
    end

    if isnothing(special_edge)
        foreach(i->linesegments!(ax, @lift([($allVertices)[Int64(F.edges[i][1])], ($allVertices)[Int64(F.edges[i][2])]]); linewidth=line_width, color=edge_color), 1:length(F.edges))
    else
        edges_here = filter(edge->!([special_edge[1],special_edge[2]]==[edge[1],edge[2]] || [special_edge[1],special_edge[2]]==[edge[2],edge[1]]), F.edges)
        foreach(i->linesegments!(ax, @lift([($allVertices)[Int64(edges_here[i][1])], ($allVertices)[Int64(edges_here[i][2])]]); linewidth=line_width, color=edge_color), 1:length(edges_here))
        linesegments!(ax, @lift([($allVertices)[Int64(special_edge[1])], ($allVertices)[Int64(special_edge[2])]]); linewidth=line_width+2.5, color=special_edge_color)
    end
    vertex_labels && foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=vertex_color), 1:(size(F.G.realization)[2]-length(F.facets)))
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:(size(F.G.realization)[2]-length(F.facets)))
    if typeof(F) <: Polytope && renderEntirePolytope
        mesh!(ax, @lift(Polyhedra.Mesh(Polyhedra.polyhedron(Polyhedra.vrep(($allVertices_asLists)), CDDLib.Library(:exact)))); shading=NoShading, color=(facet_color,alpha), transparency=true)
    elseif typeof(F) <: BodyHinge || !renderEntirePolytope
        for face in F.facets
            try
                mesh!(ax, @lift(Polyhedra.Mesh(Polyhedra.polyhedron(Polyhedra.vrep([($allVertices_asLists)[j] for j in face]), CDDLib.Library(:exact)))); shading=NoShading, color=(facet_color,alpha), transparency=true)
            catch e
                continue
            end
        end
    end
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end

    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.elevation[] = elevation
            ax.azimuth[] = azimuth + 2pi * t / rotation_frames
        end
    end
end


"""
    animate2D_diskpacking(D, F, filename)

Compute an animation for a 2-dimensional sticky disk packing.
"""
function animate2D_diskpacking(D::DeformationPath, F::SpherePacking, filename::String; alpha=0.08, framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_labels=true, disk_strokewidth::Union{Float64,Int}=8.5, line_width::Union{Float64,Int}=7, font_color=:black, sphere_color=:steelblue, markersize::Union{Float64,Int}=75, markercolor=:red3, dualgraph_color=:grey80, n_circle_segments::Int=50, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    if F.G.dimension!=2
        throw("The dimension must be 2, but is $(F.G.dimension)!")
    end
    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    limits= [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])]
    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation-maximum(F.radii), limits[2]+padding+0.5*translation+maximum(F.radii))
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation-maximum(F.radii), limits[2]+padding+0.5*translation+maximum(F.radii))
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    contacts=Observable(D._contacts[1])
    allVertices=@lift begin
        pointys = matrix_coords[$time]
        [Point2f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    linesegments!(ax, @lift(vcat([[($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]] for edge in $contacts]...)); linewidth = line_width, color=dualgraph_color)
    for index in 1:length(F.G.vertices)
        disk_vertices = @lift([Point2f(Vector($allVertices[index])+F.radii[index]*[cos(2*i*pi/n_circle_segments), sin(2*i*pi/n_circle_segments)]) for i in 1:n_circle_segments])
        diskedges = [(i,i%n_circle_segments+1) for i in 1:n_circle_segments]
        poly!(ax, @lift([($disk_vertices)[i] for i in 1:n_circle_segments]); color=(sphere_color, alpha))
        lines!(ax, @lift([($disk_vertices)[v] for v in vcat(1:n_circle_segments,1)]); linewidth = disk_strokewidth, color=sphere_color)
    end

    foreach(v->scatter!(ax, @lift([($allVertices)[v]]); markersize=markersize, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end

    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        contacts[] = D._contacts[t]
    end
end


"""
    animate3D_spherepacking(D, F, filename)

Compute an animation for a 3-dimensional sticky sphere packing.
"""
function animate3D_spherepacking(D::DeformationPath, F::SpherePacking, filename::String; alpha=0.2, framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.1, vertex_labels=true, font_color=:black, line_width::Union{Float64,Int}=7, sphere_color=:steelblue, markersize::Union{Float64,Int}=55, markercolor=:red3, dualgraph_color=:grey50, n_circle_segments::Int=50, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect = (1, 1, 1))
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    if F.G.dimension!=3
        throw("The dimension must be 3, but is $(F.G.dimension)!")
    end
    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    zlims = [minimum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...))]
    limits= [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]
    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation-maximum(F.radii), limits[2]+padding+0.5*translation+maximum(F.radii))
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation-maximum(F.radii), limits[2]+padding+0.5*translation+maximum(F.radii))
    translation = (zlims[1]-limits[1]) - (limits[2]-zlims[2])
    zlims!(ax, limits[1]-padding+0.5*translation-maximum(F.radii), limits[2]+padding+0.5*translation+maximum(F.radii))
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    contacts=Observable(D._contacts[1])
    allVertices=@lift begin
        pointys = matrix_coords[$time]
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end

    for index in 1:length(F.G.vertices)
        mesh!(ax, @lift(Sphere(($allVertices)[index], F.radii[index]));  transparency=true, color = (sphere_color,alpha))
    end
    linesegments!(ax, @lift(vcat([[($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]] for edge in $contacts]...)); linewidth = line_width, color=dualgraph_color)
    
    foreach(v->scatter!(ax, @lift([($allVertices)[v]]); markersize=markersize, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end

    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        contacts[] = D._contacts[t]
    end
end


"""
    animate3D_sphericaldiskpacking(D, F, filename)

Compute an animation for a disk packing on the 2-sphere in the Minkowski metric.
"""
function animate3D_sphericaldiskpacking(D::DeformationPath, F::SphericalDiskPacking, filename::String; alpha=0.15, framerate::Int=25, animate_rotation=false, azimuth = π / 10, elevation=pi/8, perspectiveness=0., font_color=:black, rotation_frames = 240, step::Int=1, padding=0.015, sphere_color=:lightgrey, vertex_size=60, disk_strokewidth=9, line_width=6, disk_color=:steelblue, dualgraph_color=(:red3,0.45), vertex_color=:black, vertex_labels::Bool=true, n_circle_segments=45, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]

    ax = Axis3(fig[1,1], aspect=(1,1,1), perspectiveness=perspectiveness)
    xlims!(ax,-1.5-padding, 1.5+padding)
    ylims!(ax,-1.5-padding, 1.5+padding)
    zlims!(ax,-1.5-padding, 1.5+padding)
    hidespines!(ax)
    hidedecorations!(ax)
    mesh!(ax, Sphere(Point3f(0), 1f0); transparency=true, color = (sphere_color,alpha))

    time=Observable(1)
    planePoints=@lift begin
        pointys = matrix_coords[$time]
        [Point3f(pointys[:,j]./norm(pointys[:,j])^2) for j in 1:size(pointys)[2]]
    end

    spherePoints=@lift begin
        pointys = matrix_coords[$time]
        [Point3f(pointys[:,j]./norm(pointys[:,j])) for j in 1:size(pointys)[2]]
    end

    koebePoints=@lift begin
        pointys = matrix_coords[$time]
        pointys = [Point3f(-pointys[:,j]) for j in 1:size(pointys)[2]]
        output = []
        for i in 1:length(F.G.vertices)
            rotation_axis = cross([0, 0, 1], Vector(($spherePoints)[i]))
            if isapprox(norm(rotation_axis), 0, atol=1e-6)
                angle = acos([0, 0, 1]'* ($spherePoints)[i])
                rotation_matrix = [1 0 0; 0 1 0; 0 0 cos(angle);]
            else
                rotation_axis = rotation_axis ./ norm(rotation_axis)
                angle = acos([0, 0, 1]'* Vector(($spherePoints)[i]))
                rotation_matrix = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                                    rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                                    rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
            end
            push!(output,Point3f(inv(rotation_matrix)*[0,0,norm(pointys[i])]))
        end
        output
    end
    foreach(edge->linesegments!(ax, @lift([($koebePoints)[Int64(edge[1])], ($koebePoints)[Int64(edge[2])]]); linewidth = line_width, color=dualgraph_color), F.contacts)
    disk_vertices = @lift begin
        output = []
        for i in 1:length(F.G.vertices)
            rotation_axis = cross([0, 0, 1], Vector(($spherePoints)[i]))
            if isapprox(norm(rotation_axis), 0, atol=1e-6)
                angle = acos([0, 0, 1]'* ($spherePoints)[i])
                rotation_matrix = [1 0 0; 0 1 0; 0 0 cos(angle);]
            else
                rotation_axis = rotation_axis ./ norm(rotation_axis)
                angle = acos([0, 0, 1]'* Vector(($spherePoints)[i]))
                rotation_matrix = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                                    rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                                    rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
            end
            radius = sqrt(1-norm(($planePoints)[i])^2)
            push!(output,[Point3f(inv(rotation_matrix)*([radius*cos(2*j*pi/n_circle_segments), radius*sin(2*j*pi/n_circle_segments), norm(($planePoints)[i])])) for j in 1:n_circle_segments])
        end
        output
    end
    foreach(i->lines!(ax, @lift([($disk_vertices)[i][v] for v in vcat(1:n_circle_segments,1)]); linewidth = disk_strokewidth, color=disk_color), 1:length(F.G.vertices))
    rotatedPoints = @lift begin
        output=[]
        for i in 1:length(F.G.vertices)
            rotation_axis = cross([0, 0, 1], Vector(($spherePoints)[i]))
            if isapprox(norm(rotation_axis), 0, atol=1e-6)
                angle = acos([0, 0, 1]'* ($spherePoints)[i])
                rotation_matrix = [1 0 0; 0 1 0; 0 0 cos(angle);]
            else
                rotation_axis = rotation_axis ./ norm(rotation_axis)
                angle = acos([0, 0, 1]'* Vector(($spherePoints)[i]))
                rotation_matrix = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                                    rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                                    rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
            end
            push!(output, Point3f(inv(rotation_matrix)*[0,0,1]))
        end
        output
    end
    vertex_labels && foreach(i->text!(ax, @lift([($rotatedPoints)[i]]), text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw("The chosen filetype needs to be either gif or mp4, but is $(filetype)")
    end
    
    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.elevation[] = elevation
            ax.azimuth[] = azimuth + 2pi * t / rotation_frames
        end
    end
    return fig
end


"""
    project_deformation_random(D, F, filename)

Compute a random projection of deformation paths.

This method can either take a single deformation path or a vector of deformation paths and projects it to curves in 2D or 3D.
This makes it possible to visualize high-dimensional deformation spaces. 
"""
function project_deformation_random(D::Union{DeformationPath,Vector{DeformationPath}}, projected_dimension::Int; line_width::Number=8, edge_colors=[:green3], markersize::Number=45, markercolor=:steelblue, draw_start::Bool=true)
    if !(projected_dimension in [2,3])
        throw("The projected_dimension is neither 2 nor 3.")
    end

    if D isa DeformationPath
        D = [D]
    else
        length(D) > 0 || throw("The length of the vector 'DeformationPath' is 0.")
    end

    if length(edge_colors) < length(D)
        @warn "The length of `line_colors` is $(length(edge_colors)) but needs to be at least $(length(D)). Choosing distinguishable colors instead."
        edge_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(D), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
    end
    randmats = [hcat([rand(Float64,projected_dimension) for _ in 1:length(D[i].G.variables)]...) for i in 1:length(D)]
    proj_curve = [[(pinv(randmats[i]'*randmats[i])*randmats[i]')'*entry for entry in Defo.motion_samples] for (i,Defo) in enumerate(D)]
    fig = Figure(size=(1000,1000))
    if projected_dimension==3
        ax = Axis3(fig[1,1], aspect = (1, 1, 1))
    else
        ax = Axis(fig[1,1])
    end
    hidespines!(ax)
    hidedecorations!(ax)
    if projected_dimension==3
        foreach(j->lines!(ax, [Point3f(pt) for pt in proj_curve[j]]; linewidth=line_width, color=edge_colors[j]), 1:length(proj_curve))
        draw_start && scatter!(ax, [proj_curve[1][1][1]], [proj_curve[1][1][2]], [proj_curve[1][1][3]]; markersize=markersize, color=markercolor, marker=:pentagon)
    else
        foreach(j->lines!(ax, [Point2f(pt) for pt in proj_curve[j]]; linewidth=line_width, color=edge_colors[j]), 1:length(proj_curve))
        draw_start && scatter!(ax, [proj_curve[1][1][1]], [proj_curve[1][1][2]]; markersize=markersize, color=markercolor, marker=:pentagon)
    end
    return fig
end

end 