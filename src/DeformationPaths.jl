module DeformationPaths

import HomotopyContinuation: evaluate, differentiate, newton, Expression
import LinearAlgebra: norm, pinv, nullspace, rank, qr, zeros, inv, cross, det, svd, I
import GLMakie: Sphere, mesh!, @lift, poly!, text!, Figure, record, hidespines!, hidedecorations!, lines!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Observable, Point3f, Point2f, connect, faces, Mesh, mesh
import ProgressMeter: @showprogress
import Combinatorics: powerset
import Colors: distinguishable_colors, red, green, blue, colormap, RGB

include("GeometricConstraintSystem.jl")
using .GeometricConstraintSystem: ConstraintSystem, Framework, equations!, realization!, to_Array, to_Matrix, VolumeHypergraph, plot, Polytope, DiskPacking, SphericalDiskPacking

export  ConstraintSystem, 
        Framework, 
        DeformationPath,
        VolumeHypergraph,
        animate,
        plot,
        project_deformation_random,
        Polytope,
        to_Matrix,
        to_Array,
        DiskPacking,
        SphericalDiskPacking,
        equations!,
        realization!,
        newton_correct

mutable struct DeformationPath
    G::ConstraintSystem
    step_size::Float64
    motion_samples::Vector{Vector{Float64}}
    motion_matrices::Vector{Matrix{Float64}}
    flex_mult::Vector{Float64}
    _contacts::Vector

    function DeformationPath(G::ConstraintSystem, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int, type::String; step_size::Float64=1e-2, newton_tol=1e-14)
        start_point = to_Array(G, G.realization)
        if !(type in ["framework", "hypergraph", "polytope",  "sphericaldiskpacking"])
            throw(error("The type must either be 'framework', 'diskpacking', 'sphericaldiskpacking', 'hypergraph' or 'polytope', but is $(type)."))
        end
        if type=="framework"
            K_n = Framework([[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j], G.realization; pinned_vertices=G.pinned_vertices)
        elseif type=="hypergraph"
            K_n = VolumeHypergraph(collect(powerset(G.vertices, G.dimension+1, G.dimension+1)), G.realization)
        elseif type=="polytope" || type == "diskpacking"
            K_n = ConstraintSystem(G.vertices, G.variables, vcat(G.equations, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j]]), G.realization, G.xs; pinned_vertices=G.pinned_vertices)
        elseif  type=="sphericaldiskpacking"
            minkowski_scalar_product(e1,e2) = e1'*e2-1
            inversive_distances = [minkowski_scalar_product(G.realization[:,contact[1]], G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(G.realization[:,contact[1]], G.realization[:,contact[1]]) * minkowski_scalar_product(G.realization[:,contact[2]], G.realization[:,contact[2]])) for contact in powerset(G.vertices, 2, 2)]
            K_n = ConstraintSystem(G.vertices, G.variables, [minkowski_scalar_product(G.xs[:,contact[1]], G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(G.xs[:,contact[1]], G.xs[:,contact[1]]) * minkowski_scalar_product(G.xs[:,contact[2]], G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(G.vertices, 2, 2))], G.realization, G.xs)
        end

        flex_space = compute_nontrivial_inf_flexes(G, start_point, K_n)
        size(flex_space)[2]==length(flex_mult) || throw(error("The length of 'flex_mult' match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2])."))
        prev_flex = sum(flex_mult[i] .* flex_space[:,i] for i in 1:length(flex_mult))
        prev_flex = prev_flex ./ norm(prev_flex)
        motion_samples, motion_matrices = [Float64.(start_point)], [to_Matrix(G, Float64.(start_point))]
        
        failure_to_converge = 0
        @showprogress for i in 1:num_steps
            try
                q, prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end], K_n)
                q = newton_correct(G, q; tol=newton_tol)
                failure_to_converge = 0
                if isapprox(q, motion_samples[end]; atol=1e-12)
                    throw(error("Slow Progress detected."))
                end
                push!(motion_samples, q)
                push!(motion_matrices, to_Matrix(G, Float64.(q)))                    
            catch e
                @warn e
                if failure_to_converge == 3
                    break
                else
                    # If Newton's method only diverges once and we are in a singularity,
                    # we first try to reverse the previous flex before exiting the routine.
                    failure_to_converge += 1
                    if failure_to_converge==1
                        try
                            q, prev_flex = euler_step(G, step_size/10, prev_flex, motion_samples[end], K_n)
                            q = newton_correct(G, q; tol=newton_tol)
                            push!(motion_samples, q)
                            push!(motion_matrices, to_Matrix(G, Float64.(q)))                    
                        catch
                            continue
                        end
                    elseif length(motion_samples)==1 || rank(evaluate.(G.jacobian, G.variables=>motion_samples[end]); atol=1e-8) < rank(evaluate.(G.jacobian, G.variables=>motion_samples[end-1]); atol=1e-8)
                        @warn "Direction was reversed."
                        prev_flex = -prev_flex
                    end
                end
            end
        end
        new(G, step_size, motion_samples, motion_matrices, flex_mult, [])
    end

    function DeformationPath(F::Framework, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; step_size::Float64=1e-2, newton_tol=1e-14)
        DeformationPath(F.G, flex_mult, num_steps, "framework"; step_size=step_size)
    end

    function DeformationPath(F::VolumeHypergraph, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; step_size::Float64=1e-2, newton_tol=1e-14)
        DeformationPath(F.G, flex_mult, num_steps, "hypergraph"; step_size=step_size)
    end

    function DeformationPath(F::Polytope, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; step_size::Float64=1e-2, newton_tol=1e-14)
        DeformationPath(F.G, flex_mult, num_steps, "polytope"; step_size=step_size)
    end

    function DeformationPath(F::SphericalDiskPacking, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; step_size::Float64=1e-2, newton_tol=1e-14)
        DeformationPath(F.G, flex_mult, num_steps, "sphericaldiskpacking"; step_size=step_size)
    end

    function DeformationPath(F::DiskPacking, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; motion_samples::Vector=[], _contacts::Vector=[], step_size::Float64=1e-2, prev_flex=nothing, newton_tol=1e-14)
        start_point = to_Array(F, F.G.realization)
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
        if prev_flex == nothing
            flex_space = compute_nontrivial_inf_flexes(F.G, start_point, K_n)
            size(flex_space)[2]==length(flex_mult) || throw(error("The length of 'flex_mult' must match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2])."))
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
                q, prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end], K_n)
                q = newton_correct(G, q; tol=newton_tol)
                failure_to_converge = 0
                if isapprox(q, motion_samples[end]; atol=1e-12)
                    throw(error("Slow Progress detected."))
                end

                cur_realization = to_Matrix(F,Float64.(q))
                if any(t->norm(cur_realization[:,t[1]] - cur_realization[:,t[2]]) < F.radii[t[1]] + F.radii[t[2]] - F.tolerance, powerset(F.G.vertices,2,2))
                    _F = DiskPacking(F.G.vertices, F.radii, cur_realization; pinned_vertices=F.G.pinned_vertices, tolerance=step_size)
                    DeformationPath(_F, flex_mult, num_steps-i; motion_samples=motion_samples, _contacts=_contacts, step_size=step_size, prev_flex=prev_flex, newton_tol=newton_tol)
                    break
                end

                push!(motion_samples, q)
                push!(_contacts, F.contacts)    
            catch e
                @warn e
                if failure_to_converge == 3
                    break
                else
                    # If Newton's method only diverges once and we are in a singularity,
                    # we first try to reverse the previous flex before exiting the routine.
                    failure_to_converge += 1
                    if failure_to_converge==1
                        try
                            q, prev_flex = euler_step(G, step_size/10, prev_flex, motion_samples[end], K_n)
                            q = newton_correct(G, q; tol=newton_tol)
                            push!(motion_samples, q)
                            push!(motion_matrices, to_Matrix(G, Float64.(q)))                    
                        catch
                            continue
                        end
                    elseif length(motion_samples)==1 || rank(evaluate.(G.jacobian, G.variables=>motion_samples[end]); atol=1e-8) < rank(evaluate.(G.jacobian, G.variables=>motion_samples[end-1]); atol=1e-8)
                        @warn "Direction was reversed."
                        prev_flex = -prev_flex
                    end
                end
            end
        end

        motion_matrices = [to_Matrix(F, Float64.(sample)) for sample in motion_samples]
        new(F.G, step_size, motion_samples, motion_matrices, flex_mult, _contacts)
    end

    function compute_nontrivial_inf_flexes(G::ConstraintSystem, point::Union{Vector{Float64},Vector{Int}}, K_n)
        inf_flexes = nullspace(evaluate(G.jacobian, G.variables=>point); atol=1e-8)
        trivial_inf_flexes = nullspace(evaluate(typeof(K_n)==ConstraintSystem ? K_n.jacobian : K_n.G.jacobian, (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables)=>point[1:length( (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables))]); atol=1e-8)
        s = size(trivial_inf_flexes)[2]+1
        extend_basis_matrix = trivial_inf_flexes
        for inf_flex in [inf_flexes[:,i] for i in 1:size(inf_flexes)[2]]
            tmp_matrix = hcat(trivial_inf_flexes, inf_flex)
            if !(rank(tmp_matrix; atol=1e-8) == rank(trivial_inf_flexes; atol=1e-8))
                extend_basis_matrix = hcat(extend_basis_matrix, inf_flex)
            end
        end
        Q, R = qr(extend_basis_matrix)
        Q = Q[:, s:rank(R, atol=1e-8)]
        return Q
    end

    function euler_step(G::ConstraintSystem, step_size::Float64, prev_flex::Vector{Float64}, point::Union{Vector{Int},Vector{Float64}}, K_n)
        J = evaluate(G.jacobian, G.variables=>point)
        flex_space = compute_nontrivial_inf_flexes(G, point, K_n)
        flex_coefficients = pinv(flex_space) * prev_flex
        predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
        predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
        return point+step_size*predicted_inf_flex, predicted_inf_flex
    end
end

function newton_correct(G::ConstraintSystem, point::Vector{Float64}; tol = 1e-14)
    q = Base.copy(point)
    global damping = 0.15
    start_time=Base.time()
    while(norm(evaluate(G.equations, G.variables=>q)) > tol)
        J = evaluate.(G.jacobian, G.variables=>q)
        stress_dimension = size(nullspace(J'; atol=1e-8))[2]
        if stress_dimension > 0
            rand_mat = randn(Float64, length(G.equations) - stress_dimension, length(G.equations))
            equations = rand_mat*G.equations
            J = rand_mat*J
        else
            equations = G.equations
        end

        qnew = q - damping*pinv(J)*evaluate(equations, G.variables=>q)
        if norm(evaluate(G.equations, G.variables=>qnew)) < norm(evaluate(G.equations, G.variables=>q))
            global damping = damping*1.2
        else
            global damping = damping/2
        end
        if damping < 1e-14 || Base.time()-start_time > 10
            throw(error("Newton's method did not converge in time."))
        end
        q = qnew
        if damping > 1
            global damping = 1
        end
    end
    return q
end


function animate(F, filename::String; flex_mult=nothing, num_steps::Int=100, step_size::Float64=1e-2, kwargs...)
    if flex_mult==nothing
        if typeof(F)==Framework
            K_n = Framework([[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices)
        elseif typeof(F)==VolumeHypergraph
            K_n = VolumeHypergraph(collect(powerset(F.G.vertices, F.G.dimension+1, F.G.dimension+1)), F.G.realization)
        elseif typeof(F)==DiskPacking
            K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs)
        elseif typeof(F)==Polytope
            K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs)
        elseif typeof(F)==SphericalDiskPacking
            minkowski_scalar_product(e1,e2) = e1'*e2-1
            inversive_distances = [minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[1]]) * minkowski_scalar_product(F.G.realization[:,contact[2]], F.G.realization[:,contact[2]])) for contact in powerset(F.G.vertices, 2, 2)]
            K_n = ConstraintSystem(F.G.vertices, F.G.variables, [minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[1]]) * minkowski_scalar_product(F.G.xs[:,contact[2]], F.G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(F.G.vertices, 2, 2))], F.G.realization, F.G.xs)
        else
            throw(error("Type of F is not yet supported. It is $(typeof(F))."))
        end
        flex_space = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), K_n)
        flex_mult = [1 for _ in 1:size(flex_space)[2]]
    end
    D = DeformationPath(F, filename, flex_mult, num_steps; step_size=step_size)
    animate(D, F, filename; kwargs...)
end

function animate(D::DeformationPath, F, filename::String; kwargs...)
    if typeof(F)==Framework
        if F.G.dimension==2
            animate2D_framework(D, F, filename; kwargs...)
        elseif F.G.dimension==3
            animate3D_framework(D, F, filename; kwargs...)
        else
            throw(error("The dimension of 'F' needs to be either 2 or 3, but is $(F.G.dimension)"))
        end
    elseif typeof(F)==VolumeHypergraph
        return animate2D_hypergraph(D, F, filename; kwargs...)
    elseif typeof(F)==Polytope
        return animate3D_polytope(D, F, filename; kwargs...)
    elseif typeof(F)==DiskPacking
        return animate2D_diskpacking(D, F, filename; kwargs...)
    elseif typeof(F)==SphericalDiskPacking
        return animate3D_sphericaldiskpacking(D, F, filename; kwargs...)
    else
        throw(error("The type of 'F' needs to be either Framework, DiskPacking, Polytope, SphericalDiskPacking or VolumeHypergraph, but is $(typeof(F))"))
    end


end

function animate2D_framework(D::DeformationPath, F::Framework, filename::String; recompute_deformation_samples::Bool=true, fixed_pair::Tuple{Int,Int}=(1,2), fixed_direction=[1.,0], framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.15, markercolor=:red3, pin_point_offset=0.2, vertex_size::Union{Float64,Int}=42, line_width::Union{Float64,Int}=10, edge_color=:steelblue, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
    fig = Figure(size=(800,800))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    fixed_pair[1] in D.G.vertices && fixed_pair[2] in D.G.vertices || throw(error("pinned_vertex is not a vertex of the underlying graph."))
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_pair[1]]
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
        theta = atan(matrix_coords[i][:,fixed_pair[2]][2] , matrix_coords[i][:,fixed_pair[2]][1])
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
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=edge_color), F.bars)
    foreach(v->scatter!(ax, @lift([Point2f(($allVertices)[v]-[pin_point_offset,0])]); markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=25, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))

    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw(error("The chosen filetype needs to be either gif or mp4, but is $(filetype)"))
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
    end
end

function animate3D_framework(D::DeformationPath, F::Framework, filename::String; recompute_deformation_samples::Bool=true, fixed_vertices::Union{Tuple{Int,Int}, Tuple{Int,Int,Int}}=(1,2), fixed_direction=[1.,0,0], framerate::Int=25, animate_rotation=false, rotation_start_angle = π / 4, rotation_frames = 240, markercolor=:red3, pin_point_offset=0.05, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_size::Union{Float64,Int}=42, line_width::Union{Float64,Int}=10, edge_color=:steelblue, vertex_color=:black, filetype::String="gif")
    fig = Figure(size=(800,800))
    ax = Axis3(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    length(fixed_vertices)==length(collect(Set(fixed_vertices))) && fixed_vertices[1] in D.G.vertices && fixed_vertices[2] in D.G.vertices && (length(fixed_vertices)==2 || fixed_vertices[3] in D.G.vertices) || throw(error("The elements of `fixed_vertices`` are not vertices of the underlying graph."))
    
    if isapprox(norm(fixed_direction),0;atol=1e-6)
        @warn "fixed_direction is $(norm(fixed_direction)) which is too close to 0! We thus set it to [1,0,0]"
        fixed_direction = [1.,0]
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
            rotation_matrix = [1 0 0; 0 1 0; 0 0 1;]
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
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=:steelblue), F.bars)
    foreach(v->scatter!(ax, @lift([Point3f(($allVertices)[v]-[pin_point_offset,0,0])]); markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(D.G.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw(error("The chosen filetype needs to be either gif or mp4, but is $(filetype)"))
    end

    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.azimuth[] = rotation_start_angle + 2pi * t / rotation_frames
        end
    end
end

function animate2D_hypergraph(D::DeformationPath, F::VolumeHypergraph, filename::String; recompute_deformation_samples::Bool=true, target_stretch::Union{Float64,Int}=1., fixed_triangle::Union{Tuple{Int,Int,Int},Vector{Int},Nothing}=nothing, skip_stretch::Bool=true, tip_value::Union{Float64,Int}=0.5, framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_size::Union{Float64,Int}=42, line_width::Union{Float64,Int}=6, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
    fig = Figure(size=(800,800))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    if facet_colors==nothing
        facet_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(F.facets), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
    end
    skip_scaling = false
    if fixed_triangle==nothing
        fixed_triangle=(1,2,3)
        skip_scaling=true
    else
        all(i->fixed_triangle[i] in D.G.vertices, 1:3) && (Tuple(fixed_triangle) in [Tuple(facet) for facet in F.facets]) || (Tuple([fixed_triangle[2],fixed_triangle[3],fixed_triangle[1]]) in [Tuple(facet) for facet in F.facets]) || (Tuple([fixed_triangle[3],fixed_triangle[1],fixed_triangle[2]]) in [Tuple(facet) for facet in F.facets]) || throw(error("fixed_triangle is not a vertex of the underlying graph."))
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
    foreach(i->poly!(ax, @lift([($allVertices)[Int64(v)] for v in F.facets[i]]); color=(facet_colors[i], 0.25)), 1:length(F.facets))
    foreach(i->lines!(ax, @lift([($allVertices)[Int64(v)] for v in vcat(F.facets[i], F.facets[i][1])]); linewidth=line_width, linestyle=:dash, color=facet_colors[i]), 1:length(F.facets))
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(F.G.vertices))
    foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=26, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))

    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw(error("The chosen filetype needs to be either gif or mp4, but is $(filetype)"))
    end

    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
    end
end

function animate3D_polytope(D::DeformationPath, F::Polytope, filename::String; recompute_deformation_samples::Bool=true, fixed_vertices::Union{Tuple{Int,Int}, Tuple{Int,Int,Int}}=(1,2), framerate::Int=25, animate_rotation=false, rotation_start_angle = π / 4, rotation_frames = 240, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_size::Union{Float64,Int}=42, line_width::Union{Float64,Int}=6, line_color=:steelblue, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]

    fixed_vertices[1] in D.G.vertices && fixed_vertices[2] in D.G.vertices && (length(fixed_vertices)==2 || fixed_vertices[3] in D.G.vertices) || throw(error("The elements of `fixed_vertices`` are not vertices of the underlying graph."))
    ax = Axis3(fig[1,1])
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_vertices[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
        edge_vector = Vector(matrix_coords[i][:,fixed_vertices[2]] ./ norm(matrix_coords[i][:,fixed_vertices[2]]))
        rotation_axis = cross([1,0,0], edge_vector)
        if isapprox(norm(rotation_axis), 0, atol=1e-6)
            rotation_matrix = [1 0 0; 0 1 0; 0 0 1;]
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

    foreach(i->linesegments!(ax, @lift([($allVertices)[Int64(F.edges[i][1])], ($allVertices)[Int64(F.edges[i][2])]]); linewidth=line_width, color=line_color), 1:length(F.edges))
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))
    
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw(error("The chosen filetype needs to be either gif or mp4, but is $(filetype)"))
    end

    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.azimuth[] = rotation_start_angle + 2pi * t / rotation_frames
        end
    end
end


function animate2D_diskpacking(D::DeformationPath, F::DiskPacking, filename::String; framerate::Int=25, step::Int=1, padding::Union{Float64,Int}=0.15, vertex_labels=true, disk_strokewidth::Union{Float64,Int}=8.5, line_width::Union{Float64,Int}=7, disk_color=:steelblue, markersize::Union{Float64,Int}=75, markercolor=:red3, dualgraph_color=:grey80, n_circle_segments::Int=50, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    if F.G.dimension!=2
        throw(error("The dimension must be 2, but is $(F.G.dimension)!"))
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
        poly!(ax, @lift([($disk_vertices)[i] for i in 1:n_circle_segments]); color=(disk_color, 0.08))
        lines!(ax, @lift([($disk_vertices)[v] for v in vcat(1:n_circle_segments,1)]); linewidth = disk_strokewidth, color=disk_color)
    end

    foreach(v->scatter!(ax, @lift([($allVertices)[v]]); markersize=markersize, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    vertex_labels && foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[:black]), 1:length(F.G.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw(error("The chosen filetype needs to be either gif or mp4, but is $(filetype)"))
    end

    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        contacts[] = D._contacts[t]
    end
end

function animate3D_sphericaldiskpacking(D::DeformationPath, F::SphericalDiskPacking, filename::String; framerate::Int=25, animate_rotation=false, rotation_start_angle = π / 4, rotation_frames = 240, step::Int=1, padding=0.015, sphere_color=:lightgrey, vertex_size=60, disk_strokewidth=9, line_width=6, disk_color=:steelblue, dualgraph_color=(:red3,0.45), vertex_color=:black, vertex_labels::Bool=true, n_circle_segments=50, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]

    ax = Axis3(fig[1,1], aspect=(1,1,1))
    xlims!(ax,-1.5-padding, 1.5+padding)
    ylims!(ax,-1.5-padding, 1.5+padding)
    zlims!(ax,-1.5-padding, 1.5+padding)
    hidespines!(ax)
    hidedecorations!(ax)
    mesh!(ax, Sphere(Point3f(0), 1f0); transparency=true, color = (sphere_color,0.15))

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
    vertex_labels && foreach(i->text!(ax, @lift([($rotatedPoints)[i]]), text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[:black]), 1:length(F.G.vertices))
    
    timestamps = range(1, length(D.motion_samples), step=step)
    if !(lowercase(filetype) in ["gif","mp4"])
        throw(error("The chosen filetype needs to be either gif or mp4, but is $(filetype)"))
    end
    
    if animate_rotation
        ax.viewmode = :fit # Prevent axis from resizing during animation
    end
    record(fig, "../data/$(filename).$(lowercase(filetype))", timestamps; framerate = framerate) do t
        time[] = t
        if animate_rotation
            ax.azimuth[] = rotation_start_angle + 2pi * t / rotation_frames
        end
    end
    return fig
end



function project_deformation_random(D::DeformationPath, projected_dimension::Int; line_width::Union{Float64,Int}=8, line_color=:green3, markersize::Union{Float64,Int}=45, markercolor=:steelblue, draw_start::Bool=true)
    if !(projected_dimension in [2,3])
        throw(error("The projected_dimension is neither 2 nor 3."))
    end
    randmat = hcat([rand(Float64,projected_dimension) for _ in 1:length(D.G.variables)]...)
    proj_curve = [(pinv(randmat'*randmat)*randmat')'*entry for entry in D.motion_samples]
    fig = Figure(size=(1000,1000))
    if projected_dimension==3
        ax = Axis3(fig[1,1])
    else
        ax = Axis(fig[1,1])
    end
    hidespines!(ax)
    hidedecorations!(ax)
    if projected_dimension==3
        lines!(ax, [Point3f(proj_curve[i]) for i in 1:length(proj_curve)]; linewidth=line_width, color=line_color)
        draw_start && scatter!(ax, [proj_curve[1][1]], [proj_curve[1][2]], [proj_curve[1][3]]; markersize=markersize, color=markercolor, marker=:pentagon)
    else
        lines!(ax, [Point2f(proj_curve[i]) for i in 1:length(proj_curve)]; linewidth=line_width, color=line_color)
        draw_start && scatter!(ax, [proj_curve[1][1]], [proj_curve[1][2]]; markersize=markersize, color=markercolor, marker=:pentagon)
    end
    return fig
end

end 