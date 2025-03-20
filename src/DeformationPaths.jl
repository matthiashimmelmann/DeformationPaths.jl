module DeformationPaths

import HomotopyContinuation: evaluate
import LinearAlgebra: norm, pinv, nullspace, rank, qr
import GLMakie: @lift, poly!, text!, Figure, record, hidespines!, hidedecorations!, lines!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Observable, Point3f, Point2f
import ProgressMeter: @showprogress
import Combinatorics: powerset
import Colors: distinguishable_colors, red, green, blue, colormap, RGB

include("GeometricConstraintSystem.jl")
using .GeometricConstraintSystem: ConstraintSystem, Framework, to_Array, to_Matrix, VolumeHypergraph, plot

export  ConstraintSystem, 
        Framework, 
        DeformationPath, 
        animate,
        plot,
        project_deformation_random

mutable struct DeformationPath
    G::ConstraintSystem
    step_size::Float64
    motion_samples::Vector{Vector{Float64}}
    flex_mult::Vector{Float64}

    function DeformationPath(G::ConstraintSystem, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int, type::String; step_size::Float64=1e-2)
        start_point = to_Array(G, G.realization)
        flex_space = compute_nontrivial_inf_flexes(G, start_point, type)
        size(flex_space)[2]==length(flex_mult) || throw(error("The length of 'flex_mult' must match the size of the nontrivial infinitesimal flexes, which is $(size(flex_space)[2])."))
        prev_flex = sum(flex_mult[i] .* flex_space[:,i] for i in 1:length(flex_mult))
        prev_flex = prev_flex ./ norm(prev_flex)
        motion_samples = [Float64.(start_point)]
        if !(type in ["framework", "hypergraph"])
            throw(error("The type must either be 'framework' or 'hypergraph', but is $(type)."))
        end
        @showprogress for i in 1:num_steps
            q, prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end], type)
            q = newton_correct(G, q)
            push!(motion_samples, q)
        end
        new(G, step_size, motion_samples, flex_mult)
    end

    function DeformationPath(F::Framework, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; step_size::Float64=1e-2)
        DeformationPath(F.G, flex_mult, num_steps, "framework"; step_size=step_size)
    end

    function DeformationPath(F::VolumeHypergraph, flex_mult::Union{Vector{Float64}, Vector{Int}}, num_steps::Int; step_size::Float64=1e-2)
        DeformationPath(F.G, flex_mult, num_steps, "hypergraph"; step_size=step_size)
    end

    function compute_nontrivial_inf_flexes(G::ConstraintSystem, point::Union{Vector{Float64},Vector{Int}}, type::String)
        inf_flexes = nullspace(evaluate(G.jacobian, G.variables=>point))
        realization = to_Matrix(G,point)
        if type=="framework"
            K_n = Framework([[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j], realization)
        elseif type=="hypergraph"
            K_n = VolumeHypergraph(collect(powerset(G.vertices, G.dimension+1, G.dimension+1)), realization)
        end
        trivial_inf_flexes = nullspace(evaluate(K_n.G.jacobian, G.variables=>point))
        s = size(trivial_inf_flexes)[2]+1
        extend_basis_matrix = trivial_inf_flexes
        for inf_flex in [inf_flexes[:,i] for i in 1:size(inf_flexes)[2]]
            tmp_matrix = hcat(trivial_inf_flexes, inf_flex)
            if !(rank(tmp_matrix; atol=1e-12) == rank(trivial_inf_flexes; atol=1e-12))
                extend_basis_matrix = hcat(extend_basis_matrix, inf_flex)
            end
        end
        Q, R = qr(extend_basis_matrix)
        Q = Q[:, s:rank(R, atol=1e-12)]
        return Q
    end

    function newton_correct(G::ConstraintSystem, point::Vector{Float64}; tol = 1e-14)
        q = point
        global damping = 0.15
        while(norm(evaluate(G.equations, G.variables=>q)) > tol)
            J = evaluate.(G.jacobian, G.variables=>q)
            stress_dimension = size(nullspace(J'; atol=1e-12))[2]
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
            if damping < 1e-15
                throw(error("Newton's method did not converge"))
            end
            q = qnew
            if damping > 1
                global damping = 1
            end
        end
        return q
    end

    function euler_step(G::ConstraintSystem, step_size::Float64, prev_flex::Vector{Float64}, point::Union{Vector{Int},Vector{Float64}}, type)
        J = evaluate(G.jacobian, G.variables=>point)
        flex_space = compute_nontrivial_inf_flexes(G, point, type)
        flex_coefficients = pinv(flex_space) * prev_flex
        predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
        predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
        return point+step_size*predicted_inf_flex, predicted_inf_flex
    end
end

function animate(F, filename::String; flex_mult=nothing, num_steps::Int=100, step_size::Float64=1e-2, kwargs...)
    if flex_mult==nothing
        flex_space = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), typeof(F)==Framework ? "framework" : "hypergraph")
        flex_mult = [1 for _ in 1:size(flex_space)[2]]
    end
    D = DeformationPath(F, filename, flex_mult, num_steps; step_size=step_size)
    if typeof(F)==Framework
        if F.G.dimension==2
            animate2D_framework(F, filename; kwargs...)
        elseif F.G.dimension==3
            animate3D_framework(F, filename; kwargs...)
        else
            throw(error("The dimension of 'F' needs to be either 2 or 3, but is $(F.G.dimension)"))
        end
    elseif typeof(F)==VolumeHypergraph
        return animate2D_hypergraph(F, filename; kwargs...)
    else
        throw(error("The type of 'F' needs to be either Framework or VolumeHypergraph, but is $(typeof(F))"))
    end
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
    else
        throw(error("The type of 'F' needs to be either Framework or VolumeHypergraph, but is $(typeof(F))"))
    end
end

function animate2D_framework(D::DeformationPath, F::Framework, filename::String; fixed_edge::Tuple{Int,Int}=(1,2), framerate::Int=25, step::Int=1, padding::Float64=0.15, vertex_size::Int=42, line_width::Int=10, edge_color=:steelblue, vertex_color=:black, vertex_labels::Bool=true)
    fig = Figure(size=(800,800))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    fixed_edge[1] in D.G.vertices && fixed_edge[2] in D.G.vertices || throw(error("pinned_vertex is not a vertex of the underlying graph."))
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_edge[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
    end
    fixed_direction = [1.,0]
    for i in 1:length(matrix_coords)
        theta = atan(matrix_coords[i][:,fixed_edge[2]][2], matrix_coords[i][:,fixed_edge[2]][1])
        rotation_matrix = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        # Rotate the realization to the `fixed_direction`.
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = rotation_matrix*matrix_coords[i][:,j]
        end
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
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, @lift([($allVertices)[i]]), text=["$(F.G.vertices[i])"], fontsize=25, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))

    timestamps = range(1, length(D.motion_samples), step=step)
    record(fig, "../data/$(filename).gif", timestamps; framerate = framerate) do t
        time[] = t
    end
end

function animate3D_framework(D::DeformationPath, F::Framework, filename::String; pinned_vertex::Int=1, framerate::Int=25, step::Int=1, padding::Float64=0.15, vertex_size::Int=42, line_width::Int=10, edge_color=:steelblue, vertex_color=:black)
    fig = Figure(size=(800,800))
    ax = Axis3(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    pinned_vertex in D.G.vertices || throw(error("pinned_vertex is not a vertex of the underlying graph."))
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,pinned_vertex]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
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
        pointys = matrix_coords[$time]
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=:steelblue), F.bars)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(D.G.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    record(fig, "../data/$(filename).gif", timestamps; framerate = framerate) do t
        time[] = t
    end
end

function animate2D_hypergraph(D::DeformationPath, F::VolumeHypergraph, filename::String; fixed_edge::Tuple{Int,Int}=(1,2), framerate::Int=25, step::Int=1, padding::Float64=0.15, vertex_size::Int=42, line_width::Int=6, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true)
    fig = Figure(size=(800,800))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    if facet_colors==nothing
        facet_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(F.facets), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
    end
    fixed_edge[1] in D.G.vertices && fixed_edge[2] in D.G.vertices || throw(error("pinned_vertex is not a vertex of the underlying graph."))
    for i in 1:length(matrix_coords)
        p0 = matrix_coords[i][:,fixed_edge[1]]
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = matrix_coords[i][:,j] - p0
        end
    end
    fixed_direction = [1.,0]
    for i in 1:length(matrix_coords)
        theta = atan(matrix_coords[i][:,fixed_edge[2]][2], matrix_coords[i][:,fixed_edge[2]][1])
        rotation_matrix = [cos(theta) sin(theta); -sin(theta) cos(theta)]
        # Rotate the realization to the `fixed_direction`.
        for j in 1:size(matrix_coords[i])[2]
            matrix_coords[i][:,j] = rotation_matrix*matrix_coords[i][:,j]
        end
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
    record(fig, "../data/$(filename).gif", timestamps; framerate = framerate) do t
        time[] = t
    end
end


function project_deformation_random(D::DeformationPath, projected_dimension::Int; line_width::Int=8, line_color=:green3, markersize::Int=45, markercolor=:steelblue, draw_start::Bool=true)
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