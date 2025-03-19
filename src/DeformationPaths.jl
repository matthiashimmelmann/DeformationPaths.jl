module DeformationPaths

import HomotopyContinuation: evaluate
import LinearAlgebra: norm, pinv, nullspace
import GLMakie: @lift, Figure, record, hidespines!, hidedecorations!, lines!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Observable, Point3f, Point2f
import ProgressMeter: @showprogress

include("GeometricConstraintSystem.jl")
using .GeometricConstraintSystem: ConstraintSystem, Framework, to_Array, to_Matrix, plot_framework

export ConstraintSystem, Framework, DeformationPath, animate2D_framework, animate3D_framework, plot_framework

mutable struct DeformationPath
    G::ConstraintSystem
    step_size::Float64
    motion_samples::Vector{Vector{Float64}}
    flex_choice::Union{Int, Vector{Float64}}

    function DeformationPath(G::ConstraintSystem, flex_choice::Union{Int, Vector{Float64}}, num_steps::Int; step_size::Float64=1e-2)
        start_point = to_Array(G, G.realization)
        if typeof(flex_choice)==Int || typeof(flex_choice)==Float64
            J = evaluate(G.jacobian, G.variables=>start_point)
            flex_space = nullspace(J; atol=1e-12)
            prev_flex = flex_space[:,Int(flex_choice)]
        else
            prev_flex = flex_choice
        end
        motion_samples = [Float64.(start_point)]

        @showprogress for i in 1:num_steps
            q, prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end])
            q = newton_correct(G, q)
            push!(motion_samples, q)
        end
        new(G, step_size, motion_samples, flex_choice)
    end

    function newton_correct(G::ConstraintSystem, point::Vector{Float64}; tol = 1e-13)
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

    function euler_step(G::ConstraintSystem, step_size::Float64, prev_flex::Vector{Float64}, point::Union{Vector{Int},Vector{Float64}})
        J = evaluate(G.jacobian, G.variables=>point)
        stress_dimension = size(nullspace(J'; atol=1e-10))[2]
        if stress_dimension > 0
            rand_mat = randn(Float64, length(G.equations) - stress_dimension, length(G.equations))
            J = rand_mat*J
        end
        flex_space = nullspace(J; atol=1e-10)
        flex_coefficients = pinv(flex_space) * prev_flex
        predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
        predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
        return point+step_size*predicted_inf_flex, predicted_inf_flex
    end
end

function animate2D_framework(D::DeformationPath, F::Framework, filename::String; framerate::Int=25, step::Int=1, padding::Float64=0.15, vertex_size::Int=50, line_width::Int=12, edge_color=:steelblue, vertex_color=:black)
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1,1])
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    limits= [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])]
    translation = limits[1]-xlims[1] + xlims[2]-limits[2]
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = limits[1]-ylims[1] + ylims[2]-limits[2]
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    allVertices=@lift begin
        pointys = to_Matrix(D.G, D.motion_samples[$time])
        [Point2f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=edge_color), F.bars)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=vertex_color), 1:length(F.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    record(fig, "$(filename).mp4", timestamps; framerate = framerate) do t
        time[] = t
    end
end

function animate3D_framework(D::DeformationPath, F::Framework, filename::String; framerate::Int=25, step::Int=1, padding::Float64=0.15, vertex_size::Int=50, line_width::Int=12, edge_color=:steelblue, vertex_color=:black)
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1])
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

    time=Observable(1)
    allVertices=@lift begin
        pointys = to_Matrix(D.G, D.motion_samples[$time])
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = line_width, color=:steelblue), F.bars)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = vertex_size, color=:black), 1:length(F.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    record(fig, "$(filename).mp4", timestamps; framerate = framerate) do t
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

F = Framework(vcat([[1,2],[2,3],[3,4],[1,4],[1,5],[2,6],[3,7],[4,8],[5,6],[6,7],[7,8],[5,8]],[[i,9] for i in 1:8]), Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1; 0 0 sqrt(2)]'))
display(F.G.pinned_indices)
D = DeformationPath(F.G, 1, 1500; step_size=0.02)
animate3D_framework(D,F,"coned_cube_motion.jl")
end 
