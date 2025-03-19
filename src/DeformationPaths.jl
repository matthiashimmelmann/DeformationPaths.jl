module DeformationPaths

import HomotopyContinuation: evaluate
import LinearAlgebra: norm, pinv, nullspace
import GLMakie: @lift, Figure, record, hidespines!, hidedecorations!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Observable, Point3f, Point2f

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
            flex_space = nullspace(evaluate(G.jacobian, G.variables=>start_point); atol=1e-12)
            prev_flex = flex_space[:,Int(flex_choice)]
        else
            prev_flex = flex_choice
        end
        motion_samples = [Float64.(start_point)]

        for i in 1:num_steps
            q, prev_flex = euler_step(G, step_size, prev_flex, motion_samples[end])
            q = newton_correct(G, q)
            push!(motion_samples, q)
        end
        new(G, step_size, motion_samples, flex_choice)
    end

    function newton_correct(G::ConstraintSystem, point::Vector{Float64}; tol = 1e-13)
        q = point
        global damping = 0.15
        stress_dimension = size(nullspace(evaluate.(G.jacobian, G.variables=>point)'; atol=1e-12))[2]
        while(norm(evaluate(G.equations, G.variables=>q)) > tol)
            J = evaluate.(G.jacobian, G.variables=>q)
            rand_mat = randn(Float64, length(G.equations) - stress_dimension, length(G.equations))
            if stress_dimension > 0
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
        flex_space = nullspace(evaluate(G.jacobian, G.variables=>point); atol=1e-12)
        flex_coefficients = pinv(flex_space) * prev_flex
        predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in 1:length(flex_coefficients))
        predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
        return point+step_size*predicted_inf_flex, predicted_inf_flex
    end
end

function animate2D_framework(D::DeformationPath, F::Framework, filename::String; framerate::Int=25, step::Int=1, padding=0.15)
    fig = Figure(size=(1000,1000))
    ax = Axis(fig[1,1], aspect = 1)
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    xlims!(ax, xlims[1]-padding,xlims[2]+padding)
    ylims!(ax, ylims[1]-padding,ylims[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    allVertices=@lift begin
        pointys = to_Matrix(D.G, D.motion_samples[$time])
        [Point2f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = 12, color=:steelblue), F.bars)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = 50, color=:black), 1:length(F.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    record(fig, "$(filename).mp4", timestamps; framerate = framerate) do t
        time[] = t
    end
end

function animate3D_framework(D::DeformationPath, F::Framework, filename::String; framerate::Int=25, step::Int=1, padding=0.15)
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect = 1)
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    xlims = [minimum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][1,:] for i in 1:length(matrix_coords)]...))]
    ylims = [minimum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][2,:] for i in 1:length(matrix_coords)]...))]
    zlims = [minimum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...)), maximum(vcat([matrix_coords[i][3,:] for i in 1:length(matrix_coords)]...))]
    xlims!(ax, xlims[1]-padding,xlims[2]+padding)
    ylims!(ax, ylims[1]-padding,ylims[2]+padding)
    zlims!(ax, zlims[1]-padding,zlims[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    time=Observable(1)
    allVertices=@lift begin
        pointys = to_Matrix(D.G, D.motion_samples[$time])
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end
    foreach(edge->linesegments!(ax, @lift([($allVertices)[Int64(edge[1])], ($allVertices)[Int64(edge[2])]]); linewidth = 12, color=:steelblue), F.bars)
    foreach(i->scatter!(ax, @lift([($allVertices)[i]]); markersize = 50, color=:black), 1:length(F.vertices))
    timestamps = range(1, length(D.motion_samples), step=step)
    record(fig, "$(filename).mp4", timestamps; framerate = framerate) do t
        time[] = t
    end
end

end 
