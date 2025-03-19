module GeometricConstraintSystem

import HomotopyContinuation: @var, evaluate, newton, Variable, Expression, differentiate
import GLMakie: Figure, save, hidespines!, hidedecorations!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Point3f, Point2f

export GeometricConstraintSystem, Framework, plot_framework

mutable struct ConstraintSystem
    variables::Vector{Variable}
    equations::Vector{Expression}
    realization::Union{Matrix{Int},Matrix{Float64}}
    jacobian::Matrix{Expression}
    dimension::Int
    pinned_indices::Vector{Tuple{Int,Int}}

    function ConstraintSystem(variables, equations, realization, pinned_indices)
        jacobian = hcat([differentiate(eq, variables) for eq in equations]...)'
        dimension = size(realization)[1]
        pinned_indices = [Tuple(p_i) for p_i in pinned_indices]
        size(realization)[1]==dimension && size(realization)[2]==(length(variables)+length(pinned_indices))//dimension || throw(error("The realization does not have the correct format."))
        new(variables, equations, realization, jacobian, dimension, pinned_indices)
    end
end

mutable struct Framework
    G::ConstraintSystem
    vertices::Vector{Int}
    bars::Vector{Tuple{Int,Int}}

    function Framework(vertices::Vector{Int}, bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, pinned_indices::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        all(t->length(t)==2, bars) && all(bar->bar[1] in vertices && bar[2] in vertices, bars) || throw(error("The bars don't have the correct format."))
        all(t->length(t)==2, pinned_indices) || throw(error("The pinned indices don't have the correct format."))
        dimension = size(realization)[1]
        bars = [bar[1]<=bar[2] ? Tuple(bar) : Tuple([bar[2],bar[1]]) for bar in bars]
        pinned_indices = [Tuple(p_i) for p_i in pinned_indices]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, size(realization)[1], size(realization)[2])
        xs .= x
        for pin in pinned_indices
            xs[pin[1],pin[2]] = realization[pin[1],pin[2]]
        end
        variables = [x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices))) if !((i,j) in pinned_indices)]
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        G = ConstraintSystem(variables, bar_equations, realization, pinned_indices)
        new(G, vertices, bars)
    end

    function Framework(bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = collect(Set(vcat([bar[1] for bar in bars], [bar[2] for bar in bars])))
        dimension = size(realization)[1]
        pinned_indices = vcat([(j,k) for k in 1:dimension for j in 1:1+dimension-k]...)
        Framework(vertices, bars, pinned_indices, realization)
    end
end

function to_Array(G::ConstraintSystem, p::Union{Matrix{Int},Matrix{Float64}})
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(G.realization)[1], 1:size(G.realization)[2])) if !((i,j) in G.pinned_indices)]...)
end

function to_Array(F::Framework, p::Union{Matrix{Int},Matrix{Float64}})
    return to_Array(F.G, p)
end

function to_Matrix(G::ConstraintSystem, q::Union{Vector{Float64}, Vector{Int}})
    point = zeros(typeof(q[1]),size(G.realization)[1],size(G.realization)[2])
    counts = 1
    point .= G.realization

    for i in 1:size(point)[2]
        for j in 1:size(point)[1]
            if (j,i) in G.pinned_indices
                continue
            end
            point[j,i] = q[counts]
            counts += 1
        end
    end
    return point
end

function to_Matrix(F::Framework, q::Union{Vector{Float64}, Vector{Int}})
    return to_Matrix(F.G, q)
end

function plot_framework(F::Framework, filename::String; padding=0.15, markersize=50, linewidth=12)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization
    if F.G.dimension==2
        ax = Axis(fig[1,1], aspect = 1)
    elseif F.G.dimension==3
        ax = Axis3(fig[1,1], aspect = 1)
        zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
        zlims!(ax, zlims[1]-padding, zlims[2]+padding)
    else
        throw(error("The dimension must either be 2 or 3!"))
    end
    xlims = [minimum(vcat(matrix_coords[1,:])), maximum(matrix_coords[1,:])]
    ylims = [minimum(vcat(matrix_coords[2,:])), maximum(matrix_coords[2,:])]
    xlims!(ax, xlims[1]-padding, xlims[2]+padding)
    ylims!(ax, ylims[1]-padding, ylims[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = F.G.dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = linewidth, color=:steelblue), F.bars)
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = markersize, color=:black), 1:length(F.vertices))
    save("$(filename).png", fig)
end

end
