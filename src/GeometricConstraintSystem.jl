module GeometricConstraintSystem

import HomotopyContinuation: @var, evaluate, newton, Variable, Expression, differentiate
import GLMakie: Figure, save, hidespines!, hidedecorations!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Point3f, Point2f

export GeometricConstraintSystem, Framework, plot_framework

mutable struct ConstraintSystem
    vertices::Vector{Int}
    variables::Vector{Variable}
    equations::Vector{Expression}
    realization::Union{Matrix{Int},Matrix{Float64}}
    jacobian::Matrix{Expression}
    dimension::Int

    function ConstraintSystem(vertices,variables, equations, realization)
        jacobian = hcat([differentiate(eq, variables) for eq in equations]...)'
        dimension = size(realization)[1]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        new(vertices, variables, equations, realization, jacobian, dimension)
    end
end

mutable struct Framework
    G::ConstraintSystem
    bars::Vector{Tuple{Int,Int}}

    function Framework(vertices::Vector{Int}, bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        all(t->length(t)==2, bars) && all(bar->bar[1] in vertices && bar[2] in vertices, bars) || throw(error("The bars don't have the correct format."))
        dimension = size(realization)[1]
        bars = [bar[1]<=bar[2] ? Tuple(bar) : Tuple([bar[2],bar[1]]) for bar in bars]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices)))]...)
        bar_equations = [sum( (x[:,bar[1]]-x[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        bar_equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices,variables, bar_equations, realization)
        new(G, bars)
    end

    function Framework(bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = sort(collect(Set(vcat([bar[1] for bar in bars], [bar[2] for bar in bars]))))
        dimension = size(realization)[1]
        Framework(vertices, bars, realization)
    end
end

mutable struct VolumeHypergraph
    G::ConstraintSystem
    facets::Vector{Tuple{Int,Int,Int}}

    function VolumeHypergraph(vertices::Vector{Int}, facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        all(t->length(t)==3, bars) && all(facet->facet[1] in vertices && facet[2] in vertices && facet[3] in vertices, facets) || throw(error("The facets don't have the correct format."))
        dimension = size(realization)[1]
        facets = [Tuple(sort(facet)) for facet in facets]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices)))]...)
        facet_equations = [sum( (x[:,bar[1]]-x[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        facet_equations = filter(eq->eq!=0, facet_equations)
        G = ConstraintSystem(vertices,variables, bar_equations, realization)
        new(G, facets)
    end

    function VolumeHypergraph(facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = collect(Set(vcat([bar[1] for bar in facets], [bar[2] for bar in facets], [bar[3] for bar in facets])))
        dimension = size(realization)[1]
        VolumeHypergraph(vertices, bars, realization)
    end
end


function to_Array(G::ConstraintSystem, p::Union{Matrix{Int},Matrix{Float64}})
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(G.realization)[1], 1:size(G.realization)[2]))]...)
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
            point[j,i] = q[counts]
            counts += 1
        end
    end
    return point
end

function to_Matrix(F::Framework, q::Union{Vector{Float64}, Vector{Int}})
    return to_Matrix(F.G, q)
end

function plot_framework(F::Framework, filename::String; padding::Float64=0.15, vertex_size::Int=50, line_width::Int=12, edge_color=:steelblue, vertex_color=:black)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization
    if F.G.dimension==2
        ax = Axis(fig[1,1])
    elseif F.G.dimension==3
        ax = Axis3(fig[1,1])
        zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
    else
        throw(error("The dimension must either be 2 or 3!"))
    end
    xlims = [minimum(vcat(matrix_coords[1,:])), maximum(matrix_coords[1,:])]
    ylims = [minimum(vcat(matrix_coords[2,:])), maximum(matrix_coords[2,:])]
    limits= F.G.dimension==2 ? [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])] : [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]

    xlims!(ax, limits[1]-padding, limits[2]+padding)
    ylims!(ax, limits[1]-padding, limits[2]+padding)
    F.G.dimension==2 && zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = F.G.dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=edge_color), F.bars)
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = markersize, color=vertex_color), 1:length(F.vertices))
    save("$(filename).png", fig)
end

end
