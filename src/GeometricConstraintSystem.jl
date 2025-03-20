module GeometricConstraintSystem

import HomotopyContinuation: @var, evaluate, newton, Variable, Expression, differentiate
import GLMakie: Figure, text!, poly!, lines!, save, hidespines!, hidedecorations!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Point3f, Point2f
import LinearAlgebra: det
import Colors: distinguishable_colors, red, green, blue, colormap, RGB

export GeometricConstraintSystem, Framework, plot_framework, VolumeHypergraph, plot_hypergraph

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
    facets::Vector{Vector{Int}}

    function VolumeHypergraph(vertices::Vector{Int}, facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        dimension = size(realization)[1]
        all(t->length(t)==dimension+1, facets) && all(facet->all(v->v in vertices, facet), facets) || throw(error("The facets don't have the correct format."))
        facets = [Vector(facet) for facet in facets]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices)))]...)
        facet_equations = [det(vcat([1. for _ in 1:dimension+1]', hcat([x[:,v] for v in facet]...))) - det(vcat([1. for _ in 1:dimension+1]', hcat([realization[:,v] for v in facet]...))) for facet in facets]
        facet_equations = filter(eq->eq!=0, facet_equations)
        G = ConstraintSystem(vertices,variables, facet_equations, realization)
        new(G, facets)
    end

    function VolumeHypergraph(facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in facets]...))))
        dimension = size(realization)[1]
        VolumeHypergraph(vertices, facets, realization)
    end
end


function to_Array(G::ConstraintSystem, p::Union{Matrix{Int},Matrix{Float64}})
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(G.realization)[1], 1:size(G.realization)[2]))]...)
end

function to_Array(F::Framework, p::Union{Matrix{Int},Matrix{Float64}})
    return to_Array(F.G, p)
end

function to_Array(F::VolumeHypergraph, p::Union{Matrix{Int},Matrix{Float64}})
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

function to_Matrix(F::VolumeHypergraph, q::Union{Vector{Float64}, Vector{Int}})
    return to_Matrix(F.G, q)
end

function plot(F, filename::String; kwargs...)
    if typeof(F)==Framework
        return plot_framework(F, filename; kwargs...)
    elseif typeof(F)==VolumeHypergraph
        return plot_hypergraph(F, filename; kwargs...)
    else
        throw(error("The type of 'F' needs to be either Framework or VolumeHypergraph, but is $(typeof(F))"))
    end
end

function plot_framework(F::Framework, filename::String; padding::Float64=0.15, vertex_size::Int=60, line_width::Int=12, edge_color=:steelblue, vertex_color=:black)
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
    F.G.dimension==3 && zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = F.G.dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=edge_color), F.bars)
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end

function plot_hypergraph(F::VolumeHypergraph, filename::String; padding::Float64=0.15, vertex_size::Int=60, line_width::Int=8, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization    
    if facet_colors==nothing
        facet_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(F.facets), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
    end

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
    F.G.dimension==3 && zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = F.G.dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(i->poly!(ax, [(allVertices)[Int64(v)] for v in F.facets[i]]; color=(facet_colors[i], 0.25)), 1:length(F.facets))
    foreach(i->lines!(ax, [(allVertices)[Int64(v)] for v in vcat(F.facets[i], F.facets[i][1])]; linewidth=line_width, linestyle=:dash, color=facet_colors[i]), 1:length(F.facets))
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end

end
