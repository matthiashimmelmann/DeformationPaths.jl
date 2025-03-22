module GeometricConstraintSystem

import HomotopyContinuation: @var, evaluate, newton, Variable, Expression, differentiate, System
import GLMakie: Figure, text!, poly!, lines!, save, hidespines!, hidedecorations!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Point3f, Point2f, mesh!
import LinearAlgebra: det, cross, norm
import Colors: distinguishable_colors, red, green, blue, colormap, RGB
import Combinatorics: powerset

export GeometricConstraintSystem, Framework, plot_framework, VolumeHypergraph, Polytope, plot_hypergraph, to_Matrix, to_Array, DiskPacking

mutable struct ConstraintSystem
    vertices::Vector{Int}
    variables::Vector{Variable}
    equations::Vector{Expression}
    realization::Union{Matrix{Int},Matrix{Float64}}
    jacobian::Matrix{Expression}
    dimension::Int
    xs::Union{Matrix{Variable},Matrix{Expression}}
    system::System
    pinned_vertices::Vector{Int}

    function ConstraintSystem(vertices,variables, equations, realization, xs; pinned_vertices::Vector{Int}=[])
        jacobian = hcat([differentiate(eq, variables) for eq in equations]...)'
        dimension = size(realization)[1]
        size(realization)[1]==dimension && (size(realization)[2]==length(vertices) || size(realization)[2]==length(variables)//dimension) || throw(error("The realization does not have the correct format."))
        size(xs)[1]==size(realization)[1] && size(xs)[2]==size(realization)[2] || throw(error("The matrix 'xs' does not have the correct format."))
        new(vertices, variables, equations, realization, jacobian, dimension, xs, System(equations; variables=variables), pinned_vertices)
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
        G = ConstraintSystem(vertices,variables, bar_equations, realization, x)
        new(G, bars)
    end

    function Framework(bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = sort(collect(Set(vcat([bar[1] for bar in bars], [bar[2] for bar in bars]))))
        dimension = size(realization)[1]
        Framework(vertices, bars, realization)
    end
end

mutable struct DiskPacking
    G::ConstraintSystem
    contacts::Vector{Tuple{Int,Int}}
    radii::Union{Vector{Int},Vector{Float64}}
    pinned_vertices::Vector{Int}

    function DiskPacking(vertices::Vector{Int}, radii::Union{Vector{Int},Vector{Float64}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=[])
        length(vertices)==length(radii) && length(radii)==size(realization)[2] || throw(error(("The length of the radii does not match the length of the vertices or the dimensionality of the realization.")))
        all(v->v in vertices, pinned_vertices) || throw(error("Some of the pinned_vertices are not contained in vertices."))
        dimension = size(realization)[1]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        all(t->isapprox(norm(realization[:,t[1]]-realization[:,t[2]]) >=radii[t[1]]+radii[t[2]]-1e-12), powerset(length(vertices),2,2)) || throw(error("Some of the disks are too close"))
        contacts = [Tuple([i,j]) for i in 1:length(vertices) for j in i+1:length(vertices) if isapprox(norm(realization[:,i]-realization[:,j]),radii[i]+radii[j],atol=1e-12)]
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in contacts]
        bar_equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices, variables, bar_equations, realization, xs; pinned_vertices=pinned_vertices)
        new(G, contacts, radii)
    end

    function DiskPacking(radii::Union{Vector{Int},Vector{Float64}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=[])
        vertices = [i for i in 1:length(radii)]
        DiskPacking(vertices, radii, realization; pinned_vertices=pinned_vertices)
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
        G = ConstraintSystem(vertices,variables, facet_equations, realization, x)
        new(G, facets)
    end

    function VolumeHypergraph(facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in facets]...))))
        dimension = size(realization)[1]
        VolumeHypergraph(vertices, facets, realization)
    end
end

mutable struct Polytope
    G::ConstraintSystem
    facets::Vector{Vector{Int}}
    edges::Vector{Tuple{Int,Int}}
    x_variables::Vector{Variable}
    n_variables::Vector{Variable}

    function Polytope(vertices::Vector{Int}, facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        dimension = size(realization)[1]
        dimension==3 || raise(error("The dimension needs to be 3, but is $(dimension)"))
        all(facet->all(v->v in vertices, facet), facets) && all(facet->length(facet)>=3, facets) || throw(error("The facets don't have the correct format. They need to contain at least 3 vertices each."))
        facets = [Vector(facet) for facet in facets]
        !(size(realization)[1]==dimension && size(realization)[2]==length(vertices)) && throw(error("The realization does not have the correct format."))
        normal_realization, bars = Array{Float64,2}(undef, 3, length(facets)), []
        for j in 1:length(facets)
            normal_realization[:,j] = cross(realization[:,facets[j][2]] - realization[:,facets[j][1]], realization[:,facets[j][3]] - realization[:,facets[j][2]])
            normal_realization[:,j] = normal_realization[:,j] ./ norm(normal_realization[:,j])
            for i in j+1:length(facets)
                edge = facets[i][findall(q -> q in facets[j], facets[i])]
                if length(edge) != 2
                    continue
                end
                push!(bars, Tuple(edge))
            end
        end
        _realization = hcat(realization, normal_realization)
        length(vertices)-length(bars)+length(facets)==2 || throw(error("The Euler characteristic of the Polytope needs to be 2, but is $(length(vertices)-length(bars)+length(facets))"))
                
        @var x[1:dimension, 1:length(vertices)] n[1:dimension, 1:length(facets)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices)))]...)
        normal_variables = vcat([n[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(facets)))]...)
        facet_equations = vcat([[n[:,i]'*(x[:,facets[i][j]]-x[:,(facets[i][j%length(facets[i])+1])]) for j in 1:length(facets[i])] for i in 1:length(facets)]...)
        normal_equations = [n[:,i]'*n[:,i]-1 for i in 1:length(facets)]
        bar_equations = [sum( (x[:,bar[1]]-x[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        equations = filter(eq->eq!=0, vcat(facet_equations, normal_equations, bar_equations))
        G = ConstraintSystem(vertices, vcat(variables, normal_variables), equations, _realization, hcat(x,n))
        new(G, facets, bars, variables, normal_variables)
    end

    function Polytope(facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in facets]...))))
        Polytope(vertices, facets, realization)
    end
end



function to_Array(G::ConstraintSystem, p::Union{Matrix{Int},Matrix{Float64}})
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(G.realization)[1], 1:size(G.realization)[2])) if !(j in F.pinned_vertices)]...)
end

function to_Array(F::Union{Framework,VolumeHypergraph,Polytope,DiskPacking}, p::Union{Matrix{Int},Matrix{Float64}})
    return to_Array(F.G, p)
end



function to_Matrix(G::ConstraintSystem, q::Union{Vector{Float64}, Vector{Int}})
    counts = 1
    point = Matrix{Float64}(Base.copy(G.realization))

    for i in 1:size(point)[2]
        for j in 1:size(point)[1]
            if j in G.pinned_vertices
                point[j,i] = G.realization[j,i]
                continue
            end
            point[j,i] = q[counts]
            counts += 1
        end
    end
    return point
end


function to_Matrix(F::Union{Framework,VolumeHypergraph,Polytope,DiskPacking}, q::Union{Vector{Float64}, Vector{Int}})
    return to_Matrix(F.G, q)
end


function plot(F, filename::String; kwargs...)
    if typeof(F)==Framework
        return plot_framework(F, filename; kwargs...)
    elseif typeof(F)==VolumeHypergraph
        return plot_hypergraph(F, filename; kwargs...)
    elseif typeof(F)==Polytope
        return plot_polytope(F, filename; kwargs...)
    elseif typeof(F)==DiskPacking
        return plot_diskpacking(F, filename; kwargs...)
    else
        throw(error("The type of 'F' needs to be either Framework, Polytope, DiskPacking or VolumeHypergraph, but is $(typeof(F))"))
    end
end

function plot_framework(F::Framework, filename::String; padding::Float64=0.15, vertex_size=60, line_width=12, edge_color=:steelblue, vertex_color=:black)
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

    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
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

function plot_diskpacking(F::DiskPacking, filename::String; padding::Float64=0.15, disk_strokewidth=8.5, vertex_labels::Bool=true, disk_color=:steelblue, markersize=75, markercolor=:red3, line_width=7, dualgraph_color=:grey85, n_circle_segments::Int=50)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization
    ax = Axis(fig[1,1])
    if F.G.dimension!=2
        throw(error("The dimension must be 2, but is $(F.G.dimension)!"))
    end
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=dualgraph_color), F.contacts)
    limit_vertices = allVertices
    for index in 1:length(F.G.vertices)
        disk_vertices = [Vector(allVertices[index])+F.radii[index]*Point2f([cos(2*i*pi/n_circle_segments), sin(2*i*pi/n_circle_segments)]) for i in 1:n_circle_segments]
        limit_vertices = vcat(limit_vertices, disk_vertices)
        diskedges = [(i,i%n_circle_segments+1) for i in 1:n_circle_segments]
        poly!(ax, [(disk_vertices)[i] for i in 1:n_circle_segments]; color=(disk_color, 0.1))
        lines!(ax, [(disk_vertices)[v] for v in vcat(1:n_circle_segments,1)]; linewidth = disk_strokewidth, color=disk_color)
    end
    xlims = [minimum([limit_vertices[i][1] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][1] for i in 1:length(limit_vertices)])]
    ylims = [minimum([limit_vertices[i][2] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][2] for i in 1:length(limit_vertices)])]
    limits= [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])]

    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)

    foreach(v->scatter!(ax, [(allVertices)[v]]; markersize=markersize, color=(markercolor, 0.4), marker=:utriangle), F.pinned_vertices)
    vertex_labels && foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=35, font=:bold, align = (:center, :center), color=[:black]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end


function plot_hypergraph(F::VolumeHypergraph, filename::String; padding::Float64=0.15, vertex_size=60, line_width=8, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true)
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
    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
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

function plot_polytope(F::Polytope, filename::String; padding=0.15, vertex_size=60, line_width=8, line_color=:steelblue, facet_color=:lightgrey, vertex_color=:black, vertex_labels::Bool=true)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization    

    ax = Axis3(fig[1,1])
    xlims = [minimum(vcat(matrix_coords[1,:])), maximum(matrix_coords[1,:])]
    ylims = [minimum(vcat(matrix_coords[2,:])), maximum(matrix_coords[2,:])]
    zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
    limits= [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]
    xlims!(ax, limits[1]-padding, limits[2]+padding)
    ylims!(ax, limits[1]-padding, limits[2]+padding)
    zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    #=for i in 1:length(F.facets)
        facet_matrix = vcat([[F.facets[i][1] F.facets[i][j] F.facets[i][j+1];] for j in 2:length(F.facets[i])-1]...)
        mesh!(ax, matrix_coords, facet_matrix; transparency=true, color=(facet_color, 0.25)), 1:length(F.facets)
    end=#
    foreach(i->linesegments!(ax, [(allVertices)[Int64(F.edges[i][1])], (allVertices)[Int64(F.edges[i][2])]]; linewidth=line_width, color=line_color), 1:length(F.edges))
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end


end
