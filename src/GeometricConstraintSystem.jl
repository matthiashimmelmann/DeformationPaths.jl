module GeometricConstraintSystem

import HomotopyContinuation: @var, evaluate, newton, Variable, Expression, differentiate, System
import GLMakie: Figure, text!, poly!, lines!, save, hidespines!, hidedecorations!, linesegments!, scatter!, Axis, Axis3, xlims!, ylims!, zlims!, Point3f, Point2f, mesh!, Sphere
import LinearAlgebra: det, cross, norm, inv
import Colors: distinguishable_colors, red, green, blue, colormap, RGB
import Combinatorics: powerset

export GeometricConstraintSystem, Framework, plot_framework, VolumeHypergraph, Polytope, plot_hypergraph, to_Matrix, to_Array, SpherePacking, SphericalDiskPacking, equations!, realization!

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

    function ConstraintSystem(vertices,variables, equations, realization, xs; pinned_vertices::Vector{Int}=Vector{Int}([]))
        jacobian = hcat([differentiate(eq, variables) for eq in equations]...)'
        dimension = size(realization)[1]
        size(realization)[1]==dimension && (size(realization)[2]==length(vertices) || size(realization)[2]==length(variables)//dimension) || throw(error("The realization does not have the correct format."))
        size(xs)[1]==size(realization)[1] && size(xs)[2]==size(realization)[2] || throw(error("The matrix 'xs' does not have the correct format."))
        new(vertices, variables, equations, realization, jacobian, dimension, xs, System(equations; variables=variables), pinned_vertices)
    end
end

function equations!(G::ConstraintSystem, equations::Vector{Expression})
    set(System(equations).variables)==G.variables || throw(error("The variables in `equations` do not match the original variables."))
    G.equations = equations
    G.jacobian = hcat([differentiate(eq, G.variables) for eq in equations]...)'
    G.system = System(equations; variables=G.variables)
    return nothing
end

function realization!(G::ConstraintSystem, realization::Union{Matrix{Float64},Matrix{Int}})
    size(realization)[1]==G.dimension && (size(realization)[2]==length(G.vertices) || size(realization)[2]==length(G.variables)//G.dimension) || throw(error("The realization does not have the correct format."))
    point = to_Array(G, realization)
    re_matrix = to_Matrix(G, point)
    norm(evaluate(G.equations, G.variables=>point))>1e-8 && throw(error("The point does not satisfy the constraint system's equations!"))
    G.realization = realization
    return nothing
end



mutable struct Framework
    G::ConstraintSystem
    bars::Vector{Tuple{Int,Int}}

    function Framework(vertices::Vector{Int}, bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices=Vector{Int}([]))
        all(t->length(t)==2, bars) && all(bar->bar[1] in vertices && bar[2] in vertices, bars) || throw(error("The bars don't have the correct format."))
        dimension = size(realization)[1]
        all(v->v in vertices, pinned_vertices) || throw(error("Some of the pinned_vertices are not contained in vertices."))
        bars = [bar[1]<=bar[2] ? Tuple(bar) : Tuple([bar[2],bar[1]]) for bar in bars]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        bar_equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices,variables, bar_equations, realization, xs; pinned_vertices=pinned_vertices)
        new(G, bars)
    end

    function Framework(bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([bar[1] for bar in bars], [bar[2] for bar in bars]))))
        dimension = size(realization)[1]
        Framework(vertices, bars, realization; pinned_vertices=pinned_vertices)
    end
end

function equations!(F::Framework, equations::Vector{Expression})
    F.G.equations = equations
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in equations]...)'
    F.G.system = System(equations; variables=F.G.variables)
    return nothing
end


mutable struct SpherePacking
    G::ConstraintSystem
    contacts::Vector{Tuple{Int,Int}}
    radii::Union{Vector{Int},Vector{Float64}}
    tolerance::Float64

    function SpherePacking(vertices::Vector{Int}, radii::Union{Vector{Int},Vector{Float64}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=Vector{Int}([]), tolerance::Float64=1e-8)
        length(vertices)==length(radii) && length(radii)==size(realization)[2] && all(r->r>0, radii) || throw(error(("The length of the radii does not match the length of the vertices or the dimensionality of the realization.")))
        all(v->v in vertices, pinned_vertices) || throw(error("Some of the pinned_vertices are not contained in vertices."))
        dimension = size(realization)[1]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        all(t->norm(realization[:,t[1]]-realization[:,t[2]]) >= radii[t[1]]+radii[t[2]]-tolerance, powerset(length(vertices),2,2)) || throw(error("Some of the disks are too close"))
        contacts = [Tuple([i,j]) for i in 1:length(vertices) for j in i+1:length(vertices) if isapprox(norm(realization[:,i]-realization[:,j]),radii[i]+radii[j],atol=tolerance)]
        dimension in [2,3] || raise(error("The dimension for SpherePacking must be 2."))
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - (radii[bar[1]]+radii[bar[2]])^2 for bar in contacts]
        bar_equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices, variables, bar_equations, realization, xs; pinned_vertices=pinned_vertices)
        new(G, contacts, radii, tolerance)
    end

    function SpherePacking(radii::Union{Vector{Int},Vector{Float64}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        vertices = [i for i in 1:length(radii)]
        SpherePacking(vertices, radii, realization; pinned_vertices=pinned_vertices)
    end
end

function equations!(F::SpherePacking, equations::Vector{Expression})
    F.G.equations = equations
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in equations]...)'
    F.G.system = System(equations; variables=F.G.variables)
    return nothing
end

mutable struct SphericalDiskPacking
    G::ConstraintSystem
    contacts::Vector{Tuple{Int,Int}}
    inversive_distances::Union{Vector{Int},Vector{Float64}}

    function SphericalDiskPacking(vertices::Vector{Int}, contacts::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}}, inversive_distances::Union{Vector{Int},Vector{Float64}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=Vector{Int}([]), tolerance::Float64=1e-8)
        length(contacts)==length(inversive_distances) || throw(error(("The length of the inversive distances does not match the length of the vertices or the dimensionality of the realization.")))
        all(v->v in vertices, pinned_vertices) || throw(error("Some of the pinned_vertices are not contained in vertices."))
        dimension = size(realization)[1]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension==3 || raise(error("The dimension for SphericalDiskPacking must be 3."))
        all(i->isapprox(minkowski_scalar_product(realization[:,contacts[i][1]], realization[:,contacts[i][2]])/sqrt(minkowski_scalar_product(realization[:,contacts[i][1]], realization[:,contacts[i][1]]) * minkowski_scalar_product(realization[:,contacts[i][2]], realization[:,contacts[i][2]])), inversive_distances[i], atol=tolerance), 1:length(contacts)) || throw(error("The Minkowski distances do not match the given realization."))

        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        inversive_distance_equation = [minkowski_scalar_product(xs[:,contacts[i][1]], xs[:,contacts[i][2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(xs[:,contacts[i][1]], xs[:,contacts[i][1]]) * minkowski_scalar_product(xs[:,contacts[i][2]], xs[:,contacts[i][2]]) for i in 1:length(contacts)]
        inversive_distance_equation = filter(eq->eq!=0, inversive_distance_equation)
        G = ConstraintSystem(vertices, variables, inversive_distance_equation, realization, xs; pinned_vertices=pinned_vertices)
        new(G, contacts, inversive_distances)
    end

    function SphericalDiskPacking(contacts::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}}, inversive_distances::Union{Vector{Int},Vector{Float64}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([bar[1] for bar in contacts], [bar[2] for bar in contacts]))))
        SphericalDiskPacking(vertices, contacts, inversive_distances, realization; pinned_vertices=pinned_vertices)
    end

    function SphericalDiskPacking(contacts::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}}, realization::Union{Matrix{Int},Matrix{Float64}}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        inversive_distances = [minkowski_scalar_product(realization[:,contact[1]], realization[:,contact[2]])/sqrt(minkowski_scalar_product(realization[:,contact[1]], realization[:,contact[1]]) * minkowski_scalar_product(realization[:,contact[2]], realization[:,contact[2]])) for contact in contacts]
        SphericalDiskPacking(contacts, inversive_distances, realization; pinned_vertices=pinned_vertices)
    end

    minkowski_scalar_product(e1,e2) = e1'*e2-1
end

function equations!(F::SphericalDiskPacking, equations::Vector{Expression})
    F.G.equations = equations
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in equations]...)'
    F.G.system = System(equations; variables=F.G.variables)
    return nothing
end


mutable struct VolumeHypergraph
    G::ConstraintSystem
    volumes::Vector{Vector{Int}}

    function VolumeHypergraph(vertices::Vector{Int}, volumes::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        dimension = size(realization)[1]
        all(t->length(t)==dimension+1, volumes) && all(facet->all(v->v in vertices, facet), volumes) || throw(error("The volumes don't have the correct format."))
        volumes = [Vector(facet) for facet in volumes]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw(error("The realization does not have the correct format."))
        dimension>=1 || raise(error("The dimension is not an integer bigger than 0."))
        @var x[1:dimension, 1:length(vertices)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices)))]...)
        facet_equations = [det(vcat([1. for _ in 1:dimension+1]', hcat([x[:,v] for v in facet]...))) - det(vcat([1. for _ in 1:dimension+1]', hcat([realization[:,v] for v in facet]...))) for facet in volumes]
        facet_equations = filter(eq->eq!=0, facet_equations)
        G = ConstraintSystem(vertices,variables, facet_equations, realization, x)
        new(G, volumes)
    end

    function VolumeHypergraph(volumes::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Union{Matrix{Int},Matrix{Float64}})
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in volumes]...))))
        VolumeHypergraph(vertices, volumes, realization)
    end
end

function equations!(F::VolumeHypergraph, equations::Vector{Expression})
    F.G.equations = equations
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in equations]...)'
    F.G.system = System(equations; variables=F.G.variables)
    return nothing
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

function equations!(F::Polytope, equations::Vector{Expression})
    F.G.equations = equations
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in equations]...)'
    F.G.system = System(equations; variables=F.G.variables)
    return nothing
end


function to_Array(G::ConstraintSystem, p::Union{Matrix{Int},Matrix{Float64}})
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(G.realization)[1], 1:size(G.realization)[2])) if !(j in G.pinned_vertices)]...)
end

function to_Array(F::Union{Framework,VolumeHypergraph,Polytope,SpherePacking,SphericalDiskPacking}, p::Union{Matrix{Int},Matrix{Float64}})
    return to_Array(F.G, p)
end



function to_Matrix(G::ConstraintSystem, q::Union{Vector{Float64}, Vector{Int}})
    counts = 1
    point = Matrix{Float64}(Base.copy(G.realization))

    for (i, v) in enumerate(G.vertices)
        if v in G.pinned_vertices
            point[:,i] = G.realization[:,i]
            continue
        end
        for j in 1:size(point)[1]
            point[j,i] = q[counts]
            counts += 1
        end
    end
    return point
end


function to_Matrix(F::Union{Framework,VolumeHypergraph,Polytope,SpherePacking,SphericalDiskPacking}, q::Union{Vector{Float64}, Vector{Int}})
    return to_Matrix(F.G, q)
end


function plot(F, filename::String; kwargs...)
    if typeof(F)==Framework
        return plot_framework(F, filename; kwargs...)
    elseif typeof(F)==VolumeHypergraph
        return plot_hypergraph(F, filename; kwargs...)
    elseif typeof(F)==Polytope
        return plot_polytope(F, filename; kwargs...)
    elseif typeof(F)==SpherePacking
        return plot_spherepacking(F, filename; kwargs...)
    elseif typeof(F)==SphericalDiskPacking
        return plot_sphericaldiskpacking(F, filename; kwargs...)
    else
        throw(error("The type of 'F' needs to be either Framework, Polytope, DiskPacking or VolumeHypergraph, but is $(typeof(F))"))
    end
end

function plot_framework(F::Framework, filename::String; padding::Float64=0.15, vertex_size=60, line_width=12, edge_color=:steelblue, markercolor=:red3, pin_point_offset=0.2, vertex_color=:black)
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
    foreach(v->scatter!(ax, [F.G.dimension==2 ? Point2f((allVertices)[v]-[pin_point_offset,0]) : Point3f((allVertices)[v]-[pin_point_offset,0,0])]; markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[:lightgrey]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end

function plot_spherepacking(F::SpherePacking, filename::String; padding::Float64=0.15, disk_strokewidth=8.5, vertex_labels::Bool=true, disk_color=:steelblue, markersize=75, markercolor=:red3, line_width=7, dualgraph_color=:grey80, n_circle_segments::Int=50)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization
    allVertices = F.G. dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    limit_vertices = allVertices

    if F.G.dimension==2
        ax = Axis(fig[1,1])
    elseif F.G.dimension==3
        ax = Axis3(fig[1,1])
        zlims = [minimum([limit_vertices[i][3] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][3] for i in 1:length(limit_vertices)])]
    else
        throw(error("The dimension must either be 2 or 3!"))
    end

    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=dualgraph_color), F.contacts)
    for index in 1:length(F.G.vertices)
        disk_vertices = [Vector(allVertices[index])+F.radii[index]*Point2f([cos(2*i*pi/n_circle_segments), sin(2*i*pi/n_circle_segments)]) for i in 1:n_circle_segments]
        limit_vertices = vcat(limit_vertices, disk_vertices)
        poly!(ax, [(disk_vertices)[i] for i in 1:n_circle_segments]; color=(disk_color, 0.08))
        lines!(ax, [(disk_vertices)[v] for v in vcat(1:n_circle_segments,1)]; linewidth = disk_strokewidth, color=disk_color)
    end

    xlims = [minimum([limit_vertices[i][1] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][1] for i in 1:length(limit_vertices)])]
    ylims = [minimum([limit_vertices[i][2] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][2] for i in 1:length(limit_vertices)])]
    limits= F.G.dimension==2 ? [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])] : [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]

    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    F.G.dimension==3 && zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    foreach(v->scatter!(ax, [(allVertices)[v]]; markersize=markersize, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    vertex_labels && foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[:black]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end


function plot_sphericaldiskpacking(F::SphericalDiskPacking, filename::String; padding=0.015, sphere_color=:lightgrey, vertex_size=60, disk_strokewidth=9, line_width=6, disk_color=:steelblue, dualgraph_color=(:red3,0.45), vertex_color=:black, vertex_labels::Bool=true, n_circle_segments=50)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization    

    ax = Axis3(fig[1,1], aspect=(1,1,1))
    xlims!(ax,-1.5-padding, 1.5+padding)
    ylims!(ax,-1.5-padding, 1.5+padding)
    zlims!(ax,-1.5-padding, 1.5+padding)
    hidespines!(ax)
    hidedecorations!(ax)
    mesh!(ax, Sphere(Point3f(0), 1f0); transparency=true, color = (sphere_color,0.15))

    planePoints = [Point3f(matrix_coords[:,j]./norm(matrix_coords[:,j])^2) for j in 1:size(matrix_coords)[2]]
    koebePoints =[]
    spherePoints = [Point3f(matrix_coords[:,j]./norm(matrix_coords[:,j])) for j in 1:size(matrix_coords)[2]]
    rotatedPoints=[]
    #foreach(edge->linesegments!(ax, [spherePoints[Int64(edge[1])], spherePoints[Int64(edge[2])]]; linewidth = line_width, color=dualgraph_color), F.contacts)
    for i in 1:length(F.G.vertices)
        rotation_axis = cross([0, 0, 1], spherePoints[i])
        if isapprox(norm(rotation_axis), 0, atol=1e-5)
            angle = acos([0, 0, 1]'* spherePoints[i])
            rotation_matrix = [1 0 0; 0 1 0; 0 0 cos(angle);]
        else
            rotation_axis = rotation_axis ./ norm(rotation_axis)
            angle = acos([0, 0, 1]'* spherePoints[i])
            rotation_matrix = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                                rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                                rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
        end
        radius = sqrt(1-norm(planePoints[i])^2)
        disk_vertices = [Point3f(inv(rotation_matrix)*([radius*cos(2*j*pi/n_circle_segments), radius*sin(2*j*pi/n_circle_segments), norm(planePoints[i])])) for j in 1:n_circle_segments]
        push!(koebePoints, Point3f(inv(rotation_matrix)*([0,0,norm(matrix_coords[:,i])])))
        push!(rotatedPoints, Point3f(inv(rotation_matrix)*[0,0,1]))
        lines!(ax, [(disk_vertices)[v] for v in vcat(1:n_circle_segments,1)]; linewidth = disk_strokewidth, color=disk_color)
    end
    foreach(edge->linesegments!(ax, [koebePoints[Int64(edge[1])], koebePoints[Int64(edge[2])]]; linewidth = line_width, color=dualgraph_color), F.contacts)
    vertex_labels && foreach(i->text!(ax, [(rotatedPoints)[i]], text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[:black]), 1:length(F.G.vertices))
    save("../data/$(filename).png", fig)
    return fig
end
    

function plot_hypergraph(F::VolumeHypergraph, filename::String; padding::Float64=0.15, vertex_size=60, line_width=8, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization    
    if facet_colors==nothing
        facet_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(F.volumes), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
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
    foreach(i->poly!(ax, [(allVertices)[Int64(v)] for v in F.volumes[i]]; color=(facet_colors[i], 0.25)), 1:length(F.volumes))
    foreach(i->lines!(ax, [(allVertices)[Int64(v)] for v in vcat(F.volumes[i], F.volumes[i][1])]; linewidth=line_width, linestyle=:dash, color=facet_colors[i]), 1:length(F.volumes))
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
