
"""
    plot(F, filename)

Plot the geometric constraint system `F`.

This is a wrapper method for plotting of the individual geometric constraint systems. 
It calls one of the following: [`plot_flexes!`](@ref), [`plot_framework`](@ref), [`plot_frameworkonsurface`](@ref), [`plot_hypergraph`](@ref), [`plot_polytope`](@ref), [`plot_spherepacking`](@ref) and [`plot_sphericaldiskpacking`](@ref).
"""
function plot(F::AllTypes, filename::Union{String, Nothing}; kwargs...)
    if typeof(F)==Framework || typeof(F)==AngularFramework
        global fig = plot_framework(F, filename; kwargs...)
    elseif typeof(F)==FrameworkOnSurface
        global fig = plot_frameworkonsurface(F, filename; kwargs...)    
    elseif typeof(F)==VolumeHypergraph
        global fig = plot_hypergraph(F, filename; kwargs...)
    elseif typeof(F)==Polytope || typeof(F)==BodyHinge
        global fig = plot_polytope(F, filename; kwargs...)
    elseif typeof(F)==SpherePacking
        global fig = plot_spherepacking(F, filename; kwargs...)
    elseif typeof(F)==SphericalDiskPacking
        global fig = plot_sphericaldiskpacking(F, filename; kwargs...)
    else
        throw("The type of 'F' needs to be either Framework, AngularFramework, Polytope, BodyHinge, SpherePacking, SphericalDiskPacking, FrameworkOnSurface or VolumeHypergraph, but is $(typeof(F))")
    end
    
    return fig
end

"""
    plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, linewidth, arrowsize)

Add infinitesimal flexes to a plot of a geometric constraint system
"""
function plot_flexes!(ax, F::AllTypes, flex_Real::Int, flex_color, flex_scale, linewidth, arrowsize)
    if F isa Framework
        K_n = Framework([[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices)
    elseif F isa AngularFramework
        K_n = AngularFramework([[i,j,k] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) for k in 1:length(F.G.vertices) if (i<j && j<k) || (i<k && k<j) || (j<i && i<k)], F.G.realization; pinned_vertices=F.G.pinned_vertices)
    elseif F isa FrameworkOnSurface
        K_n = deepcopy(G)
        add_equations!(K_n, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(G.vertices) for j in 1:length(F.G.vertices) if i<j]])
    elseif F isa VolumeHypergraph
        K_n = VolumeHypergraph(collect(powerset(F.G.vertices, F.G.dimension+1, F.G.dimension+1)), F.G.realization)
    elseif F isa Polytope || F isa SpherePacking || F isa BodyHinge
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    elseif  F isa SphericalDiskPacking
        minkowski_scalar_product(e1,e2) = e1'*e2-1
        inversive_distances = [minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[1]]) * minkowski_scalar_product(F.G.realization[:,contact[2]], F.G.realization[:,contact[2]])) for contact in powerset(F.G.vertices, 2, 2)]
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, [minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[1]]) * minkowski_scalar_product(F.G.xs[:,contact[2]], F.G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(F.G.vertices, 2, 2))], F.G.realization, F.G.xs)
    else
        throw("The type must either be 'framework', 'frameworkonsurface', 'diskpacking', 'sphericaldiskpacking', 'hypergraph', 'bodyhinge' or 'polytope', but is '$(type)'.")
    end

    flex = to_Matrix(F,compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), K_n)[:,flex_Real])
    if F.G.dimension==2
        pts = [Point2f(F.G.realization[:,i]) for i in 1:length(F.G.vertices)]
        dirs = [Vec2f(flex[:,i]) for i in 1:length(F.G.vertices)]
        arrows!(ax, pts, dirs; lengthscale=flex_scale, arrowcolor = flex_color, linecolor = flex_color, linewidth=linewidth, arrowsize=arrowsize)
    else
        pts = [Point3f(F.G.realization[:,i]) for i in 1:length(F.G.vertices)]
        dirs = [Vec3f(flex[:,i]) for i in 1:length(F.G.vertices)]
        arrows!(ax, pts, dirs; lengthscale=flex_scale*8, arrowcolor = flex_color, linecolor = flex_color, arrowsize=0.135)
    end
end

"""
    plot_framework(F, filename)

Plot a bar-joint or angular framework.
"""
function plot_framework(F::Union{Framework,AngularFramework}, filename::Union{String, Nothing}; padding::Real=0.15, vertex_size=55, azimuth=π / 10, elevation=pi/10, perspectiveness=0., vertex_labels=true, line_width=12, edge_color=:steelblue, angle_color=:lightgrey, font_color=:lightgrey, angle_size=0.3, markercolor=:red3, pin_point_offset=0.1, vertex_color=:black, plot_flexes=false, flex_Real=1, flex_color=:green3, flex_scale=0.35, arrowsize=40)
    fig = Figure(size=(1000,1000))
    matrix_coords = Base.copy(F.G.realization)
    centroid = sum(matrix_coords[:,i] for i in 1:size(matrix_coords)[2]) ./ size(matrix_coords)[2]
    for i in 1:size(matrix_coords)[2]
        matrix_coords[:,i] = matrix_coords[:,i] - centroid
    end

    if F.G.dimension==2
        ax = Axis(fig[1,1])
    elseif F.G.dimension==3
        ax = Axis3(fig[1,1], aspect = (1, 1, 1), azimuth=azimuth, elevation=elevation, perspectiveness=perspectiveness)
        zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
    else
        throw("The dimension must either be 2 or 3!")
    end
    xlims = [minimum(vcat(matrix_coords[1,:])), maximum(matrix_coords[1,:])]
    ylims = [minimum(vcat(matrix_coords[2,:])), maximum(matrix_coords[2,:])]
    limits= F.G.dimension==2 ? [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])] : [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]

    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    if F.G.dimension==3
        translation = (zlims[1]-limits[1]) - (limits[2]-zlims[2])
        zlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    end
    hidespines!(ax)
    hidedecorations!(ax)

    allVertices = F.G.dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    if typeof(F)==AngularFramework
        L = angle_size*maximum([minimum([norm(matrix_coords[:,bar[1]]-matrix_coords[:,bar[2]]) for bar in F.bars]),0])
        for angle in F.angles
            poly_points=[F.G.dimension==2 ? Point2f(matrix_coords[:,angle[2]]) : Point3f(matrix_coords[:,angle[2]])]
            if norm(matrix_coords[:,angle[2]]-matrix_coords[:,angle[1]])==0 || norm(matrix_coords[:,angle[2]]-matrix_coords[:,angle[3]])==0
                continue
            end
            for t in 0:0.05:1
                curpt = matrix_coords[:,angle[2]] + L*((1-t)*(matrix_coords[:,angle[1]]-matrix_coords[:,angle[2]]) + t*(matrix_coords[:,angle[3]]-matrix_coords[:,angle[2]])) ./ norm((1-t)*(matrix_coords[:,angle[1]]-matrix_coords[:,angle[2]]) + t*(matrix_coords[:,angle[3]]-matrix_coords[:,angle[2]]))
                push!(poly_points, F.G.dimension==2 ? Point2f(curpt) : Point3f(curpt))
            end
            poly!(ax, poly_points, color=(angle_color,0.5), strokewidth=0)
            lines!(ax, poly_points[2:end], color=angle_color, linewidth=line_width/2)
        end
    end
    
    if plot_flexes
        plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, line_width-2, arrowsize)
    end

    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=edge_color), F.bars)
    foreach(v->scatter!(ax, [F.G.dimension==2 ? Point2f((allVertices)[v]-[pin_point_offset,0]) : Point3f((allVertices)[v]-[pin_point_offset,0,0])]; markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    vertex_labels && foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    if !isnothing(filename)
        save("../data/$(filename).png", fig)
    end
    return fig
end


"""
    plot_frameworkonsurface(F, filename)

Plot a bar-joint framework constrained to a surface.
"""
function plot_frameworkonsurface(F::FrameworkOnSurface, filename::Union{String, Nothing}; padding::Real=0.15, azimuth=pi/10, alpha=0.45, elevation=pi/8, perspectiveness=0., vertex_size=55, line_width=10, edge_color=:steelblue, markercolor=:red3, pin_point_offset=0.2, vertex_color=:black, vertex_labels=true, font_color=:lightgrey, surface_color=:grey80, surface_samples=150, plot_flexes=false, flex_Real=1, flex_color=:green3, flex_scale=0.35, arrowsize=40)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization

    if F.G.dimension==3
        ax = Axis3(fig[1,1], aspect = (1, 1, 1), azimuth=azimuth, elevation=elevation, perspectiveness=perspectiveness)
    else
        throw("The dimension must either be 2 or 3!")
    end
    xlims = [minimum(vcat(matrix_coords[1,:])), maximum(matrix_coords[1,:])]
    ylims = [minimum(vcat(matrix_coords[2,:])), maximum(matrix_coords[2,:])]
    zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
    limits= [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]

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

    mesh!(ax, msh; color=(surface_color,alpha), transparency=true)

    if plot_flexes
        plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, line_width-2, arrowsize)
    end

    allVertices = [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=edge_color), F.bars)
    foreach(v->scatter!(ax, [Point3f((allVertices)[v]-[pin_point_offset,0,0])]; markersize=vertex_size, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    vertex_labels && foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    if !isnothing(filename)
        save("../data/$(filename).png", fig)
    end
    return fig
end


"""
    plot_spherepacking(F, filename)

Plot a sphere packing.
"""
function plot_spherepacking(F::SpherePacking, filename::Union{String, Nothing}; padding::Real=0.15, azimuth=pi/10, alpha=0.1, elevation=pi/8, perspectiveness=0., disk_strokewidth=8.5, vertex_labels::Bool=true, font_color=:black, sphere_color=:steelblue, D2_markersize=75, D3_markersize=55, markercolor=:red3, line_width=7, D2_dualgraph_color=:grey80, D3_dualgraph_color=:grey50, n_circle_segments::Int=50, plot_flexes=false, flex_Real=1, flex_color=:green3, flex_scale=0.35, arrowsize=40, kwargs...)
    fig = Figure(size=(1000,1000))
    matrix_coords = Base.copy(F.G.realization)
    centroid = sum(matrix_coords[:,i] for i in 1:size(matrix_coords)[2]) ./ size(matrix_coords)[2]
    for i in 1:size(matrix_coords)[2]
        matrix_coords[:,i] = matrix_coords[:,i] - centroid
    end

    allVertices = F.G. dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    limit_vertices = allVertices

    if F.G.dimension==2
        ax = Axis(fig[1,1])
    elseif F.G.dimension==3
        ax = Axis3(fig[1,1], aspect = (1, 1, 1), azimuth=azimuth, elevation=elevation, perspectiveness=perspectiveness)
    else
        throw("The dimension must either be 2 or 3!")
    end

    if !haskey(kwargs, "dualgraph_color")
        dualgraph_color = F.G.dimension==2 ? D2_dualgraph_color : D3_dualgraph_color
    end

    foreach(edge->linesegments!(ax, [(allVertices)[Int64(edge[1])], (allVertices)[Int64(edge[2])]]; linewidth = line_width, color=dualgraph_color), F.contacts)
    for index in 1:length(F.G.vertices)
        if F.G.dimension==2
            disk_vertices = [Vector(allVertices[index])+F.radii[index]*Point2f([cos(2*i*pi/n_circle_segments), sin(2*i*pi/n_circle_segments)]) for i in 1:n_circle_segments]
            limit_vertices = vcat(limit_vertices, disk_vertices)
            poly!(ax, [(disk_vertices)[i] for i in 1:n_circle_segments]; color=(sphere_color, alpha))
            lines!(ax, [(disk_vertices)[v] for v in vcat(1:n_circle_segments,1)]; linewidth = disk_strokewidth, color=sphere_color)
        else
            mesh!(ax, Sphere(allVertices[index], F.radii[index]); transparency=true, color = (sphere_color,alpha))
            append!(limit_vertices, vcat([[allVertices[index]+F.radii[index]*Vector{Float64}(I[1:3, i]), allVertices[index]-F.radii[index]*Vector{Float64}(I[1:3, i])] for i in 1:3]...))
        end
    end

    xlims = [minimum([limit_vertices[i][1] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][1] for i in 1:length(limit_vertices)])]
    ylims = [minimum([limit_vertices[i][2] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][2] for i in 1:length(limit_vertices)])]
    if F.G.dimension==3
        zlims = [minimum([limit_vertices[i][3] for i in 1:length(limit_vertices)]), maximum([limit_vertices[i][3] for i in 1:length(limit_vertices)])]
    end
    limits= F.G.dimension==2 ? [minimum([xlims[1], ylims[1]]), maximum([xlims[2], ylims[2]])] : [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]

    translation = (xlims[1]-limits[1]) - (limits[2]-xlims[2])
    xlims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    translation = (ylims[1]-limits[1]) - (limits[2]-ylims[2])
    ylims!(ax, limits[1]-padding+0.5*translation, limits[2]+padding+0.5*translation)
    F.G.dimension==3 && zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)

    if plot_flexes
        plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, line_width-2, arrowsize)
    end

    if !haskey(kwargs, "markersize")
        markersize = F.G.dimension==2 ? D2_markersize : D3_markersize
    end
    
    foreach(v->scatter!(ax, [(allVertices)[v]]; markersize=markersize, color=(markercolor, 0.4), marker=:rtriangle), F.G.pinned_vertices)
    vertex_labels && foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    if !isnothing(filename)
        save("../data/$(filename).png", fig)
    end
    return fig
end


"""
    plot_sphericaldiskpacking(F, filename)

Plot a spherical disk packing.
"""
function plot_sphericaldiskpacking(F::SphericalDiskPacking, filename::Union{String, Nothing}; azimuth=pi/10, elevation=pi/8, perspectiveness=0., padding=0.015, alpha=0.15, sphere_color=:lightgrey, font_color=:black, vertex_size=60, disk_strokewidth=9, line_width=6, disk_color=:steelblue, dualgraph_color=(:red3,0.45), vertex_color=:black, vertex_labels::Bool=true, n_circle_segments=50, plot_flexes=false, flex_Real=1, flex_color=:green3, flex_scale=0.35, arrowsize=40)
    fig = Figure(size=(1000,1000))
    matrix_coords = F.G.realization    

    ax = Axis3(fig[1,1], aspect=(1,1,1), azimuth=azimuth, elevation=elevation, perspectiveness=perspectiveness)
    xlims!(ax,-1.5-padding, 1.5+padding)
    ylims!(ax,-1.5-padding, 1.5+padding)
    zlims!(ax,-1.5-padding, 1.5+padding)
    hidespines!(ax)
    hidedecorations!(ax)
    mesh!(ax, Sphere(Point3f(0), 1f0); transparency=true, color = (sphere_color,alpha))

    if plot_flexes
        plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, line_width-2, arrowsize)
    end

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
    vertex_labels && foreach(i->text!(ax, [(rotatedPoints)[i]], text=["$(F.G.vertices[i])"], fontsize=32, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    if !isnothing(filename)
        save("../data/$(filename).png", fig)
    end
    return fig
end
    

"""
    plot_hypergraph(F, filename)

Plot a volume-constrained hypergraph.
"""
function plot_hypergraph(F::VolumeHypergraph, filename::Union{String, Nothing}; padding::Real=0.15, alpha=0.25, azimuth=pi/10, elevation=pi/8, perspectiveness=0., vertex_size=60, line_width=8, facet_colors=nothing, vertex_color=:black, font_color=:lightgrey, vertex_labels::Bool=true, plot_flexes=false, flex_Real=1, flex_color=:green3, flex_scale=0.35, arrowsize=40)
    fig = Figure(size=(1000,1000))
    matrix_coords = Base.copy(F.G.realization)
    centroid = sum(matrix_coords[:,i] for i in 1:size(matrix_coords)[2]) ./ size(matrix_coords)[2]
    for i in 1:size(matrix_coords)[2]
        matrix_coords[:,i] = matrix_coords[:,i] - centroid
    end

    if facet_colors==nothing
        facet_colors = map(col -> (red(col), green(col), blue(col)), distinguishable_colors(length(F.volumes), [RGB(1,1,1), RGB(0,0,0)], dropseed=true, lchoices = range(20, stop=70, length=15), hchoices = range(0, stop=360, length=30)))
    end

    if F.G.dimension==2
        ax = Axis(fig[1,1])
    elseif F.G.dimension==3
        ax = Axis3(fig[1,1], aspect = (1, 1, 1), azimuth=azimuth, elevation=elevation, perspectiveness=perspectiveness)
        zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
    else
        throw("The dimension must either be 2 or 3!")
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

    if plot_flexes
        plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, line_width-2, arrowsize)
    end

    allVertices = F.G.dimension==2 ? [Point2f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]] : [Point3f(matrix_coords[:,j]) for j in 1:size(matrix_coords)[2]]
    foreach(i->poly!(ax, [(allVertices)[Int64(v)] for v in F.volumes[i]]; color=(facet_colors[i], alpha)), 1:length(F.volumes))
    foreach(i->lines!(ax, [(allVertices)[Int64(v)] for v in vcat(F.volumes[i], F.volumes[i][1])]; linewidth=line_width, linestyle=:dash, color=facet_colors[i]), 1:length(F.volumes))
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:length(F.G.vertices))
    foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:length(F.G.vertices))
    if !isnothing(filename)
        save("../data/$(filename).png", fig)
    end
    return fig
end


"""
    plot_polytope(F, filename)

Plot a polytope.
"""
function plot_polytope(F::Union{Polytope,BodyHinge}, filename::Union{String, Nothing}; padding=0.1, vertex_size=60, alpha=0.55, line_width=12, edge_color=:steelblue, perspectiveness=0., azimuth=π/10, elevation=pi/8, facet_color=:grey98, font_color=:lightgrey, vertex_color=:black, vertex_labels::Bool=true, plot_flexes=false, flex_Real=1, flex_color=:green3, flex_scale=0.35, arrowsize=40)
    fig = Figure(size=(1000,1000))
    ax = Axis3(fig[1,1], aspect = (1, 1, 1), azimuth=azimuth, elevation=elevation, perspectiveness=perspectiveness)

    matrix_coords = F isa Polytope ? Base.copy(F.G.realization)[:,1:(size(F.G.realization)[2]-length(F.facets))] : Base.copy(F.G.realization)
    centroid = F isa Polytope ? sum([matrix_coords[:,i] for i in 1:(size(F.G.realization)[2]-length(F.facets))]) ./ (size(F.G.realization)[2]-length(F.facets)) : sum([matrix_coords[:,i] for i in 1:(size(F.G.realization)[2])]) ./ (size(F.G.realization)[2])
    for i in 1:(size(F.G.realization)[2]-length(F.facets))
        matrix_coords[:,i] = matrix_coords[:,i] - centroid
    end

    xlims = [minimum(vcat(matrix_coords[1,:])), maximum(matrix_coords[1,:])]
    ylims = [minimum(vcat(matrix_coords[2,:])), maximum(matrix_coords[2,:])]
    zlims = [minimum(vcat(matrix_coords[3,:])), maximum(matrix_coords[3,:])]
    limits = [minimum([xlims[1], ylims[1], zlims[1]]), maximum([xlims[2], ylims[2], zlims[2]])]
    xlims!(ax, limits[1]-padding, limits[2]+padding)
    ylims!(ax, limits[1]-padding, limits[2]+padding)
    zlims!(ax, limits[1]-padding, limits[2]+padding)
    hidespines!(ax)
    hidedecorations!(ax)
    allVertices = F isa Polytope ? [Point3f(matrix_coords[:,j]) for j in 1:(size(F.G.realization)[2]-length(F.facets))] : [Point3f(matrix_coords[:,j]) for j in 1:(size(F.G.realization)[2])]

    if typeof(F) <: Polytope
        P = Polyhedra.polyhedron(Polyhedra.vrep([matrix_coords[:,j] for j in 1:(size(F.G.realization)[2]-length(F.facets))]))
        m = Polyhedra.Mesh(P)
        mesh!(ax, m; color=(facet_color,alpha), shading=NoShading, transparency=true)
    else
        for face in F.facets
            P = Polyhedra.polyhedron(Polyhedra.vrep([matrix_coords[:,j] for j in face]))
            m = Polyhedra.Mesh(P)
            mesh!(ax, m; color=(facet_color,alpha), shading=NoShading, transparency=true)
        end
    end

    if plot_flexes
        plot_flexes!(ax, F, flex_Real, flex_color, flex_scale, line_width-2, arrowsize)
    end

    foreach(i->linesegments!(ax, [(allVertices)[Int64(F.edges[i][1])], (allVertices)[Int64(F.edges[i][2])]]; linewidth=line_width, color=edge_color), 1:length(F.edges))
    foreach(i->scatter!(ax, [(allVertices)[i]]; markersize = vertex_size, color=vertex_color), 1:(size(F.G.realization)[2]-length(F.facets)))
    vertex_labels && foreach(i->text!(ax, [(allVertices)[i]], text=["$(F.G.vertices[i])"], fontsize=28, font=:bold, align = (:center, :center), color=[font_color]), 1:(size(F.G.realization)[2]-length(F.facets)))
    if !isnothing(filename)
        save("../data/$(filename).png", fig)
    end
    return fig
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
function animate2D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::String; recompute_deformation_samples::Bool=true, fixed_vertices::Tuple{Int,Int}=(1,2), fixed_direction::Vector{<:Real}=[1.,0], framerate::Int=25, step::Int=1, padding::Real=0.15, markercolor=:red3, pin_point_offset=0.1, vertex_size::Real=55, line_width::Real=12, angle_color=:lightgrey, font_color=:lightgrey, angle_size=0.3, edge_color=:steelblue, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
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
function animate3D_framework(D::DeformationPath, F::Union{Framework,AngularFramework}, filename::String; recompute_deformation_samples::Bool=true, fixed_vertices::Union{Tuple{Int,Int}, Tuple{Int,Int,Int}}=(1,2), fixed_direction=[1.,0,0], framerate::Int=25, animate_rotation=false, azimuth = π / 4, elevation=pi/8, perspectiveness=0., rotation_frames = 240, markercolor=:red3, pin_point_offset=0.05, step::Int=1, padding::Real=0.15, vertex_size::Real=55, vertex_labels=false, font_color=:lightgrey, line_width::Real=12, angle_color=:lightgrey, angle_size=0.3, edge_color=:steelblue, vertex_color=:black, filetype::String="gif")
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
function animate3D_frameworkonsurface(D::DeformationPath, F::FrameworkOnSurface, filename::String; alpha=0.45, framerate::Int=25, animate_rotation=false, azimuth = pi/4, elevation=pi/8, perspectiveness=0., rotation_frames = 480, markercolor=:red3, pin_point_offset=0.05, step::Int=1, padding::Real=0.15, vertex_size::Real=55, line_width::Real=10, edge_color=:steelblue, vertex_labels=true, font_color=:lightgrey, vertex_color=:black, filetype::String="gif", surface_color=:grey80, surface_samples=150)
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
function animate2D_hypergraph(D::DeformationPath, F::VolumeHypergraph, filename::String; alpha=0.2, recompute_deformation_samples::Bool=true, target_stretch::Real=1., fixed_triangle::Union{Tuple{Int,Int,Int},Vector{Int},Nothing}=nothing, font_color=:black, skip_stretch::Bool=true, tip_value::Real=0.5, framerate::Int=25, step::Int=1, padding::Real=0.15, vertex_size::Real=42, line_width::Real=6, facet_colors=nothing, vertex_color=:black, vertex_labels::Bool=true, filetype::String="gif")
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
function animate3D_polytope(D::DeformationPath, F::Union{Polytope,BodyHinge}, filename::String; renderEntirePolytope::Bool=true, recompute_deformation_samples::Bool=true, fixed_vertices::Union{Tuple{Int,Int}, Tuple{Int,Int,Int}}=(1,2), alpha=0.6, font_color=:lightgrey, facet_color=:grey98, framerate::Int=25, animate_rotation=false, azimuth = π / 10, elevation=pi/8, perspectiveness=0., rotation_frames = 240, step::Int=1, padding::Real=0.1, vertex_size::Real=45, line_width::Real=8.5, edge_color=:steelblue, special_edge=nothing, special_edge_color=:red3, vertex_color=:black, vertex_labels::Bool=false, filetype::String="gif")
    fig = Figure(size=(1000,1000))
    matrix_coords = [to_Matrix(F, D.motion_samples[i]) for i in 1:length(D.motion_samples)]
    F isa BodyHinge || (fixed_vertices[1] in 1:(size(F.G.realization)[2]-length(F.facets)) && fixed_vertices[2] in 1:(size(F.G.realization)[2]-length(F.facets)) && (length(fixed_vertices)==2 || fixed_vertices[3] in 1:(size(F.G.realization)[2]-length(F.facets)))) || throw("The elements of `fixed_vertices` are not vertices of the underlying graph.")
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
        pointys = F isa Polytope ? matrix_coords[$time][:,1:(size(F.G.realization)[2]-length(F.facets))] : matrix_coords[$time][:,1:(size(F.G.realization)[2])]
        [Point3f(pointys[:,j]) for j in 1:size(pointys)[2]]
    end

    allVertices_asLists = @lift begin
        pointys = F isa Polytope ? matrix_coords[$time][:,1:(size(F.G.realization)[2]-length(F.facets))] : matrix_coords[$time][:,1:(size(F.G.realization)[2])]
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
function animate2D_diskpacking(D::DeformationPath, F::SpherePacking, filename::String; alpha=0.08, framerate::Int=25, step::Int=1, padding::Real=0.15, vertex_labels=true, disk_strokewidth::Real=8.5, line_width::Real=7, font_color=:black, sphere_color=:steelblue, markersize::Real=75, markercolor=:red3, dualgraph_color=:grey80, n_circle_segments::Int=50, filetype::String="gif")
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
function animate3D_spherepacking(D::DeformationPath, F::SpherePacking, filename::String; alpha=0.2, framerate::Int=25, step::Int=1, padding::Real=0.1, vertex_labels=true, font_color=:black, line_width::Real=7, sphere_color=:steelblue, markersize::Real=55, markercolor=:red3, dualgraph_color=:grey50, n_circle_segments::Int=50, filetype::String="gif")
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
function project_deformation_random(D::Union{DeformationPath,Vector{DeformationPath}}, projected_dimension::Int; line_width::Real=8, edge_colors=[:green3], markersize::Real=45, markercolor=:steelblue, draw_start::Bool=true)
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