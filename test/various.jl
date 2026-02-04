@testset "24cell" begin
    realization = Matrix([1 1 0 0; 1 -1 0 0; -1 1 0 0; -1 -1 0 0; 1 0 1 0; 1 0 -1 0; -1 0 1 0; -1 0 -1 0; 1 0 0 1; 1 0 0 -1; -1 0 0 1; -1 0 0 -1; 0 1 1 0; 0 1 -1 0; 0 -1 1 0; 0 -1 -1 0; 0 1 0 1; 0 1 0 -1; 0 -1 0 1; 0 -1 0 -1; 0 0 1 1; 0 0 1 -1; 0 0 -1 1; 0 0 -1 -1;]')
    #=edges, facets = [], []
    for i in 1:size(realization)[2]
        new_edges = [(i,j) for j in i+1:size(realization)[2] if 2==sum( (realization[:,i]-realization[:,j]) .^ 2 )]
        append!(edges, new_edges)
    end
    for i in 1:size(realization)[2]
        _edges = edges[findall(edge->i in collect(edge), edges)]
        edge_subsets = combinations(_edges, 4)
        for subset in edge_subsets
            vs = sort([edge[1]==i ? edge[2] : edge[1] for edge in subset])
            _incidences_at_1 = Set([edge[1]==vs[1] ? edge[2] : edge[1] for edge in edges if vs[1] in collect(edge) && !(i in edge)])
            _incidences_at_2 = Set([edge[1]==vs[2] ? edge[2] : edge[1] for edge in edges if vs[2] in collect(edge) && !(i in edge)])
            _incidences_at_3 = Set([edge[1]==vs[3] ? edge[2] : edge[1] for edge in edges if vs[3] in collect(edge) && !(i in edge)])
            _incidences_at_4 = Set([edge[1]==vs[4] ? edge[2] : edge[1] for edge in edges if vs[4] in collect(edge) && !(i in edge)])
            inters = collect(intersect(intersect(_incidences_at_1,_incidences_at_2), intersect(_incidences_at_3,_incidences_at_4)))
            if !isempty(inters) && !(Tuple(sort([i,vs[1],vs[2],vs[3],vs[4],inters[1]])) in facets)
                push!(facets, Tuple(sort([i,vs[1],vs[2],vs[3],vs[4],inters[1]])))
            end
        end
    end=#
    F = FacetPolytope([(1, 2, 5, 6, 9, 10), (1, 5, 9, 13, 17, 21), (1, 5, 10, 13, 18, 22), (1, 6, 9, 14, 17, 23), (1, 6, 10, 14, 18, 24), (1, 3, 13, 14, 17, 18), (2, 5, 9, 15, 19, 21), (2, 5, 10, 15, 20, 22), (2, 6, 9, 16, 19, 23), (2, 6, 10, 16, 20, 24), (2, 4, 15, 16, 19, 20), (3, 4, 7, 8, 11, 12), (3, 7, 11, 13, 17, 21), (3, 7, 12, 13, 18, 22), (3, 8, 11, 14, 17, 23), (3, 8, 12, 14, 18, 24), (4, 7, 11, 15, 19, 21), (4, 7, 12, 15, 20, 22), (4, 8, 11, 16, 19, 23), (4, 8, 12, 16, 20, 24), (5, 7, 13, 15, 21, 22), (6, 8, 14, 16, 23, 24), (9, 11, 17, 19, 21, 23), (10, 12, 18, 20, 22, 24)], realization)
    inf_flexes = compute_inf_flexes(F.G, to_Array(F.G, F.G.realization))
    @test size(inf_flexes)[2] == 52
    stresses = compute_equilibrium_stresses(F.G, to_Array(F.G, F.G.realization))
    @test size(stresses)[2] == 4
end


@testset "sphericaldiskpacking" begin
    F = SphericalDiskPacking([(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(3,5),(4,5),(2,6),(3,6),(4,6),(5,6)], Matrix([sqrt(2) 0 0; 0 sqrt(2) 0; 0 0 sqrt(2); 0 -sqrt(2) 0; 0 0 -sqrt(2); -sqrt(2) 0 0]'); pinned_vertices=[1])
    plot(F,"sphericaldiskpacking"; disk_color=teal, sphere_color=soft_teal, dualgraph_color=(coral,0.75))
    D = DeformationPath(F, [1], 100; step_size=0.01, show_progress=false)
    if is_no_ci
        animate(D,F; filetype="mp4")
    end
end


@testset "bodybar_prism" begin
    F = BodyBar([[1,4],[2,5],[3,6]], [[1,2,3],[4,5,6]], Matrix([cos(0*2*pi/3) sin(0*2*pi/3) 0; cos(2*pi/3) sin(2*pi/3) 0; cos(2*2*pi/3) sin(2*2*pi/3) 0; cos(0*2*pi/3+pi/3) sin(0*2*pi/3+pi/3) 1.3; cos(2*pi/3+pi/3) sin(2*pi/3+pi/3) 1.3; cos(2*2*pi/3+pi/3) sin(2*2*pi/3+pi/3) 1.3;]'); pinned_vertices=[1,2,3])
    plot(F; edge_color=logocolors.blue, special_edges=[[1,4],[2,5],[3,6]], special_edge_color=logocolors.purple, vertex_color=logocolors.red, facet_color=logocolors.green, flex_color=logocolors.purple, plot_flexes=false, azimuth=0.185*pi, flex_scale=0.1, elevation=pi/12, alpha=0.5, line_width=20, vertex_size=70, padding=0.5)
end


@testset "bricard_octahedron" begin
    F = BodyHinge([[1,2,4], [1,2,3], [6,2,4], [6,2,3], [1,3,5], [6,3,5], [1,3,5], [6,3,5]], Matrix([0 0 -1; -1 -1 0; 1 -1 0; 1 1 0; -1 1 0; 0 0 1]'))
    plot(F; edge_color=teal, vertex_color=teal, facet_color=soft_teal, alpha=0.45, azimuth=0.15*pi, elevation=pi/13)
    D = DeformationPath(F, [], 400; step_size=0.025, show_progress=false)

    if is_no_ci
        #=for i in 25:25:length(D.motion_samples)
            _F = BodyHinge([[1,2,4], [1,2,3], [6,2,4], [6,2,3], [1,3,5], [6,3,5], [1,3,5], [6,3,5]], D.motion_matrices[i])
            plot(_F, "bricard_octahedron$i"; edge_color=teal, vertex_color=teal, facet_color=soft_teal, alpha=0.45, azimuth=0.15*pi, elevation=pi/13)
        end=#
        animate(D,F,"bricard_octahedron"; filetype="mp4", edge_color=teal, vertex_color=teal, facet_color=soft_teal, alpha=0.45, azimuth=0.15*pi, elevation=pi/13, padding=0.01)
    end
end

@testset "bodybar_cube" begin
    F = BodyBar([[1,5],[2,6],[3,7],[4,8]], [[1,2,3,4],[5,6,7,8]], Matrix([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1;]'))
    plot(F; edge_color=teal, vertex_color=teal, facet_color=soft_teal, alpha=0.5)
    D = DeformationPath(F, [], 100; step_size=0.01)
    animate(D,F; filetype="mp4", edge_color=teal, vertex_color=teal, facet_color=soft_teal, alpha=0.5)
end


@testset "bodyhinge_pyramid" begin
    F = BodyHinge([[1,2,3],[1,3,4],[1,4,5],[1,5,6],[1,6,2]], Matrix([0 0 1; cos(2*pi/5) sin(2*pi/5) 0; cos(4*pi/5) sin(4*pi/5) 0; cos(6*pi/5) sin(6*pi/5) 0; cos(8*pi/5) sin(8*pi/5) 0; cos(10*pi/5) sin(10*pi/5) 0;]'))
    plot(F; edge_color=teal, vertex_color=teal, facet_color=soft_teal, alpha=0.5)
    D = DeformationPath(F, [], 100; step_size=0.0075)
    animate(D,F; filetype="mp4", animate_rotation=true, rotation_frames=450, padding=0.01)
end

@testset "bodyhinge" begin
    F = BodyHinge([[1,2,3,4],[3,4,5,6]], Matrix([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 1 1 1]'))
    plot(F)
    D = DeformationPath(F, [], 100; step_size=0.025)
    animate(D,F; filetype="mp4")
end


@testset "thales" begin
    F = AngularFramework([[1,3,2]], Matrix([-1 0; 1 0; -sqrt(1/2) sqrt(1/2);]'); pinned_vertices=[1,2])
    plot(F; padding=0.1, pin_point_offset=0.075, edge_color=teal, vertex_size=75, fontsize=36, line_width=14)
    D = DeformationPath(F, [1], 200; step_size=0.025)
    if is_no_ci
        animate(D,F; padding=0.075, pin_point_offset=0.075, filetype="mp4")
    end
end


@testset "squareonhyperboloid" begin
    F = FrameworkOnSurface([[1,2],[2,3],[3,4],[1,4]], Matrix([-sqrt(1/2) -sqrt(1/2) -1; -1 0 0; 0 1 0; sqrt(1/2) sqrt(1/2) 1]'), x->x[1]^2+x[2]^2-x[3]^2-1)
    plot(F; vertex_labels=false, azimuth=3*pi/4, padding=0.35)
    D = DeformationPath(F, [1,1], 200; step_size=0.035, show_progress=false)
    animate(D,F; vertex_labels=false, animate_rotation=true, filetype="mp4", edge_color=teal)
end