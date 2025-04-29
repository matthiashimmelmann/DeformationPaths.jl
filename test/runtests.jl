import DeformationPaths:    Framework,
                            AngularFramework,
                            DeformationPath,
                            Polytope,
                            animate, 
                            VolumeHypergraph,
                            SpherePacking,
                            plot,
                            SphericalDiskPacking,
                            to_Array,
                            equations!,
                            to_Matrix,
                            newton_correct,
                            realization!,
                            project_deformation_random,
                            FrameworkOnSurface,
                            is_rigid,
                            is_inf_rigid
using Test
using HomotopyContinuation
using LinearAlgebra



@testset "is_rigid_test" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[1,5],[3,5],[4,5]], Matrix([0. 0; 1 0; 2 0; 1 1; 1 2]'))
    @test !is_inf_rigid(F)
    @test is_rigid(F)
end


@testset "square" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 1 1; 0 1]'))
    plot(F,"square")
    D = DeformationPath(F, [1], 200; step_size=0.025)
    animate(D,F,"square_motion")
end


@testset "prestress_stable_component" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[1,5],[3,5],[4,5],[1,6]], Matrix([0. 0; 1 0; 2 0; 1 1; 1 2; 0 -1]'))
    plot(F,"prestress_stable_component")
    D = DeformationPath(F, [1,1], 200; step_size=0.025)
    animate(D,F,"prestress_stable_component_motion")
end


@testset "K_4" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[2,4],[1,3],[1,5]], Matrix([0. 0; 1 0; 1 1; 0 1; 0 -1]'))
    plot(F,"K_4")
    D = DeformationPath(F, [1], 200; step_size=0.025)
    animate(D,F,"K_4_motion")
end


@testset "twoprism" begin
    F = Framework([[1,2],[2,3],[3,1],[1,4],[2,5],[3,6],[4,5],[5,6],[6,4]], Matrix([0 0; 2 0; 1 1; 0 2; 2 2; 1 3]'))
    plot(F,"twoprism")
    D = DeformationPath(F, [1], 220; step_size=0.025)
    animate(D,F,"twoprism_motion")
end


@testset "completebipartite" begin
    F = Framework([[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6]], Matrix([0 0; 0 1; 1 -1; 1 0; 1 1; 1 2]'))
    plot(F,"completebipartite")
    D = DeformationPath(F, [1], 500; step_size=0.025)
    animate(D,F,"completebipartite_motion")
end


@testset "coned_cube" begin
    F = Framework(vcat([[1,2],[2,3],[3,4],[1,4],[1,5],[2,6],[3,7],[4,8],[5,6],[6,7],[7,8],[5,8]],[[i,9] for i in 1:8]), Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1; 0 0 sqrt(2)]'))
    plot(F,"coned_cube")
    D = DeformationPath(F, [0.5,0.5], 500; step_size=0.02)
    animate(D,F,"coned_cube_motion")
end


@testset "double_watt" begin
    F = Framework([[1,2],[2,3],[2,4],[3,9],[3,4],[3,5],[4,5],[5,6],[6,7],[7,8],[7,9],[8,9],[8,10],[9,10],[10,11]], Matrix([0 0; 1 0; 2 1; 1 2; 3 2; 4 2; 5 2; 7 2; 6 1; 7 0; 8 0;]'); pinned_vertices=[1,6,11])
    plot(F,"double_watt"; padding=0.35, pin_point_offset=0.2)
    D = DeformationPath(F, [0.5,0.5], 500; step_size=0.05)
    animate(D,F,"double_watt_motion"; padding=0.35, fixed_vertices=(1,6), fixed_direction=[4,2], pin_point_offset=0.2)
end


@testset "thales" begin
    F = AngularFramework([[1,3,2]], Matrix([-1 0; 1 0; -sqrt(1/2) sqrt(1/2);]'); pinned_vertices=[1,2])
    plot(F,"thales"; padding=0.1, pin_point_offset=0.075)
    D = DeformationPath(F, [1], 250; step_size=0.025)
    animate(D,F,"thales_motion"; padding=0.075, pin_point_offset=0.075)
end


@testset "squareonhyperboloid" begin
    F = FrameworkOnSurface([[1,2],[2,3],[3,4],[1,4]], Matrix([-sqrt(1/2) -sqrt(1/2) -1; -1 0 0; 0 1 0; sqrt(1/2) sqrt(1/2) 1]'), x->x[1]^2+x[2]^2-x[3]^2-1)
    plot(F,"squareonhyperboloid")
    D = DeformationPath(F, [1,1], 350; step_size=0.035)
    animate(D,F,"squareonhyperboloid_motion"; animate_rotation=true, filetype="mp4")
end


@testset "sphericaldiskpacking" begin
    F = SphericalDiskPacking([(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(3,5),(4,5),(2,6),(3,6),(4,6),(5,6)], Matrix([sqrt(2) 0 0; 0 sqrt(2) 0; 0 0 sqrt(2); 0 -sqrt(2) 0; 0 0 -sqrt(2); -sqrt(2) 0 0]'); pinned_vertices=[1])
    plot(F,"sphericaldiskpacking")
    D = DeformationPath(F, [1], 250; step_size=0.01)
    animate(D,F,"sphericaldiskpacking_motion")
end


@testset "diskpacking" begin
    F = SpherePacking([1.,1.,1.,1.,1.], Matrix([0 0; 2 0; 3 sqrt(3); 4 0; 6 0]'); pinned_vertices=[1,5])
    plot(F,"diskpacking")
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0; 1.75 -sqrt(2^2-(1.75)^2); 3.5 0; 4.5 sqrt(3)]'); pinned_vertices=[1])
    D = DeformationPath(F, [1,1], 250; step_size=0.025)
    animate(D,F,"diskpacking_motion")
end


@testset "squarediskpacking" begin
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0; 2 0; 0 2; 2 2]'); pinned_vertices=[1])
    plot(F,"squarediskpacking")
    D = DeformationPath(F, [1], 250; step_size=0.01)
    animate(D,F,"squarediskpacking_motion")
end


@testset "spherepacking" begin
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0 0; 2 0 0; 0 2 0; 0 0 2]'), pinned_vertices = [1,2])
    plot(F,"spherepacking")
    D = DeformationPath(F, [1,1,1], 500; step_size=0.04)
    animate(D,F,"spherepacking_motion")
end


@testset "cube" begin
    F = Polytope([[1,2,3,4],[5,6,7,8],[1,2,5,6],[2,3,6,7],[3,4,7,8],[1,4,5,8]], Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1]'))
    plot(F,"cube")
    D = DeformationPath(F, [1,1,1], 100; step_size=0.025)
    animate(D,F,"cube_motion")
end


@testset "cuboctahedron" begin
    F = Polytope([[1,5,9],[1,5,3,7],[1,7,11],[1,9,2,11],[2,9,6],[2,11,8],[3,5,10],[3,7,12],[3,10,4,12],[4,10,6],[4,12,8],[6,4,8,2],[5,9,6,10],[7,11,8,12]], Matrix([1 1 0; -1 1 0; 1 -1 0; -1 -1 0; 1 0 1; -1 0 1; 1 0 -1; -1 0 -1; 0 1 1; 0 -1 1; 0 1 -1; 0 -1 -1;]'))
    plot(F,"cuboctahedron")
    D = DeformationPath(F, [], 200; step_size=0.01, newton_tol=1e-2)
    animate(D,F,"cuboctahedron_motion")
end


@testset "two_triangles" begin
    F = VolumeHypergraph([[1,2,3],[2,3,4]], Matrix([0 0; 1 0; 0 1; 1 1]'))
    plot(F,"two_triangles")
    D = DeformationPath(F, [1], 100; step_size=0.01)
    animate(D, F,"two_triangles_motion"; fixed_triangle=(1,2,3),tip_value=0,skip_stretch=false)
end


@testset "two_triangles" begin
    F = VolumeHypergraph([[1,2,3],[2,3,4]], Matrix([0 0; 1 0; 0 1; 1 1]'))
    plot(F,"two_triangles")
    D = DeformationPath(F, [1], 100; step_size=0.01)
    animate(D, F,"two_triangles_motion"; fixed_triangle=(1,2,3),tip_value=0,skip_stretch=false)
end


@testset "octehedral_decomposition" begin
    F = VolumeHypergraph([[1,3,6],[1,2,5],[2,3,4],[1,5,6],[6,4,5]], Matrix([0 0; 3 0; 0 3; 1 1; 1 0.5; 0.5 1]'))
    plot(F,"octahedral_decomposition")
    D = DeformationPath(F, [0.333, 1], 350; step_size=0.002)
    animate(D, F,"octahedral_decomposition_motion"; fixed_triangle=(6,4,5), skip_stretch=true, target_stretch=0.5, tip_value=0.5)
end


@testset "volume_tetrahedron" begin
    p = rand(Float64,4,3)
    p[1,:] = [0,0,0]
    p[2,2:3] = [0,0]
    p[3,3] = 0 
    F = Framework([[1,2], [1,3], [1,4], [2,4],[2,3], [3,4]], Matrix(p'))
    triangles = [(1,2,3),(1,2,4),(1,3,4),(2,3,4)]
    area_equations = Vector{Expression}([])
    for (i,triang) in enumerate(triangles)
        parallelogram_area = sum(cross(F.G.xs[:,triang[2]]-F.G.xs[:,triang[1]], F.G.xs[:,triang[3]]-F.G.xs[:,triang[1]]).^2)
        push!(area_equations, 1/4*parallelogram_area)
    end

    volume_constraint = 1/6*det(hcat([1 for _ in 1:4], F.G.xs[1,:], F.G.xs[2,:], F.G.xs[3,:]))
    push!(area_equations, volume_constraint)
    area_equations = area_equations .- evaluate.(area_equations, F.G.variables=>vcat([p[i,:] for i in 1:size(p)[1]]...))
    equations!(F, area_equations)
    point = newton_correct(F.G, to_Array(F, F.G.realization))
    realization!(F.G, to_Matrix(F, point))
    D = DeformationPath(F, [1], 500; step_size=0.01, newton_tol=1e-15)
    animate(D,F,"volume_tetrahedron"; padding=0.05, fixed_vertices=(1,2,3), animate_rotation=false, rotation_frames=1500, fixed_direction=[1,0,0], filetype="mp4")
    project_deformation_random(D, 3)
end