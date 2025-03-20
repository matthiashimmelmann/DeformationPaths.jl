import DeformationPaths:    Framework, 
                            DeformationPath,
                            Polytope,
                            animate, 
                            VolumeHypergraph, 
                            plot
using Test

@testset "cube" begin
    F = Polytope([[1,2,3,4],[5,6,7,8],[1,2,5,6],[2,3,6,7],[3,4,7,8],[1,4,5,8]], Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1]'))
    plot(F,"cube")
    D = DeformationPath(F, [1,1,1], 100; step_size=0.01)
    animate(D,F,"cube_motion")
end

@testset "square" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 1 1; 0 1]'))
    plot(F,"square")
    D = DeformationPath(F, [1], 200; step_size=0.025)
    animate(D,F,"square_motion")
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

@testset "two_triangles" begin
    F = VolumeHypergraph([[1,2,3],[2,3,4]], Matrix([0 0; 1 0; 0 1; 1 1]'))
    plot(F,"two_triangles")
    D = DeformationPath(F, [1], 100; step_size=0.01)
    animate(D, F,"two_triangles_motion"; fixed_triangle=(1,2,3),tip_value=0)
end

@testset "octehedral_decomposition" begin
    F = VolumeHypergraph([[1,3,6],[1,2,5],[2,3,4],[1,2,3],[6,4,5]], Matrix([0 0; 3 0; 0 3; 1 1; 1 0.5; 0.5 1]'))
    plot(F,"octahedral_decomposition")
    D = DeformationPath(F, [1, 1, 1], 150; step_size=0.01)
    animate(D, F,"octahedral_decomposition_motion"; fixed_triangle=(6,4,5), tip_value=1.)
end
