@testset "square" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 1 1; 0 1]'))
    plot(F)
    @test !is_rigid(F)
    D = DeformationPath(F, [1], 200; step_size=0.025)
    animate(D,F; filetype="mp4")
end


@testset "rigid_prestress_stable" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[1,5],[3,5],[4,5]], Matrix([0. 0; 1 0; 2 0; 1 1; 1 2]'))
    @test !is_inf_rigid(F)
    @test is_second_order_rigid(F)
    #@test is_rigid(F)
end


@testset "flexible_prestress_stable_component" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[1,5],[3,5],[4,5],[1,6]], Matrix([0. 0; 1 0; 2 0; 1 1; 1 2; 0 -1]'))
    plot(F)
    D = DeformationPath(F, [1,1], 200; step_size=0.025, show_progress=false)
    animate(D,F; filetype="mp4")
end


@testset "K_4" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[2,4],[1,3],[1,5]], Matrix([0. 0; 1 0; 1 1; 0 1; 0 -1]'))
    plot(F)
    D = DeformationPath(F, [1], 200; step_size=0.025, show_progress=false)
    animate(D,F)
end


@testset "twoprism" begin
    F = Framework([[1,2],[2,3],[3,1],[1,4],[2,5],[3,6],[4,5],[5,6],[6,4]], Matrix([0 0; 2 0; 1 1; 0 2; 2 2; 1 3]'))
    plot(F)
    D = DeformationPath(F, [1], 220; step_size=0.025, show_progress=false)
    animate(D,F)
end


@testset "completebipartite" begin
    F = Framework([[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6]], Matrix([0 0; 0 1; 1 -1; 1 0; 1 1; 1 2]'))
    plot(F)
    D = DeformationPath(F, [1], 500; step_size=0.025, show_progress=false)
    animate(D,F; filetype="mp4")
end


@testset "coned_cube" begin
    F = Framework(vcat([[1,2],[2,3],[3,4],[1,4],[1,5],[2,6],[3,7],[4,8],[5,6],[6,7],[7,8],[5,8]],[[i,9] for i in 1:8]), Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1; 0 0 sqrt(2)]'))
    plot(F)
    D = DeformationPath(F, [0.5,0.5], 500; step_size=0.02, show_progress=false)
    animate(D,F; filetype="mp4")
end


@testset "double_watt" begin
    F = Framework([[1,2],[2,3],[2,4],[3,9],[3,4],[3,5],[4,5],[5,6],[6,7],[7,8],[7,9],[8,9],[8,10],[9,10],[10,11]], Matrix([0 0; 1 0; 2 1; 1 2; 3 2; 4 2; 5 2; 7 2; 6 1; 7 0; 8 0;]'); pinned_vertices=[1,6,11])
    plot(F; padding=0.35, pin_point_offset=0.2)
    D = DeformationPath(F, [0.5,0.5], 500; step_size=0.05, show_progress=false)
    animate(D,F; padding=0.35, fixed_vertices=(1,6), fixed_direction=[4,2], pin_point_offset=0.2, filetype="mp4")
end
