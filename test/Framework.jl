@testset "3Prism" begin
    F = Framework([(1,2), (1,3), (2,3), (4,5), (5,6), (4,6), (1,4), (2,5), (3,6)], Matrix( [0 0; 0 1; sqrt(3)/2 0.5; 1.05 0; 1.05 1; 1.05+sqrt(3)/2 0.5]'))
    @test !is_prestress_stable(F)
    @test !is_second_order_rigid(F)
    inf_flexes = compute_inf_flexes(F.G, to_Array(F.G, F.G.realization))
    @test size(inf_flexes)[2] == 1+3
    stresses = compute_equilibrium_stresses(F.G, to_Array(F.G, F.G.realization))
    @test size(stresses)[2] == 1

    F = Framework([(1,2), (1,3), (2,3), (4,5), (5,6), (4,6), (1,4), (2,5), (3,6)], Matrix( [0 0; 0 1; sqrt(3)/2 0.5; 1.05 0; 1.05 1; 1.05+sqrt(3)/2 0.5]'); pinned_vertices=[1,4])
    D = DeformationPath(F, [-1], 27; step_size=0.025)
    F2 = Framework([(1,2), (1,3), (2,3), (4,5), (5,6), (4,6), (1,4), (2,5), (3,6)], D.motion_matrices[end]) 
    fig, ax = plot(F2; edge_color=:lightgrey, flex_color=coral, show_pins=false, vertex_labels=false, padding=0.1, vertex_color=:lightgrey, vertex_size=7)
    plot!(ax, F; edge_color=teal, flex_color=coral, plot_flexes=false, show_pins=false, vertex_labels=false, padding=nothing)
    add_shadow!(ax, F, D; flex_color=coral)
    points = [Point2f(D.motion_matrices[1][:,j]) for j in 1:size(D.motion_matrices[end])[2]]
    scatter!(ax, points; color=:black, markersize=55)
end


@testset "square" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 1 1; 0 1]'))
    plot(F)
    @test !is_rigid(F)
    D = DeformationPath(F, [1], 350; step_size=0.025)
    @test !is_prestress_stable(F)

    F1 = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 1 1; 0 1]'); pinned_vertices=[1,2])
    D1 = DeformationPath(F1, [1], 357; step_size=0.025)
    if is_no_ci
        animate(D1,F1,"square"; edge_color=teal, padding=0.1, vertex_size=15, vertex_color=teal, vertex_labels=false, show_pins=false, filetype="mp4")
    end

    F2 = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 1 1; 1 0]'); pinned_vertices=[1,2])
    D2 = DeformationPath(F2, [1], 350; step_size=0.025)

    F3 = Framework([[1,2],[2,3],[3,4],[1,4]], Matrix([0. 0; 1 0; 0 0; 0 1]'); pinned_vertices=[1,2])
    D3 = DeformationPath(F3, [1], 350; step_size=0.025)

    if is_no_ci
        project_deformation_random([D1,D2,D3], F, 2, "square_realizations"; animate=true, padding=nothing, vertex_size=65, line_width=8)
    end
end


@testset "triangle" begin
    F = Framework([[1,2],[2,3],[1,3]], Matrix([0. 0; 1 0; 0.5 sqrt(3)/2]'))
    plot(F; edge_color=teal, padding=1, vertex_size=15, vertex_color=teal, vertex_labels=false)
    @test is_rigid(F)
    @test is_inf_rigid(F)
end


@testset "rigid_prestress_stable" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[1,5],[3,5],[4,5]], Matrix([0. 0; 1 0; 2 0; 1 1; 1 2]'); pinned_vertices=[1,4])
    plot(F; edge_color=teal, flex_color=coral, padding=0.25, plot_flexes=true, flex_Real=[1], show_pins=false, flex_scale=0.5, vertex_labels=false)
    @test !is_inf_rigid(F)
    @test is_prestress_stable(F)
    @test is_second_order_rigid(F)
    @test is_rigid(F)
end


@testset "coned_cube" begin
    F = Framework(vcat([[1,2],[2,3],[3,4],[1,4],[1,5],[2,6],[3,7],[4,8],[5,6],[6,7],[7,8],[5,8]],[[i,9] for i in 1:8]), Matrix([-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1; -1 -1 1; 1 -1 1; 1 1 1; -1 1 1; 0 0 1.65]'))
    plot(F; edge_color=teal, flex_color=coral, padding=0.5, plot_flexes=true, flex_Real=[1,0], show_pins=false, flex_scale=0.2, vertex_labels=false)
    D = DeformationPath(F, [0.5,0.5], 500; step_size=0.02, show_progress=false)
    animate(D,F; filetype="mp4")
end


@testset "flexible_prestress_stable_component" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[1,5],[3,5],[4,5],[1,6]], Matrix([0. 0; 1 0; 2 0; 1 1; 1 2; -sqrt(1/2) sqrt(1/2)]'); pinned_vertices=[1,4])
    plot(F; edge_color=teal, flex_color=coral, show_pins=false, flex_Real=[-1,-1], padding=0.25, plot_flexes=true, flex_scale=0.5, vertex_labels=false)
    D = DeformationPath(F, [1,1], 200; step_size=0.025, show_progress=false)
    if is_no_ci
        animate(D,F; filetype="mp4")
    end
end


@testset "K_4_edge" begin
    F = Framework([[1,2],[2,3],[3,4],[1,4],[2,4],[1,3],[2,5]], Matrix([0. 0; 1 0; 1 1; 0 1; 2 0]'); pinned_vertices=[1,2])
    D = DeformationPath(F, [-1], 30; step_size=0.025)
    F2 = Framework([[1,2],[2,3],[3,4],[1,4],[2,4],[1,3],[2,5]], D.motion_matrices[end]; pinned_vertices=[1,2]) 
    fig, ax = plot(F2; edge_color=:lightgrey, flex_color=coral, show_pins=false, vertex_labels=false, vertex_color=:lightgrey, vertex_size=7)
    plot!(ax, F; edge_color=teal, flex_color=coral, plot_flexes=false, show_pins=false, vertex_labels=false, padding=0.15)
    add_shadow!(ax, F, D; flex_color=coral)
    points = [Point2f(D.motion_matrices[1][:,j]) for j in 1:size(D.motion_matrices[end])[2]]
    scatter!(ax, points; color=:black, markersize=55)
end


@testset "rigid_test" begin
    F = Framework([[1,2],[1,4],[1,5],[4,5],[4,3],[5,3],[2,6],[2,7],[6,7],[3,6],[3,7]], Matrix([0 0; 2 0; 1 1; 0.5-1/6 1/2+1/6; 0.5+1/6 1/2-1/6; 1.5+1/6 1/2+1/6; 1.5-1/6 1/2-1/6;]'))
    plot(F; edge_color=teal, vertex_labels=false)
    @test is_prestress_stable(F)
end

@testset "3Frustum" begin
    F = Framework([[1,2],[2,3],[1,3],[4,5],[5,6],[4,6],[1,4],[2,5],[3,6]], Matrix([cos(2*pi/3) sin(2*pi/3); cos(4*pi/3) sin(4*pi/3); cos(6*pi/3) sin(6*pi/3); 2*cos(2*pi/3) 2*sin(2*pi/3); 2*cos(4*pi/3) 2*sin(4*pi/3); 2*cos(6*pi/3) 2*sin(6*pi/3);]'), pinned_vertices=[1,2,3])
    plot(F; edge_color=teal, flex_color=coral, plot_flexes=true, show_pins=false, flex_scale=0.85, vertex_labels=false)
    @test is_prestress_stable(F)
end


@testset "double_watt" begin
    F = Framework([[1,2],[2,3],[2,4],[3,9],[3,4],[3,5],[4,5],[5,6],[6,7],[7,8],[7,9],[8,9],[8,10],[9,10],[10,11]], Matrix([0 0; 1 0; 2 1; 1 2; 3 2; 4 2; 5 2; 7 2; 6 1; 7 0; 8 0;]'); pinned_vertices=[1,6,11])
    plot(F; padding=0.35, pin_point_offset=0.2, edge_color=teal)
    D = DeformationPath(F, [0.5,0.5], 500; step_size=0.05)
    if is_no_ci
        animate(D,F; padding=0.352, edge_color=teal, fixed_vertices=(1,6), fixed_direction=[4,2], pin_point_offset=0.2, filetype="mp4")
        for i in 1:50
            project_deformation_random(D, F, 2, "Double_Watt$i"; padding=nothing, vertex_size=85, line_width=11)
        end
    end
end



@testset "twoprism" begin
    F = Framework([[1,2],[2,3],[3,1],[1,4],[2,5],[3,6],[4,5],[5,6],[6,4]], Matrix([0 0; 2 0; 1 1; 0 2; 2 2; 1 3]'))
    plot(F)
    D = DeformationPath(F, [1], 220; step_size=0.025, show_progress=false)
    if is_no_ci
        animate(D,F)
    end
end


@testset "completebipartite" begin
    F = Framework([[1,3],[1,4],[1,5],[1,6],[2,3],[2,4],[2,5],[2,6]], Matrix([0 0; 0 1; 1 -1; 1 0; 1 1; 1 2]'))
    plot(F)
    D = DeformationPath(F, [1], 500; step_size=0.025, show_progress=false)
    animate(D,F; filetype="mp4")
    @test !is_prestress_stable(F)
end


@testset "bricard_octahedron" begin
    F = Framework([(1,2),(1,3),(1,4),(1,5),(2,6),(3,6),(4,6),(5,6),(2,4),(4,3),(3,5),(5,2)], Matrix([0 0 -1; -1 -1 0; 1 -1 0; 1 1 0; -1 1 0; 0 0 1]'))
    plot(F)
    D = DeformationPath(F, [1], 200; step_size=0.02, show_progress=false)
end