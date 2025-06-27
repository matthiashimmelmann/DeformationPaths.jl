@testset "two_triangles" begin
    F = VolumeHypergraph([[1,2,3],[2,3,4]], Matrix([0 0; 1 0; 0 1; 1 1]'))
    plot(F,"two_triangles")
    D = DeformationPath(F, [1], 100; step_size=0.01)
    animate(D, F,"two_triangles_motion"; fixed_triangle=(1,2,3),tip_value=0,skip_stretch=false, filetype="mp4")
end


@testset "two_triangles" begin
    F = VolumeHypergraph([[1,2,3],[2,3,4]], Matrix([0 0; 1 0; 0 1; 1 1]'))
    plot(F,"two_triangles")
    D = DeformationPath(F, [1], 100; step_size=0.01)
    animate(D, F,"two_triangles_motion"; fixed_triangle=(1,2,3),tip_value=0,skip_stretch=false, filetype="mp4")
end


@testset "octehedral_decomposition" begin
    F = VolumeHypergraph([[1,3,6],[1,2,5],[2,3,4],[1,5,6],[6,4,5]], Matrix([0 0; 3 0; 0 3; 1 1; 1 0.5; 0.5 1]'))
    plot(F,"octahedral_decomposition")
    D = DeformationPath(F, [0.333, 1], 350; step_size=0.002)
    animate(D, F,"octahedral_decomposition_motion"; fixed_triangle=(6,4,5), skip_stretch=true, target_stretch=0.5, tip_value=0.5, filetype="mp4")
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