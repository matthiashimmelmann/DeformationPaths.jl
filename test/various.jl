@testset "bodyhinge_pyramid" begin
    F = BodyHinge([[1,2,3],[1,3,4],[1,4,5],[1,5,6],[1,6,2]], Matrix([0 0 1; cos(2*pi/5) sin(2*pi/5) 0; cos(4*pi/5) sin(4*pi/5) 0; cos(6*pi/5) sin(6*pi/5) 0; cos(8*pi/5) sin(8*pi/5) 0; cos(10*pi/5) sin(10*pi/5) 0;]'))
    plot(F)
    D = DeformationPath(F, [], 100; step_size=0.0075)
    animate(D,F; filetype="mp4", animate_rotation=true, rotation_frames=450, padding=0.01)
end

@testset "bodyhinge" begin
    F = BodyHinge([[1,2,3,4],[3,4,5,6]], Matrix([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 1 1 1]'))
    plot(F)
    D = DeformationPath(F, [], 100; step_size=0.025)
    animate(D,F; filetype="mp4")
end


@testset "bricard_octahedron" begin
    F = BodyHinge([[1,2,4], [1,2,3], [6,2,4], [6,2,3], [1,3,5], [6,3,5], [1,3,5], [6,3,5]], Matrix([0 0 -1; -1 -1 0; 1 -1 0; 1 1 0; -1 1 0; 0 0 1]'), pinned_vertices=[1,2,3])
    plot(F)
    D = DeformationPath(F, [], 400; step_size=0.025, show_progress=false)
    if is_no_ci
        animate(D,F,"../animations/bricard_octahedron"; filetype="mp4", renderEntirePolytope=true, azimuth=0.15*pi, elevation=pi/2, animate_rotation=false,  padding=0.01)
    end
end



@testset "thales" begin
    F = AngularFramework([[1,3,2]], Matrix([-1 0; 1 0; -sqrt(1/2) sqrt(1/2);]'); pinned_vertices=[1,2])
    plot(F; padding=0.1, pin_point_offset=0.075)
    D = DeformationPath(F, [1], 200; step_size=0.025)
    if is_no_ci
        animate(D,F; padding=0.075, pin_point_offset=0.075, filetype="mp4")
    end
end


@testset "squareonhyperboloid" begin
    F = FrameworkOnSurface([[1,2],[2,3],[3,4],[1,4]], Matrix([-sqrt(1/2) -sqrt(1/2) -1; -1 0 0; 0 1 0; sqrt(1/2) sqrt(1/2) 1]'), x->x[1]^2+x[2]^2-x[3]^2-1)
    plot(F; vertex_labels=false, azimuth=3*pi/4, padding=0.35)
    D = DeformationPath(F, [1,1], 200; step_size=0.035, show_progress=false)
    if is_no_ci
        animate(D,F; vertex_labels=false, animate_rotation=true, filetype="mp4")
    end
end


@testset "sphericaldiskpacking" begin
    F = SphericalDiskPacking([(1,2),(1,3),(1,4),(1,5),(2,3),(2,4),(3,5),(4,5),(2,6),(3,6),(4,6),(5,6)], Matrix([sqrt(2) 0 0; 0 sqrt(2) 0; 0 0 sqrt(2); 0 -sqrt(2) 0; 0 0 -sqrt(2); -sqrt(2) 0 0]'); pinned_vertices=[1])
    plot(F)
    D = DeformationPath(F, [1], 100; step_size=0.01, show_progress=false)
    if is_no_ci
        animate(D,F; filetype="mp4")
    end
end
