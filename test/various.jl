@testset "bodyhinge_pyramid" begin
    F = BodyHinge([[1,2,3],[1,3,4],[1,4,5],[1,5,6],[1,6,2]], Matrix([0 0 1; cos(2*pi/5) sin(2*pi/5) 0; cos(4*pi/5) sin(4*pi/5) 0; cos(6*pi/5) sin(6*pi/5) 0; cos(8*pi/5) sin(8*pi/5) 0; cos(10*pi/5) sin(10*pi/5) 0;]'))
    plot(F,"bodyhinge_pyramid")
    D = DeformationPath(F, [], 200; step_size=0.0075)
    animate(D,F,"bodyhinge_pyramid_motion"; filetype="mp4", animate_rotation=true, rotation_frames=450, padding=0.01)
end

@testset "bodyhinge" begin
    F = BodyHinge([[1,2,3,4],[3,4,5,6]], Matrix([0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 1 1; 1 1 1]'))
    plot(F,"bodyhinge")
    D = DeformationPath(F, [], 200; step_size=0.025)
    animate(D,F,"bodyhinge_motion"; filetype="mp4")
end


@testset "thales" begin
    F = AngularFramework([[1,3,2]], Matrix([-1 0; 1 0; -sqrt(1/2) sqrt(1/2);]'); pinned_vertices=[1,2])
    plot(F,"thales"; padding=0.1, pin_point_offset=0.075)
    D = DeformationPath(F, [1], 250; step_size=0.025)
    animate(D,F,"thales_motion"; padding=0.075, pin_point_offset=0.075, filetype="mp4")
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
    animate(D,F,"sphericaldiskpacking_motion"; filetype="mp4")
end
