@testset "diskpacking" begin
    F = SpherePacking([1.,1.,1.,1.,1.], Matrix([0 0; 2 0; 3 sqrt(3); 4 0; 6 0]'); pinned_vertices=[1,5])
    plot(F,"diskpacking")
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0; 1.75 -sqrt(2^2-(1.75)^2); 3.5 0; 4.5 sqrt(3)]'); pinned_vertices=[1])
    D = DeformationPath(F, [1,1], 250; step_size=0.025)
    animate(D,F,"diskpacking_motion"; filetype="mp4")
end


@testset "squarediskpacking" begin
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0; 2 0; 0 2; 2 2]'); pinned_vertices=[1])
    plot(F,"squarediskpacking")
    D = DeformationPath(F, [1], 250; step_size=0.01)
    animate(D,F,"squarediskpacking_motion"; filetype="mp4")
end


@testset "spherepacking" begin
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0 0; 2 0 0; 0 2 0; 0 0 2]'), pinned_vertices = [1,2])
    plot(F,"spherepacking")
    D = DeformationPath(F, [1,1,1], 500; step_size=0.04)
    animate(D,F,"spherepacking_motion"; filetype="mp4")
end
