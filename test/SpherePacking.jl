@testset "diskpacking" begin
    F = SpherePacking([1.,1.,1.,1.,1.], Matrix([0 0; 2 0; 3 sqrt(3); 4 0; 6 0]'); pinned_vertices=[1,5])
    plot(F)
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0; 1.75 -sqrt(2^2-(1.75)^2); 3.5 0; 4.5 sqrt(3)]'); pinned_vertices=[1])
    D = DeformationPath(F, [1,1], 250; step_size=0.025, show_progress=false)
    animate(D,F; filetype="mp4")
end


@testset "squarediskpacking" begin
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0; 2 0; 0 2; 2 2]'); pinned_vertices=[1])
    plot(F)
    D = DeformationPath(F, [1], 250; step_size=0.01, show_progress=false)
    animate(D,F; filetype="mp4")
end


@testset "spherepacking" begin
    F = SpherePacking([1.,1.,1.,1.], Matrix([0 0 0; 2 0 0; 0 2 0; 0 0 2]'), pinned_vertices = [1,2])
    plot(F)
    D = DeformationPath(F, [1,1,1], 500; step_size=0.04, show_progress=false)
    animate(D,F; filetype="mp4")
end

@testset "hypostaticspherepacking" begin
    #INFO cf. Miranda Holmes-Cerfon
    realization = Matrix([0 0 0; 1.000000000000000 0.000000000000000 0.000000000000000; -0.500000000000000 0.866025403784439 0.000000007289415; 1.000000003306545 1.603750749657996 0.453609204877056; 0.999999994048218 0.577350265753363 -0.816496583357531; -0.000000003967855 1.539600726278313 -0.544331039863996; 0.000000003306546 1.603750740716201 0.453609226708918; 0.999999996032145 1.539600715548160 -0.544331060431297; 1.500000000000000 0.866025403784439 -0.000000007289415; 0.500000000000000 0.866025403784439 -0.000000000000000]')
    F = SpherePacking([1. for _ in 1:10], realization)
    plot(F, "hypostaticspherepacking"; sphere_radius=0.25,alpha=0.45,sphere_color=soft_teal, vertex_labels=false)
end