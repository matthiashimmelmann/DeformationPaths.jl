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
                            is_inf_rigid,
                            is_second_order_rigid,
                            BodyHinge
using Test
using HomotopyContinuation
using LinearAlgebra

@testset "SemialgebraicOpt" begin
    include("Framework.jl")
    include("Polytope.jl")
    include("SpherePacking.jl")
    include("VolumeHypergraph.jl")
    include("various.jl")
end