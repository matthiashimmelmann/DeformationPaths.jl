import DeformationPaths:    Framework,
                            AngularFramework,
                            DeformationPath,
                            DeformationPath_EdgeContraction,
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
                            BodyHinge,
                            triangle_shrinking,
                            fix_antipodals!,
                            tetrahedral_symmetry!,
                            add_equations!,
                            compute_nontrivial_inf_flexes,
                            ConstraintSystem,
                            compute_nonblocked_flex,
                            stich_deformation_paths
using Test
using HomotopyContinuation
using LinearAlgebra

is_no_ci = !(get(ENV, "GITHUB_ACTIONS", "false") == "true")
@testset "DeformationPaths.jl" begin
    include("Framework.jl")
    include("Polytope.jl")
    include("SpherePacking.jl")
    include("VolumeHypergraph.jl")
    include("various.jl")
end