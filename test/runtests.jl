import DeformationPaths:    Framework,
                            AngularFramework,
                            DeformationPath,
                            DeformationPath_EdgeContraction,
                            Polytope,
                            animate, 
                            VolumeHypergraph,
                            SpherePacking,
                            plot,
                            plot!,
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
                            BodyBar,
                            triangle_shrinking,
                            fix_antipodals!,
                            tetrahedral_symmetry!,
                            add_equations!,
                            compute_nontrivial_inf_flexes,
                            ConstraintSystem,
                            compute_nonblocked_flex,
                            stich_deformation_paths,
                            add_shadow!
using Test
using HomotopyContinuation
using LinearAlgebra
using IterTools
using Colors
import GLMakie: save, scatter!, Point2f, MultiLightShading

is_no_ci = !(get(ENV, "GITHUB_ACTIONS", "false") == "true")

teal = RGB(0/255, 128/255, 128/255)
soft_teal = RGB(160/255,218/255,218/255)
coral=RGB(255/255, 127/255, 80/255)
logocolors = Colors.JULIA_LOGO_COLORS

@testset "DeformationPaths.jl" begin
    include("SpherePacking.jl")
    include("various.jl")
    include("VolumeHypergraph.jl")
    include("Framework.jl") 
    include("Polytope.jl")
end