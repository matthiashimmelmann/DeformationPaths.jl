
"""
    is_rigid(F[; tol, tol, tested_random_flexes, symmetric_newton])

Heuristically checks if a geometric constraint system `F` is (continuously) rigid. 
"""
function is_rigid(F::AllTypes; tol_rank_drop::Real=1e-7, tol::Real=1e-13, tested_random_flexes::Int=8, symmetric_newton::Bool=false)::Bool
    if is_inf_rigid(F; tol_rank_drop=tol_rank_drop)
        return true
    end
    for _ in 1:tested_random_flexes
        D = DeformationPath(F, [], 5; step_size=sqrt(tol_rank_drop), tol=tol, random_flex=true, symmetric_newton=symmetric_newton)
        if any(sample->norm(sample-D.motion_samples[1], Inf)>tol_rank_drop, D.motion_samples)
            return false
        end
    end
    return true
end

"""
    is_inf_rigid(F[; tol_rank_drop])

Checks if a geometric constraint system `F` is infinitesimally rigid.
"""
function is_inf_rigid(F::AllTypes; tol_rank_drop::Real=1e-8)::Bool
    if typeof(F)==Framework
        K_n = Framework([[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices).G
    elseif typeof(F)==AngularFramework
        K_n = AngularFramework([[i,j,k] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) for k in 1:length(F.G.vertices) if (i<j && j<k) || (i<k && k<j) || (j<i && i<k)], F.G.realization; pinned_vertices=F.G.pinned_vertices).G
    elseif typeof(F)==FrameworkOnSurface
        K_n = deepcopy(G)
        add_equations!(K_n, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(G.vertices) for j in 1:length(G.vertices) if i<j]])
    elseif typeof(F)==VolumeHypergraph
        K_n = VolumeHypergraph(collect(powerset(F.G.vertices, F.G.dimension+1, F.G.dimension+1)), F.G.realization).G
    elseif typeof(F)==Polytope || typeof(F)==SpherePacking || typeof(F)==BodyHinge
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    elseif typeof(F)==SphericalDiskPacking
        minkowski_scalar_product(e1,e2) = e1'*e2-1
        inversive_distances = [minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[2]])/sqrt(minkowski_scalar_product(F.G.realization[:,contact[1]], F.G.realization[:,contact[1]]) * minkowski_scalar_product(F.G.realization[:,contact[2]], F.G.realization[:,contact[2]])) for contact in powerset(F.G.vertices, 2, 2)]
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, [minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(F.G.xs[:,contact[1]], F.G.xs[:,contact[1]]) * minkowski_scalar_product(F.G.xs[:,contact[2]], F.G.xs[:,contact[2]]) for (i,contact) in enumerate(powerset(F.G.vertices, 2, 2))], F.G.realization, F.G.xs)
    else
        throw("Type of F is not yet supported. It is $(typeof(F)).")
    end
    inf_flexes = nullspace(evaluate(F.G.jacobian, F.G.variables=>to_Array(F, F.G.realization)); atol=tol_rank_drop)
    trivial_inf_flexes = nullspace(evaluate(typeof(K_n)==ConstraintSystem ? K_n.jacobian : K_n.G.jacobian, (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables)=>to_Array(F, F.G.realization)[1:length( (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables))]); atol=tol_rank_drop)
    #println("flexes: $(size(inf_flexes)[2]), nontrivial: $(size(inf_flexes)[2]-size(trivial_inf_flexes)[2])")
    return length(inf_flexes) == length(trivial_inf_flexes)
end


"""
    is_second_order_rigid(F)

Checks if a geometric constraint system `F` is second-order rigid.

See also [`compute_nonblocked_flex`](@ref)
"""
function is_second_order_rigid(F::AllTypes; kwargs...)::Bool
    flex_mult = compute_nonblocked_flex(F; kwargs...)
    if length(flex_mult)==0
        return true
    else
        return false
    end
end