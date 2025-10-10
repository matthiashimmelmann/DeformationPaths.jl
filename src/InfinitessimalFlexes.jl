"""
    compute_nontrivial_inf_flexes(G, point, K_n[; tol])

Compute the nontrivial infinitesimal flexes of a geometric constraint system `G` in `point`.
"""
function compute_nontrivial_inf_flexes(G::ConstraintSystem, point::Vector{<:Real}, K_n::ConstraintSystem; tol::Real=1e-8)::Matrix{<:Real}
    inf_flexes = nullspace(evaluate(G.jacobian, G.variables=>point); atol=tol)
    trivial_inf_flexes = nullspace(evaluate(typeof(K_n)==ConstraintSystem ? K_n.jacobian : K_n.G.jacobian, (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables)=>point[1:length( (typeof(K_n)==ConstraintSystem ? K_n.variables : K_n.G.variables))]); atol=tol)
    s = size(trivial_inf_flexes)[2]+1
    extend_basis_matrix = trivial_inf_flexes
    for inf_flex in [inf_flexes[:,i] for i in 1:size(inf_flexes)[2]]
        tmp_matrix = hcat(trivial_inf_flexes, inf_flex)
        if !(rank(tmp_matrix; atol=tol) == rank(trivial_inf_flexes; atol=tol))
            extend_basis_matrix = hcat(extend_basis_matrix, inf_flex)
        end
    end
    Q, R = qr(extend_basis_matrix)
    Q = Q[:, s:rank(R, atol=tol)]
    return Q
end


"""
    compute_nonblocked_flex(F[; tol_rank_drop, tol])

Compute an infinitesimal flex of `F` that is not blocked by an equilibrium stress.
"""
function compute_nonblocked_flex(F::AllTypes; fast_search::Bool=true, tol_rank_drop::Real=1e-6, tol::Real=1e-12)::Vector
    if typeof(F)==Framework
        K_n = Framework([[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices).G
    elseif typeof(F)==Polytope || typeof(F)==SpherePacking || typeof(F)==BodyHinge
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in 1:length(F.G.vertices) for j in 1:length(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    else
        throw("Type of F is not yet supported. It is $(typeof(F)).")
    end
    flexes = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), K_n; tol=tol_rank_drop)
    if size(flexes)[2]==0
        return []
    end
    rigidity_matrix = evaluate.(F.G.jacobian, F.G.variables=>to_Array(F, F.G.realization))
    stresses = nullspace(rigidity_matrix'; atol=tol_rank_drop)
        
    @var λ[1:size(flexes)[2]] ω[1:size(stresses)[2]]
    parametrized_flex = flexes*λ
    parametrized_stress = stresses*ω
    stress_energy = parametrized_stress'*evaluate.(F.G.jacobian, F.G.variables=>Vector{Expression}(parametrized_flex))*parametrized_flex
    stress_poly_system = differentiate(stress_energy, ω)
    projective_stress_system = vcat(stress_poly_system, sum(λ .^ 2) - 1)

    if fast_search
        J_stress_energy = Matrix{Expression}(hcat([differentiate(eq, λ) for eq in projective_stress_system]...)')
        for _ in 1:size(flexes)[2]*size(stresses)[2]*3
            rand_flex_parameter = randn(Float64, size(flexes)[2])
            rand_flex_parameter = rand_flex_parameter ./ norm(rand_flex_parameter)
            try
                q = newton_correct(projective_stress_system, λ, J_stress_energy, rand_flex_parameter; tol = tol, time_penalty=2)
                return q
            catch e
                display(e)
                continue
            end
        end
    else
        try 
            res = solve(projective_stress_system; show_progress = true)
            return real_solutions(res)[1]
        catch
            nothing
        end

        N = size(flexes)[2]*size(stresses)[2]*3 # ambient dimension
        #TODO add Paul's test
        for i in 1:length(λ)
            try
                L₀ = rand_subspace(length(λ); codim = length(λ)-i)
                R_L₀ = solve(projective_stress_system; target_subspace = L₀)
                Ω = solve(
                    f,
                    solutions(R_L₀);
                    start_subspace = L₀,
                    target_subspaces = [rand_subspace(N; codim = length(λ)-i, real = true) for _ in 1:N],
                    transform_result = (R,p) -> real_solutions(R),
                    flatten = true
                )
                return Ω[1]
            catch e
                continue
            end
        end
    end
    return []
end
