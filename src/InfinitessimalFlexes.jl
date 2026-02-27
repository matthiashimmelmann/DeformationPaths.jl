"""
    compute_inf_flexes(G, point[; tol])

Compute all infinitesimal flexes of a geometric constraint system `G` in `point`.
"""
function compute_inf_flexes(G::ConstraintSystem, point::Vector{<:Real}; tol::Real=1e-8)::Matrix{<:Real}
    inf_flexes = nullspace(evaluate(G.jacobian, G.variables=>point); atol=tol)
    return inf_flexes
end


"""
    compute_equilibrium_stresses(G, point[; tol])

Compute all equilibrium stresses of a geometric constraint system `G` in `point`.
"""
function compute_equilibrium_stresses(G::ConstraintSystem, point::Vector{<:Real}; tol::Real=1e-8)::Matrix{<:Real}
    stresses = nullspace(evaluate(G.jacobian, G.variables=>point)'; atol=tol)
    return stresses
end

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
    compute_jk_flexes(G, j, k, point)

"""
function compute_jk_flexes(G::ConstraintSystem, j::Int, k::Int, point::Vector{<:Real}; tol::Real=1e-8)::Vector{Expression}
    rig_mat_pinv = pinv(evaluate(G.jacobian, G.variables=>point))
    Q = compute_inf_flexes(G, point; tol=tol)
    first_order = [Q[:,i] for i in 1:size(Q)[2]]
    second_order = [rig_mat_pinv * (-evaluate(G.jacobian, G.variables=>first_order[i])*first_order[i]) for i in 1:length(first_order)]
    third_order = [rig_mat_pinv * (-2*evaluate(G.jacobian, G.variables=>first_order[i])*second_order[i] - evaluate(G.jacobian, G.variables=>second_order[i])*first_order[i]) for i in 1:length(first_order)]
    fourth_order = [rig_mat_pinv * (-3*evaluate(G.jacobian, G.variables=>first_order[i])*third_order[i] - 3*evaluate(G.jacobian, G.variables=>second_order[i])*second_order[i] - evaluate(G.jacobian, G.variables=>third_order[i])*first_order[i]) for i in 1:length(first_order)]
    @var _a[1:4,1:length(first_order)] t
    return point + t*_a[1,:]'*first_order + t^2*_a[2,:]'*second_order + t^3*_a[3,:]'*third_order + t^4*_a[4,:]'*fourth_order
end

"""
    compute_nonblocked_flex(F[; tol_rank_drop, tol])

Compute an infinitesimal flex of `F` that is not blocked by an equilibrium stress.
"""
function compute_nonblocked_flex(F::AllTypes; fast_search::Bool=false, tol_rank_drop::Real=1e-6, tol::Real=1e-12)::Vector
    if typeof(F)==Framework
        K_n = Framework([[i,j] for i in eachindex(F.G.vertices) for j in eachindex(F.G.vertices) if i<j], F.G.realization; pinned_vertices=F.G.pinned_vertices).G
    elseif typeof(F)==Polytope || typeof(F)==SpherePacking || typeof(F)==BodyHinge || typeof(F)==FacetPolytope
        K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(F.G.vertices) for j in eachindex(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    else
        throw("Type of F is not yet supported. It is $(typeof(F)).")
    end
    flexes = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), K_n; tol=tol_rank_drop)
    if size(flexes)[2]==0
        return []
    end
    rigidity_matrix = evaluate.(F.G.jacobian, F.G.variables=>to_Array(F, F.G.realization))
    stresses = nullspace(rigidity_matrix'; atol=tol_rank_drop)
    if size(stresses)[2]==0
        return flexes[1,:]
    end
    @var λ[1:size(flexes)[2]] ω[1:size(stresses)[2]]
    parametrized_flex = flexes*λ
    parametrized_stress = stresses*ω
    stress_energy = parametrized_stress'*evaluate.(F.G.jacobian, F.G.variables=>Vector{Expression}(parametrized_flex))*parametrized_flex
    stress_poly_system = differentiate(stress_energy, ω)
    projective_stress_system = vcat(stress_poly_system, sum(λ .^ 2) - 1)

    codim = rank(evaluate.(differentiate(projective_stress_system, λ), λ=>randn(ComplexF64, length(λ))); atol=1e-10)
    if codim == length(λ)
        ED_stress_system = projective_stress_system
    else
        rand_pt = randn(Float64, length(λ))
        ED_matrix = hcat(length(stress_poly_system)==1 ? differentiate(stress_poly_system, λ)' : differentiate(stress_poly_system, λ), λ - rand_pt)
        ED_stress_system = vcat(projective_stress_system, minors(ED_matrix, codim+1))
    end
    sols = real_solutions(solve(ED_stress_system))
    return isempty(sols) ? [] : sols[1]
end


"""
    minors(A, k)

Compute the (k+1)x(k+1) minors of the matrix `A`
"""
function minors(A, k)
    n, m = size(A)
    rowsets = collect(combinations(1:n, k))
    colsets = collect(combinations(1:m, k))

    result = Dict()
    for r in rowsets, c in colsets
        result[(r, c)] = det(A[r, c])
    end
    return values(result)
end
