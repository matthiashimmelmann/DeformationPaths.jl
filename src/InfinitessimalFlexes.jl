export  compute_inf_flexes,
        compute_trivial_inf_flexes,
        compute_equilibrium_stresses,
        compute_nontrivial_inf_flexes,
        compute_nonblocked_flex

"""
    compute_inf_flexes(G, point[; tol])

Compute all infinitesimal flexes of a geometric constraint system `G` in `point`.
"""
function compute_inf_flexes(G::ConstraintSystem, point::Vector{<:Real}; tol::Real=1e-8)::Matrix{<:Real}
    inf_flexes = nullspace(evaluate(G.jacobian, G.variables=>point); atol=tol)
    return inf_flexes
end
function compute_inf_flexes(G::ConstraintSystem; tol::Real=1e-8)::Matrix{<:Real}
    return compute_inf_flexes(G, to_Array(G, G.realization); tol=tol)
end
function compute_inf_flexes(F::AllTypes; tol::Real=1e-8)::Matrix{<:Real}
    return compute_inf_flexes(F.G; tol=tol)
end


"""
    compute_trivial_inf_flexes(G, point[; tol])

Compute all trivial infinitesimal flexes of a geometric constraint system `G` in `point`
"""
function compute_trivial_inf_flexes(G::ConstraintSystem, point::Vector{<:Real}; tol::Real=1e-8)::Matrix{<:Real}
    dim = G.dimension
    if length(G.vertices)<size(G.realization)[2]
        K_n = ConstraintSystem(G.vertices, G.variables, vcat(G.equations, [sum( (G.xs[:,bar[1]]-G.xs[:,bar[2]]) .^2) - sum( (G.realization[:,bar[1]]-G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(G.vertices) for j in eachindex(G.vertices) if i<j]]), G.realization, G.xs; pinned_GCS=G.pinned_GCS, pinned_vertices=G.pinned_vertices)
        jacobian = evaluate.(K_n.jacobian, K_n.variables=>to_Array(K_n, K_n.realization))
        basis = nullspace(jacobian; atol=tol)
    else
        translations = [
            vcat([[i!=j ? 0 : 1 for i in 1:dim] for _ in 1:length(G.vertices)]...) for j in 1:dim
        ]
        basis_skew_symmetric = []
        for i in 2:dim
            for j in 1:i-1
                A = zeros(dim, dim)
                A[i, j] = 1
                A[j, i] = -1
                push!(basis_skew_symmetric, A)
            end
        end
        inf_rot = [
            vcat([A * G.realization[:,v] for v in 1:length(G.vertices)]...)
            for A in basis_skew_symmetric
        ]
        matrix_inf_flexes = hcat(vcat(translations, inf_rot)...)
        projector = Matrix(I, length(G.vertices)*dim, length(G.vertices)*dim)
        projector = projector[:, vcat([[dim*(i-1)+j for j in 1:dim] for i in G.vertices if !(i in G.pinned_vertices)]...)]
        M = hcat(matrix_inf_flexes, -projector)
        Zerospace = nullspace(M; atol=tol)
        r = size(matrix_inf_flexes,2)
        X = Zerospace[1:r, :]
        basis = matrix_inf_flexes*X
        basis = basis[vcat([[dim*(i-1)+j for j in 1:dim] for i in G.vertices if !(i in G.pinned_vertices)]...),:]
    end
    return basis
end
function compute_trivial_inf_flexes(G::ConstraintSystem; tol::Real=1e-8)::Matrix{<:Real}
    return compute_trivial_inf_flexes(G, to_Array(G, G.realization); tol=tol)
end
function compute_trivial_inf_flexes(F::AllTypes; tol::Real=1e-8)::Matrix{<:Real}
    return compute_trivial_inf_flexes(F.G; tol=tol)
end


"""
    compute_equilibrium_stresses(G, point[; tol])

Compute all equilibrium stresses of a geometric constraint system `G` in `point`.
"""
function compute_equilibrium_stresses(G::ConstraintSystem, point::Vector{<:Real}; tol::Real=1e-8)::Matrix{<:Real}
    stresses = nullspace(evaluate(G.jacobian, G.variables=>point)'; atol=tol)
    return stresses
end
function compute_equilibrium_stresses(G::ConstraintSystem; tol::Real=1e-8)::Matrix{<:Real}
    return compute_equilibrium_stresses(G, to_Array(G, G.realization); tol=tol)
end
function compute_equilibrium_stresses(F::AllTypes; tol::Real=1e-8)::Matrix{<:Real}
    return compute_equilibrium_stresses(F.G; tol=tol)
end


"""
    compute_nontrivial_inf_flexes(G, point[; tol])

Compute the nontrivial infinitesimal flexes of a geometric constraint system `G` in `point`.
"""
function compute_nontrivial_inf_flexes(G::ConstraintSystem, point::Vector{<:Real}; tol::Real=1e-8)::Matrix{<:Real}
    inf_flexes = compute_inf_flexes(G, point; tol)
    trivial_inf_flexes = compute_trivial_inf_flexes(G, point; tol)
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
function compute_nontrivial_inf_flexes(G::ConstraintSystem; tol::Real=1e-8)::Matrix{<:Real}
    return compute_nontrivial_inf_flexes(G, to_Array(G, G.realization); tol=tol)
end
function compute_nontrivial_inf_flexes(F::AllTypes; tol::Real=1e-8)::Matrix{<:Real}
    return compute_nontrivial_inf_flexes(F.G; tol=tol)
end


"""
    compute_nonblocked_flex(F[; tol_rank_drop, tol])

Compute an infinitesimal flex of `F` that is not blocked by an equilibrium stress.
"""
function compute_nonblocked_flex(F::AllTypes; fast_search::Bool=false, tol_rank_drop::Real=1e-6, tol::Real=1e-12)::Vector
    flexes = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization); tol=tol_rank_drop)
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
    sols = real_solutions(solve(Vector{Expression}(ED_stress_system)))    
    return isempty(sols) ? [] : sols[1]
end