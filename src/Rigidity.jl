export  is_rigid,
        is_inf_rigid,
        is_second_order_rigid,
        is_prestress_stable,
        coned_rigidity_phase_space

"""
    is_rigid(F[; tol, tol, tested_random_flexes, symmetric_newton])

Heuristically checks if a geometric constraint system `F` is (continuously) rigid. 

See [`DeformationPath(G::DeformationPaths.ConstraintSystem, motion_samples::Vector{<:Vector{<:Real}})`](@ref) for a description of the possible parameters.
"""
function is_rigid(F::AllTypes; tol_rank_drop::Real=1e-3, tol::Real=1e-10, tested_random_flexes::Int=8, symmetric_newton::Bool=false, time_penalty::Union{Real,Nothing}=2)::Bool
    #TODO needs work
    if is_inf_rigid(F; tol_rank_drop=tol_rank_drop)
        return true
    end
    for _ in 1:tested_random_flexes
        D = DeformationPath(F, [], 5; show_progress=false, step_size=tol_rank_drop, tol=tol, random_flex=true, symmetric_newton=symmetric_newton, time_penalty=time_penalty)
        if any(sample->norm(sample-D.motion_samples[1], Inf)>tol_rank_drop*0.1, D.motion_samples)
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
    inf_flexes = compute_inf_flexes(F; tol=tol_rank_drop)
    trivial_inf_flexes = compute_trivial_inf_flexes(F; tol=tol_rank_drop)
    return length(inf_flexes) == length(trivial_inf_flexes)
end


"""
    is_second_order_rigid(F)

Checks if a geometric constraint system `F` is second-order rigid.

See also [`compute_nonblocked_flex`](@ref) for the possible keywords.
"""
function is_second_order_rigid(F::AllTypes; kwargs...)::Bool
    nonblocked_flex = compute_nonblocked_flex(F; fast_search=false, kwargs...)
    if isempty(nonblocked_flex)
        return true
    else
        return false
    end
end


"""
    is_prestress_stable(F)

Checks if a geometric constraint system `F` is prestress stable.
"""
function is_prestress_stable(F::AllTypes; tol_rank_drop::Real=1e-6, tol::Real=1e-10)::Bool
    flexes = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization); tol=tol_rank_drop)
    if size(flexes)[2]==0
        return true
    end
    rigidity_matrix = evaluate.(F.G.jacobian, F.G.variables=>to_Array(F, F.G.realization))
    stresses = nullspace(rigidity_matrix'; atol=tol_rank_drop)
    if size(stresses)[2]==0
        return false
    end
    @var λ[1:size(flexes)[2]] ω[1:size(stresses)[2]] μ[1:size(flexes)[2]]
    parametrized_flex = flexes*λ
    parametrized_stress = stresses*ω
    if size(flexes)[2]==1
        if all(coef -> isapprox(coef,0,atol=tol), coefficients(parametrized_stress'*(evaluate.(F.G.jacobian, F.G.variables=>flexes[:,1])*flexes[:,1]), ω))
            return false
        else
            return true
        end
    end
    if size(stresses)[2]==1
        Hessian = evaluate(differentiate(differentiate(stresses[:,1]'*evaluate.(F.G.jacobian, F.G.variables=>parametrized_flex)*parametrized_flex, λ), λ), λ=>[0 for _ in 1:length(λ)])
        if all(ev->ev>tol, real.(eigvals(Hessian))) || all(ev->ev<-tol, real.(eigvals(Hessian)))
            return true
        else
            return false
        end
    end

    stress_energy = parametrized_stress'*evaluate.(F.G.jacobian, F.G.variables=>Vector{Expression}(parametrized_flex))*parametrized_flex
    Hessian = differentiate(differentiate(stress_energy, λ), λ)
    matrices = [[evaluate(differentiate(Hessian[j,k], [ω[i]])[1], [ω[i]]=>[0.]) for j in axes(Hessian,1), k in axes(Hessian,2)] for i in eachindex(ω)]
    #=for matrix in matrices
        string_output = "{"
        for i in axes(matrix,1)
            string_output *= "{"
            for j in axes(matrix,2)
                string_output *= "$(matrix[i,j])"
                if j != size(matrix)[2]
                    string_output *= ","
                end
            end
            string_output *= "}"
            if i != size(matrix)[1]
                string_output *= ","
            end
        end
        string_output *= "}"
    end=#

    #INFO Test needs work
    return any(matrix->all(ev->ev>tol, real.(eigvals(matrix))), matrices) || any(matrix->all(ev->ev<-tol, real.(eigvals(matrix))), matrices)
end


"""
    coned_rigidity_phase_space(P, start_points, end_points[; check, discretization_size, show_progress])

Compute the rigidity phase space of the coned 3-dimensional polytope `P` (or alternatively a bar-joint framework `P`).
Depending on the `check`` keyword, we check the rigidity for all points contained in the cuboid given by the corner points `start_points` and `end_points`.
The `check` keyword is a Symbol that can take the values of 
- `:FOR`, in which case the method [`is_inf_rigid`](@ref) is called,
- `:PSS`, in which case the method [`is_prestress_stable`](@ref) is called,
- `:SOR`, in which case [`is_second_order_rigid`](@ref) is called and
- `:RIG`, in which case [`is_rigid`](@ref) is called.
The keyword's default value is `:PSS`.
"""
function coned_rigidity_phase_space(P::Union{Polytope,Framework}, start_points::Vector{<:Real}, end_points::Vector{<:Real}; check::Symbol=:PSS, discretization_size::Real=0.1, show_progress::Bool=true)
    cone_point = maximum(P.G.vertices)+1
    rigid_points = Vector{Vector{Float64}}([])
    if length(start_points)!=P.G.dimension || length(end_points)!=P.G.dimension || length(start_points)!=length(end_points) || !(P.G.dimension in [2,3])
        throw(error("The length of `start_points` and `end_points` both need to be equal to 2 or 3."))
    end
    if !(check in [:FOR, :PSS, :SOR, :RIG])
        throw(error("The `check` keyword needs to be either `:FOR` for infinitesimal rigidity, `:PSS` for prestress stability, `:SOR` for second-order rigidity or `:RIG` for (continuous) rigidity."))
    end

    F = Framework(vcat(P.edges,[(i, cone_point) for i in P.G.vertices]), hcat(P.G.realization[:,1:length(P.G.vertices)], P.G.dimension==2 ? [0,0] : [0,0,0]))
    discretization = P.G.dimension==2 ? [[x,y] for x in start_points[1]:discretization_size:end_points[1], y in start_points[2]:discretization_size:end_points[2]] : [[x,y,z] for x in start_points[1]:discretization_size:end_points[1], y in start_points[2]:discretization_size:end_points[2], z in start_points[3]:discretization_size:end_points[3]]
    @showprogress enabled=show_progress for pt in discretization
        F.G.realization[:,end] .= pt
        if (check==:FOR && is_inf_rigid(F)) || (check==:PSS && is_prestress_stable(F)) || (check==:SOR && is_second_order_rigid(F)) || (check==:RIG && is_rigid(F))
            push!(rigid_points, pt)
        end
    end
    return rigid_points
end