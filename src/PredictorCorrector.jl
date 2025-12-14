"""
    euler_step(G, step_size, prev_flex, point, K_n)

Euler step predicting the next point along the approximate motion.

# Returns
- `predicted_point::Vector{<:Real}`: The next point predicted by Euler's method.
- `predicted_inf_flex::Vector{<:Real}`: The tangent vector predicted by Euler's method.
"""
function euler_step(G::ConstraintSystem, step_size::Real, prev_flex::Vector{<:Real}, point::Vector{<:Real}, K_n::ConstraintSystem; tol=1e-5)::Tuple{Vector{<:Real}, Vector{<:Real}}
    flex_space = compute_nontrivial_inf_flexes(G, point, K_n; tol=tol)
    if size(flex_space)[2]==0
        throw("The space of nontrivial infinitesimal motions is empty.")
    end
    flex_coefficients = pinv(flex_space) * prev_flex
    predicted_inf_flex = sum(flex_space[:,i] .* flex_coefficients[i] for i in eachindex(flex_coefficients))
    predicted_inf_flex = predicted_inf_flex ./ norm(predicted_inf_flex)
    return point+step_size*predicted_inf_flex, predicted_inf_flex
end


"""
    newton_correct(G, point)

Apply Newton's method to correct `point` back to the constraints in `G`.

# Arguments
- `G::ConstraintSystem`: The underlying geometric constraint system.
- `point::Vector{<:Real}`: The initial point that Newton's method is applied to.

For further keywords, see [`newton_correct`](@ref).

# Returns
- `q::Vector{<:Real}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`
"""
function newton_correct(G::ConstraintSystem, point::Vector{<:Real}; kwargs...)::Vector{<:Real}
    return newton_correct(G.equations, G.variables, G.jacobian, point; kwargs...)
end

"""
    newton_correct(equations, variables, jac, point[; tol, time_penalty])

Apply Newton's method to correct `point` back to the constraints in `equations`.

# Arguments
- `equations::Vector{Expression}`: Equations to correct `point` to.
- `variables::Vector{Variable}`: Variables from the affine coordinate ring.
- `jac::Matrix{Expression}`: Jacobian matrix corresponding to `equations` and `variables`.
- `point::Vector{<:Real}`: The initial point that Newton's method is applied to.
- `tol::Real` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Union{Real,Nothing}` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.
- `armijo_linesearch::Bool` (optional): Flag for activating the Armijo backtracking line search procedure. Default value: `true`.

# Returns
- `q::Vector{<:Real}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`
"""
function newton_correct(equations::Vector{Expression}, variables::Vector{Variable}, jac::Matrix{Expression}, point::Vector{<:Real}; tol::Real = 1e-13, armijo_linesearch::Bool=true, time_penalty::Union{Real,Nothing}=2)::Vector{<:Real}
    #TODO needs work
    q = Base.copy(point)
    start_time=Base.time()
    if armijo_linesearch
        while(norm(evaluate(equations, variables=>q)) > tol)
            J = evaluate.(jac, variables=>q)
            stress_dimension = size(nullspace(J'; atol=1e-8))[2]
            if stress_dimension > 0
                rand_mat = randn(Float64, length(equations) - stress_dimension, length(equations))
                new_equations = rand_mat*equations
                J = rand_mat*J
            else
                new_equations = equations
            end

            #Armijo Line Search
            r_val = evaluate(new_equations, variables=>q)
            v = -(J \ r_val)
            global damping, damping_too_small = 0.9, 0
            qnew = q + damping*v
            while norm(evaluate(equations, variables=>qnew)) > norm(evaluate(equations, variables=>q)) - 0.5 * damping * v'*v
                global damping = damping*0.7
                qnew = q + damping*v
                if damping < 1e-12
                    global damping_too_small += 1
                    qnew = q + 0.025*v
                    break
                end
                if damping_too_small >= 3 || (!isnothing(time_penalty) && Base.time()-start_time > length(point)/time_penalty)
                    throw("Newton's method did not converge in time. damping=$damping and time=$(Base.time()-start_time)")
                end
            end
            if damping >= 1e-1
                global damping_too_small = maximum([0, damping_too_small-1])
                q = qnew
            end
        end
    else
        global damping = 0.15
        while(norm([eq(variables=>q) for eq in equations]) > tol)
            J = evaluate.(jac, variables=>q)
            stress_dimension = size(nullspace(J'; atol=1e-8))[2]
            if stress_dimension > 0
                rand_mat = randn(Float64, length(equations) - stress_dimension, length(equations))
                new_equations = rand_mat*equations
                J = rand_mat*J
            else
                new_equations = equations
            end

            qnew = q - damping * (J \ evaluate.(new_equations, variables=>q))
            if norm(evaluate(equations, variables=>qnew)) < norm(evaluate(equations, variables=>q))
                global damping = damping*1.2
            else
                global damping = damping/2
            end
            if damping > 1
                global damping = 1
            end
            if damping < 1e-12 || (!isnothing(time_penalty) && Base.time()-start_time > length(point)/time_penalty)
                throw("Newton's method did not converge in time. damping=$damping and time=$(Base.time()-start_time)")
            end
            q = qnew
        end

    end
    return q
end

"""
    symmetric_newton_correct(G, point[; tol, time_penalty])

Apply symmetric Newton's method to correct `point` back to the constraints in `equations`.

The symmetric Newton corrector evaluates the Jacobian matrix less often.

# Arguments
- `G::ConstraintSystem`: The underlying geometric constraint system.
- `point::Vector{<:Real}`: The initial point that Newton's method is applied to.
- `tol::Real` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Union{Real,Nothing}` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.

# Returns
- `q::Vector{<:Real}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`

See also [`symmetric_newton_correct`](@ref)
"""
function symmetric_newton_correct(G::ConstraintSystem, point::Vector{<:Real}; tol::Real = 1e-13, time_penalty::Union{Real,Nothing}=2)::Vector{<:Real}
    return symmetric_newton_correct(G.equations, G.variables, G.jacobian, point; tol = tol, time_penalty=time_penalty)
end


"""
    symmetric_newton_correct(equations, variables, jac, point[; tol, time_penalty])

Apply symmetric Newton's method to correct `point` back to the constraints in `equations`.

The symmetric Newton corrector evaluates the Jacobian matrix less often.

# Arguments
- `equations::Vector{Expression}`: Equations to correct `point` to.
- `variables::Vector{Variable}`: Variables from the affine coordinate ring.
- `jac::Matrix{Expression}`: Jacobian matrix corresponding to `equations` and `variables`.
- `point::Vector{<:Real}`: The initial point that Newton's method is applied to.
- `tol::Real` (optional): Numerical tolerance that is used as a stopping criterion for Newton's method. Default value: `1e-13`.
- `time_penalty::Union{Real,Nothing}` (optional): If Newton's method takes too long, we stop the iteration and throw an error. Here, "too long" is measured in terms of `length(point)/time_penalty` seconds. Default value: `2`.

# Returns
- `q::Vector{<:Real}`: A point `q` such that the Euclidean norm of the evaluated equations is at most `tol`
"""
function symmetric_newton_correct(equations::Vector{Expression}, variables::Vector{Variable}, jacobian::Matrix{Expression}, p::Vector{<:Real}; tol::Real = 1e-13, time_penalty::Union{Real,Nothing}=2)::Vector{<:Real}
    global _q = Base.copy(p)
    global qnew = _q
    global damping = 0.2
    J = Matrix{Float64}(evaluate.(jacobian, variables=>_q))
    new_equations, J_new = Base.copy(equations), Base.copy(J)
    index = 0
    # Randomize the linear system of equations
    while length(equations)>0 && norm(evaluate.(equations, variables=>_q), Inf) > tol
        #println("sym: $damping\t $(norm(evaluate(eqns, vars=>q)))")
        if index%6 == 0
            stress_dimension = size(nullspace(J'; atol=1e-8))[2]
            if stress_dimension > 0
                rand_mat = randn(Float64, length(equations) - stress_dimension, length(equations))
                new_equations = rand_mat*equations
                J_new = rand_mat*J
            else
                new_equations = equations
                J_new = J
            end
        end
        #println("$damping\t $(norm(evaluate(equations, variables=>_q)))")
        #qnew = _q - damping * (J_new \ evaluate(new_equations, variables=>_q))
        r_val = evaluate(new_equations, variables=>_q)
        v = -(J_new \ r_val)
        qnew = _q + damping*v
        if norm(evaluate(equations, variables=>qnew), Inf) < norm(evaluate(equations, variables=>_q), Inf)
			global damping = damping*1.2
		else
			global damping = damping/2
		end
		if damping > 1
			global damping = 1
        elseif damping < 1e-12  || (!isnothing(time_penalty) && Base.time()-start_time > length(point)/time_penalty)
            throw("Newton's method did not converge in time. damping=$damping and time=$(Base.time()-start_time)")
        end

        _q = qnew
        index=index+1
    end
    return _q
end

