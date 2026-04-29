export  read_realizations,
        save_realizations,
        save_to_Houdini,
        minors

"""
    save_to_Houdini(points)

Save the `points` to a `.poly` file with the name `filename` which can later be opened in SideFX Houdini. 
"""
function save_to_Houdini(points::Vector, filename::String)
    all(pt->pt isa <:Vector{<:Real} && length(pt)==3, points) || throw(error("All points need to be vectors of length 3."))
    open("$(filename).poly","w") do file
        write(file, "POINTS\n")
        for i in 1:length(points)
            write(file, "$(i): $(points[i][1]) $(points[i][2]) $(points[i][3])\n")
        end
        write(file, "POLYS\nEND")
    end
end


"""
    save_realizations(D, filename)

Saves the `motion_matrices` from the `DeformationPath` `D` in a `.txt` file called `<filename>.txt`.
"""
function save_realizations(D::DeformationPath, filename::String)
    open("$(filename).txt","w") do file
        for mat in D.motion_matrices
            write(file, "[")
            for row in 1:size(mat)[1]
                for val in mat[row,1:end-1]
                    write(file, "$(val) ")
                end
                if row<size(mat)[1]
                    write(file, "$(mat[row,end]); ")
                else
                    write(file, "$(mat[row,end])]\n")
                end
            end
        end
    end
end


"""
    read_realizations(G, filename)

Reads the `motion_matrices` for the `ConstraintSystem` `G` from a `.txt` file called `<filename>.txt` and returns a `DeformationPath`.
"""
function read_realizations(G::ConstraintSystem, filename::String; kwargs...)::DeformationPath
    !isfile("$(filename).txt") && throw(error("A file with the name `$(filename).txt` does not exist."))
    motion_matrices = Vector{Matrix{Float64}}([])
    open("$(filename).txt","r") do file
        while !eof(file)  
            s = readline(file)[2:end-1]
            rows = split(s, "; ")
            row_array=[]
            _realization = Base.copy(G.realization)
            for row in rows
                col = split(row, " ")
                col_array = [parse(Float64, col[i]) for i in eachindex(col)]
                push!(row_array, col_array)
            end
            if length(row_array)!=G.dimension || !all(row->length(row)==size(G.realization)[2] && length(row)==length(row_array[1]), row_array)
                throw(error("The dimensions of the `row_array` do not match the given geometric constraint system!"))
            end
            for row in eachindex(row_array)
                _realization[row,:] .= row_array[row]
            end
            push!(motion_matrices, _realization)
        end
    end
    return DeformationPath(G, motion_matrices; kwargs...)
end


"""
    read_realizations(F, filename)

Reads the `motion_matrices` for the geometric constraint system `F` from a `.txt` file called `<filename>.txt` and returns a `DeformationPath`.
"""
function read_realizations(F::AllTypes, filename::String; kwargs...)::DeformationPath
    read_realizations(F.G, filename; kwargs...)
end


"""
    minors(A, k)

Compute the (k+1)x(k+1) minors of the matrix `A`.
"""
function minors(A, k)::Vector
    n, m = size(A)
    rowsets = collect(combinations(1:n, k))
    colsets = collect(combinations(1:m, k))

    result = Dict()
    for r in rowsets, c in colsets
        result[(r, c)] = det(A[r, c])
    end
    return collect(values(result))
end
