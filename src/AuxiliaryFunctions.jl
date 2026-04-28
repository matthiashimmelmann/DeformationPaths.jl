"""
    save_to_Houdini(points)

Save the `points` to a `.poly` file with the name `filename` which can later be opened in SideFX Houdini. 
"""
function save_to_Houdini(points::Vector{Vector}, filename::String)
    all(pt->length(pt)==3 && pt isa Vector{<:Real}, points) || throw(error("All points need to be vectors of length 3."))
    open("$(filename).poly","w") do file
        write(file, "POINTS\n")
        for i in 1:length(points)
            write(file, "$(i): $(points[i][1]) $(points[i][2]) $(points[i][3])\n")
        end
        write(file, "POLYS\nEND")
    end
end