"""
Class for Constructing a general constraint system.

# Attributes
- `vertices::Vector{Int}`: Vertices describing the geometric constraint system. They are given as a list of integers.
- `variables::Vector{Variable}`: Coordinate variables representing the positions of the `vertices`.
- `equations::Vector{Expression}`: Equations corresponding to the constraint system.
- `realization::Matrix{<:Real}`: Realization of the geometric constraint system satisfying the `equations`.
- `jacobian::Matrix{Expression}`: Jacobian matrix for the `equations` and `variables`.
- `dimension::Int`: Dimension in which the geometric constraint system lives.
- `xs::Union{Matrix{Variable}, Matrix{Expression}}`: Matrix representing the possible realizations of a geometric constraint system.
- `system::System`: Polynomial system consisting of the `equations` and `variables` in the form of `HomotopyContinuation`.
- `pinned_vertices::Vector{Int}`: Pinned vertices of the system. These vertices remain unchanged and are added as constraints. Pinning can, for instance, be used to factor out rigid motions.
"""
mutable struct ConstraintSystem
    vertices::Vector{Int}
    variables::Vector{Variable}
    equations::Vector{Expression}
    realization::Matrix{<:Real}
    jacobian::Matrix{Expression}
    dimension::Int
    xs::Union{Matrix{Variable},Matrix{Expression}}
    system::System
    pinned_vertices::Vector{Int}
end

"""
    ConstraintSystem(vertices::Vector{Int}, variables::Vector{Variable}, equations::Vector{Expression}, realization::Matrix{<:Real}, xs)

Constructor of a `ConstraintSystem` object.
"""
function ConstraintSystem(vertices::Vector{Int}, variables::Vector{Variable}, equations::Vector{Expression}, realization::Matrix{<:Real}, xs; pinned_vertices::Vector{Int}=Vector{Int}([]))::ConstraintSystem
    jacobian = hcat([differentiate(eq, variables) for eq in equations]...)'
    dimension = size(realization)[1]
    size(realization)[1]==dimension && (size(realization)[2]==length(vertices) || size(realization)[2]==length(variables)//dimension+length(pinned_vertices)) || size(realization)[2]==size(xs)[2] || throw("The realization does not have the correct format.")
    size(xs)[1]==size(realization)[1] && size(xs)[2]==size(realization)[2] || throw("The matrix 'xs' does not have the correct format.")
    ConstraintSystem(vertices, variables, equations, realization, jacobian, dimension, xs, System(equations; variables=variables), pinned_vertices)
end


function Base.:(==)(G1::ConstraintSystem, G2::ConstraintSystem)
    """
    Overloads the equality operator for `ConstraintSystem`.
    """
    return G1.vertices==G2.vertices && length(G1.variables)==length(G2.variables) &&  all(i->G1.variables[i]==G2.variables[i], eachindex(G1.variables)) && length(G1.equations)==length(G2.equations) && all(i->G1.equations[i]==G2.equations[i], eachindex(G1.equations))
end


function Base.show(io::IO, G::ConstraintSystem)::Nothing
    """
    Display method for the base clase `ConstraintSystem`.
    """
    print(io,"Constraint System:\n")
    print(io,"\t\tVertices:\t$(G.vertices)\n")
    print(io,"\t\tEquations:\t[$(join(G.equations[1:min(2,length(G.equations))], ", "))$(length(G.equations)<=2 ? "" : ", ...")]\n")
    print(io, "\t\tRealization:\t")
    for i in 1:min(4,size(G.realization)[2])
        if size(G.realization)[1]==1
            print(io, (i==1 ? "" : "\t\t\t\t")*"$(G.realization[1,i])\n")
        elseif size(G.realization)[1]==2
            print(io, (i==1 ? "" : "\t\t\t\t")*"$(G.realization[1,i])\t$(G.realization[2,i])\n")
        elseif size(G.realization)[1]==3
            print(io, (i==1 ? "" : "\t\t\t\t")*"$(G.realization[1,i])\t$(G.realization[2,i])\t$(G.realization[3,i])\n")
        else
            print(io, (i==1 ? "" : "\t\t\t\t")*"$(G.realization[1,i])\t$(G.realization[2,i])\t$(G.realization[3,i]) ...\n")
        end
    end
    if !(isempty(G.pinned_vertices))
        print(io, "\tPinned Vertices: $(G.pinned_vertices)")
    end
    return nothing
end

#TODO Allow random realizations

"""
    Framework([vertices,] bars, realization[; pinned_vertices])

Class for bar-joint frameworks.
"""
mutable struct Framework
    G::ConstraintSystem
    bars::Vector{Tuple{Int,Int}}

    function Framework(vertices::Vector{Int}, bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        all(t->length(t)==2, bars) && all(bar->bar[1] in vertices && bar[2] in vertices, bars) || throw("The bars don't have the correct format.")
        realization = Float64.(realization)
        dimension = size(realization)[1]
        all(v->v in vertices, pinned_vertices) || throw("Some of the pinned_vertices are not contained in vertices.")
        bars = [bar[1]<=bar[2] ? Tuple(bar) : Tuple([bar[2],bar[1]]) for bar in bars]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw("The realization does not have the correct format. Length(vertices): $(length(vertices)); realization: $(size(realization)[2])")
        dimension>=1 || throw("The dimension is not an integer bigger than 0.")
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        bar_equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices,variables, bar_equations, realization, xs; pinned_vertices=pinned_vertices)
        new(G, bars)
    end

    function Framework(bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([bar[1] for bar in bars], [bar[2] for bar in bars]))))
        Framework(vertices, bars, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    AngularFramework([vertices,] angles, realization[; pinned_vertices])

Class for angular bar-joint frameworks.
"""
mutable struct AngularFramework
    G::ConstraintSystem
    bars::Vector{Tuple{Int,Int}}
    angles::Vector{Tuple{Int,Int,Int}}

    function AngularFramework(vertices::Vector{Int}, angles::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        all(t->length(t)==3, angles) && all(angle->angle[1] in vertices && angle[2] in vertices && angle[3] in vertices, angles) || throw("The angles don't have the correct format.")
        realization = Float64.(realization)
        dimension = size(realization)[1]
        all(v->v in vertices, pinned_vertices) || throw("Some of the pinned_vertices are not contained in vertices.")
        angles = [(angle[1],angle[2],angle[3]) for angle in angles]
        bars = vcat([[angle[1], angle[2]] for angle in angles], [[angle[2], angle[3]] for angle in angles])
        bars = [bar[1]<=bar[2] ? Tuple(bar) : Tuple([bar[2],bar[1]]) for bar in bars]
        bars = collect(Set(bars))
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw("The realization does not have the correct format.")
        dimension>=1 || throw("The dimension is not an integer bigger than 0.")
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        angle_constraints = [((xs[:,angle[1]]-xs[:,angle[2]])'*(xs[:,angle[3]]-xs[:,angle[2]])) * sqrt(sum((realization[:,angle[1]]-realization[:,angle[2]]).^2)*sum((realization[:,angle[3]]-realization[:,angle[2]]).^2)) - sqrt(sum((xs[:,angle[1]]-xs[:,angle[2]]).^2)*sum((xs[:,angle[3]]-xs[:,angle[2]]).^2)) * ((realization[:,angle[1]]-realization[:,angle[2]])'*(realization[:,angle[3]]-realization[:,angle[2]])) for angle in angles]
        G = ConstraintSystem(vertices,variables, angle_constraints, realization, xs; pinned_vertices=pinned_vertices)
        new(G, bars, angles)
    end

    function AngularFramework(angles::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([angle[1] for angle in angles], [angle[2] for angle in angles], [angle[3] for angle in angles]))))
        AngularFramework(vertices, angles, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    FrameworkOnSurface([vertices,] bars, realization, surface[; pinned_vertices])

Class for bar-joint frameworks constrained to a surface.
"""
mutable struct FrameworkOnSurface
    G::ConstraintSystem
    bars::Vector{Tuple{Int,Int}}
    surface::Function

    function FrameworkOnSurface(vertices::Vector{Int}, bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Matrix{<:Real}, surface::Function; pinned_vertices=Vector{Int}([]))
        try
            surface([1,1,1])
        catch
            throw("The specified implicit surface function does not take 3 arguments.")
        end
        realization = Float64.(realization)
        dimension = size(realization)[1]
        dimension==3 || throw("The dimension for FrameworkOnSurface needs to be equal to 3.")
        
        F = Framework(vertices, bars, realization; pinned_vertices=pinned_vertices)
        G = F.G
        add_equations!(G, [surface(G.xs[:,i]) for i in eachindex(G.vertices)])
        new(G, F.bars, surface)
    end

    function FrameworkOnSurface(bars::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, realization::Matrix{<:Real}, surface::Function; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([bar[1] for bar in bars], [bar[2] for bar in bars]))))
        FrameworkOnSurface(vertices, bars, realization, surface; pinned_vertices=pinned_vertices)
    end
end


"""
    SpherePacking([vertices,] radii, realization[; pinned_vertices])

Class for sticky sphere packings.
"""
mutable struct SpherePacking
    G::ConstraintSystem
    contacts::Vector{Tuple{Int,Int}}
    radii::Vector{<:Real}
    tolerance::Real

    function SpherePacking(vertices::Vector{Int}, radii::Vector{<:Real}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]), tolerance::Real=1e-6)
        length(vertices)==length(radii) && length(radii)==size(realization)[2] && all(r->r>0, radii) || throw("The length of the radii does not match the length of the vertices or the dimensionality of the realization.")
        all(v->v in vertices, pinned_vertices) || throw("Some of the pinned_vertices are not contained in vertices.")
        realization = Float64.(realization)
        dimension = size(realization)[1]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw("The realization does not have the correct format.")
        all(t->norm(realization[:,t[1]]-realization[:,t[2]]) >= radii[t[1]]+radii[t[2]]-tolerance, powerset(length(vertices),2,2)) || throw("Some of the disks are too close")
        contacts = [Tuple([i,j]) for i in eachindex(vertices) for j in i+1:length(vertices) if isapprox(norm(realization[:,i]-realization[:,j]),radii[i]+radii[j],atol=tolerance)]
        dimension in [2,3] || throw("The dimension for SpherePacking must be 2.")
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - (radii[bar[1]]+radii[bar[2]])^2 for bar in contacts]
        bar_equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices, variables, bar_equations, realization, xs; pinned_vertices=pinned_vertices)
        new(G, contacts, radii, tolerance)
    end

    function SpherePacking(radii::Vector{<:Real}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        vertices = [i for i in eachindex(radii)]
        SpherePacking(vertices, radii, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    SphericalDiskPacking([vertices,] contacts, [inversive_distances,] realization[; pinned_vertices])

Class for spherical disk packings in Minkowski space.
"""
mutable struct SphericalDiskPacking
    G::ConstraintSystem
    contacts::Vector{Tuple{Int,Int}}
    inversive_distances::Vector{<:Real}

    function SphericalDiskPacking(vertices::Vector{Int}, contacts::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}}, inversive_distances::Vector{<:Real}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]), tolerance::Real=1e-8)
        length(contacts)==length(inversive_distances) || throw("The length of the inversive distances does not match the length of the vertices or the dimensionality of the realization.")
        all(v->v in vertices, pinned_vertices) || throw("Some of the pinned_vertices are not contained in vertices.")
        realization = Float64.(realization)
        dimension = size(realization)[1]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw("The realization does not have the correct format.")
        dimension==3 || throw("The dimension for SphericalDiskPacking must be 3.")
        all(i->isapprox(minkowski_scalar_product(realization[:,contacts[i][1]], realization[:,contacts[i][2]])/sqrt(minkowski_scalar_product(realization[:,contacts[i][1]], realization[:,contacts[i][1]]) * minkowski_scalar_product(realization[:,contacts[i][2]], realization[:,contacts[i][2]])), inversive_distances[i], atol=tolerance), 1:length(contacts)) || throw("The Minkowski distances do not match the given realization.")

        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        inversive_distance_equation = [minkowski_scalar_product(xs[:,contacts[i][1]], xs[:,contacts[i][2]])^2 - inversive_distances[i]^2 * minkowski_scalar_product(xs[:,contacts[i][1]], xs[:,contacts[i][1]]) * minkowski_scalar_product(xs[:,contacts[i][2]], xs[:,contacts[i][2]]) for i in eachindex(contacts)]
        inversive_distance_equation = filter(eq->eq!=0, inversive_distance_equation)
        G = ConstraintSystem(vertices, variables, inversive_distance_equation, realization, xs; pinned_vertices=pinned_vertices)
        new(G, contacts, inversive_distances)
    end

    function SphericalDiskPacking(contacts::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}}, inversive_distances::Vector{<:Real}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([bar[1] for bar in contacts], [bar[2] for bar in contacts]))))
        SphericalDiskPacking(vertices, contacts, inversive_distances, realization; pinned_vertices=pinned_vertices)
    end

    function SphericalDiskPacking(contacts::Union{Vector{Tuple{Int,Int}},Vector{Vector{Int}}}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        inversive_distances = [minkowski_scalar_product(realization[:,contact[1]], realization[:,contact[2]])/sqrt(minkowski_scalar_product(realization[:,contact[1]], realization[:,contact[1]]) * minkowski_scalar_product(realization[:,contact[2]], realization[:,contact[2]])) for contact in contacts]
        SphericalDiskPacking(contacts, inversive_distances, realization; pinned_vertices=pinned_vertices)
    end

    minkowski_scalar_product(e1,e2) = e1'*e2-1
end


"""
    VolumeHypergraph([vertices,] volumes, realization[; pinned_vertices])

Class for volume-constrained hypergraphs.
"""
mutable struct VolumeHypergraph
    G::ConstraintSystem
    volumes::Vector{Vector{Int}}

    function VolumeHypergraph(vertices::Vector{Int}, volumes::Union{Vector{Vector{Int}}, Vector{<:Tuple{Int,Int,Vararg{Int}}}}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        realization = Float64.(realization)
        dimension = size(realization)[1]
        all(t->length(t)==dimension+1, volumes) && all(facet->all(v->v in vertices, facet), volumes) || throw("The volumes don't have the correct format.")
        volumes = [Vector(facet) for facet in volumes]
        size(realization)[1]==dimension && size(realization)[2]==length(vertices) || throw("The realization does not have the correct format.")
        dimension>=1 || throw("The dimension is not an integer bigger than 0.")
        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] = realization[:,v]
        end
        variables = vcat([x[t[1],t[2]] for t in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(t[2] in pinned_vertices)]...)
        facet_equations = [det(vcat([1. for _ in 1:dimension+1]', hcat([xs[:,v] for v in facet]...))) - det(vcat([1. for _ in 1:dimension+1]', hcat([realization[:,v] for v in facet]...))) for facet in volumes]
        facet_equations = filter(eq->eq!=0, facet_equations)
        G = ConstraintSystem(vertices, variables, facet_equations, realization, xs; pinned_vertices=Vector{Int64}(pinned_vertices))
        new(G, volumes)
    end

    function VolumeHypergraph(volumes::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices::Vector{Int}=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in volumes]...))))
        VolumeHypergraph(vertices, volumes, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    Polytope([vertices,] facets, realization[; pinned_vertices])

Class for 3-dimensional polytopes with edge-length and facet planarity constraints.

For the computation of the normals, we assume that the origin lies in the barycenter. Otherwise, we translate the polytope accordingly.
"""
mutable struct Polytope
    G::ConstraintSystem
    facets::Vector{Vector{Int}}
    edges::Vector{Tuple{Int,Int}}
    x_variables::Vector{Variable}
    n_variables::Vector{Variable}

    function Polytope(vertices::Vector{Int}, facets::Union{Vector{Vector{Int}}, Vector{<:Tuple{Int, Int, Int, Vararg{Int}}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        realization = Float64.(realization)
        dimension = size(realization)[1]
        dimension==3 || throw("The dimension needs to be 3, but is $(dimension)")
        all(facet->all(v->v in vertices, facet), facets) && all(facet->length(facet)>=3, facets) || throw("The facets don't have the correct format. They need to contain at least 3 vertices each.")
        facets = [[f for f in facet] for facet in facets]
        all(v->v in vertices, pinned_vertices) || throw("`pinned_vertices` does not have the correct format.")
        centered_realization = hcat([realization[:,j] - sum(realization[:,i] for i in eachindex(vertices)) ./ length(vertices)  for j in eachindex(vertices)]...)
        if size(centered_realization)[1]==dimension && size(centered_realization)[2]==length(vertices)
            normal_realization, bars = Array{Float64,2}(undef, 3, length(facets)), []
            for j in eachindex(facets)
                normal_realization[:,j] = cross(centered_realization[:,facets[j][2]] - centered_realization[:,facets[j][1]], centered_realization[:,facets[j][3]] - centered_realization[:,facets[j][2]])
                normal_realization[:,j] = normal_realization[:,j] ./ norm(normal_realization[:,j])
                normal_realization[:,j] = normal_realization[:,j] ./ (centered_realization[:,facets[j][1]]'*normal_realization[:,j])
                for i in j+1:length(facets)
                    edge = facets[i][findall(q -> q in facets[j], facets[i])]
                    if length(edge) != 2
                        continue
                    end
                    push!(bars, Tuple(edge))
                end
            end
            _realization = hcat(centered_realization, normal_realization)
        elseif size(centered_realization)[1]==dimension && size(centered_realization)[2]==length(vertices)+length(facets)
            _realization, bars = Base.copy(centered_realization), []
            for j in eachindex(facets)
                _realization[:, length(vertices)+j] = _realization[:, length(vertices)+j] ./ (_realization[:, length(vertices)+j]'*_realization[:, facets[j][1]])
                for i in j+1:length(facets)
                    edge = facets[i][findall(q -> q in facets[j], facets[i])]
                    if length(edge) != 2
                        continue
                    end
                    push!(bars, Tuple(edge))
                end
            end
        else
            throw("The realization does not have the correct format. The size of the realization is $(size(realization)), while the vertices suggest a size of $((dimension,length(vertices)))!")
        end
        bars = collect(Set(bars))
        length(vertices)-length(bars)+length(facets)==2 || throw("The Euler characteristic of the Polytope needs to be 2, but is $(length(vertices)-length(bars)+length(facets))")
                
        @var x[1:dimension, 1:length(vertices)] n[1:dimension, 1:length(facets)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(j in pinned_vertices)]...)
        normal_variables = vcat([n[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(facets)))]...)
        xs = Expression.(hcat(x,n))
        for v in pinned_vertices
            xs[:,v] .= _realization[:,v]
        end
        facet_equations = vcat([[n[:,i]'*xs[:,facets[i][j]] - 1 for j in eachindex(facets[i])] for i in eachindex(facets)]...)
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        equations = filter(eq->eq!=0, vcat(facet_equations, bar_equations))
        all(eq->isapprox(evaluate(eq, vcat(variables, normal_variables)=>vcat([_realization[i,j] for (i,j) in collect(Iterators.product(1:size(_realization)[1], 1:size(_realization)[2])) if !(j in pinned_vertices)]...)), 0; atol=1e-4), equations) || throw(error("The given realization does not satisfy the constraints."))
        G = ConstraintSystem(vertices, vcat(variables, normal_variables), equations, _realization, xs; pinned_vertices=Vector{Int64}(pinned_vertices))
        new(G, facets, bars, variables, normal_variables)
    end

    function Polytope(facets::Union{Vector{Vector{Int}}, Vector{<:Tuple{Int, Int, Int, Vararg{Int}}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in facets]...))))
        Polytope(vertices, facets, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    FacetPolytope([vertices,] facets, realization[; pinned_vertices])

Class for 3-dimensional polytopes with facet planarity constraints.

For the computation of the normals, we assume that the origin lies in the barycenter. Otherwise, we translate the polytope accordingly.
"""
mutable struct FacetPolytope
    G::ConstraintSystem
    facets::Vector{Vector{Int}}
    x_variables::Vector{Variable}
    n_variables::Vector{Variable}

    function FacetPolytope(vertices::Vector{Int}, facets::Union{Vector{Vector{Int}}, Vector{<:Tuple{Int, Int, Int, Vararg{Int}}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        realization = Float64.(realization)
        dimension = size(realization)[1]
        all(facet->all(v->v in vertices, facet), facets) && all(facet->length(facet)>=dimension+1, facets) || throw("The facets don't have the correct format. They need to contain at least 3 vertices each.")
        facets = [[f for f in facet] for facet in facets]
        all(v->v in vertices, pinned_vertices) || throw("`pinned_vertices` does not have the correct format.")
        centered_realization = hcat([realization[:,j] - sum(realization[:,i] for i in eachindex(vertices)) ./ length(vertices)  for j in eachindex(vertices)]...)
        if size(centered_realization)[1]==dimension && size(centered_realization)[2]==length(vertices)
            normal_realization = Array{Float64,2}(undef, dimension, length(facets))
            for j in eachindex(facets)
                edge_mat = Matrix(hcat([centered_realization[:,fac] for fac in facets[j]]...)')
                normal_realization[:,j] = pinv(edge_mat) * [1. for _ in 1:length(facets[j])]
            end
            _realization = hcat(centered_realization, normal_realization)
        elseif size(centered_realization)[1]==dimension && size(centered_realization)[2]==length(vertices)+length(facets)
            _realization = Base.copy(centered_realization)
            for j in eachindex(facets)
                _realization[:, length(vertices)+j] = _realization[:, length(vertices)+j] ./ (_realization[:, length(vertices)+j]'*_realization[:, facets[j][1]])
            end
        else
            throw("The realization does not have the correct format. The size of the realization is $(size(realization)), while the vertices suggest a size of $((dimension,length(vertices)))!")
        end
                
        @var x[1:dimension, 1:length(vertices)] n[1:dimension, 1:length(facets)]
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(j in pinned_vertices)]...)
        normal_variables = vcat([n[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(facets)))]...)
        xs = Expression.(hcat(x,n))
        for v in pinned_vertices
            xs[:,v] .= _realization[:,v]
        end
        facet_equations = vcat([[n[:,i]'*xs[:,facets[i][j]] - 1 for j in eachindex(facets[i])] for i in eachindex(facets)]...)
        equations = filter(eq->eq!=0, facet_equations)
        all(eq->isapprox(evaluate(eq, vcat(variables, normal_variables)=>vcat([_realization[i,j] for (i,j) in collect(Iterators.product(1:size(_realization)[1], 1:size(_realization)[2])) if !(j in pinned_vertices)]...)), 0; atol=1e-4), equations) || throw(error("The given realization does not satisfy the constraints."))
        G = ConstraintSystem(vertices, vcat(variables, normal_variables), equations, _realization, xs; pinned_vertices=Vector{Int64}(pinned_vertices))
        new(G, facets, variables, normal_variables)
    end

    function FacetPolytope(facets::Union{Vector{Vector{Int}}, Vector{<:Tuple{Int, Int, Int, Vararg{Int}}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in facets]...))))
        FacetPolytope(vertices, facets, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    fix_antipodals!(F)

Entangles the antipodal points in a polytope `F` so that their position is constrained to antipodal points on a sphere. 
"""
function fix_antipodals!(F::Polytope)
    vertex_list = [i for i in F.G.vertices]
    while !isempty(vertex_list)
        v = pop!(vertex_list)
        for w in vertex_list
            if isapprox(norm(F.G.realization[:,v]+F.G.realization[:,w]), 0)
                index = findfirst(t->w==t, vertex_list)
                deleteat!(vertex_list, index)
                F.G.vertices = filter(t->t!=v, F.G.vertices)
                F.G.variables = filter(t->!(t in Variable.(F.G.xs[:,v])), F.G.variables)
                F.G.equations = evaluate(F.G.equations, Variable.(F.G.xs[:,v])=>-F.G.xs[:,w])
                F.G.xs[:,v] .= -F.G.xs[:,w]
                break
            end
        end
    end
    F.G.equations = filter(eq->eq!=0, F.G.equations)
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in F.G.equations]...)'
end


"""
    tetrahedral_symmetry!(F)

Add constraints to a polytope `F` representing a tetrahedral symmetry. 
"""
function tetrahedral_symmetry!(F::Polytope)
    vertex_list, antipodals = [i for i in F.G.vertices], []
    symmetry_list = []
    triang = F.facets[findfirst(facet->length(facet)==3&&(vertex_list[1] in facet), F.facets)]
    n = cross(F.G.realization[:,triang[3]]-F.G.realization[:,triang[1]], F.G.realization[:,triang[2]]-F.G.realization[:,triang[1]])
    rotation_axis = n ./ norm(n)
    angle = 2*pi/3
    R1 = [ cos(angle)+rotation_axis[1]^2*(1-cos(angle)) rotation_axis[1]*rotation_axis[2]*(1-cos(angle))-rotation_axis[3]*sin(angle) rotation_axis[1]*rotation_axis[3]*(1-cos(angle))+rotation_axis[2]*sin(angle); 
                    rotation_axis[1]*rotation_axis[2]*(1-cos(angle))+rotation_axis[3]*sin(angle) cos(angle)+rotation_axis[2]^2*(1-cos(angle)) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))-rotation_axis[1]*sin(angle); 
                    rotation_axis[1]*rotation_axis[3]*(1-cos(angle))-rotation_axis[2]*sin(angle) rotation_axis[2]*rotation_axis[3]*(1-cos(angle))+rotation_axis[1]*sin(angle) cos(angle)+rotation_axis[3]^2*(1-cos(angle));]
    
    edge = F.G.realization[:,triang[2]]-F.G.realization[:,triang[1]]
    mirror_axis = edge ./ norm(edge)
    R2 = Matrix(I, 3, 3) - 2 .* mirror_axis*mirror_axis'
    while !isempty(vertex_list)
        v = pop!(vertex_list)
        helper = [v]
        _vertex_list = Base.copy(vertex_list)
        for w in _vertex_list
            try
                if isapprox(norm(F.G.realization[:,w]-R1*F.G.realization[:,v]), 0; atol=1e-4)
                    push!(helper,w)
                    index = findfirst(t->w==t, vertex_list)
                    deleteat!(vertex_list, index)
                    F.G.vertices = filter(t->t!=w, F.G.vertices)
                    F.G.variables = filter(t->!(t in Variable.(F.G.xs[:,w])), F.G.variables)
                    F.G.equations = evaluate(F.G.equations, Variable.(F.G.xs[:,w])=>R1*F.G.xs[:,v])
                    F.G.xs[:,w] .= R1*F.G.xs[:,v]
                elseif isapprox(norm(F.G.realization[:,w]-(R1*R1)*F.G.realization[:,v]), 0; atol=1e-4)
                    push!(helper,w)
                    index = findfirst(t->w==t, vertex_list)
                    deleteat!(vertex_list, index)
                    F.G.vertices = filter(t->t!=w, F.G.vertices)
                    F.G.variables = filter(t->!(t in Variable.(F.G.xs[:,w])), F.G.variables)
                    F.G.equations = evaluate(F.G.equations, Variable.(F.G.xs[:,w])=>(R1*R1)*F.G.xs[:,v])
                    F.G.xs[:,w] .= (R1*R1)*F.G.xs[:,v]
                end
            catch e
                continue
            end
        end
        push!(symmetry_list, helper)
    end
    F.G.equations = filter(eq->eq!=0, F.G.equations)
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in F.G.equations]...)'
end

"""
    find_isolated_points(F)
"""
function find_isolated_points(F::Polytope)
    return solve(System(F.G.equations, F.G.variables))
end

"""
    triangle_shrinking(F)

Evenly shrink the triangular facets of a given polytope and compute the nontrivial infinitesimal flexes in each step.
"""
function triangle_shrinking(F::Polytope; filename="truncatedTetrahedron")
    K_n = ConstraintSystem(F.G.vertices, F.G.variables, vcat(F.G.equations, [sum( (F.G.xs[:,bar[1]]-F.G.xs[:,bar[2]]) .^2) - sum( (F.G.realization[:,bar[1]]-F.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(F.G.vertices) for j in eachindex(F.G.vertices) if i<j]]), F.G.realization, F.G.xs; pinned_vertices=F.G.pinned_vertices)
    initial_flexes = compute_nontrivial_inf_flexes(F.G, to_Array(F, F.G.realization), K_n)
    triangles = filter(facet->length(facet)==3, F.facets)
    triangle_centers = [sum(F.G.realization[:,k] for k in triang) ./ 3 for triang in triangles]
    
    for t in 0:-0.1:-1
        _realization = Base.copy(F.G.realization)
        for (i,triang) in enumerate(triangles)
            for k in triang
                _realization[:,k] .= t .* triangle_centers[i] .+ (1-t) .* F.G.realization[:,k]
            end
        end
        #println([_realization[:,k] for k in triangles[1]])
        P = Polytope(F.facets, _realization)
        K_n = ConstraintSystem(P.G.vertices, P.G.variables, vcat(P.G.equations, [sum( (P.G.xs[:,bar[1]]-P.G.xs[:,bar[2]]) .^2) - sum( (P.G.realization[:,bar[1]]-P.G.realization[:,bar[2]]) .^2) for bar in [[i,j] for i in eachindex(P.G.vertices) for j in eachindex(P.G.vertices) if i<j]]), P.G.realization, P.G.xs; pinned_vertices=P.G.pinned_vertices)
        final_flexes = compute_nontrivial_inf_flexes(P.G, to_Array(P, P.G.realization), K_n)
        display(final_flexes)
        plot(P, "$(filename)$(t)"; vertex_labels=false, vertex_size=16, vertex_color=:steelblue, padding=0.01, azimuth=0., elevation=0.035*pi, alpha=0.65)
    end
end


"""
    BodyHinge([vertices,] facets, realization[; pinned_vertices])

Class for body-hinge frameworks, which are essentially polytopes with rigid facets.
"""
mutable struct BodyHinge
    G::ConstraintSystem
    facets::Vector{Vector{Int}}
    edges::Vector{Tuple{Int,Int}}

    function BodyHinge(vertices::Vector{Int}, facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        realization = Float64.(realization)
        dimension = size(realization)[1]
        dimension==3 || throw("The dimension needs to be 3, but is $(dimension)")
        all(v->v in vertices, pinned_vertices) || throw("pinned_vertices does not have the correct format.")
        all(facet->all(v->v in vertices, facet), facets) && all(facet->length(facet)>=3, facets) || throw("The facets don't have the correct format. They need to contain at least 3 vertices each.")
        facets = [Vector(facet) for facet in facets]
        !(size(realization)[1]==dimension && size(realization)[2]==length(vertices)) && throw("The realization does not have the correct format.")
        edges = []
        for facet in facets
            foreach(i->push!(edges, (facet[i],facet[i%length(facet)+1])), 1:length(facet))
        end

        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] .= realization[:,v]
        end

        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(j in pinned_vertices)]...)
        bars = [(i,j) for facet in facets for i in facet for j in facet if i<j]
        bars = collect(Set(bars))
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in bars]
        equations = filter(eq->eq!=0, vcat(bar_equations))
        G = ConstraintSystem(vertices, variables, equations, realization, xs; pinned_vertices=Vector{Int64}(pinned_vertices))
        new(G, facets, edges)
    end

    function BodyHinge(facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in facets]...))))
        BodyHinge(vertices, facets, realization; pinned_vertices=pinned_vertices)
    end
end


"""
    BodyBar([vertices,] facets, realization[; pinned_vertices])

Class for body-bar frameworks, which are essentially polytopes with rigid facets and connecting .
"""
mutable struct BodyBar
    G::ConstraintSystem
    facets::Vector{Vector{Int}}
    edges::Vector{Tuple{Int,Int}}

    function BodyBar(vertices::Vector{Int}, edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        realization = Float64.(realization)
        dimension = size(realization)[1]
        dimension==3 || throw("The dimension needs to be 3, but is $(dimension)")
        all(v->v in vertices, pinned_vertices) || throw("pinned_vertices does not have the correct format.")
        all(facet->all(v->v in vertices, facet), facets) && all(facet->length(facet)>=3, facets) || throw("The facets don't have the correct format. They need to contain at least 3 vertices each.")
        all(edge->all(v->v in vertices, edge), edges) && all(edge->length(edge)==2, edges) || throw("The edges don't have the correct format. They need to contain 2 vertices each.")

        facets = [Vector(facet) for facet in facets]
        edges = [(edge[1],edge[2]) for edge in edges]
        !(size(realization)[1]==dimension && size(realization)[2]==length(vertices)) && throw("The realization does not have the correct format.")

        @var x[1:dimension, 1:length(vertices)]
        xs = Array{Expression,2}(undef, dimension, length(vertices))
        xs .= x
        for v in pinned_vertices
            xs[:,v] .= realization[:,v]
        end
        variables = vcat([x[i,j] for (i,j) in collect(Iterators.product(1:dimension, 1:length(vertices))) if !(j in pinned_vertices)]...)
        bars = [(i,j) for facet in facets for i in facet for j in facet if i<j]
        bars = collect(Set(bars))
        bar_equations = [sum( (xs[:,bar[1]]-xs[:,bar[2]]) .^2) - sum( (realization[:,bar[1]]-realization[:,bar[2]]) .^2) for bar in vcat(edges,bars)]
        equations = filter(eq->eq!=0, bar_equations)
        G = ConstraintSystem(vertices, variables, equations, realization, xs; pinned_vertices=Vector{Int64}(pinned_vertices))
        new(G, facets, edges)
    end

    function BodyBar(edges::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int}}}, facets::Union{Vector{Vector{Int}}, Vector{Tuple{Int,Int,Int}}}, realization::Matrix{<:Real}; pinned_vertices=Vector{Int}([]))
        vertices = sort(collect(Set(vcat([[v for v in facet] for facet in vcat(edges,facets)]...))))
        BodyBar(vertices, edges, facets, realization; pinned_vertices=pinned_vertices)
    end
end




AllTypes = Union{SpherePacking,Framework,AngularFramework,FrameworkOnSurface,SphericalDiskPacking,VolumeHypergraph,Polytope,BodyHinge,BodyBar}
AllTypesWithoutSpherePacking = Union{Framework,AngularFramework,FrameworkOnSurface,SphericalDiskPacking,VolumeHypergraph,Polytope,BodyHinge,BodyBar}
AllTypesWithoutPolytope = Union{SpherePacking,Framework,AngularFramework,FrameworkOnSurface,SphericalDiskPacking,VolumeHypergraph,BodyHinge,BodyBar}

function Base.show(io::IO, F::AllTypes)
    """
    Unifying display method for all types of geometric constraint systems.
    """
    print(io,"$(string(nameof(typeof(F)))):\n")
    print(io,"\t$(F.G)")
    if typeof(F) in [Framework, AngularFramework, FrameworkOnSurface]
        print(io,"\tBars:\t\t\t$(F.bars)")
        if F isa AngularFramework
            print(io,"\n\tAngles:\t\t$(F.angles)")
        elseif F isa FrameworkOnSurface
            @var x y z
            print(io,"\n\tSurface:\t\t$(F.surface([x,y,z])) = 0")
        end
    end
    if typeof(F) in [Polytope, BodyHinge, BodyBar]
        print(io,"\tFacets:\t\t$(F.facets[1])")
        for i in 2:min(3,length(F.facets))
            print(io,", $(F.facets[i])")
        end
        print(io,", ...")
        print(io,"\n\tEdges:\t\t$(F.edges[1])")
        for i in 2:min(5,length(F.edges))
            print(io,", $(F.edges[i])")
        end
        print(io,", ...")
    end
    if typeof(F) in [SpherePacking, SphericalDiskPacking]
        print(io,"\tContacts:\t\t$(F.contacts)\n")
        if F isa SpherePacking
            print(io,"\tRadii:\t\t$(F.radii)")
        else
            print(io,"\tInversive Distances:\t$(F.inversive_distances)")
        end
    end
    if F isa VolumeHypergraph
        print(io,"\tVolumes:\t\t$(F.volumes)")
    end
end




"""
    equations!(G, equations)

Set the equations of `G` to `equations`.
"""
function equations!(G::ConstraintSystem, equations::Vector{Expression})::Nothing
    Set(System(equations).variables)==Set(G.variables) && length(System(equations).variables)==length(G.variables) || throw("The variables in `equations` do not match the original variables.")
    G.equations = equations
    G.jacobian = hcat([differentiate(eq, G.variables) for eq in equations]...)'
    G.system = System(equations; variables=G.variables)
    return nothing
end


"""
    add_equations!(G, equations)

Add the `equations` to those of `G`.
"""
function add_equations!(G::ConstraintSystem, equations::Vector{Expression})::Nothing
    equations!(G, vcat(G.equations, equations))
    return nothing
end


"""
    equations!(F, equations)

Set the equations of the geometric constraint system `F`.
"""
function equations!(F::AllTypes, equations::Vector{Expression})
    Set(System(equations).variables)==Set(F.G.variables) && length(System(equations).variables)==length(F.G.variables) || throw("The variables in `equations` do not match the original variables.")
    F.G.equations = equations
    F.G.jacobian = hcat([differentiate(eq, F.G.variables) for eq in equations]...)'
    F.G.system = System(equations; variables=F.G.variables)
    return nothing
end


"""
    add_equations!(F, equations)

Add equations to the geometric constraint system `F`.
"""
function add_equations!(F::AllTypes, equations::Vector{Expression})
    equations!(F, vcat(F.G.equations, equations))
end

"""
    realization!(G, realization)

Set the realization of the constraint system `G` to `realization`.
"""
function realization!(G::ConstraintSystem, realization::Matrix{<:Real}; check_accuracy=false)::Nothing
    point = to_Array(G, realization)
    !check_accuracy || norm(evaluate(G.equations, G.variables=>point))>1e-8 && throw("The point does not satisfy the constraint system's equations!")
    size(realization)[1]==G.dimension && (size(realization)[2]==length(G.vertices) || size(realization)[2]==length(G.variables)//G.dimension+length(G.pinned_vertices)) || size(realization)[2]==size(G.xs)[2] || throw("The given realization does not have the correct format.")
    size(G.xs)[1]==size(realization)[1] && size(G.xs)[2]==size(realization)[2] || throw("The given realization does not have the correct format.")
    G.realization = realization
    return nothing
end

"""
    realization!(F, realization)

Set the realization of the geometric constraint system `F` to `realization`.
"""
function realization!(F::AllTypesWithoutPolytope, realization::Matrix{<:Real})::Nothing
    realization!(F.G, realization)
    return nothing
end


"""
    realization!(F, realization)

Set the realization of the geometric constraint system `F` to `realization`.
"""
function realization!(F::Polytope, realization::Matrix{<:Real})::Nothing
    if size(realization)[1]==F.G.dimension && size(realization)[2]==length(F.G.vertices)
        normal_realization = Array{Float64,2}(undef, 3, length(F.facets))
        for j in eachindex(F.facets)
            normal_realization[:,j] = cross(realization[:,F.facets[j][2]] - realization[:,F.facets[j][1]], realization[:,F.facets[j][3]] - realization[:,F.facets[j][2]])
            normal_realization[:,j] = normal_realization[:,j] ./ norm(normal_realization[:,j])
            normal_realization[:,j] = normal_realization[:,j] ./ (realization[:,F.facets[j][1]]'*normal_realization[:,j])
        end
        _realization = hcat(realization, normal_realization)
    elseif size(realization)[1]==F.G.dimension && size(realization)[2]==length(F.G.vertices)+length(F.facets)
        _realization = Base.copy(realization)
        for j in eachindex(F.facets)
            _realization[:, length(F.G.vertices)+j] = _realization[:, length(F.G.vertices)+j] ./ (_realization[:, length(F.G.vertices)+j]'*_realization[:, F.facets[j][1]])
        end
    else
        throw("The realization does not have the correct format. The size of the realization is $(size(realization)), while the vertices suggest a size of $((dimension,length(F.G.vertices)))!")
    end
    realization!(F.G, _realization)
    return nothing
end


"""
    to_Array(G, p)

Transform a realization `p` to a vector of coordinates.
"""
function to_Array(G::ConstraintSystem, p::Matrix{<:Real})::Vector{<:Real}
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(G.realization)[1], 1:size(G.realization)[2])) if !(j in G.pinned_vertices)]...)
end


"""
    to_Array(F, p)

Transform a realization `p` to a vector of coordinates.
"""
function to_Array(F::Union{SpherePacking,Framework,AngularFramework,FrameworkOnSurface,SphericalDiskPacking,VolumeHypergraph,BodyHinge,BodyBar}, p::Matrix{<:Real})::Vector{<:Real}
    return to_Array(F.G, p)
end


"""
    to_Array(F, p)

Transform a realization `p` to a vector of coordinates.
"""
function to_Array(F::Polytope, p::Matrix{<:Real})::Vector{<:Real}
    return vcat([p[i,j] for (i,j) in collect(Iterators.product(1:size(F.G.realization)[1], vcat(F.G.vertices, size(F.G.realization)[2]-length(F.facets)+1:size(F.G.realization)[2]))) if !(j in F.G.pinned_vertices)]...)
end


"""
    to_Matrix(G, q)

Transform a vector of coordinates `q` to a realization matrix.
"""
function to_Matrix(G::ConstraintSystem, q::Vector{<:Real}; flexes=false)::Matrix{<:Real}
    point = Matrix{Float64}(Base.copy(G.realization))
    point = evaluate.(G.xs, G.variables=>q)
    if flexes
        for i in eachindex(G.pinned_vertices)
            point[:,G.pinned_vertices[i]] = [0 for _ in 1:size(point)[1]]
        end
    end
    return point
end


"""
    to_Matrix(F, q)

Transform a vector of coordinates `q` to a realization matrix.
"""
function to_Matrix(F::AllTypes, q::Vector{<:Real}; flexes=false)::Matrix{<:Real}
    return to_Matrix(F.G, q; flexes=flexes)
end