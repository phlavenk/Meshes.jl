
function point_on_plane(point::Point, plane::Plane; atol = atol(coordtype(point)))
    n = normal(plane)
    p = origin(plane)
    A = point
    As = As =  (A - p) ⋅ n 
    return abs(As) < atol
end

function precise_partition(vert, method::BisectPointPartition)
    _eps = 1e-15
    n = method.normal
    p = method.point       
    left = Int[]
    right = Int[]
    intersect = Int[]        
    for (i, A) in enumerate(vert)
        As =  (A - p) ⋅ n 
        if As < (Z + _eps)
            push!(left, i)
        end
        if As > (Z - _eps)
            push!(left, i)
        end
    end

end

"""
    precise_partition(mesh::Meshes.SimpleMesh, method::BisectPointPartition) -> left[]], [right[]], intersect[]

Divides indices of triangles of the mesh in similar fashion the partition function does, but not based on centroids but looks on each point in the triangle. Returns also explicitly triangles that are crossing the plane.
"""
function precise_partition(mesh::Meshes.SimpleMesh, method::BisectPointPartition)
    _eps = 1e-15
    n = method.normal
    p = method.point       
    left = Int[]
    right = Int[]
    intersect = Int[]        
    for i in 1:nelements(mesh)        
        A, B, C = Meshes.vertices(mesh[i])  
        As =  (A - p) ⋅ n 
        Bs =  (B - p) ⋅ n 
        Cs =  (C - p) ⋅ n 
        Z = zero(coordtype(mesh))
        @show "all", As, Bs, Cs    
        isright = isleft = false    
        if all((As, Bs, Cs) .< (Z + _eps))
            push!(left, i)
            isleft = true
    #        @show "left", As, Bs, Cs
        end
        if all((As, Bs, Cs) .> (Z - _eps))
            push!(right, i)
            isright = true 
     #       @show "right", As, Bs, Cs        
        end
        if !(isright | isleft)
      #      @show "inter", As, Bs, Cs
            push!(intersect, i)
        end
    end
    return left, right, intersect
end

function remove_unused_vertices(points, connects)
    b = fill(false, size(points))
    for c in connects
        for i in indices(c)
            b[i] = true
        end
    end    
    _points = points[b]
    _ind = fill(0, size(points))
    c = 0
    for i in 1:size(points,1)
        if b[i]
            c += 1
            _ind[i] = c 
        end
    end
    @show [getindex.(Ref(_ind), indices(c)) for c in connects]
    _connects = connect.([getindex.(Ref(_ind), indices(c)) for c in connects])
    return _points, _connects    
end

function generate_elementlist(vertices, elements)
    element_list = [ Array{Int32}(undef, 0) for i = 1:length(vertices)]
    for i = 1:length(elements)
        for j = 1:3
            id = indices(elements[i])[j]            
            push!(element_list[id], i)
        end
    end
    element_list
end

# generate nearest neighbor list
# for each node it records ids of nodes that are connected to the current node
function generate_nnlist(nodes, elements, element_list)
    nnlist = [Array{Int32}(undef, 0) for i = 1:length(nodes)]
    for i = 1:length(nodes)
        for element in element_list[i]
            for j = 1:3
                node_id = indices(elements[element])[j]
                if !in(node_id, nnlist[i])
                    push!(nnlist[i], node_id)
                end
            end
        end
    end
    nnlist
end


function split_healing(mesh::Meshes.SimpleMesh, method::BisectPointPartition)    
    plane = Plane(method.point, method.normal)      
    
    points = vertices(mesh)
    elems  = elements(mesh)
    topo   = topology(mesh)
    connec = elements(topo)    
    
    left, right, cross_triangle = precise_partition(mesh, method)
    @show left, right, cross_triangle
    if length(cross_triangle) == 0
        _points = points
        l_connec = collect(connec)[left]
        r_connec = collect(connec)[right]
    else
        new_points = Point3[]    
        new_ngons = Connectivity[]
        for i in cross_triangle
            vv = vertices(mesh[i])
            segs = collect(segments(Chain(vv..., vv[begin])))        
            element_indices = first(mesh[i].vertices.indices)
            _inters = [intersect(s, plane) for s in segs]
            @show i, vv, _inters                        
            ib = _inters .!== nothing
            local tri
            n = nvertices(mesh) + length(new_points)
            Ai, Bi, Ci = element_indices
            Mi, Ni = n+1, n+2
            append!(new_points, _inters[ib])
            if ib == [true, true, false]            
                new_triangles = connect.([(Ai, Mi, Ni), (Ai, Ni, Ci), (Mi, Bi, Ni)])
            elseif ib == [true, false, true]
                new_triangles = connect.([(Ai, Mi, Ni), (Mi, Bi, Ni), (Ni, Bi, Ci)])
            elseif ib == [false, true, true]
                new_triangles = connect.([(Ai, Bi, Mi), (Ai, Mi, Ni), (Mi, Ni, Ci)])
            end
            append!(new_ngons, new_triangles)       
        end        
        _points = vcat(points, new_points)
        _mesh = SimpleMesh(_points, new_ngons)
        _connec = elements(topology(_mesh))
        L, R = partition(_mesh, bp)        
        l_connec = vcat(collect(connec)[left], collect(_connec)[L.inds])    
        r_connec = vcat(collect(connec)[right], collect(_connec)[R.inds])
    end
    l_points, l_connec =  remove_unused_vertices(_points, l_connec)
    l_mesh = SimpleMesh(l_points, l_connec)
    r_points, r_connec =  remove_unused_vertices(_points, r_connec)
    r_mesh = SimpleMesh(r_points, r_connec)
    return l_mesh, r_mesh
end

function heal_plane_cut(_mesh, pl::Plane)
    _vertices = vertices(_mesh)
    _elements = collect(elements(topology(_mesh)))
    _elist = generate_elementlist(_vertices, _elements)
    _nnlist = generate_nnlist(_vertices, _elements, _elist)
    set_i = Set{Int}()
    for e in _elements
        for i in indices(e)
            push!(set_i, i)
        end
    end
    ii = filter(i->point_on_plane(_vertices[i], pl), set_i)
    last_point = pop!(ii)
    shape = [last_point]
    i = 1
    while !isempty(ii)
        @show last_point, _nnlist[last_point]
        for p in _nnlist[last_point]
            @show i, p, last_point, ii, shape
            if p ∈ ii                
                last_point = p
                @show p, last_point
                push!(shape, p)
                pop!(ii, p)
                break
            end
        end
        yield()
        i=+1
    end        
    ngon = connect(tuple(shape...))
    return ngon
end
    

