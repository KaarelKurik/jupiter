struct CatmullClarkRes
    refined_mesh::Mesh
    quarter_edge_map::Dict{NTuple{2, Int}, NTuple{2, Int}}
end

function catmullclark(m::Mesh) # returns a new mesh, assumes closed input
    face_vertices = map(m.faces) do face
        verts = m.vertices[face]
        mean(verts)
    end
    mesh_edges = edges(m)
    he_to_face = half_edge_name_to_face_index(m)
    edge_vertices = Dict(
        begin
            a,b = edge
            fa = he_to_face[(a,b)]
            fb = he_to_face[(b,a)]
            points = [face_vertices[fa], face_vertices[fb], m.vertices[a], m.vertices[b]]
            edge => mean(points)
        end
        for edge in mesh_edges
    )
    v_to_faces = vertex_to_faces(m)
    v_to_vertices = vertex_to_vertices(m)
    v_vertices = map(enumerate(m.vertices)) do (ix,v)
        face_neighbors = v_to_faces[ix]
        vertex_neighbors = v_to_vertices[ix]
        face_average = mean(map(nix -> face_vertices[nix], face_neighbors))
        edge_average = mean(map(nix -> (v + m.vertices[nix])/2, vertex_neighbors))
        n = length(face_neighbors)
        ((n-3)*v + 2*edge_average + face_average)/n
    end

    fixed_ev_pairs = collect(edge_vertices)
    all_vertices = [face_vertices ; map(x -> x[2], fixed_ev_pairs) ; v_vertices]
    block_lengths = [length(face_vertices), length(fixed_ev_pairs), length(v_vertices)]
    block_ends = accumulate(+, block_lengths)
    face_vixes = [i for i in 1:block_ends[1]]
    edge_vixes = Dict(zip(map(x -> x[1], fixed_ev_pairs), (block_ends[1]+1):block_ends[2]))
    v_vixes = [i for i in (block_ends[2]+1):block_ends[3]]

    old_v_to_new_v = v_vixes
    he_to_next = half_edge_name_to_next_name(m)

    new_faces = map(collect(he_to_face)) do (he,f)
        nhe = he_to_next[he]
        edge = sort(he)
        nedge = sort(nhe)
        point_indices = [face_vixes[f], edge_vixes[edge], v_vixes[he[2]], edge_vixes[nedge]]
        point_indices
    end

    # for every half-edge in the original,
    # produce a half-edge in the new that is a shortening of the original,
    # hence "quarter edge"
    quarter_edges = Dict()
    for face in m.faces
        rim = facerim(face)
        for he in rim
            edge = sort(he)
            quarter_edges[he] = (old_v_to_new_v[he[1]], edge_vixes[edge])
        end
    end

    refined_mesh = Mesh(all_vertices, new_faces)
    CatmullClarkRes(refined_mesh, quarter_edges)
end

"""
requires vertex to sit at center of quad star
"""
function limit_position(m::Mesh, vertex_index::Int)

    fan = handlefan(m, vertex_index)
    near_vert_indices = [next(handle).name[1] for handle in fan] 
    diag_vert_indices = [next(handle).name[2] for handle in fan]
    n = length(fan)
    central_vert = m.vertices[vertex_index]
    near_vert_avg = sum(m.vertices[ix] for ix in near_vert_indices) / n
    diag_vert_avg = sum(m.vertices[ix] for ix in diag_vert_indices) / n
    (n / (n+5)) * central_vert + (4/(n+5)) * near_vert_avg + (1/(n+5)) * diag_vert_avg
end

function limit_positions(m::Mesh)
    [limit_position(m,i) for i in 1:length(m.vertices)]
end