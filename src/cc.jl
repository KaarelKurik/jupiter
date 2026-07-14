struct CatmullClarkRes{M}
    refined_mesh::M
    quarter_edge_map::Dict{NTuple{2, Int}, NTuple{2, Int}}
end

function catmullclark(m::Mesh) # returns a new mesh, assumes closed input
    face_vertices = map(m.faces) do face
        verts = m.vertices[face]
        mean(verts)
    end
    mesh_edges = edges(m) # first-encounter order fixes edge-vertex numbering
    he_to_face = Dict{NTuple{2, Int}, Int}((f[i], f[rotoindex(1, length(f), i+1)]) => ix
                                           for (ix, f) in pairs(m.faces) for i in 1:length(f)) # lookup only
    edge_vertices = map(mesh_edges) do (a, b)
        fa = he_to_face[(a,b)]
        fb = he_to_face[(b,a)]
        mean([face_vertices[fa], face_vertices[fb], m.vertices[a], m.vertices[b]])
    end
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

    all_vertices = [face_vertices ; edge_vertices ; v_vertices]
    block_lengths = [length(face_vertices), length(edge_vertices), length(v_vertices)]
    block_ends = accumulate(+, block_lengths)
    face_vixes = [i for i in 1:block_ends[1]]
    edge_vixes = Dict(zip(mesh_edges, (block_ends[1]+1):block_ends[2]))
    v_vixes = [i for i in (block_ends[2]+1):block_ends[3]]

    old_v_to_new_v = v_vixes

    # one refined face per half-edge, in face-slot order — numbering is
    # deterministic by construction, not inherited from a Dict collect
    new_faces = [
        begin
            a, b, c = f[i], f[rotoindex(1, length(f), i+1)], f[rotoindex(1, length(f), i+2)]
            [face_vixes[ix], edge_vixes[sort((a, b))], v_vixes[b], edge_vixes[sort((b, c))]]
        end
        for (ix, f) in pairs(m.faces) for i in 1:length(f)
    ]

    # for every half-edge in the original,
    # produce a half-edge in the new that is a shortening of the original,
    # hence "quarter edge"
    quarter_edges = Dict{NTuple{2, Int}, NTuple{2, Int}}()
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
    near_vert_indices = [head_index(handle) for handle in fan]
    diag_vert_indices = [head_index(next(handle)) for handle in fan]
    n = length(fan)
    central_vert = m.vertices[vertex_index]
    near_vert_avg = sum(m.vertices[ix] for ix in near_vert_indices) / n
    diag_vert_avg = sum(m.vertices[ix] for ix in diag_vert_indices) / n
    (n / (n+5)) * central_vert + (4/(n+5)) * near_vert_avg + (1/(n+5)) * diag_vert_avg
end

function limit_positions(m::Mesh)
    [limit_position(m,i) for i in 1:length(m.vertices)]
end