# convention: vertex induced neighbor is the first in its neighbor list

# both vertices and faces consecutively 1-indexed. Connectivity the hot loop
# touches is flat, id-indexed arrays (GPU-shaped); the one construction-time
# Dict is lookup-only, so no hash iteration order is ever observed — all
# numbering and canonical-frame choices are deterministic by construction.
# Half-edge id order: faces in order, slots ccw within each face.
struct Mesh
    vertices::Vector{Vector{Float64}}
    faces::Vector{Vector{Int}} # ccw
    vertex_neighbors::Vector{Vector{Int}}
    he_tail::Vector{Int}
    he_head::Vector{Int}
    he_next::Vector{Int}
    he_prev::Vector{Int}
    he_twin::Vector{Int}
    he_offset::Vector{Int} # ccw count from the tail vertex's canonical half-edge
    vertex_he::Vector{Int} # canonical outgoing half-edge: (vertex, first neighbor)
    vertex_valence::Vector{Int}
end

struct HalfEdgeHandle
    mesh::Mesh
    id::Int
end

function rotoindex(l,r,i)
    (i-1)%(r-l+1) + l
end

function mean(x)
    sum(x)/length(x)
end


# encounter order, not a Set collect: the first entry becomes the vertex's
# canonical frame, and that gauge choice must be deterministic by
# construction, not an artifact of hash iteration
function vertex_neighbors(faces::Vector{Vector{Int}}, vertex_total)
    neighbors::Vector{Vector{Int}} = [Int[] for i = 1:vertex_total]
    for f in faces
        for i = 1:length(f)
            cv = f[i]
            nv = f[rotoindex(1,length(f),i+1)]
            nv in neighbors[cv] || push!(neighbors[cv], nv)
            cv in neighbors[nv] || push!(neighbors[nv], cv)
        end
    end
    neighbors
end

function next(h::HalfEdgeHandle)
    HalfEdgeHandle(h.mesh, h.mesh.he_next[h.id])
end

function twin(h::HalfEdgeHandle)
    HalfEdgeHandle(h.mesh, h.mesh.he_twin[h.id])
end

function prev(h::HalfEdgeHandle)
    HalfEdgeHandle(h.mesh, h.mesh.he_prev[h.id])
end

function ccw(h::HalfEdgeHandle)
    twin(prev(h))
end

function ccw(h::HalfEdgeHandle, n::Int)
    cur = h
    for i=1:n
        cur = ccw(cur)
    end
    cur
end

function cw(h::HalfEdgeHandle)
    next(twin(h))
end

function valence(mesh::Mesh, vertex_index)
    mesh.vertex_valence[vertex_index]
end

function vertex_index(h::HalfEdgeHandle)
    h.mesh.he_tail[h.id]
end

function head_index(h::HalfEdgeHandle)
    h.mesh.he_head[h.id]
end

function he_name(h::HalfEdgeHandle)
    (vertex_index(h), head_index(h))
end

function valence(h::HalfEdgeHandle)
    valence(h.mesh, vertex_index(h))
end

"""
name -> id by walking the tail vertex's fan; Dict-free so the same lookup
works against device-side copies of the flat tables
"""
function half_edge_id(m::Mesh, name::NTuple{2, Int})
    e = m.vertex_he[name[1]]
    for _ in 1:m.vertex_valence[name[1]]
        m.he_head[e] == name[2] && return e
        e = m.he_twin[m.he_prev[e]] # ccw
    end
    error("no half-edge ", name)
end

HalfEdgeHandle(m::Mesh, name::NTuple{2, Int}) = HalfEdgeHandle(m, half_edge_id(m, name))

Mesh(vertices::Vector{Vector{Float64}}, faces::Vector{Vector{Int}}) = begin
    vn = vertex_neighbors(faces, length(vertices))
    total = sum(length, faces)
    ids = Dict{NTuple{2, Int}, Int}() # transient; ids themselves are face-slot order
    id = 0
    for f in faces, i in 1:length(f)
        ids[(f[i], f[mod1(i + 1, length(f))])] = (id += 1)
    end
    he_tail = zeros(Int, total); he_head = zeros(Int, total)
    he_next = zeros(Int, total); he_prev = zeros(Int, total)
    he_twin = zeros(Int, total); he_offset = zeros(Int, total)
    id = 0
    for f in faces
        L = length(f)
        for i in 1:L
            id += 1
            he_tail[id] = f[i]
            he_head[id] = f[mod1(i + 1, L)]
            he_next[id] = ids[(f[mod1(i + 1, L)], f[mod1(i + 2, L)])]
            he_prev[id] = ids[(f[mod1(i - 1, L)], f[i])]
            he_twin[id] = ids[(f[mod1(i + 1, L)], f[i])]
        end
    end
    vertex_he = [ids[(v, first(vn[v]))] for v in 1:length(vertices)]
    vertex_valence = length.(vn)
    for v in 1:length(vertices)
        e = vertex_he[v]
        for j in 0:(vertex_valence[v] - 1)
            he_offset[e] = j
            e = he_twin[he_prev[e]] # ccw, matching the fan walk above
        end
    end
    Mesh(vertices, faces, vn, he_tail, he_head, he_next, he_prev, he_twin,
         he_offset, vertex_he, vertex_valence)
end

Mesh(metamesh::MetaMesh) = begin
    vertices = [convert(Vector{Float64}, c) for c in coordinates(metamesh)]
    faces = [convert(Vector{Int}, f) for f in GeometryBasics.faces(metamesh)]
    Mesh(vertices, faces)
end

function save_obj(path, vertices, faces)
    open(path, "w") do io
        for v in vertices
            println(io, "v ", v[1], " ", v[2], " ", v[3])
        end
        for f in faces
            println(io, "f ", join(f, " "))
        end
    end
end

save_obj(path, m::Mesh) = save_obj(path, m.vertices, m.faces)

function load_obj(path) # MeshIO triangulates on load; YZ fitting needs the quads kept intact
    vertices = Vector{Float64}[]
    faces = Vector{Int}[]
    for line in eachline(path)
        parts = split(line)
        isempty(parts) && continue
        if parts[1] == "v"
            push!(vertices, parse.(Float64, parts[2:4]))
        elseif parts[1] == "f"
            push!(faces, [parse(Int, first(split(p, '/'))) for p in parts[2:end]])
        end
    end
    Mesh(vertices, faces)
end

function edges(m::Mesh) # sorted vertex pairs, in first-encounter order over faces
    seen = Set{NTuple{2, Int}}() # membership only, order never observed
    out = NTuple{2, Int}[]
    for face in m.faces
        for j in 1:length(face)
            e = sort((face[j], face[(j%length(face))+1]))
            e in seen || (push!(seen, e); push!(out, e))
        end
    end
    out
end

function vertex_to_vertices(m::Mesh)
    m.vertex_neighbors
end

function vertex_to_faces(m::Mesh)
    out = Dict()
    for (ix, face) in pairs(m.faces)
        for v in face
            neighbors = get!(out, v) do
                []
            end
            push!(neighbors, ix)
        end
    end
    out
end

function facerim(f)
    out = []
    for j in 1:length(f)
        push!(out, (f[j], f[(j%length(f))+1]))
    end
    out
end

handlefan(m::Mesh, vertex_index::Int) = handlefan(vertex_induced_handle(m, vertex_index))

function vertex_induced_handle(m::Mesh, vertex_index::Int)
    HalfEdgeHandle(m, m.vertex_he[vertex_index])
end

function half_edge_offset(h::HalfEdgeHandle)
    h.mesh.he_offset[h.id]
end

function handlefan(h::HalfEdgeHandle)
    [h ; accumulate((x,_)->ccw(x), 1:(valence(h)-1), init=h)]
end

function cubemesh()
    vertices::Vector{Vector{Float64}} = vec([Float64[i,j,k] for i in (-1,1), j in (-1,1), k in (-1,1)])
    faces::Vector{Vector{Int}} = [ # no particular order
        [1,3,4,2],
        [5,1,2,6],
        [7,3,1,5],
        [7,5,6,8],
        [4,8,6,2],
        [3,7,8,4]
    ]
    Mesh(vertices, faces)
end
