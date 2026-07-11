# convention: vertex induced neighbor is the first in its neighbor list

struct HalfEdge
    vertices::NTuple{4, Int}
    face_index::Int
end

struct Mesh # both vertices and faces consecutively 1-indexed; concrete fields keep the mesh-walking hot loop type-stable
    vertices::Vector{Vector{Float64}}
    faces::Vector{Vector{Int}} # ccw
    vertex_neighbors::Vector{Vector{Int}}
    half_edges::Dict{NTuple{2,Int}, HalfEdge}
    half_edge_offsets::Dict{NTuple{2,Int}, Int}
end

struct HalfEdgeHandle
    mesh::Mesh
    name::NTuple{2, Int}
end

function rotoindex(l,r,i)
    (i-1)%(r-l+1) + l
end

function mean(x)
    sum(x)/length(x)
end


function vertex_neighbors(faces::Vector{Vector{Int}}, vertex_total)
    neighbors::Vector{Set{Int}} = [Set{Int}() for i = 1:vertex_total]
    for f in faces
        for i = 1:length(f)
            cv = f[i]
            nv = f[rotoindex(1,length(f),i+1)]
            push!(neighbors[cv], nv)
            push!(neighbors[nv], cv)
        end
    end
    [collect(ns) for ns in neighbors]
end

function half_edges(faces::Vector{Vector{Int}})
    out = Dict{NTuple{2,Int}, HalfEdge}()
    for (faceix,f) in enumerate(faces)
        for i = 1:length(f)
            vertices = ntuple(ix -> f[rotoindex(1,length(f), ix + i - 1)], 4)
            he = HalfEdge(vertices, faceix)
            out[vertices[2:3]] = he
        end
    end
    out
end

function half_edge_offsets(vertex_neighbors::Vector{Vector{Int}}, half_edges::Dict{NTuple{2,Int}, HalfEdge})
    out = Dict{NTuple{2,Int}, Int}()
    for vix=1:length(vertex_neighbors)
        canonical_neighbor = first(vertex_neighbors[vix])
        cur_he_name = (vix, canonical_neighbor)
        cur_he = half_edges[cur_he_name]
        out[cur_he_name] = 0
        for j=2:length(vertex_neighbors[vix])
            cur_he_name = cur_he.vertices[2:-1:1]
            cur_he = half_edges[cur_he_name]
            out[cur_he_name] = j-1
        end
    end
    out
end

function handle(m::Mesh, h::HalfEdge)
    HalfEdgeHandle(m, h.vertices[2:3])
end

function next(h::HalfEdgeHandle)
    HalfEdgeHandle(h.mesh, h.mesh.half_edges[h.name].vertices[3:4])
end

function twin(h::HalfEdgeHandle)
    HalfEdgeHandle(h.mesh, h.mesh.half_edges[h.name].vertices[3:-1:2])
end

function prev(h::HalfEdgeHandle)
    HalfEdgeHandle(h.mesh, h.mesh.half_edges[h.name].vertices[1:2])
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
    length(mesh.vertex_neighbors[vertex_index])
end

function vertex_index(h::HalfEdgeHandle)
    h.name[1]
end

function valence(h::HalfEdgeHandle)
    valence(h.mesh, vertex_index(h))
end

Mesh(vertices::Vector{Vector{Float64}}, faces::Vector{Vector{Int}}) = begin
    vn = vertex_neighbors(faces, length(vertices))
    he = half_edges(faces)
    Mesh(vertices, faces, vn, he, half_edge_offsets(vn, he))
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

function edges(m::Mesh)
    out = Set()
    for face in m.faces
        for j in 1:length(face)
            half_edge = (face[j], face[(j%length(face))+1])
            push!(out, sort(half_edge))
        end
    end
    out
end

function half_edge_name_to_face_index(m::Mesh)
    Dict(name => he.face_index for (name,he) in pairs(m.half_edges))
end

function half_edge_name_to_next_name(m::Mesh)
    Dict(name => he.vertices[3:4] for (name,he) in m.half_edges)
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

function handlefan(m::Mesh, vertex_index::Int)
    fan = [vertex_induced_handle(m, vertex_index)]
    for i = 2:length(m.vertex_neighbors[vertex_index])
        push!(fan, ccw(last(fan)))
    end
    fan
end

function vertex_induced_handle(m::Mesh, vertex_index::Int)
    neighbor = first(m.vertex_neighbors[vertex_index])
    handle = HalfEdgeHandle(m, (vertex_index, neighbor))
    handle
end

function half_edge_offset(h::HalfEdgeHandle)
    h.mesh.half_edge_offsets[h.name]
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