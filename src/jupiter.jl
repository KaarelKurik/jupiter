module jupiter

using FileIO
using GeometryBasics
using LinearAlgebra
using TensorOperations
using TaylorDiff

greet() = print("Hello World!")

# Let's just make half-edges a thick concept, fuck it

struct Mesh{T} # both vertices and faces consecutively 1-indexed
    vertices::Vector{Vector{T}}
    vertex_valences::Vector{UInt32}
    faces::Vector{Vector{UInt32}} # ccw
end

struct Refinement
    quarter_edge::Dict{HalfEdge, HalfEdge}
end

struct Chart
    half_throat
    half_edge
end

struct SituatedPhase
    chart
    pos
    vel
end

struct SituatedPos
    chart
    pos
end

struct SituatedPhaseVel
    chart
    posvel
    velvel
end

struct SituatedTangent
    chart
    vel
end

struct HalfEdge
    mesh::Mesh
    next::HalfEdge
    prev::HalfEdge
    twin::HalfEdge
    vertex_index
end

struct HalfThroat
    mesh
end

function mean(x)
    sum(x)/length(x)
end

function catmullclark(m::Mesh) # returns a new mesh, assumes closed input
    face_vertices = map(m.faces) do face
        verts = m.vertices[face]
        mean(verts)
    end
    mesh_edges = edges(m)
    he_to_face = half_edge_to_face(m)
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
    v_vertices = map(pairs(m.vertices)) do (ix,v)
        face_neighbors = v_to_faces[ix]
        vertex_neighbors = v_to_vertices[ix]
        face_average = mean(map(nix -> face_vertices[nix], face_neighbors))
        edge_average = mean(map(nix -> (v + m.vertices[nix])/2, vertex_neighbors))
        n = length(face_neighbors)
        ((n-3)*v + 2*edge_average + face_average)/n
    end

    fixed_ev_pairs = pairs(edge_vertices)
    all_vertices = [face_vertices ; map(x -> x[2], fixed_ev_pairs) ; v_vertices]
    block_lengths = [length(face_vertices), length(fixed_ev_pairs), length(v_vertices)]
    block_ends = accumulate(+, block_lengths)
    face_vixes = [i for i in 1:block_ends[1]]
    edge_vixes = Dict(zip(map(x -> x[1], fixed_ev_pairs), (block_ends[1]+1):block_ends[2]))
    v_vixes = [i for i in (block_ends[2]+1):block_ends[3]]

    he_to_next = half_edge_to_next(m)

    new_faces = map(pairs(he_to_face)) do (he,f)
        nhe = he_to_next[he]
        edge = sort(he)
        nedge = sort(nhe)
        point_indices = [face_vixes[f], edge_vixes[edge], v_vixes[he[2]], edge_vixes[nedge]]
        point_indices
    end
    Mesh(all_vertices, new_faces)
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

function half_edge_to_face(m::Mesh)
    out = Dict()
    for (ix,face) in pairs(m.faces)
        for j in 1:length(face)
            half_edge = (face[j], face[(j%length(face))+1])
            out[half_edge] = ix
        end
    end
    out
end

function half_edge_to_next(m::Mesh)
    out = Dict()
    for face in m.faces
        for j in 1:length(face)
            a,b = (face[j], face[(j%length(face))+1])
            c = face[((j+1)%length(face))+1]
            out[(a,b)] = (b,c)
        end
    end
    out
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

function vertex_to_vertices(m::Mesh)
    out = Dict()
    for face in m.faces
        for j in 1:length(face)
            a, b = (face[j], face[(j%length(face))+1])
            a_neighbors = get!(out, a) do
                Set()
            end
            b_neighbors = get!(out, b) do
                Set()
            end
            push!(a_neighbors, b)
            push!(b_neighbors, a)
        end
    end
    out
end

function wedge_index_and_angle(n_wedges, pos)
    α = atan(pos[2], pos[1])
    ix = floor((α * n_wedges)/(2pi))
    remainder = α - ix * (2pi / n_wedges)
    (ix, remainder)
end

function surface(env, chart, uv)
    chart.poly(u => uv[1], v => uv[2])
end

function surface_normal_out(env, chart, uv)
    s = Base.Fix1(Base.Fix1(surface, env), chart)
    du = derivative(s, uv, [1., 0.], Val(1))
    dv = derivative(s, uv, [0., 1.], Val(1))
    cross(du, dv)
end

function collar(env, chart, pos)
    uv = pos[1:2]
    surface(env, chart, uv) - pos[3] * surface_normal_out(env, chart, uv)
end

function outer_metric(env, chart, pos)
    cl = Base.Fix1(Base.Fix1(collar, env), chart)
    collar_jac = reduce(hcat, [derivative(cl, pos, [Float64(i == j) for i in 1:3], Val(1)) for j in 1:3])
    collar_jac' * collar_jac
end

function inner_metric(env, chart, pos)
    sf = Base.Fix1(Base.Fix1(surface, env), chart)
    params = inner_metric_params(env, chart)
    out = zeros(3,3)
    sf_jac = reduce(hcat, [derivative(sf, pos[1:2], [Float64(i == j) for i in 1:2], Val(1))] for j in 1:2)
    out[1:2, 1:2] = params.cross_scale * sf_jac
    out[3,3] = params.depth_scale
    out
end

function metric(env, chart, pos)
    om = outer_metric(env, chart, pos)
    im = inner_metric(env, chart, pos)
    depth_interpolate(pos[3], om, im)
end

function christoffel(env, v::SituatedPhase)
    mf = Base.Fix1(Base.Fix1(metric, env), v.chart)
    metric_derivs = reduce(hcat, [derivative(mf, v.pos, [Float64(i == j) for i in 1:3], Val(1)) for j in 1:3])
    inv_m = inv(mf(v.pos))
    @tensor begin
        cs[k,i,j] := 0.5 * inv_m[k,u] * (metric_derivs[i,j,u] + metric_derivs[j,u,i] - metric_derivs[u,i,j])
    end
    cs
end

function wvel_along_v(env, v::SituatedPhase, w::SituatedPhase)
    c = christoffel(env, v)
    @tensor begin
        wvel[k] := -v.vel[i] * w.vel[j] * c[k,i,j]
    end
    wvel
end

function wedge_map(source_n, target_n, pos) # need to check if this works as expected with TaylorDiff
    a0 = complex(pos[1], pos[2])
    a1 = a0^(source_n / 4)
    a2 = 1/sqrt(2) - a1
    a3 = a2^(4 / target_n)
    [real(a3), imag(a3)]
end

function valence(mesh::Mesh, vertex_index)
    mesh.vertex_valences[vertex_index]
end

function vertex_index(chart::Chart)
    half_edge(chart).vertex_index
end

function valence(chart::Chart)
    valence(mesh(chart), vertex_index(chart))
end

function twin(half_edge::HalfEdge)
    half_edge.twin
end

function prev(half_edge::HalfEdge)
    half_edge.prev
end

function half_throat(chart::Chart)
    chart.half_throat
end

function mesh(half_throat::HalfThroat)
    half_throat.mesh
end

function mesh(chart::Chart)
    mesh(half_throat(chart))
end

function half_edge(chart::Chart)
    chart.half_edge
end

function half_edge_ccw(half_edge::HalfEdge, rotcount)
    cur = half_edge
    for i = 1:rotcount
        cur = twin(prev(half_edge))
    end
    cur
end

function chart_by_half_edge(half_throat::HalfThroat, half_edge::HalfEdge)
    Chart(half_throat, half_edge)
end

function neighbor_chart(chart::Chart, target_index)
    ce = half_edge_ccw(half_edge(chart), target_index)
    chart_by_half_edge(half_throat, twin(ce))
end

function view_phase_at_target(v::SituatedPhase, chart_valence, target_index)
    central_angle = 2pi / chart_valence
    s,c = sincos(-central_angle * target_index)
    rotation = [
        c -s 0
        s c 0
        0 0 1
    ]
    rpos = rotation * v.pos
    rvel = rotation * v.vel
    (pos=rpos, vel=rvel)
end

function chart_transition(v::SituatedPhase, target_index) # no check on input validity
    # index local and ccw, starting with 0
    source_n = valence(v.chart)
    neighbor = neighbor_chart(v.chart, target_index)
    target_n = valence(neighbor)
    v_adj = view_phase_at_target(v, chart_valence, target_index)
    wm = Base.Fix1(Base.Fix1(wedge_map, source_n), target_n)
    newpos = wm(v_adj.pos)
    newvel = [derivative(wm, v_adj.pos, v_adj.vel, Val(1)) ; v_adj.vel[3]]
    SituatedPhase(neighbor, newpos, newvel)
end


end # module jupiter
