# how ambient rays get into the throat: the Mouth interface, its tessellated
# implementation with a BVH, and Newton-refined entry against the exact surface.

"""
a way for ambient rays to enter a half-throat's d=0 surface; the contract is
enter_mouth(env, mouth, origin, dir) -> SituatedPhase or nothing on a miss,
plus half_throat(mouth). Concrete strategies (tessellation, distance fields,
...) are interchangeable behind this.
"""
abstract type Mouth end

function half_throat(mouth::Mouth)
    mouth.half_throat
end

struct MouthTriangle # one tessellation triangle of a mouth, with its chart provenance; isbits — the mesh lives on the mouth's half_throat
    corners::NTuple{3, SVector{3, Float64}} # ambient positions
    sts::NTuple{3, SVector{2, Float64}} # square coords of the corners in face_he's frame
    face_he::Int # face-inducing half-edge id
end

struct BVHNode
    lo::SVector{3, Float64}
    hi::SVector{3, Float64}
    left::Int # 0 for a leaf
    right::Int
    first::Int # leaves: range into the ordering
    count::Int
end

struct BVH
    nodes::Vector{BVHNode} # root is the last node
    order::Vector{Int}
    skip::Vector{Int} # threaded DFS: where traversal resumes after a node's subtree (0 = done)
end

function build_bvh(los, his, leaf_size=4)
    centroids = [(lo + hi) / 2 for (lo, hi) in zip(los, his)]
    nodes = BVHNode[]
    order = collect(1:length(los))
    function build(a, b)
        ids = view(order, a:b)
        lo = reduce((x, y) -> min.(x, y), los[ids])
        hi = reduce((x, y) -> max.(x, y), his[ids])
        if length(ids) <= leaf_size
            push!(nodes, BVHNode(lo, hi, 0, 0, a, length(ids)))
        else
            axis = argmax(hi - lo)
            sort!(ids, by = i -> centroids[i][axis])
            mid = (a + b) ÷ 2
            l = build(a, mid)
            r = build(mid + 1, b)
            push!(nodes, BVHNode(lo, hi, l, r, 0, 0))
        end
        length(nodes)
    end
    build(1, length(order))
    # thread the right-first DFS: children precede parents in push order, so a
    # reverse sweep sees each parent before assigning its children's skips
    skip = zeros(Int, length(nodes))
    for ix in length(nodes):-1:1
        node = nodes[ix]
        node.left == 0 && continue
        skip[node.right] = node.left
        skip[node.left] = skip[ix]
    end
    BVH(nodes, order, skip)
end

function slab_test(node::BVHNode, origin, inv_dir, t_bound)
    t0 = (node.lo - origin) .* inv_dir
    t1 = (node.hi - origin) .* inv_dir
    tnear = maximum(min.(t0, t1))
    tfar = minimum(max.(t0, t1))
    tfar >= max(tnear, 0) && tnear <= t_bound
end

struct TessellatedMouth <: Mouth # tessellated d=0 surface of a half-throat
    half_throat::HalfThroat
    triangles::Vector{MouthTriangle}
    bvh::BVH
end

function TessellatedMouth(half_throat::HalfThroat, triangles::Vector{MouthTriangle})
    los = [reduce((x, y) -> min.(x, y), tri.corners) for tri in triangles]
    his = [reduce((x, y) -> max.(x, y), tri.corners) for tri in triangles]
    TessellatedMouth(half_throat, triangles, build_bvh(los, his))
end

function TessellatedMouth(env, half_throat::HalfThroat, samples_per_edge::Int)
    pl = placement(half_throat)
    m = mesh(half_throat)
    n = samples_per_edge
    triangles = MouthTriangle[]
    for face in m.faces
        h = HalfEdgeHandle(m, (face[1], face[2]))
        chart = induced_chart(half_throat, vertex_index(h))
        wedge = half_edge_offset(h)
        grid = [pl.linear * surface(env, chart, square_coords_to_chart(valence(h), wedge, [i / n, j / n])) + pl.translation
                for i = 0:n, j = 0:n]
        at(ij) = grid[ij[1] + 1, ij[2] + 1]
        st(ij) = SVector(ij[1] / n, ij[2] / n)
        for j = 0:(n - 1), i = 0:(n - 1)
            a = (i, j); b = (i + 1, j); c = (i + 1, j + 1); d = (i, j + 1)
            push!(triangles, MouthTriangle((at(a), at(b), at(c)), (st(a), st(b), st(c)), h.id))
            push!(triangles, MouthTriangle((at(a), at(c), at(d)), (st(a), st(c), st(d)), h.id))
        end
    end
    TessellatedMouth(half_throat, triangles)
end

function ray_triangle_intersection(origin, dir, a, b, c) # Moeller-Trumbore; loose edge tolerances, Newton cleans up
    e1 = b - a
    e2 = c - a
    p = cross(dir, e2)
    det = e1' * p
    abs(det) < 1e-12 && return nothing
    s = (origin - a) / det
    u = s' * p
    (u < -1e-6 || u > 1 + 1e-6) && return nothing
    q = cross(s, e1)
    v = dir' * q
    (v < -1e-6 || u + v > 1 + 1e-6) && return nothing
    t = e2' * q
    t <= 1e-9 && return nothing
    (t, u, v)
end

# threaded (stackless) traversal, the GPU-native shape: device recursion is a
# wall, and per-thread traversal stacks are register/local-memory pressure —
# on the CPU side both also fight the allocator (MArray setindex! defeats
# escape analysis via pointer_from_objref; a wide NTuple stack exceeds vararg
# specialization and boxes per push). Descending into the right child first
# makes the visit sequence — and with it the best-hit evolution — identical
# to the former recursive version (which processed right subtrees before
# left). The running best stays the fixed-type (t, (u, v), triangle index;
# 0 = no hit yet) — a `nothing`-seeded accumulator would re-trip the
# self-call inference widening described at the dual-pass primitives (ad.jl).
function nearest_mouth_hit(mouth::TessellatedMouth, origin, dir)
    origin = SVector{3, Float64}(origin)
    dir = SVector{3, Float64}(dir)
    inv_dir = map(d -> 1 / (d == 0 ? 1e-300 : d), dir) # dodge 0*Inf NaNs in the slab test
    nodes = mouth.bvh.nodes
    best_t, best_uv, best_ix = Inf, (0.0, 0.0), 0
    ix = length(nodes) # root is the last node
    while ix != 0
        node = nodes[ix]
        if !slab_test(node, origin, inv_dir, best_t)
            ix = mouth.bvh.skip[ix]
        elseif node.left == 0
            for k in node.first:(node.first + node.count - 1)
                tri_ix = mouth.bvh.order[k]
                hit = ray_triangle_intersection(origin, dir, mouth.triangles[tri_ix].corners...)
                (hit === nothing || hit[1] >= best_t) && continue
                best_t = hit[1]
                best_uv = (hit[2], hit[3])
                best_ix = tri_ix
            end
            ix = mouth.bvh.skip[ix]
        else
            ix = node.right
        end
    end
    best_ix == 0 && return nothing
    ((best_t, best_uv[1], best_uv[2]), mouth.triangles[best_ix])
end

"""
shared entry solve: BVH hit + Newton against the exact surface; returns
(chart, d = 0 position, ambient collar differential at it, ray parameter of
the entry point) or nothing on a miss
"""
function mouth_entry(env, mouth::TessellatedMouth, origin, dir)
    origin = SVector{3, Float64}(origin)
    dir = SVector{3, Float64}(dir)
    best = nearest_mouth_hit(mouth, origin, dir)
    best === nothing && return nothing
    (τ, bu, bv), tri = best
    pl = placement(mouth.half_throat)
    h = HalfEdgeHandle(mesh(mouth.half_throat), tri.face_he)
    chart = induced_chart(mouth.half_throat, vertex_index(h))
    n = valence(h)
    wedge = half_edge_offset(h)
    amb(st) = pl.linear * surface(env, chart, square_coords_to_chart(n, wedge, st)) + pl.translation
    bary = (1 - bu - bv) * tri.sts[1] + bu * tri.sts[2] + bv * tri.sts[3]
    x = SVector(bary[1], bary[2], τ)
    for _ in 1:10 # Newton against the exact surface, in (s, t, ray parameter)
        st = SVector(x[1], x[2])
        residual = amb(st) - origin - x[3] * dir
        maximum(abs, residual) < 1e-12 && break
        jac = hcat(directional(amb, st, SVector(1.0, 0.0)), directional(amb, st, SVector(0.0, 1.0)), -dir)
        x = x - jac \ residual
    end
    w = square_coords_to_chart(n, wedge, SVector(x[1], x[2]))
    pos = SVector(w[1], w[2], 0.0)
    cl = p -> pl.linear * collar(env, chart, p) + pl.translation
    (chart, pos, jacobian_columns(cl, pos, Val(3)), x[3])
end

"""
ambient ray -> the SituatedPhase entering this mouth at d = 0, or nothing on a
miss; origin and dir live in the ambient space of the half_throat's placement
"""
function enter_mouth(env, mouth::TessellatedMouth, origin, dir)
    entry = mouth_entry(env, mouth, origin, dir)
    entry === nothing && return nothing
    chart, pos, collar_jac, _ = entry
    settle_phase(env, SituatedPhase(chart, pos, collar_jac \ SVector{3, Float64}(dir)))
end

"""
enter carrying a frame: the frame's columns pull back through the collar
differential — the exact inverse of to_ambient's export, so cross-and-return
is identity. Returns (settled phase, chart frame, entry ray parameter) or
nothing on a miss; the parameter is in units of vel, for path bookkeeping.
"""
function enter_transport(env, mouth::TessellatedMouth, origin, vel, E::SMatrix{3, N}) where {N}
    entry = mouth_entry(env, mouth, origin, vel)
    entry === nothing && return nothing
    chart, pos, collar_jac, τ = entry
    v, E2 = settle_transport(env, SituatedPhase(chart, pos, collar_jac \ SVector{3, Float64}(vel)),
                             hcat(ntuple(col -> collar_jac \ E[:, col], Val(N))...))
    (v, E2, τ)
end
