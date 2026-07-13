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

struct MouthTriangle # one tessellation triangle of a mouth, with its chart provenance
    corners::NTuple{3, SVector{3, Float64}} # ambient positions
    sts::NTuple{3, SVector{2, Float64}} # square coords of the corners in face_handle's frame
    face_handle::HalfEdgeHandle
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
    BVH(nodes, order)
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
            push!(triangles, MouthTriangle((at(a), at(b), at(c)), (st(a), st(b), st(c)), h))
            push!(triangles, MouthTriangle((at(a), at(c), at(d)), (st(a), st(c), st(d)), h))
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

function nearest_mouth_hit(mouth::TessellatedMouth, origin, dir)
    origin = SVector{3, Float64}(origin)
    dir = SVector{3, Float64}(dir)
    inv_dir = map(d -> 1 / (d == 0 ? 1e-300 : d), dir) # dodge 0*Inf NaNs in the slab test
    best_t, best_uv, best_ix = visit_bvh_node(mouth, length(mouth.bvh.nodes), origin, dir, inv_dir, Inf, (0.0, 0.0), 0)
    best_ix == 0 && return nothing
    ((best_t, best_uv[1], best_uv[2]), mouth.triangles[best_ix])
end

# recursion instead of an explicit stack keeps traversal allocation-free. The
# running best is threaded as (t, (u, v), triangle index; 0 = no hit yet) so
# the recursive signature stays fixed — a Union-typed accumulator that starts
# as `nothing` and becomes a tuple re-trips the same self-call inference
# widening described at the dual-pass primitives (ad.jl). Right child first
# reproduces the former explicit stack's pop order exactly.
function visit_bvh_node(mouth::TessellatedMouth, ix, origin, dir, inv_dir, best_t, best_uv, best_ix)
    node = mouth.bvh.nodes[ix]
    slab_test(node, origin, inv_dir, best_t) || return (best_t, best_uv, best_ix)
    if node.left == 0
        for k in node.first:(node.first + node.count - 1)
            tri_ix = mouth.bvh.order[k]
            hit = ray_triangle_intersection(origin, dir, mouth.triangles[tri_ix].corners...)
            (hit === nothing || hit[1] >= best_t) && continue
            best_t = hit[1]
            best_uv = (hit[2], hit[3])
            best_ix = tri_ix
        end
    else
        best_t, best_uv, best_ix = visit_bvh_node(mouth, node.right, origin, dir, inv_dir, best_t, best_uv, best_ix)
        best_t, best_uv, best_ix = visit_bvh_node(mouth, node.left, origin, dir, inv_dir, best_t, best_uv, best_ix)
    end
    (best_t, best_uv, best_ix)
end

"""
ambient ray -> the SituatedPhase entering this mouth at d = 0, or nothing on a
miss; origin and dir live in the ambient space of the half_throat's placement
"""
function enter_mouth(env, mouth::TessellatedMouth, origin, dir)
    origin = SVector{3, Float64}(origin)
    dir = SVector{3, Float64}(dir)
    best = nearest_mouth_hit(mouth, origin, dir)
    best === nothing && return nothing
    (τ, bu, bv), tri = best
    pl = placement(mouth.half_throat)
    chart = induced_chart(mouth.half_throat, vertex_index(tri.face_handle))
    n = valence(tri.face_handle)
    wedge = half_edge_offset(tri.face_handle)
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
    collar_jac = jacobian_columns(cl, pos, Val(3))
    settle_phase(env, SituatedPhase(chart, pos, collar_jac \ dir))
end
