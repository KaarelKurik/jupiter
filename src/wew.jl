
# Let's just make half-edges a thick concept, fuck it

struct Surface
    mesh::Mesh
    chart_polys
    packed_polys::Vector{Array{Float64, 3}} # fast evaluation form of chart_polys; [component, u-power+1, v-power+1]
end

Surface(mesh, chart_polys) = Surface(mesh, chart_polys, [pack_polys(ps) for ps in chart_polys])

"""
dense coefficient table of a chart's 3 polynomials, for allocation-free Horner
evaluation; chart_polys stays the mathematical source of truth
"""
function pack_polys(polys)
    td = maximum(maxdegree.(polys))
    packed = zeros(3, td + 1, td + 1)
    for (component, p) in enumerate(polys), t in terms(p)
        eu, ev = exponents(monomial(t))
        packed[component, eu + 1, ev + 1] = coefficient(t)
    end
    packed
end

function eval_packed(packed::Array{Float64, 3}, uv)
    u, v = uv[1], uv[2]
    SVector(ntuple(Val(3)) do component
        acc_u = zero(u) * zero(v)
        for i in size(packed, 2):-1:1
            acc_v = zero(acc_u)
            for j in size(packed, 3):-1:1
                acc_v = acc_v * v + packed[component, i, j]
            end
            acc_u = acc_u * u + acc_v
        end
        acc_u
    end)
end

struct ThroatParams
    cross_scale::Float64
    depth_scale::Float64
    cylinder_depth::Float64 # depth by which the metric is fully cylindrical
    transition_depth::Float64 # depth of the handover to the other half
    function ThroatParams(cross_scale, depth_scale, cylinder_depth, transition_depth)
        @assert transition_depth >= cylinder_depth
        new(cross_scale, depth_scale, cylinder_depth, transition_depth)
    end
end

struct Placement # local-to-global; orientation bookkeeping lives here, not in the half-to-half gluing
    linear::SMatrix{3, 3, Float64, 9}
    translation::SVector{3, Float64}
end

identity_placement() = Placement([1. 0 0; 0 1 0; 0 0 1], zeros(3))

struct Throat
    surface::Surface
    params::ThroatParams
    placements::NTuple{2, Placement}
end

Throat(surface, params) = Throat(surface, params, (identity_placement(), identity_placement()))

struct HalfThroat
    throat::Throat
    side::Int # 1 or 2
end

struct Chart
    half_throat::HalfThroat
    half_edge_handle::HalfEdgeHandle
end

struct SituatedPhase{P, V} # pos/vel parametric so the integrator's phases are fully typed
    chart::Chart
    pos::P
    vel::V
end

struct AmbientRay{P, V} # left through a mouth; pos/vel in the ambient flat space of half_throat's placement
    half_throat::HalfThroat
    pos::P
    vel::V
end

function wedge_index(n_wedges, pos)
    α = atan(pos[2], pos[1])
    Int(mod(floor((α * n_wedges)/(2pi)), n_wedges))
end

function polynomial_surface(chart, uv)
    eval_packed(packed_polys(geometry(half_throat(chart)))[vertex_index(chart)], uv)
end

primal(x) = x
primal(d::ForwardDiff.Dual) = primal(ForwardDiff.value(d)) # recurse: christoffel nests duals

# the one point of contact with the AD library; see taylordiff-bugs.md for why we left TaylorDiff
directional(f, x, v) = ForwardDiff.derivative(t -> f(x .+ t .* v), 0.0)

# defined twice under two names: when one jacobian pass runs inside another
# (jacobian_columns → collar → surface_normal_out → jacobian_columns),
# inference meets the same method on its own stack and its termination
# heuristic widens the inner call to Any, boxing the entire geodesic step
# loop. A structurally identical twin per nesting level keeps every level
# inferred; inner passes use `nested_jacobian_columns`.
for jc in (:jacobian_columns, :nested_jacobian_columns)
    @eval function $jc(f, x, ::Val{N}) where {N}
        tag = typeof(ForwardDiff.Tag(f, eltype(x)))
        seeded = SVector{N}(ntuple(i -> ForwardDiff.Dual{tag}(x[i], ntuple(j -> eltype(x)(i == j), Val(N))), Val(N)))
        fd = f(seeded)
        hcat(ntuple(j -> ForwardDiff.partials.(fd, j), Val(N))...)
    end
end

"""
all coordinate directional derivatives of f at x in one dual pass; returns an
SMatrix of columns [df/dx_1 ... df/dx_n] like the hcat-of-directionals it
replaces. N must equal length(x) (passed as Val so everything stack-allocates).
`nested_jacobian_columns` is its mandatory twin for passes running inside
another pass; see the comment on the definition loop.
"""
jacobian_columns

flat_bump(x) = primal(x) > 0 ? exp(-1 / x) : zero(x)

"""
placeholder pending a deliberate choice; the contract is
f(0)=1, f(1)=0, C^inf and flat at both ends
(flat at 0: smoothness across wedge seams and chart centers;
flat at 1: smoothness where a chart's support ends).
f(x)+f(1-x)=1 additionally makes corner weights sum to exactly 1.
"""
function blend_scalar(x)
    flat_bump(1 - x) / (flat_bump(x) + flat_bump(1 - x))
end

"""
chart coords -> [0,1]^2 coords of the face at wedge_index,
with u along the wedge's starting edge
"""
function wedge_square_coords(n_wedges, wedge_index, pos)
    s, c = sincos(-2pi * wedge_index / n_wedges)
    rpos = SVector(c * pos[1] - s * pos[2], s * pos[1] + c * pos[2])
    sqrt(2) * fake_complex_pow(rpos, n_wedges / 4)
end

function square_coords_to_chart(n_wedges, wedge_index, st)
    z = fake_complex_pow(st / sqrt(2), 4 / n_wedges)
    s, c = sincos(2pi * wedge_index / n_wedges)
    SVector(c * z[1] - s * z[2], s * z[1] + c * z[2])
end

"""
same face point as seen from the next corner ccw
"""
next_corner_coords(st) = SVector(st[2], 1 - st[1])

function corner_contribution(chart, corner, st)
    w = blend_scalar(st[1]) * blend_scalar(st[2])
    c = induced_chart(half_throat(chart), vertex_index(corner))
    (w, polynomial_surface(c, square_coords_to_chart(valence(corner), half_edge_offset(corner), st)))
end

function surface(env, chart, uv)
    n = valence(chart)
    k = wedge_index(n, primal.(uv))
    corner = ccw(half_edge_handle(chart), k)
    st = wedge_square_coords(n, k, uv)
    weight, acc = corner_contribution(chart, corner, st)
    acc = weight * acc
    total_weight = weight
    for _ in 1:3
        corner = next(corner)
        st = next_corner_coords(st)
        weight, val = corner_contribution(chart, corner, st)
        acc = acc + weight * val
        total_weight += weight
    end
    acc / total_weight
end

function surface_polynomials(chart::Chart)
    chart_polys(geometry(half_throat(chart)))[vertex_index(chart)]
end

"""
samples each face on a regular grid; boundary vertices are duplicated
between faces, so the result is raw (vertices, faces) for export, not a Mesh
"""
function sample_surface(half_throat::HalfThroat, samples_per_edge::Int)
    m = mesh(half_throat)
    n = samples_per_edge
    vertices = Vector{Float64}[]
    faces = Vector{Int}[]
    for face in m.faces
        h = HalfEdgeHandle(m, (face[1], face[2]))
        chart = induced_chart(half_throat, vertex_index(h))
        wedge = half_edge_offset(h)
        base = length(vertices)
        for j = 0:n, i = 0:n
            w = square_coords_to_chart(valence(h), wedge, [i/n, j/n])
            push!(vertices, surface(nothing, chart, w))
        end
        corner(i, j) = base + j * (n + 1) + i + 1
        for j = 0:(n-1), i = 0:(n-1)
            push!(faces, [corner(i,j), corner(i+1,j), corner(i+1,j+1), corner(i,j+1)])
        end
    end
    (vertices, faces)
end

# f(pos)-shaped view of f(env, chart, pos); Fix1 (not a closure) so the callable's
# type stays concrete for AD tags and inference in the hot paths below
situate(f, env, chart) = Base.Fix1(Base.Fix1(f, env), chart)

function surface_normal_out(env, chart, uv)
    s = situate(surface, env, chart)
    jac = nested_jacobian_columns(s, uv, Val(2)) # inner pass: runs inside collar's own jacobian, see the twin definition
    generic_normalize(cross(jac[:, 1], jac[:, 2]))
end

generic_normalize(x) = x ./ sqrt(sum(abs2, x)) # LinearAlgebra.normalize has scaling branches unfriendly to dual numbers

function collar(env, chart, pos)
    uv = SVector(pos[1], pos[2])
    surface(env, chart, uv) - pos[3] * surface_normal_out(env, chart, uv)
end

function outer_metric(env, chart, pos)
    cl = situate(collar, env, chart)
    collar_jac = jacobian_columns(cl, pos, Val(3))
    collar_jac' * collar_jac
end

function inner_metric_params(env, chart)
    params(half_throat(chart))
end

function inner_metric(env, chart, pos)
    sf = situate(surface, env, chart)
    params = inner_metric_params(env, chart)
    sf_jac = jacobian_columns(sf, SVector(pos[1], pos[2]), Val(2))
    g = params.cross_scale * sf_jac' * sf_jac
    z = zero(eltype(g))
    SMatrix{3, 3}(g[1, 1], g[2, 1], z, g[1, 2], g[2, 2], z, z, z, eltype(g)(params.depth_scale))
end

function depth_interpolate(t, om, im) # t = depth / cylinder_depth
    w = blend_scalar(t)
    w * om + (1 - w) * im
end

function metric(env, chart, pos)
    t = pos[3] / inner_metric_params(env, chart).cylinder_depth
    # the blend is exactly outer for d ≤ 0 and exactly inner past cylinder_depth
    # (C^∞-flat ends make the skipped half's weight an exact dual zero there)
    primal(t) <= 0 && return outer_metric(env, chart, pos)
    primal(t) >= 1 && return inner_metric(env, chart, pos)
    depth_interpolate(t, outer_metric(env, chart, pos), inner_metric(env, chart, pos))
end

function christoffel(env, v::SituatedPhase)
    mf = situate(metric, env, v.chart)
    # one dual pass with a 3-partial seed: metric value and all three derivatives together
    tag = typeof(ForwardDiff.Tag(mf, eltype(v.pos)))
    seeded = SVector{3}(ntuple(i -> ForwardDiff.Dual{tag}(v.pos[i], ntuple(j -> eltype(v.pos)(i == j), Val(3))), Val(3)))
    md = mf(seeded)
    dg = ntuple(j -> ForwardDiff.partials.(md, j), Val(3)) # dg[c][a,b] = ∂g_ab/∂x_c
    inv_m = inv(ForwardDiff.value.(md))
    SArray{Tuple{3, 3, 3}}(ntuple(Val(27)) do n
        k = (n - 1) % 3 + 1
        i = ((n - 1) ÷ 3) % 3 + 1
        j = (n - 1) ÷ 9 + 1
        0.5 * sum(inv_m[k, u] * (dg[j][u, i] + dg[i][j, u] - dg[u][i, j]) for u in 1:3)
    end)
end

function wvel_along_v(env, v::SituatedPhase, w::SituatedPhase)
    c = christoffel(env, v)
    SVector{3}(ntuple(k -> -sum(v.vel[i] * w.vel[j] * c[k, i, j] for i in 1:3, j in 1:3), Val(3)))
end

function to_ambient(env, v::SituatedPhase)
    pl = placement(half_throat(v.chart))
    cl = situate(collar, env, v.chart)
    AmbientRay(half_throat(v.chart), pl.linear * cl(v.pos) + pl.translation, pl.linear * directional(cl, v.pos, v.vel))
end

function half_transition(v::SituatedPhase) # the gluing is natural: identity in (u,v), reversal in d
    td = params(half_throat(v.chart)).transition_depth
    c2 = Chart(other_half(half_throat(v.chart)), half_edge_handle(v.chart))
    SituatedPhase(c2, SVector(v.pos[1], v.pos[2], 2 * td - v.pos[3]), SVector(v.vel[1], v.vel[2], -v.vel[3]))
end

exits_mouth(v::SituatedPhase) = v.pos[3] < 0 && v.vel[3] <= 0 # once here, the flat outer region owns the ray

"""
re-express v in a chart that contains it comfortably: hop laterally while the
containing face's square coords leave [0, 1/2]^2, hand over to the other half
past transition_depth; closed on phases, so mouth exit is the caller's concern
"""
function settle_phase(env, v::SituatedPhase, max_hops=8)
    for _ in 1:max_hops
        if v.pos[3] > params(half_throat(v.chart)).transition_depth
            v = half_transition(v)
            continue
        end
        n = valence(v.chart)
        k = wedge_index(n, v.pos)
        st = wedge_square_coords(n, k, v.pos)
        maximum(st) <= 0.5 && return v
        v = chart_transition(v, st[1] >= st[2] ? k : Int(mod(k + 1, n)))
    end
    v # hop budget exhausted; shouldn't happen for step sizes small next to a face
end

function geodesic_flow(env, v::SituatedPhase)
    (v.vel, wvel_along_v(env, v, v))
end

function geodesic_step(env, v::SituatedPhase, h) # one RK4 step, staying in v's chart
    stage(k, s) = SituatedPhase(v.chart, v.pos + s * k[1], v.vel + s * k[2])
    k1 = geodesic_flow(env, v)
    k2 = geodesic_flow(env, stage(k1, h / 2))
    k3 = geodesic_flow(env, stage(k2, h / 2))
    k4 = geodesic_flow(env, stage(k3, h))
    dpos = (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    dvel = (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
    SituatedPhase(v.chart, v.pos + h * dpos, v.vel + h * dvel)
end

function trace_geodesic(env, v::SituatedPhase, h, max_steps)
    for _ in 1:max_steps
        v = settle_phase(env, geodesic_step(env, v, h))
        exits_mouth(v) && return to_ambient(env, v)
    end
    v # still inside; caller decides whether that's a problem
end

function flow6(env, chart, u) # phase packed as [pos; vel]
    vel = SVector(u[4], u[5], u[6])
    ph = SituatedPhase(chart, SVector(u[1], u[2], u[3]), vel)
    vcat(vel, wvel_along_v(env, ph, ph))
end

"""
one Dormand–Prince 5(4) attempt in a fixed chart: the 5th-order state and the
embedded 4th-order error estimate (sup norm)
"""
function dopri_step(env, chart, u, h)
    f(x) = flow6(env, chart, x)
    k1 = f(u)
    k2 = f(u + h * (1/5)k1)
    k3 = f(u + h * ((3/40)k1 + (9/40)k2))
    k4 = f(u + h * ((44/45)k1 - (56/15)k2 + (32/9)k3))
    k5 = f(u + h * ((19372/6561)k1 - (25360/2187)k2 + (64448/6561)k3 - (212/729)k4))
    k6 = f(u + h * ((9017/3168)k1 - (355/33)k2 + (46732/5247)k3 + (49/176)k4 - (5103/18656)k5))
    u5 = u + h * ((35/384)k1 + (500/1113)k3 + (125/192)k4 - (2187/6784)k5 + (11/84)k6)
    k7 = f(u5)
    err = h * ((71/57600)k1 - (71/16695)k3 + (71/1920)k4 - (17253/339200)k5 + (22/525)k6 - (1/40)k7)
    (u5, maximum(abs, err))
end

"""
error-controlled trace. h0 seeds the controller; each step is additionally
capped so it cannot outrun its chart (stages evaluate in one fixed chart,
trustworthy only out to ~face scale). max_attempts counts accepted AND
rejected steps: it is the work budget, not the arc length.
"""
function trace_geodesic(env, v::SituatedPhase, h0, max_attempts, tol)
    h = float(h0)
    for _ in 1:max_attempts
        h = min(h, 0.25 / (maximum(abs, v.vel) + 1e-12))
        u = vcat(SVector{3}(v.pos), SVector{3}(v.vel))
        u5, err = dopri_step(env, v.chart, u, h)
        scale = tol * (1 + maximum(abs, u))
        if err <= scale
            v = settle_phase(env, SituatedPhase(v.chart, SVector(u5[1], u5[2], u5[3]), SVector(u5[4], u5[5], u5[6])))
            exits_mouth(v) && return to_ambient(env, v)
        end
        h = h * clamp(0.9 * (scale / (err + floatmin()))^(1/5), 0.2, 5.0)
    end
    v
end

function half_throat(ray::AmbientRay)
    ray.half_throat
end

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
# as `nothing` and becomes a tuple re-trips the same inference widening that
# jacobian_columns needed its twin for. Right child first reproduces the
# former explicit stack's pop order exactly.
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

function reference_wedge_map(source_n, target_n, pos) # complex powers don't differentiate; test oracle only
    a0 = complex(pos[1], pos[2])
    a1 = a0^(source_n / 4)
    a2 = 1/sqrt(2) - a1
    a3 = a2^(4 / target_n)
    [real(a3), imag(a3)]
end

# complementary half-angle form for x<0, where the first form is 0/0 near the negative real axis
safe_atan2(y, x) = primal(x) >= 0 ? 2 * atan(y / (sqrt(x^2 + y^2) + x)) : 2 * atan((sqrt(x^2 + y^2) - x) / y)

function fake_complex_pow(z, power)
    iszero(primal(z[1])) && iszero(primal(z[2])) && return SVector(zero(z[1]), zero(z[2])) # safe_atan2 is 0/0 at the origin
    α = safe_atan2(z[2], z[1])
    n = sqrt(z[1]^2 + z[2]^2)
    β = power * α
    s,c = sincos(β)
    n^power * SVector(c, s)
end

function wedge_map(source_n, target_n, pos)
    a1 = fake_complex_pow(pos, source_n/4)
    a2 = SVector(1/sqrt(2), 0.0) - a1
    a3 = fake_complex_pow(a2, 4/target_n)
    a3
end

function vertex_index(chart::Chart)
    vertex_index(half_edge_handle(chart))
end

function valence(chart::Chart)
    valence(mesh(chart), vertex_index(chart))
end

function half_throat(chart::Chart)
    chart.half_throat
end

function throat(half_throat::HalfThroat)
    half_throat.throat
end

function side(half_throat::HalfThroat)
    half_throat.side
end

function other_half(half_throat::HalfThroat)
    HalfThroat(throat(half_throat), 3 - side(half_throat))
end

function geometry(throat::Throat) # "surface" the accessor would collide with surface the evaluator
    throat.surface
end

function geometry(half_throat::HalfThroat)
    geometry(throat(half_throat))
end

function params(throat::Throat)
    throat.params
end

function params(half_throat::HalfThroat)
    params(throat(half_throat))
end

function placement(half_throat::HalfThroat)
    throat(half_throat).placements[side(half_throat)]
end

function chart_polys(surface::Surface)
    surface.chart_polys
end

function packed_polys(surface::Surface)
    surface.packed_polys
end

function mesh(surface::Surface)
    surface.mesh
end

function mesh(half_throat::HalfThroat)
    mesh(geometry(half_throat))
end

function mesh(chart::Chart)
    mesh(half_throat(chart))
end

function half_edge_handle(chart::Chart)
    chart.half_edge_handle
end

function induced_chart(half_throat::HalfThroat, vertex_index::Int)
    Chart(half_throat, vertex_induced_handle(mesh(half_throat), vertex_index))
end

function neighbor_chart(chart::Chart, target_offset)
    ce = ccw(half_edge_handle(chart), target_offset)
    induced_chart(half_throat(chart), vertex_index(twin(ce)))
end

function view_phase_at_target(v, chart_valence, target_offset)
    central_angle = 2pi / chart_valence
    s,c = sincos(-central_angle * target_offset)
    rotation = @SMatrix [
        c -s 0
        s c 0
        0 0 1
    ]
    rpos = rotation * v.pos
    rvel = rotation * v.vel
    (pos=rpos, vel=rvel)
end

function chart_transition(v::SituatedPhase, target_offset) # no check on input validity
    # index local and ccw, starting with 0
    shared_edge_near = ccw(half_edge_handle(v.chart), target_offset)
    shared_edge_far = twin(shared_edge_near)
    neighbor = induced_chart(half_throat(v.chart), vertex_index(shared_edge_far))
    source_n = valence(v.chart)
    target_n = valence(neighbor)
    v_near_adj = view_phase_at_target(v, source_n, target_offset)
    wm = Base.Fix1(Base.Fix1(wedge_map, source_n), target_n)
    wp = wm(v_near_adj.pos)
    wv = directional(wm, v_near_adj.pos, v_near_adj.vel)
    v_far_adj_pos = SVector(wp[1], wp[2], v_near_adj.pos[3])
    v_far_adj_vel = SVector(wv[1], wv[2], v_near_adj.vel[3])
    far_offset = half_edge_offset(shared_edge_far)
    v_far = view_phase_at_target((pos=v_far_adj_pos, vel=v_far_adj_vel), valence(neighbor), -far_offset)
    SituatedPhase(neighbor, v_far.pos, v_far.vel)
end

"""
x:3, y:4, wedge:valence
"""
function nonzero_fitting_points(valence)
    trunc_square = [i/4 + (j/4) * 1im for i=1:3, j=0:3]/sqrt(2)
    wedge = map(z -> z^(4/valence), trunc_square)
    wedges =  map(x -> [real(x), imag(x)], stack(wedge * cispi(2 * i / valence) for i = 0:(valence-1)))
    wedges
end

function fitting_points(valence)
    [[[0,0]] ; vec(nonzero_fitting_points(valence))]
end

function monomial_run(total_degree)
    @polyvar u v
    run = sort(reduce(vcat, [u^k * v^(s-k) for k = 0:s] for s=0:total_degree))
    run
end

"""
points by monos
"""
function monomial_value_matrix(valence, total_degree)
    @polyvar u v
    run = monomial_run(total_degree)
    fp = fitting_points(valence)
    mvm = [m((u,v)=>p) for p in fp, m in run]
    mvm
end

function yz_degree_bound(valence)
    min(14, valence+1)
end

"""
shape is monos by points
"""
const yz_fitting_matrix_cache = Dict{Int, Matrix{Float64}}()
function yz_fitting_matrix(valence)
    get!(yz_fitting_matrix_cache, valence) do
        pinv(monomial_value_matrix(valence, yz_degree_bound(valence)))
    end
end

"""
yields indices for nonzero part of one wedge with
x:3, y:4
"""
function grid_fitting_indices(grid_handle::HalfEdgeHandle)
    base = accumulate((x,_)->cw(next(x)), 1:3, init=grid_handle)
    grid_arrangement = map(base) do bh
        [bh ; accumulate((x,_)->twin(next(next(x))), 1:3, init=bh)]
    end
    grid = [grid_arrangement[i][j] for i=1:3, j=1:4]
    vertex_indices = [h.name[1] for h in grid]
    vertex_indices
end

"""
x:3, y:4, wedge:valence(h)
"""
function chart_nonzero_fitting_indices(h::HalfEdgeHandle)
    fan = handlefan(h)
    stack(map(grid_fitting_indices, fan))
end

"""
call on refined mesh handle
"""
function chart_fitting_indices(h::HalfEdgeHandle)
    [[vertex_index(h)] ; vec(chart_nonzero_fitting_indices(h))]
end

function chart_fit_polynomials(h::HalfEdgeHandle, mesh_limit_positions)
    fm = yz_fitting_matrix(valence(h)) # monos by points
    indices = chart_fitting_indices(h)
    limits = stack(mesh_limit_positions[indices])' # points by 3
    run = monomial_run(yz_degree_bound(valence(h))) # monos
    coefs = fm * limits # monos by 3
    vec(sum(run .* coefs, dims=1))
end

chart_fit_polynomials(h::HalfEdgeHandle) = chart_fit_polynomials(h, limit_positions(h.mesh))

"""
returns chart polynomials in vertex order
"""
function fit_geometry(m::Mesh) # this is classic YZ, so assumes quad mesh to start with
    r1 = catmullclark(m)
    r2 = catmullclark(r1.refined_mesh)

    base_handles = [vertex_induced_handle(m, ix) for ix in 1:length(m.vertices)]
    grid_handle_names = [r2.quarter_edge_map[r1.quarter_edge_map[bh.name]] for bh in base_handles]
    grid_handles = [HalfEdgeHandle(r2.refined_mesh, shn) for shn in grid_handle_names]

    limits = limit_positions(r2.refined_mesh)
    chart_polys = [chart_fit_polynomials(h, limits) for h in grid_handles]
    chart_polys
end