
greet() = print("Hello World!")

# Let's just make half-edges a thick concept, fuck it

struct Surface
    mesh
    chart_polys
end

struct ThroatParams
    cross_scale
    depth_scale
    cylinder_depth # depth by which the metric is fully cylindrical
    transition_depth # depth of the handover to the other half
    function ThroatParams(cross_scale, depth_scale, cylinder_depth, transition_depth)
        @assert transition_depth >= cylinder_depth
        new(cross_scale, depth_scale, cylinder_depth, transition_depth)
    end
end

struct Placement # local-to-global; orientation bookkeeping lives here, not in the half-to-half gluing
    linear
    translation
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

struct SituatedPhase
    chart::Chart
    pos
    vel
end

struct SituatedPos
    chart::Chart
    pos
end

struct AmbientRay # left through a mouth; pos/vel in the ambient flat space of half_throat's placement
    half_throat::HalfThroat
    pos
    vel
end

function wedge_index_and_angle(n_wedges, pos)
    α = atan(pos[2], pos[1])
    ix = floor((α * n_wedges)/(2pi))
    remainder = α - ix * (2pi / n_wedges)
    (ix, remainder)
end

function polynomial_surface(chart, uv)
    @polyvar u v
    [p(u => uv[1], v => uv[2]) for p in surface_polynomials(chart)]
end

primal(x) = x
primal(d::ForwardDiff.Dual) = primal(ForwardDiff.value(d)) # recurse: christoffel nests duals

basis_direction(x, j) = [Float64(i == j) for i in eachindex(x)]

# the one point of contact with the AD library; see taylordiff-bugs.md for why we left TaylorDiff
directional(f, x, v) = ForwardDiff.derivative(t -> f(x .+ t .* v), 0.0)

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
    rpos = [c * pos[1] - s * pos[2], s * pos[1] + c * pos[2]]
    sqrt(2) * fake_complex_pow(rpos, n_wedges / 4)
end

function square_coords_to_chart(n_wedges, wedge_index, st)
    z = fake_complex_pow(st / sqrt(2), 4 / n_wedges)
    s, c = sincos(2pi * wedge_index / n_wedges)
    [c * z[1] - s * z[2], s * z[1] + c * z[2]]
end

"""
same face point as seen from the next corner ccw
"""
next_corner_coords(st) = [st[2], 1 - st[1]]

function surface(env, chart, uv)
    n = valence(chart)
    wix, _ = wedge_index_and_angle(n, primal.(uv))
    k = Int(mod(wix, n))
    h = ccw(half_edge_handle(chart), k)
    corners = [h ; accumulate((x,_)->next(x), 1:3, init=h)]
    st = wedge_square_coords(n, k, uv)
    sts = [[st] ; accumulate((x,_)->next_corner_coords(x), 1:3, init=st)]
    weights = [blend_scalar(p[1]) * blend_scalar(p[2]) for p in sts]
    vals = map(corners, sts) do corner, p
        c = induced_chart(half_throat(chart), vertex_index(corner))
        w = square_coords_to_chart(valence(corner), half_edge_offset(corner), p)
        polynomial_surface(c, w)
    end
    sum(weights .* vals) / sum(weights)
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

function surface_normal_out(env, chart, uv)
    s = Base.Fix1(Base.Fix1(surface, env), chart)
    du = directional(s, uv, basis_direction(uv, 1))
    dv = directional(s, uv, basis_direction(uv, 2))
    generic_normalize(cross(du, dv))
end

generic_normalize(x) = x ./ sqrt(sum(abs2, x)) # LinearAlgebra.normalize has scaling branches unfriendly to dual numbers

function collar(env, chart, pos)
    uv = pos[1:2]
    surface(env, chart, uv) - pos[3] * surface_normal_out(env, chart, uv)
end

function outer_metric(env, chart, pos)
    cl = Base.Fix1(Base.Fix1(collar, env), chart)
    collar_jac = reduce(hcat, [directional(cl, pos, basis_direction(pos, j)) for j in 1:3])
    collar_jac' * collar_jac
end

function inner_metric_params(env, chart)
    params(half_throat(chart))
end

function inner_metric(env, chart, pos)
    sf = Base.Fix1(Base.Fix1(surface, env), chart)
    params = inner_metric_params(env, chart)
    sf_jac = reduce(hcat, [directional(sf, pos[1:2], basis_direction(pos[1:2], j)) for j in 1:2])
    g = params.cross_scale * sf_jac' * sf_jac
    out = zeros(promote_type(eltype(g), typeof(params.depth_scale)), 3, 3)
    out[1:2, 1:2] = g
    out[3,3] = params.depth_scale
    out
end

function depth_interpolate(t, om, im) # t = depth / cylinder_depth
    w = blend_scalar(t)
    w * om + (1 - w) * im
end

function metric(env, chart, pos)
    om = outer_metric(env, chart, pos)
    im = inner_metric(env, chart, pos)
    depth_interpolate(pos[3] / inner_metric_params(env, chart).cylinder_depth, om, im)
end

function christoffel(env, v::SituatedPhase)
    mf = Base.Fix1(Base.Fix1(metric, env), v.chart)
    metric_derivs = stack([directional(mf, v.pos, basis_direction(v.pos, j)) for j in 1:3])
    inv_m = inv(mf(v.pos))
    @tensor begin
        cs[k,i,j] := 0.5 * inv_m[k,u] * (metric_derivs[u,i,j] + metric_derivs[j,u,i] - metric_derivs[i,j,u])
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

function to_ambient(env, v::SituatedPhase)
    pl = placement(half_throat(v.chart))
    cl = Base.Fix1(Base.Fix1(collar, env), v.chart)
    AmbientRay(half_throat(v.chart), pl.linear * cl(v.pos) + pl.translation, pl.linear * directional(cl, v.pos, v.vel))
end

function half_transition(v::SituatedPhase) # the gluing is natural: identity in (u,v), reversal in d
    td = params(half_throat(v.chart)).transition_depth
    c2 = Chart(other_half(half_throat(v.chart)), half_edge_handle(v.chart))
    SituatedPhase(c2, [v.pos[1], v.pos[2], 2 * td - v.pos[3]], [v.vel[1], v.vel[2], -v.vel[3]])
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
        wix, _ = wedge_index_and_angle(n, v.pos)
        k = Int(mod(wix, n))
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

function half_throat(ray::AmbientRay)
    ray.half_throat
end

struct MouthTriangle # one tessellation triangle of a mouth, with its chart provenance
    corners # ambient positions
    sts # square coords of the corners in face_handle's frame
    face_handle::HalfEdgeHandle
end

struct Mouth # tessellated d=0 surface of a half-throat, for ambient rays to enter through
    half_throat::HalfThroat
    triangles::Vector{MouthTriangle}
end

function Mouth(env, half_throat::HalfThroat, samples_per_edge::Int)
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
        st(ij) = [ij[1] / n, ij[2] / n]
        for j = 0:(n - 1), i = 0:(n - 1)
            a = (i, j); b = (i + 1, j); c = (i + 1, j + 1); d = (i, j + 1)
            push!(triangles, MouthTriangle((at(a), at(b), at(c)), (st(a), st(b), st(c)), h))
            push!(triangles, MouthTriangle((at(a), at(c), at(d)), (st(a), st(c), st(d)), h))
        end
    end
    Mouth(half_throat, triangles)
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

"""
ambient ray -> the SituatedPhase entering this mouth at d = 0, or nothing on a
miss; origin and dir live in the ambient space of the half_throat's placement
"""
function enter_mouth(env, mouth::Mouth, origin, dir)
    best_t = Inf
    best = nothing
    for tri in mouth.triangles
        hit = ray_triangle_intersection(origin, dir, tri.corners...)
        (hit === nothing || hit[1] >= best_t) && continue
        best_t = hit[1]
        best = (hit, tri)
    end
    best === nothing && return nothing
    (τ, bu, bv), tri = best
    pl = placement(mouth.half_throat)
    chart = induced_chart(mouth.half_throat, vertex_index(tri.face_handle))
    n = valence(tri.face_handle)
    wedge = half_edge_offset(tri.face_handle)
    amb(st) = pl.linear * surface(env, chart, square_coords_to_chart(n, wedge, st)) + pl.translation
    x = [(1 - bu - bv) * tri.sts[1] + bu * tri.sts[2] + bv * tri.sts[3] ; τ]
    for _ in 1:10 # Newton against the exact surface, in (s, t, ray parameter)
        residual = amb(x[1:2]) - origin - x[3] * dir
        maximum(abs, residual) < 1e-12 && break
        jac = [directional(amb, x[1:2], [1.0, 0.0]) directional(amb, x[1:2], [0.0, 1.0]) -dir]
        x = x - jac \ residual
    end
    pos = [square_coords_to_chart(n, wedge, x[1:2]) ; 0.0]
    cl = p -> pl.linear * collar(env, chart, p) + pl.translation
    collar_jac = reduce(hcat, [directional(cl, pos, basis_direction(pos, j)) for j in 1:3])
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
    iszero(primal(z[1])) && iszero(primal(z[2])) && return [zero(z[1]), zero(z[2])] # safe_atan2 is 0/0 at the origin
    α = safe_atan2(z[2], z[1])
    n = sqrt(z[1]^2 + z[2]^2)
    β = power * α
    s,c = sincos(β)
    n^power * [c, s]
end

function wedge_map(source_n, target_n, pos)
    a1 = fake_complex_pow(pos, source_n/4)
    a2 = [1/sqrt(2), 0] - a1
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
    rotation = [
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
    v_far_adj_pos = [wm(v_near_adj.pos) ; v_near_adj.pos[3]]
    v_far_adj_vel = [directional(wm, v_near_adj.pos, v_near_adj.vel) ; v_near_adj.vel[3]]
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
we should precompute and cache these

shape is monos by points
"""
function yz_fitting_matrix(valence)
    td = yz_degree_bound(valence)
    pinv(monomial_value_matrix(valence, td))
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

function chart_fit_polynomials(h::HalfEdgeHandle)
    fm = yz_fitting_matrix(valence(h)) # monos by points
    indices = chart_fitting_indices(h)
    # this'll work but the limit positions ought to be precomputed over the whole mesh, I figure
    limit_positions = stack(map(ix -> limit_position(h.mesh, ix), indices))' # points by 3
    run = monomial_run(yz_degree_bound(valence(h))) # monos
    coefs = fm * limit_positions # monos by 3
    vec(sum(run .* coefs, dims=1))
end

"""
returns chart polynomials in vertex order
"""
function fit_geometry(m::Mesh) # this is classic YZ, so assumes quad mesh to start with
    r1 = catmullclark(m)
    r2 = catmullclark(r1.refined_mesh)

    base_handles = [vertex_induced_handle(m, ix) for ix in 1:length(m.vertices)]
    grid_handle_names = [r2.quarter_edge_map[r1.quarter_edge_map[bh.name]] for bh in base_handles]
    grid_handles = [HalfEdgeHandle(r2.refined_mesh, shn) for shn in grid_handle_names]

    chart_polys = [chart_fit_polynomials(h) for h in grid_handles]
    chart_polys
end