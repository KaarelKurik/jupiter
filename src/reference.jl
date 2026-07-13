# The maintained ground truth: the semantic core in its simplest correct form,
# on plain arrays, with no performance concessions. The optimized production
# code (chart.jl/geodesic.jl) is an elaboration of this file; the "reference
# equivalence" testset holds the two together. When semantics change
# deliberately, change this file first, then make the optimized code agree.
#
# Deliberately shared with the parent (not part of the elaboration): types,
# mesh/handle machinery and accessors, and the scalar chart primitives
# (wedge coordinates, blend_scalar, wedge_map and its AD-friendly complex
# arithmetic). Those are already in their simplest form; if production ever
# elaborates one, move its simple body here first.

module Reference

using LinearAlgebra
using ForwardDiff
using TypedPolynomials

using ..jupiter: Chart, SituatedPhase, AmbientRay,
    valence, vertex_index, half_edge_handle, half_throat, induced_chart,
    other_half, params, placement, surface_polynomials, ccw, twin, next,
    half_edge_offset, wedge_index, wedge_square_coords,
    square_coords_to_chart, next_corner_coords, blend_scalar, wedge_map,
    primal, generic_normalize

directional(f, x, v) = ForwardDiff.derivative(t -> f(x .+ t .* v), 0.0)

jacobian(f, x) = reduce(hcat,
    [directional(f, x, [Float64(i == j) for i in eachindex(x)]) for j in eachindex(x)])

"""
a chart polynomial is nothing but its sum of monomials
"""
function polynomial_surface(chart, uv)
    [sum(coefficient(t) * prod(uv .^ exponents(monomial(t))) for t in terms(p))
     for p in surface_polynomials(chart)]
end

"""
blend of the containing face's 4 corner charts, weighted per corner by
blend_scalar in that corner's unit-square coordinates
"""
function surface(env, chart, uv)
    n = valence(chart)
    k = wedge_index(n, primal.(uv))
    corners = [ccw(half_edge_handle(chart), k)]
    sts = [wedge_square_coords(n, k, uv)]
    for _ in 1:3
        push!(corners, next(corners[end]))
        push!(sts, next_corner_coords(sts[end]))
    end
    weights = [blend_scalar(st[1]) * blend_scalar(st[2]) for st in sts]
    values = map(corners, sts) do corner, st
        corner_chart = induced_chart(half_throat(chart), vertex_index(corner))
        polynomial_surface(corner_chart,
                           square_coords_to_chart(valence(corner), half_edge_offset(corner), st))
    end
    sum(weights .* values) / sum(weights)
end

function surface_normal_out(env, chart, uv)
    jac = jacobian(w -> surface(env, chart, w), uv)
    generic_normalize(cross(jac[:, 1], jac[:, 2]))
end

function collar(env, chart, pos)
    uv = pos[1:2]
    surface(env, chart, uv) - pos[3] * surface_normal_out(env, chart, uv)
end

function outer_metric(env, chart, pos)
    jac = jacobian(p -> collar(env, chart, p), pos)
    jac' * jac
end

function inner_metric(env, chart, pos)
    p = params(half_throat(chart))
    jac = jacobian(uv -> surface(env, chart, uv), pos[1:2])
    g = p.cross_scale * jac' * jac
    out = zeros(eltype(g), 3, 3)
    out[1:2, 1:2] = g
    out[3, 3] = p.depth_scale
    out
end

function depth_interpolate(t, om, im) # t = depth / cylinder_depth
    w = blend_scalar(t)
    w * om + (1 - w) * im
end

function metric(env, chart, pos)
    p = params(half_throat(chart))
    depth_interpolate(pos[3] / p.cylinder_depth,
                      outer_metric(env, chart, pos), inner_metric(env, chart, pos))
end

"""
Γ^k_ij = ½ g^{ku} (∂_j g_ui + ∂_i g_ju − ∂_u g_ij), each ∂ its own pass
"""
function christoffel(env, v::SituatedPhase)
    mf = p -> metric(env, v.chart, p)
    dg = [directional(mf, v.pos, [Float64(i == j) for i in 1:3]) for j in 1:3]
    ginv = inv(mf(v.pos))
    [0.5 * sum(ginv[k, u] * (dg[j][u, i] + dg[i][j, u] - dg[u][i, j]) for u in 1:3)
     for k in 1:3, i in 1:3, j in 1:3]
end

function wvel_along_v(env, v::SituatedPhase, w::SituatedPhase)
    c = christoffel(env, v)
    [-sum(v.vel[i] * w.vel[j] * c[k, i, j] for i in 1:3, j in 1:3) for k in 1:3]
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

function view_phase_at_target(v, chart_valence, target_offset)
    central_angle = 2pi / chart_valence
    s, c = sincos(-central_angle * target_offset)
    rotation = [
        c -s 0
        s c 0
        0 0 1
    ]
    (pos=rotation * v.pos, vel=rotation * v.vel)
end

function chart_transition(v::SituatedPhase, target_offset)
    shared_edge_near = ccw(half_edge_handle(v.chart), target_offset)
    shared_edge_far = twin(shared_edge_near)
    neighbor = induced_chart(half_throat(v.chart), vertex_index(shared_edge_far))
    source_n = valence(v.chart)
    target_n = valence(neighbor)
    v_near_adj = view_phase_at_target(v, source_n, target_offset)
    wm = pos -> wedge_map(source_n, target_n, pos)
    v_far_adj_pos = [wm(v_near_adj.pos) ; v_near_adj.pos[3]]
    v_far_adj_vel = [directional(wm, v_near_adj.pos, v_near_adj.vel) ; v_near_adj.vel[3]]
    far_offset = half_edge_offset(shared_edge_far)
    v_far = view_phase_at_target((pos=v_far_adj_pos, vel=v_far_adj_vel), target_n, -far_offset)
    SituatedPhase(neighbor, v_far.pos, v_far.vel)
end

function half_transition(v::SituatedPhase) # identity in (u,v), reversal in d
    td = params(half_throat(v.chart)).transition_depth
    c2 = Chart(other_half(half_throat(v.chart)), half_edge_handle(v.chart))
    SituatedPhase(c2, [v.pos[1], v.pos[2], 2 * td - v.pos[3]], [v.vel[1], v.vel[2], -v.vel[3]])
end

exits_mouth(v::SituatedPhase) = v.pos[3] < 0 && v.vel[3] <= 0

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
    v
end

function to_ambient(env, v::SituatedPhase)
    pl = placement(half_throat(v.chart))
    cl = p -> collar(env, v.chart, p)
    AmbientRay(half_throat(v.chart), pl.linear * cl(v.pos) + pl.translation,
               pl.linear * directional(cl, v.pos, v.vel))
end

function trace_geodesic(env, v::SituatedPhase, h, max_steps)
    for _ in 1:max_steps
        v = settle_phase(env, geodesic_step(env, v, h))
        exits_mouth(v) && return to_ambient(env, v)
    end
    v
end

# ---- adaptive tracing: Dormand–Prince 5(4), error-controlled step size ----

function flow6(env, chart, u) # phase packed as [pos; vel]
    ph = SituatedPhase(chart, u[1:3], u[4:6])
    [u[4:6]; wvel_along_v(env, ph, ph)]
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
        u = [v.pos; v.vel]
        u5, err = dopri_step(env, v.chart, u, h)
        scale = tol * (1 + maximum(abs, u))
        if err <= scale
            v = settle_phase(env, SituatedPhase(v.chart, u5[1:3], u5[4:6]))
            exits_mouth(v) && return to_ambient(env, v)
        end
        h = h * clamp(0.9 * (scale / (err + floatmin()))^(1/5), 0.2, 5.0)
    end
    v
end

# ---- camera transport: parallel transport of a frame along a geodesic ----
#
# A camera is a point plus an arbitrary frame (columns of E). Transport is
# connection transport, not metric transport: with non-isometric placements
# there is no global metric, so a frame returning from a loop may come back
# rescaled or sheared — and should. Locally, though, the connection is
# metric-compatible: the Gram matrix E' g E is invariant along the transport,
# which is what the tests pin down.

"""
right-hand sides for jointly integrating a phase and a frame parallel along
it; E is a matrix of tangent-vector columns, dE^k_c = -Γ^k_ij v^i E^j_c
(the velocity obeys the same law with E = v — that is what geodesic means)
"""
function transport_flow(env, v::SituatedPhase, E)
    Γ = christoffel(env, v)
    dE = [-sum(v.vel[i] * E[j, col] * Γ[k, i, j] for i in 1:3, j in 1:3)
          for k in 1:3, col in 1:size(E, 2)]
    (v.vel, wvel_along_v(env, v, v), dE)
end

function transport_step(env, v::SituatedPhase, E, h) # one RK4 step of (pos, vel, frame), staying in v's chart
    stage(k, s) = (SituatedPhase(v.chart, v.pos + s * k[1], v.vel + s * k[2]), E + s * k[3])
    k1 = transport_flow(env, v, E)
    k2 = transport_flow(env, stage(k1, h / 2)...)
    k3 = transport_flow(env, stage(k2, h / 2)...)
    k4 = transport_flow(env, stage(k3, h)...)
    d = map((a, b, c, e) -> (a + 2b + 2c + e) / 6, k1, k2, k3, k4)
    (SituatedPhase(v.chart, v.pos + h * d[1], v.vel + h * d[2]), E + h * d[3])
end

"""
apply a phase transition f to every column of a frame at v: transitions act on
tangent vectors by their differential, i.e. exactly as they act on velocities
"""
function map_frame(v::SituatedPhase, E, f)
    reduce(hcat, [f(SituatedPhase(v.chart, v.pos, E[:, col])).vel for col in 1:size(E, 2)])
end

function chart_transition(v::SituatedPhase, E, target_offset)
    (chart_transition(v, target_offset), map_frame(v, E, w -> chart_transition(w, target_offset)))
end

function half_transition(v::SituatedPhase, E)
    (half_transition(v), map_frame(v, E, half_transition))
end

function settle_transport(env, v::SituatedPhase, E, max_hops=8) # settle_phase with the frame riding along
    for _ in 1:max_hops
        if v.pos[3] > params(half_throat(v.chart)).transition_depth
            v, E = half_transition(v, E)
            continue
        end
        n = valence(v.chart)
        k = wedge_index(n, v.pos)
        st = wedge_square_coords(n, k, v.pos)
        maximum(st) <= 0.5 && return (v, E)
        v, E = chart_transition(v, E, st[1] >= st[2] ? k : Int(mod(k + 1, n)))
    end
    (v, E)
end

function to_ambient(env, v::SituatedPhase, E) # frame columns out via the collar differential, like the velocity
    cl = p -> collar(env, v.chart, p)
    pl = placement(half_throat(v.chart))
    (to_ambient(env, v),
     reduce(hcat, [pl.linear * directional(cl, v.pos, E[:, col]) for col in 1:size(E, 2)]))
end

"""
trace a geodesic while parallel-transporting the frame E along it; returns
(AmbientRay, ambient frame) on mouth exit, (SituatedPhase, E) when out of steps
"""
function trace_transport(env, v::SituatedPhase, E, h, max_steps)
    for _ in 1:max_steps
        v, E = settle_transport(env, transport_step(env, v, E, h)...)
        exits_mouth(v) && return to_ambient(env, v, E)
    end
    (v, E)
end

"""
ray emission from a camera inside the throat: camera = chart point + frame
(columns right, up, forward — any basis, like render.jl's Camera); x, y in
[-1, 1] pick the pixel direction, and the result is ready for trace_geodesic
"""
function emit_ray(chart, pos, frame, tan_half_fov, x, y)
    SituatedPhase(chart, pos, frame * [x * tan_half_fov, y * tan_half_fov, 1.0])
end

end
