# The maintained ground truth: the semantic core in its simplest correct form,
# on plain arrays, with no performance concessions. The optimized code in
# wew.jl is an elaboration of this file; the "reference equivalence" testset
# holds the two together. When semantics change deliberately, change this file
# first, then make the optimized code agree.
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
    half_edge_offset, wedge_index_and_angle, wedge_square_coords,
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
    wix, _ = wedge_index_and_angle(n, primal.(uv))
    k = Int(mod(wix, n))
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
        wix, _ = wedge_index_and_angle(n, v.pos)
        k = Int(mod(wix, n))
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

end
