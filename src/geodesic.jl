# the metric and its geodesics: collar embedding, outer/inner metric blend,
# christoffels, phase settling across chart/half transitions, and the two
# tracers (fixed-step RK4, error-controlled DP5(4)).

function collar(env, chart, pos) # surface point pushed inward: s(u,v) − d·n̂(u,v)
    uv = SVector(pos[1], pos[2])
    val, jac = value_and_jacobian_columns(situate(surface, env, chart), uv, Val(2))
    val - pos[3] * normal_from_columns(jac)
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
    g = carrier(pos[1])(params.cross_scale) * sf_jac' * sf_jac
    z = zero(eltype(g))
    SMatrix{3, 3}(g[1, 1], g[2, 1], z, g[1, 2], g[2, 2], z, z, z, eltype(g)(params.depth_scale))
end

function depth_interpolate(t, om, im) # t = depth / cylinder_depth
    w = blend_scalar(t)
    w * om + (1 - w) * im
end

function metric(env, chart, pos)
    t = pos[3] / carrier(pos[3])(inner_metric_params(env, chart).cylinder_depth)
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
    half = carrier(v.pos[1])(0.5)
    SArray{Tuple{3, 3, 3}}(ntuple(Val(27)) do n
        k = (n - 1) % 3 + 1
        i = ((n - 1) ÷ 3) % 3 + 1
        j = (n - 1) ÷ 9 + 1
        half * sum(inv_m[k, u] * (dg[j][u, i] + dg[i][j, u] - dg[u][i, j]) for u in 1:3)
    end)
end

# the transport law's right-hand side: covariant rate of w carried along vel
christoffel_pull(Γ, vel, w) = SVector{3}(ntuple(k -> -sum(vel[i] * w[j] * Γ[k, i, j] for i in 1:3, j in 1:3), Val(3)))

function wvel_along_v(env, v::SituatedPhase, w::SituatedPhase)
    christoffel_pull(christoffel(env, v), v.vel, w.vel)
end

function to_ambient(env, v::SituatedPhase)
    pl = placement(half_throat(v.chart))
    C = carrier(v.pos[1])
    lin = SMatrix{3, 3, C}(pl.linear)
    cl = situate(collar, env, v.chart)
    AmbientRay(half_throat(v.chart), lin * cl(v.pos) + SVector{3, C}(pl.translation), lin * directional(cl, v.pos, v.vel))
end

function half_transition(v::SituatedPhase) # the gluing is natural: identity in (u,v), reversal in d
    td = carrier(v.pos[3])(params(half_throat(v.chart)).transition_depth)
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
    h = carrier(v.pos[1])(h)
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
    c = carrier(u[1]) # tableau ratios computed in Float64, one rounding into the carrier
    k1 = f(u)
    k2 = f(u + h * c(1/5)*k1)
    k3 = f(u + h * (c(3/40)*k1 + c(9/40)*k2))
    k4 = f(u + h * (c(44/45)*k1 - c(56/15)*k2 + c(32/9)*k3))
    k5 = f(u + h * (c(19372/6561)*k1 - c(25360/2187)*k2 + c(64448/6561)*k3 - c(212/729)*k4))
    k6 = f(u + h * (c(9017/3168)*k1 - c(355/33)*k2 + c(46732/5247)*k3 + c(49/176)*k4 - c(5103/18656)*k5))
    u5 = u + h * (c(35/384)*k1 + c(500/1113)*k3 + c(125/192)*k4 - c(2187/6784)*k5 + c(11/84)*k6)
    k7 = f(u5)
    err = h * (c(71/57600)*k1 - c(71/16695)*k3 + c(71/1920)*k4 - c(17253/339200)*k5 + c(22/525)*k6 - c(1/40)*k7)
    (u5, maximum(abs, err))
end

"""
error-controlled trace. h0 seeds the controller; each step is additionally
capped so it cannot outrun its chart (stages evaluate in one fixed chart,
trustworthy only out to ~face scale). max_attempts counts accepted AND
rejected steps: it is the work budget, not the arc length.
"""
function trace_geodesic(env, v::SituatedPhase, h0, max_attempts, tol)
    C = carrier(v.pos[1])
    h = C(h0)
    tol = C(tol)
    for _ in 1:max_attempts
        h = min(h, C(0.25) / (maximum(abs, v.vel) + C(1e-12)))
        u = vcat(SVector{3}(v.pos), SVector{3}(v.vel))
        u5, err = dopri_step(env, v.chart, u, h)
        scale = tol * (1 + maximum(abs, u))
        if err <= scale
            v = settle_phase(env, SituatedPhase(v.chart, SVector(u5[1], u5[2], u5[3]), SVector(u5[4], u5[5], u5[6])))
            exits_mouth(v) && return to_ambient(env, v)
        end
        h = h * clamp(C(0.9) * (scale / (err + floatmin(C)))^C(1/5), C(0.2), C(5.0))
    end
    v
end

# ---- camera transport: parallel transport of a frame along a geodesic ----
# elaboration of the reference (see reference.jl for the law and the comment on
# connection vs metric transport); the conserved object is the Gram matrix
# E' g E. One Γ evaluation feeds the velocity and every frame column.

function transport_flow(env, v::SituatedPhase, E::SMatrix{3, N}) where {N}
    Γ = christoffel(env, v)
    (v.vel, christoffel_pull(Γ, v.vel, v.vel),
     hcat(ntuple(col -> christoffel_pull(Γ, v.vel, E[:, col]), Val(N))...))
end

function transport_step(env, v::SituatedPhase, E::SMatrix, h) # one RK4 step of (pos, vel, frame), staying in v's chart
    stage(k, s) = (SituatedPhase(v.chart, v.pos + s * k[1], v.vel + s * k[2]), E + s * k[3])
    k1 = transport_flow(env, v, E)
    k2 = transport_flow(env, stage(k1, h / 2)...)
    k3 = transport_flow(env, stage(k2, h / 2)...)
    k4 = transport_flow(env, stage(k3, h)...)
    dpos = (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]) / 6
    dvel = (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]) / 6
    dE = (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]) / 6
    (SituatedPhase(v.chart, v.pos + h * dpos, v.vel + h * dvel), E + h * dE)
end

"""
apply a phase transition f to every column of a frame at v: transitions act on
tangent vectors by their differential, i.e. exactly as they act on velocities
"""
function map_frame(v::SituatedPhase, E::SMatrix{3, N}, f) where {N}
    hcat(ntuple(col -> f(SituatedPhase(v.chart, v.pos, E[:, col])).vel, Val(N))...)
end

function chart_transition(v::SituatedPhase, E::SMatrix, target_offset)
    (chart_transition(v, target_offset), map_frame(v, E, w -> chart_transition(w, target_offset)))
end

function half_transition(v::SituatedPhase, E::SMatrix)
    (half_transition(v), map_frame(v, E, half_transition))
end

# mirrors settle_phase; keep the hop decisions in sync (they must see the same v)
function settle_transport(env, v::SituatedPhase, E::SMatrix, max_hops=8)
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

function to_ambient(env, v::SituatedPhase, E::SMatrix{3, N}) where {N} # frame columns out via the collar differential
    cl = situate(collar, env, v.chart)
    lin = SMatrix{3, 3, carrier(v.pos[1])}(placement(half_throat(v.chart)).linear)
    (to_ambient(env, v), hcat(ntuple(col -> lin * directional(cl, v.pos, E[:, col]), Val(N))...))
end

"""
trace a geodesic while parallel-transporting the frame E (SMatrix of
tangent-vector columns) along it; returns (AmbientRay, ambient frame) on mouth
exit, (SituatedPhase, E) when out of steps
"""
function trace_transport(env, v::SituatedPhase, E::SMatrix, h, max_steps)
    h = carrier(v.pos[1])(h)
    for _ in 1:max_steps
        v, E = settle_transport(env, transport_step(env, v, E, h)...)
        exits_mouth(v) && return to_ambient(env, v, E)
    end
    (v, E)
end

"""
ray emission from a camera inside the throat: camera = chart point + frame
(columns right, up, forward — any basis, like Camera); x, y in [-1, 1] pick
the pixel direction, and the result is ready for trace_geodesic. Deliberately
NOT normalized: the raw frame scale gives step budgets "computational effort"
semantics (and shows holonomy rescale honestly); pipe through metric_normalize
for geometry-pegged budgets, which is what rendering wants.
"""
function emit_ray(chart, pos, frame, tan_half_fov, x, y)
    C = carrier(pos[1])
    SituatedPhase(chart, pos, frame * SVector(C(x * tan_half_fov), C(y * tan_half_fov), one(C)))
end

function metric_normalize(env, v::SituatedPhase) # unit metric speed: g(vel, vel) = 1
    g = metric(env, v.chart, v.pos)
    SituatedPhase(v.chart, v.pos, v.vel / sqrt(v.vel' * g * v.vel))
end
