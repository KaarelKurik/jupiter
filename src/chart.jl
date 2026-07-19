# chart-level geometry: the AD-friendly wedge complex arithmetic, wedge/square
# coordinates, corner blending, the blended surface and its normal, and chart
# transitions.

# complementary half-angle form for x<0, where the first form is 0/0 near the negative real axis
@inline safe_atan2(y, x) = primal(x) >= 0 ? 2 * atan(y / (sqrt(x^2 + y^2) + x)) : 2 * atan((sqrt(x^2 + y^2) - x) / y)

@inline function fake_complex_pow(z, power)
    iszero(primal(z[1])) && iszero(primal(z[2])) && return SVector(zero(z[1]), zero(z[2])) # safe_atan2 is 0/0 at the origin
    p = carrier(z[1])(power)
    α = safe_atan2(z[2], z[1])
    n = sqrt(z[1]^2 + z[2]^2)
    β = p * α
    s,c = sincos(β)
    n^p * SVector(c, s)
end

function wedge_map(source_n, target_n, pos)
    C = carrier(pos[1])
    a1 = fake_complex_pow(pos, source_n/4)
    a2 = SVector(C(1/sqrt(2)), zero(C)) - a1
    a3 = fake_complex_pow(a2, 4/target_n)
    a3
end

function reference_wedge_map(source_n, target_n, pos) # complex powers don't differentiate; test oracle only
    a0 = complex(pos[1], pos[2])
    a1 = a0^(source_n / 4)
    a2 = 1/sqrt(2) - a1
    a3 = a2^(4 / target_n)
    [real(a3), imag(a3)]
end

@inline function wedge_index(n_wedges, pos)
    α = atan(pos[2], pos[1])
    Int(mod(floor((α * n_wedges)/(2pi)), n_wedges))
end

"""
chart coords -> [0,1]^2 coords of the face at wedge_index,
with u along the wedge's starting edge
"""
@inline function wedge_square_coords(n_wedges, wedge_index, pos)
    C = carrier(pos[1])
    s, c = sincos(C(-2pi * wedge_index / n_wedges)) # angle exact in Float64, one rounding into C
    rpos = SVector(c * pos[1] - s * pos[2], s * pos[1] + c * pos[2])
    C(sqrt(2)) * fake_complex_pow(rpos, n_wedges / 4)
end

@inline function square_coords_to_chart(n_wedges, wedge_index, st)
    C = carrier(st[1])
    z = fake_complex_pow(st / C(sqrt(2)), 4 / n_wedges)
    s, c = sincos(C(2pi * wedge_index / n_wedges))
    SVector(c * z[1] - s * z[2], s * z[1] + c * z[2])
end

"""
same face point as seen from the next corner ccw
"""
next_corner_coords(st) = SVector(st[2], 1 - st[1])

@inline flat_bump(x) = primal(x) > 0 ? exp(-1 / x) : zero(x)

"""
placeholder pending a deliberate choice; the contract is
f(0)=1, f(1)=0, C^inf and flat at both ends
(flat at 0: smoothness across wedge seams and chart centers;
flat at 1: smoothness where a chart's support ends).
f(x)+f(1-x)=1 additionally makes corner weights sum to exactly 1 —
convenience, not load-bearing: surface() divides by the weight total
regardless, so blend candidates need not satisfy it.
"""
@inline function blend_scalar(x)
    flat_bump(1 - x) / (flat_bump(x) + flat_bump(1 - x))
end

@inline function polynomial_surface(chart, uv)
    eval_packed(packed_polys(geometry(half_throat(chart)))[vertex_index(chart)], uv)
end

@inline function corner_contribution(chart, corner, st)
    w = blend_scalar(st[1]) * blend_scalar(st[2])
    c = induced_chart(half_throat(chart), vertex_index(corner))
    (w, polynomial_surface(c, square_coords_to_chart(valence(corner), half_edge_offset(corner), st)))
end

@inline function surface(env, chart, uv)
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

# ---- order-3 jet of the blended surface (hand-rolled derivatives; jets.jl).
# Each function mirrors its plain twin above op-for-op in the value lane.

@inline function wedge_square_coords_cjet(n_wedges, wedge_index, uv)
    C = carrier(uv[1])
    s, c = sincos(C(-2pi * wedge_index / n_wedges)) # same one-rounding angle as wedge_square_coords
    rz = cjet_scale(complex(c, s), cjet_identity(uv)) # complex mult = the same rotation arithmetic
    cjet_scale(C(sqrt(2)), cpow_jet(rz, n_wedges / 4))
end

# (st, 1−s from the next corner) = i + (−i)·w: holomorphic, coefficients exact
@inline next_corner_cjet(j::CJet) = CJet(complex(imag(j.w0), 1 - real(j.w0)),
                                 complex(imag(j.w1), -real(j.w1)),
                                 complex(imag(j.w2), -real(j.w2)),
                                 complex(imag(j.w3), -real(j.w3)))

@inline function square_coords_to_chart_cjet(n_wedges, wedge_index, stj::CJet)
    C = carrier(real(stj.w0))
    z = cpow_jet(cjet_rdiv(stj, C(sqrt(2))), 4 / n_wedges)
    s, c = sincos(C(2pi * wedge_index / n_wedges))
    cjet_scale(complex(c, s), z)
end

@inline function corner_contribution_jet(chart, corner, stj::CJet)
    sj = re_jet(stj)
    tj = im_jet(stj)
    bs = compose1(blend_jet(sj.f)..., sj)
    bt = compose1(blend_jet(tj.f)..., tj)
    w = leibniz(*, bs, bt)
    c = induced_chart(half_throat(chart), vertex_index(corner))
    ζj = square_coords_to_chart_cjet(valence(corner), half_edge_offset(corner), stj)
    pd = eval_packed_partials(packed_polys(geometry(half_throat(chart)))[vertex_index(c)],
                              real(ζj.w0), imag(ζj.w0))
    val = vjet3(wirtinger_compose(pd[1], ζj.w1, ζj.w2, ζj.w3),
                wirtinger_compose(pd[2], ζj.w1, ζj.w2, ζj.w3),
                wirtinger_compose(pd[3], ζj.w1, ζj.w2, ζj.w3))
    (w, val)
end

"""
order-3 (u,v) jet of surface(): value lane bit-identical to surface(), the ten
derivative lanes closed-form. christoffel assembles g and ∂g from this.
"""
@inline function surface_jet(env, chart, uv)
    n = valence(chart)
    k = wedge_index(n, primal.(uv))
    corner = ccw(half_edge_handle(chart), k)
    stj = wedge_square_coords_cjet(n, k, uv)
    w, val = corner_contribution_jet(chart, corner, stj)
    acc = leibniz(*, w, val)
    tw = w
    for _ in 1:3
        corner = next(corner)
        stj = next_corner_cjet(stj)
        w, val = corner_contribution_jet(chart, corner, stj)
        acc = jadd(acc, leibniz(*, w, val))
        tw = jadd(tw, w)
    end
    jdiv(acc, tw)
end

normal_from_columns(jac) = generic_normalize(cross(jac[:, 1], jac[:, 2]))

function surface_normal_out(env, chart, uv)
    _, jac = value_and_jacobian_columns(situate(surface, env, chart), uv, Val(2))
    normal_from_columns(jac)
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

function view_phase_at_target(v, chart_valence, target_offset)
    central_angle = 2pi / chart_valence
    s,c = sincos(carrier(v.pos[1])(-central_angle * target_offset)) # angle exact in Float64, one rounding in
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
