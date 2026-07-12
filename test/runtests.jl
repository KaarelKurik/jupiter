using jupiter
using Test
using LinearAlgebra

const J = jupiter

# Shared fixture: YZ surface fitted to the cube, wrapped in a throat.
m = J.cubemesh()
surf = J.Surface(m, J.fit_geometry(m))
th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75))
ht = J.HalfThroat(th, 1)
c = J.induced_chart(ht, 1)
mouth = J.TessellatedMouth(nothing, ht, 8)

@testset "jupiter" begin

@testset "wedge map" begin
    for (sn, tn) in [(3, 3), (3, 4), (4, 3), (5, 4)]
        @test maximum(abs, J.wedge_map(sn, tn, [0.3, 0.05]) -
                           J.reference_wedge_map(sn, tn, [0.3, 0.05])) < 1e-12
    end
end

@testset "chart transitions round-trip" begin
    pos = [0.3, 0.02, 0.1]
    vel = [0.1, 0.2, 0.05]
    for offset in 0:(J.valence(c) - 1)
        v2 = J.chart_transition(J.SituatedPhase(c, pos, vel), offset)
        back = J.half_edge_offset(J.twin(J.ccw(J.half_edge_handle(c), offset)))
        v3 = J.chart_transition(v2, back)
        @test maximum(abs, v3.pos - pos) < 1e-12
        @test maximum(abs, v3.vel - vel) < 1e-12
    end
end

@testset "blended surface" begin
    # the same face point, evaluated through each of the face's 4 corner charts
    h0 = J.vertex_induced_handle(m, 1)
    corners = [h0; accumulate((x, _) -> J.next(x), 1:3, init=h0)]
    st0 = [0.37, 0.21]
    sts = [[st0]; accumulate((x, _) -> J.next_corner_coords(x), 1:3, init=st0)]
    vals = map(corners, sts) do corner, st
        cc = J.induced_chart(ht, J.vertex_index(corner))
        w = J.square_coords_to_chart(J.valence(corner), J.half_edge_offset(corner), st)
        J.surface(nothing, cc, w)
    end
    @test maximum(maximum(abs, a - b) for a in vals, b in vals) < 1e-12
    @test sum(J.blend_scalar(p[1]) * J.blend_scalar(p[2]) for p in sts) ≈ 1

    @test maximum(abs, J.surface(nothing, c, [0.0, 0.0]) -
                       J.polynomial_surface(c, [0.0, 0.0])) < 1e-12

    nrm = J.surface_normal_out(nothing, c, [0.15, 0.1])
    @test sum(nrm .* J.generic_normalize(m.vertices[1])) > 0.5
end

@testset "packed polynomial evaluation" begin
    # eval_packed against the reference sum of monomials
    for vix in 1:length(m.vertices), uv in ([0.0, 0.0], [0.3, -0.2], [-0.45, 0.4])
        cc = J.induced_chart(ht, vix)
        @test maximum(abs, J.polynomial_surface(cc, uv) -
                           J.Reference.polynomial_surface(cc, uv)) < 1e-13
    end
end

@testset "seam smoothness" begin
    # finite differences up to 4th order along an arc crossing a wedge seam,
    # compared against a no-seam baseline arc: no discontinuity signature
    function diffmax(chart, r, ca, hw, npts, order)
        angs = range(ca - hw, ca + hw, length=npts)
        d = [J.surface(nothing, chart, r .* [cos(a), sin(a)]) for a in angs]
        for _ in 1:order
            d = [d[i+1] - d[i] for i in 1:length(d)-1]
        end
        maximum(maximum(abs, x) for x in d)
    end
    n1 = J.valence(c)
    for k in 1:4
        @test diffmax(c, 0.3, 2pi / n1, 0.05, 201, k) <=
              2 * diffmax(c, 0.3, pi / n1, 0.05, 201, k)
    end
end

@testset "surface sampling" begin
    verts, _ = J.sample_surface(ht, 6)
    @test all(v -> all(isfinite, v), verts)
end

@testset "metric" begin
    pos0 = [0.2, 0.15, 0.0]
    posb = [0.2, 0.15, 0.25]
    posc = [0.2, 0.15, 0.6]
    @test maximum(abs, J.metric(nothing, c, pos0) - J.outer_metric(nothing, c, pos0)) < 1e-14
    @test maximum(abs, J.metric(nothing, c, posc) - J.inner_metric(nothing, c, posc)) < 1e-14
    @test eigmin(Symmetric(J.metric(nothing, c, posb))) > 0
end

@testset "christoffel" begin
    posb = [0.2, 0.15, 0.25]
    v1 = J.SituatedPhase(c, posb, [0.1, 0.2, 0.3])
    cs = J.christoffel(nothing, v1)
    @test maximum(abs(cs[k, i, j] - cs[k, j, i]) for k = 1:3, i = 1:3, j = 1:3) < 1e-12

    h = 1e-5
    mf = p -> J.metric(nothing, c, p)
    dg = [(mf(posb + h * Float64.(1:3 .== j)) - mf(posb - h * Float64.(1:3 .== j))) / (2h)
          for j in 1:3]
    ginv = inv(mf(posb))
    cs_fd = [0.5 * sum(ginv[k, u] * (dg[j][u, i] + dg[i][j, u] - dg[u][i, j]) for u in 1:3)
             for k = 1:3, i = 1:3, j = 1:3]
    @test maximum(abs, cs - cs_fd) < 1e-6
end

@testset "geodesics" begin
    energy(v) = v.vel' * J.metric(nothing, v.chart, v.pos) * v.vel

    # lateral trace crossing into a neighbor chart conserves energy
    v0 = J.settle_phase(nothing, J.SituatedPhase(c, [0.3, 0.02, 0.25], [0.5, 0.2, 0.0]))
    e0 = energy(v0)
    v = v0
    hopped = false
    for _ in 1:200
        v = J.settle_phase(nothing, J.geodesic_step(nothing, v, 0.01))
        hopped |= J.vertex_index(v.chart) != J.vertex_index(v0.chart)
    end
    @test hopped
    @test abs(energy(v) - e0) / e0 < 1e-6

    # settling is a change of description, not of physical point
    vfar = J.SituatedPhase(c, [0.5, 0.05, 0.2], [0.1, 0.2, 0.3])
    vset = J.settle_phase(nothing, vfar)
    @test J.vertex_index(vset.chart) != J.vertex_index(vfar.chart)
    @test maximum(abs, J.surface(nothing, vfar.chart, vfar.pos[1:2]) -
                       J.surface(nothing, vset.chart, vset.pos[1:2])) < 1e-12
    @test vset.pos[3] == vfar.pos[3]
    @test abs(energy(vset) - energy(vfar)) / energy(vfar) < 1e-12

    # half-to-half handover: involution, side flip, exact energy in the cylinder region
    vd = J.SituatedPhase(c, [0.2, 0.15, 0.7], [0.1, 0.05, 0.5])
    vh = J.half_transition(vd)
    vhh = J.half_transition(vh)
    @test J.side(J.half_throat(vh.chart)) == 2
    @test maximum(abs, [vhh.pos - vd.pos; vhh.vel - vd.vel]) == 0
    @test abs(energy(vh) - energy(vd)) < 1e-14

    # mouth exit: outward geodesic becomes an AmbientRay whose flat speed
    # matches the conserved energy (collar pullback is isometric)
    vout = J.SituatedPhase(c, [0.2, 0.15, 0.2], [0.05, 0.05, -0.6])
    r = J.trace_geodesic(nothing, vout, 0.01, 200)
    @test r isa J.AmbientRay
    @test J.side(J.half_throat(r)) == 1
    @test abs(sum(abs2, r.vel) - energy(vout)) / energy(vout) < 1e-6

    # full traversal: in through mouth 1, out through mouth 2
    vin = J.SituatedPhase(c, [0.2, 0.15, 0.01], [0.05, 0.05, 0.8])
    r2 = J.trace_geodesic(nothing, vin, 0.02, 400)
    @test r2 isa J.AmbientRay
    @test J.side(J.half_throat(r2)) == 2
    @test abs(sum(abs2, r2.vel) - energy(vin)) / energy(vin) < 1e-6

    # reversibility over a short arc that includes a chart hop
    va = J.settle_phase(nothing, J.SituatedPhase(c, [0.3, 0.02, 0.25], [0.5, 0.2, 0.1]))
    vb = va
    for _ in 1:20
        vb = J.settle_phase(nothing, J.geodesic_step(nothing, vb, 0.01))
    end
    vr = J.SituatedPhase(vb.chart, vb.pos, -vb.vel)
    for _ in 1:20
        vr = J.settle_phase(nothing, J.geodesic_step(nothing, vr, 0.01))
    end
    @test J.vertex_index(vr.chart) == J.vertex_index(va.chart)
    @test maximum(abs, vr.pos - va.pos) < 1e-9
    @test maximum(abs, vr.vel + va.vel) < 1e-9
end

@testset "bvh" begin
    function brute_hit(origin, dir)
        best_t = Inf
        best = nothing
        for tri in mouth.triangles
            hit = J.ray_triangle_intersection(origin, dir, tri.corners...)
            (hit === nothing || hit[1] >= best_t) && continue
            best_t = hit[1]
            best = (hit, tri)
        end
        best
    end
    for k in 1:40
        θ = 2pi * k / 40
        φ = pi * mod(k, 7) / 7 - pi / 2
        origin = 3 * [cos(θ) * cos(φ), sin(θ) * cos(φ), sin(φ)]
        target = 0.4 * [sin(3k), cos(2k), sin(5k)]
        dir = J.generic_normalize(target - origin)
        bf = brute_hit(origin, dir)
        bvh = J.nearest_mouth_hit(mouth, origin, dir)
        @test (bf === nothing) == (bvh === nothing)
        bf === nothing || @test bf[1][1] == bvh[1][1]
    end
end

@testset "mouth entry" begin
    uv0 = [0.1, 0.05]
    pos0 = [uv0; 0.0]
    vel0 = [0.02, 0.01, 0.6]
    # build the entering ray by exiting in reverse, then check entry inverts it
    r = J.to_ambient(nothing, J.SituatedPhase(c, pos0, -vel0))
    dir = -r.vel
    origin = r.pos + 2 * r.vel
    vin = J.enter_mouth(nothing, mouth, origin, dir)
    @test vin isa J.SituatedPhase
    @test maximum(abs, J.surface(nothing, vin.chart, vin.pos[1:2]) -
                       J.surface(nothing, c, uv0)) < 1e-9
    @test abs(vin.pos[3]) < 1e-9
    rr = J.to_ambient(nothing, J.SituatedPhase(vin.chart, vin.pos, -vin.vel))
    @test maximum(abs, rr.pos - r.pos) < 1e-9
    @test maximum(abs, rr.vel - r.vel) < 1e-9
    # flat speed^2 matches chart-metric energy (collar isometry)
    @test abs(vin.vel' * J.metric(nothing, vin.chart, vin.pos) * vin.vel -
              sum(abs2, dir)) < 1e-12
    # ray pointing away from the surface misses
    @test J.enter_mouth(nothing, mouth, origin, -dir) === nothing
    # full pipeline: ambient ray enters side 1, exits as an ambient ray on side 2
    r2 = J.trace_geodesic(nothing, vin, 0.02, 400)
    @test r2 isa J.AmbientRay
    @test J.side(J.half_throat(r2)) == 2
end

@testset "renderer" begin
    mouth2 = J.TessellatedMouth(nothing, J.HalfThroat(th, 2), 8)
    scene = J.Scene(th, (mouth, mouth2), J.checker_sky)
    cam = J.look_at_camera([0.0, -3.0, 0.4], zeros(3), [0.0, 0.0, 1.0], pi / 4)

    # a ray aimed at the mouth escapes to the other side's sky within budget
    hit_ray = J.camera_ray(cam, ht, 0.0, 0.0)
    out = J.trace_ray(nothing, scene, J.RayBudget(0.05, 300, 4), hit_ray)
    @test out isa J.AmbientRay
    @test J.side(J.half_throat(out)) == 2

    # the same ray under a starved step budget is unresolved
    @test J.trace_ray(nothing, scene, J.RayBudget(0.05, 5, 4), hit_ray) === nothing

    # a ray aimed well away never enters and comes back unchanged
    miss_ray = J.AmbientRay(ht, [0.0, -3.0, 0.4], [0.0, -1.0, 0.3])
    back = J.trace_ray(nothing, scene, J.RayBudget(0.05, 300, 4), miss_ray)
    @test back === miss_ray

    # tiny render: raymap + shade agrees with direct render, all pixels valid
    raymap = J.render_raymap(nothing, scene, J.RayBudget(0.05, 300, 4), cam, 1, 8, 6)
    img = J.shade(raymap, J.checker_sky)
    @test img == J.render(nothing, scene, J.RayBudget(0.05, 300, 4), cam, 1, 8, 6)
    @test all(isfinite, img)
    @test all(0 .<= img .<= 1)

    # raymap survives a disk round trip, and re-shading is trace-free
    path = tempname()
    J.save_raymap(path, raymap)
    raymap2 = J.load_raymap(path)
    @test raymap2.side == raymap.side
    @test raymap2.pos == raymap.pos
    @test raymap2.vel == raymap.vel
    flat_sky(side, dir) = side == 1 ? [1.0, 1.0, 1.0] : [0.5, 0.5, 0.5]
    img2 = J.shade(raymap2, flat_sky)
    @test size(img2) == size(img)
    @test all(0 .<= img2 .<= 1)
end

@testset "reference equivalence" begin
    # the optimized core against src/reference.jl, pointwise and along traces;
    # any deliberate semantic change must move both in lockstep
    R = J.Reference

    uvs = ([0.15, 0.1], [0.3, 0.3], [-0.25, 0.2], [0.1, -0.3], [-0.2, -0.15])
    for uv in uvs
        @test maximum(abs, R.surface(nothing, c, uv) - J.surface(nothing, c, uv)) < 1e-12
    end
    for uv in uvs[1:3], d in (-0.05, 0.1, 0.3, 0.6)
        pos = [uv; d]
        @test maximum(abs, R.metric(nothing, c, pos) - J.metric(nothing, c, pos)) < 1e-11
    end

    for pos in ([0.2, 0.15, 0.25], [0.3, 0.02, 0.45], [-0.2, 0.1, 0.1])
        v = J.SituatedPhase(c, pos, [0.1, 0.2, 0.3])
        @test maximum(abs, R.christoffel(nothing, v) - J.christoffel(nothing, v)) < 1e-9
        sa = J.geodesic_step(nothing, v, 0.02)
        sb = R.geodesic_step(nothing, v, 0.02)
        @test maximum(abs, [sa.pos - sb.pos; sa.vel - sb.vel]) < 1e-11
    end

    # transitions and settling agree
    vfar = J.SituatedPhase(c, [0.5, 0.05, 0.2], [0.1, 0.2, 0.3])
    for offset in 0:(J.valence(c) - 1)
        ta = J.chart_transition(vfar, offset)
        tb = R.chart_transition(vfar, offset)
        @test J.vertex_index(ta.chart) == J.vertex_index(tb.chart)
        @test maximum(abs, [ta.pos - tb.pos; ta.vel - tb.vel]) < 1e-12
    end
    sa = J.settle_phase(nothing, vfar)
    sb = R.settle_phase(nothing, vfar)
    @test J.vertex_index(sa.chart) == J.vertex_index(sb.chart)
    @test maximum(abs, [sa.pos - sb.pos; sa.vel - sb.vel]) < 1e-12

    va = J.to_ambient(nothing, J.SituatedPhase(c, [0.2, 0.15, 0.0], [0.1, 0.2, -0.3]))
    vb = R.to_ambient(nothing, J.SituatedPhase(c, [0.2, 0.15, 0.0], [0.1, 0.2, -0.3]))
    @test maximum(abs, [va.pos - vb.pos; va.vel - vb.vel]) < 1e-12

    # side-by-side traces: lateral with chart hops, and across the half handover
    function both_trace(v0, n, h)
        va = vb = v0
        for _ in 1:n
            va = J.settle_phase(nothing, J.geodesic_step(nothing, va, h))
            vb = R.settle_phase(nothing, R.geodesic_step(nothing, vb, h))
        end
        (va, vb)
    end
    for (v0, n) in ((J.SituatedPhase(c, [0.3, 0.02, 0.25], [0.5, 0.2, 0.0]), 60),
                    (J.SituatedPhase(c, [0.2, 0.15, 0.6], [0.05, 0.05, 0.5]), 60))
        ta, tb = both_trace(J.settle_phase(nothing, v0), n, 0.01)
        @test J.vertex_index(ta.chart) == J.vertex_index(tb.chart)
        @test J.side(J.half_throat(ta.chart)) == J.side(J.half_throat(tb.chart))
        @test maximum(abs, [ta.pos - tb.pos; ta.vel - tb.vel]) < 1e-9
    end

    # full drivers agree through a mouth exit
    vout = J.SituatedPhase(c, [0.2, 0.15, 0.2], [0.05, 0.05, -0.6])
    ra = J.trace_geodesic(nothing, vout, 0.01, 200)
    rb = R.trace_geodesic(nothing, vout, 0.01, 200)
    @test rb isa J.AmbientRay
    @test J.side(J.half_throat(ra)) == J.side(J.half_throat(rb))
    @test maximum(abs, [ra.pos - rb.pos; ra.vel - rb.vel]) < 1e-9
end

@testset "metric tensoriality across transition" begin
    pe = [0.3, 0.02, 0.25]
    vp(k) = J.SituatedPhase(c, pe, Float64.(1:3 .== k))
    t1 = J.chart_transition(vp(1), 0)
    Jm = reduce(hcat, [J.chart_transition(vp(k), 0).vel for k in 1:3])
    @test maximum(abs, Jm' * J.metric(nothing, t1.chart, t1.pos) * Jm -
                       J.metric(nothing, c, pe)) < 1e-12
end

end
