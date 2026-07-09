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

@testset "metric tensoriality across transition" begin
    pe = [0.3, 0.02, 0.25]
    vp(k) = J.SituatedPhase(c, pe, Float64.(1:3 .== k))
    t1 = J.chart_transition(vp(1), 0)
    Jm = reduce(hcat, [J.chart_transition(vp(k), 0).vel for k in 1:3])
    @test maximum(abs, Jm' * J.metric(nothing, t1.chart, t1.pos) * Jm -
                       J.metric(nothing, c, pe)) < 1e-12
end

end
