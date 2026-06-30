
greet() = print("Hello World!")

# Let's just make half-edges a thick concept, fuck it

struct Surface
    mesh
    chart_polys
end

struct HalfThroat
    surface
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

function wedge_index_and_angle(n_wedges, pos)
    α = atan(pos[2], pos[1])
    ix = floor((α * n_wedges)/(2pi))
    remainder = α - ix * (2pi / n_wedges)
    (ix, remainder)
end

function polynomial_surface(chart, uv)
    @polyvar u v
    [p(u => uv[1], v => v[2]) for p in surface_polynomials(chart)]
end

function surface(env, chart, uv)
    #todo by blending
end

function surface_polynomials(chart::Chart)
    half_throat(chart).chart_polys[vertex_index(chart)]
end

function surface_normal_out(env, chart, uv)
    s = Base.Fix1(Base.Fix1(surface, env), chart)
    du = derivative(s, uv, [1., 0.], Val(1))
    dv = derivative(s, uv, [0., 1.], Val(1))
    normalize(cross(du, dv)) # normalize probably doesn't work with TaylorDiff
end

function collar(env, chart, pos)
    uv = pos[1:2]
    surface(env, chart, uv) - pos[3] * surface_normal_out(env, chart, uv)
end

function outer_metric(env, chart, pos)
    cl = Base.Fix1(Base.Fix1(collar, env), chart)
    collar_jac = reduce(hcat, [derivative(cl, pos, [Float64(i == j) for i in 1:3], Val(1)) for j in 1:3])
    collar_jac' * collar_jac
end

function inner_metric(env, chart, pos)
    sf = Base.Fix1(Base.Fix1(surface, env), chart)
    params = inner_metric_params(env, chart)
    sf_jac = reduce(hcat, [derivative(sf, pos[1:2], [Float64(i == j) for i in 1:2], Val(1)) for j in 1:2])
    g = params.cross_scale * sf_jac' * sf_jac
    out = zeros(promote_type(eltype(g), typeof(params.depth_scale)), 3, 3)
    out[1:2, 1:2] = g
    out[3,3] = params.depth_scale
    out
end

function metric(env, chart, pos)
    om = outer_metric(env, chart, pos)
    im = inner_metric(env, chart, pos)
    depth_interpolate(pos[3], om, im)
end

function christoffel(env, v::SituatedPhase)
    mf = Base.Fix1(Base.Fix1(metric, env), v.chart)
    metric_derivs = stack([derivative(mf, v.pos, [Float64(i == j) for i in 1:3], Val(1)) for j in 1:3])
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

function reference_wedge_map(source_n, target_n, pos) # this doesn't work with TaylorDiff
    a0 = complex(pos[1], pos[2])
    a1 = a0^(source_n / 4)
    a2 = 1/sqrt(2) - a1
    a3 = a2^(4 / target_n)
    [real(a3), imag(a3)]
end

safe_atan2(y, x) = 2 * atan(y / (sqrt(x^2 + y^2) + x))

function fake_complex_pow(z, power)
    α = safe_atan2(z[2], z[1])
    n = sqrt(z[1]^2 + z[2]^2)
    β = power * α
    s,c = sincos(β)
    n^power * [c, s]
end

function wedge_map(source_n, target_n, pos) # ideally this would work with TaylorDiff
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

function mesh(half_throat::HalfThroat)
    half_throat.mesh
end

function mesh(chart::Chart)
    mesh(half_throat(chart))
end

function half_edge_handle(chart::Chart)
    chart.half_edge_handle
end

function induced_chart(half_throat::HalfThroat, vertex_index::Int)
    half_throat.charts[vertex_index]
end

function induced_chart(half_edge_handle::HalfEdgeHandle)
    #todo
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
    neighbor = induced_chart(shared_edge_far)
    source_n = valence(v.chart)
    target_n = valence(neighbor)
    v_near_adj = view_phase_at_target(v, source_n, target_offset)
    wm = Base.Fix1(Base.Fix1(wedge_map, source_n), target_n)
    v_far_adj_pos = [wm(v_near_adj.pos) ; v_near_adj.pos[3]]
    v_far_adj_vel = [derivative(wm, v_near_adj.pos, v_near_adj.vel, Val(1)) ; v_near_adj.vel[3]]
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