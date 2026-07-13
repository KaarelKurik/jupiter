# YZ chart fitting: least-squares polynomial charts per vertex from two rounds
# of Catmull-Clark refinement, sampled at wedge grid points.

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
