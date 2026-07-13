# the type layer: Surface/Throat/chart/phase types and their thin accessors.
# Thick data lives on Surface; HalfThroat and Chart are handles.

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

function eval_packed(packed::AbstractArray{Float64, 3}, uv) # Abstract so device-side array types dispatch here too
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

function half_throat(ray::AmbientRay)
    ray.half_throat
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

function surface_polynomials(chart::Chart)
    chart_polys(geometry(half_throat(chart)))[vertex_index(chart)]
end
