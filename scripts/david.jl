# David as the throat surface: a Michelangelo-shaped hole between two universes.
#
#   julia --project --threads=auto scripts/david.jl [scale]     # scale = pixel scale, default 1
#
# Traces one raymap per view into out/, then shades each under the checker sky
# and the space equirects (the cheap pass — trace once, shade many). For the
# textured shading, convert the tracked JPGs once per checkout:
#   magick res/textures/space1.jpg out/space1.ppm
#   magick res/textures/space2.jpg out/space2.ppm
#
# Mesh prep lives here rather than in src/ on purpose: res/models/david_small.stl
# is an unwelded triangle soup (4248 STL vertices, 700 distinct positions) and
# fit_geometry is classic YZ, which wants quads. Welding by exact position plus
# one Catmull-Clark pre-subdivision to all-quads is enough to make it a Surface —
# that is the "David mesh path" plan item, done at scene-authoring level.
#
# Throat params differ from trefoil's: the collar depth must stay inside the
# local cross-section curvature, and David's thinnest features (wrist, ankle,
# the sling cord) are much finer than the trefoil's ~0.5 tube. cylinder_depth
# 0.08 at scale 2.0 is the deepest that traces clean. The tracer is the
# error-controlled DP5 branch (tolerance > 0) — the fixed-step RK4 branch
# produces NaN velocities on a few hundredths of a percent of rays here, where
# the fitted charts are worst-conditioned over the mesh's sliver quads. Under
# DP5 those rays are honestly reported as unresolved instead.
using jupiter
using FileIO
using GeometryBasics
using LinearAlgebra
const J = jupiter

const root = joinpath(@__DIR__, "..")

"""
welded, quad-remeshed, centered David: STL triangle soup -> half-edge Mesh
"""
function david_mesh(; scale = 2.0)
    raw = load(joinpath(root, "res", "models", "david_small.stl"))
    verts = coordinates(raw)
    key(p) = (Float64(p[1]), Float64(p[2]), Float64(p[3])) # STL stores exact duplicates, so exact keys weld it
    ids = Dict{NTuple{3, Float64}, Int}()
    welded = Vector{Float64}[]
    for p in verts
        k = key(p)
        haskey(ids, k) || (push!(welded, collect(k)); ids[k] = length(welded))
    end
    faces = [[ids[key(verts[GeometryBasics.value(f[i])])] for i in 1:3]
             for f in GeometryBasics.faces(raw)]
    lo = reduce((a, b) -> min.(a, b), welded)
    hi = reduce((a, b) -> max.(a, b), welded)
    center = (lo .+ hi) ./ 2
    welded = [scale .* (v .- center) for v in welded]
    J.catmullclark(J.Mesh(welded, faces)).refined_mesh # triangles -> all quads
end

px = parse(Float64, get(ARGS, 1, "1"))
res(w, h) = (round(Int, px * w), round(Int, px * h))

# camera pos, look-at target, fov, resolution — David spans z in [-4.74, 4.74]
views = [("full", [0.0, -16.0, 0.0], [0.0, 0.0, 0.0], pi / 3.2, res(900, 1200)),
         ("bust", [1.5, -5.0, 3.2],  [0.0, 0.0, 3.4], pi / 3.2, res(1000, 1000)),
         ("head", [1.0, -2.6, 4.3],  [0.0, 0.0, 4.1], pi / 3.5, res(900, 900))]

mesh = david_mesh()
t0 = time()
surface = J.Surface(mesh, J.fit_geometry(mesh))
println("fitted ", length(mesh.vertices), " charts in ", round(time() - t0, digits=1), "s")

throat = J.Throat(surface, J.ThroatParams(0.4, 1.0, 0.08, 0.12))
mouths = (J.TessellatedMouth(nothing, J.HalfThroat(throat, 1), 2),
          J.TessellatedMouth(nothing, J.HalfThroat(throat, 2), 2))
scene = J.Scene(throat, mouths, J.checker_sky)
budget = J.RayBudget(0.05, 4000, 4, 1e-8)

skies = Pair{String, Any}["checker" => J.checker_sky]
if isfile(joinpath(root, "out", "space1.ppm")) && isfile(joinpath(root, "out", "space2.ppm"))
    push!(skies, "space" => J.TexturedSky(J.load_ppm(joinpath(root, "out", "space1.ppm")),
                                          J.load_ppm(joinpath(root, "out", "space2.ppm"))))
else
    println("out/space[12].ppm missing — checker only (see header for the magick lines)")
end

mkpath(joinpath(root, "out"))
for (name, pos, target, fov, (w, h)) in views
    camera = J.look_at_camera(pos, target, [0.0, 0.0, 1.0], fov)
    local started = time()
    raymap = J.render_raymap(nothing, scene, budget, camera, 1, w, h)
    println(name, " ", w, "x", h, ": traced in ", round(time() - started, digits=1), "s;  ",
            "unresolved pixels: ", count(==(0), raymap.side))
    J.save_raymap(joinpath(root, "out", "david_$name.raymap"), raymap)
    for (sky_name, sky) in skies
        J.save_ppm(joinpath(root, "out", "david_$(name)_$(sky_name).ppm"), J.shade(raymap, sky))
    end
    flush(stdout)
end

# undeflected outline of the full view, to hold the lensed renders against
J.save_ppm(joinpath(root, "out", "david_full_flat.ppm"),
           J.render_flat(nothing, scene, J.look_at_camera(views[1][2], views[1][3],
                                                          [0.0, 0.0, 1.0], views[1][4]),
                         1, views[1][5]...))
