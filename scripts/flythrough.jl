# fly-through demo: coast a camera from side 1 of the cube throat straight
# through to side 2 under the textured skies, rendering a keyframe after each
# leg; the last frame turns around to look back at the mouth just left behind.
#   julia --project --threads=auto scripts/flythrough.jl [WxH]
using jupiter
using LinearAlgebra
const J = jupiter

root = joinpath(@__DIR__, "..")
w, h = parse.(Int, get(ARGS, 1, "192x144") |> s -> split(s, "x"))

m = J.cubemesh()
surf = J.Surface(m, J.fit_geometry(m))
th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75))
mouths = (J.TessellatedMouth(nothing, J.HalfThroat(th, 1), 10),
          J.TessellatedMouth(nothing, J.HalfThroat(th, 2), 10))
sky = J.TexturedSky(J.load_ppm(joinpath(root, "res", "textures", "sky1.ppm")),
                    J.load_ppm(joinpath(root, "res", "textures", "sky2.ppm")))
scene = J.Scene(th, mouths, sky)
budget = J.RayBudget(0.05, 400, 4)
thf = tan(pi / 3.2 / 2)

campos = J.SVector(0.9, -2.8, 0.9)
fwd = J.generic_normalize(-campos)
right = J.generic_normalize(cross(fwd, [0.0, 0.0, 1.0]))
cam = J.FlyingCamera(J.AmbientRay(J.HalfThroat(th, 1), campos, fwd),
                     J.SMatrix{3, 3, Float64}([right cross(right, fwd) fwd]))

where(c) = c.state isa J.AmbientRay ?
    "ambient side $(J.side(J.half_throat(c.state)))" :
    "in-chart d=$(round(c.state.pos[3], digits=2))"

mkpath(joinpath(root, "out"))
legs = fill(0.55, 8)
for (k, Δτ) in enumerate(legs)
    global cam = J.coast(nothing, scene, cam, Δτ, 0.02)
    t0 = time()
    raymap = J.keyframe_raymap(nothing, scene, budget, cam, thf, w, h)
    println("frame ", k, " (", where(cam), "): traced in ", round(time() - t0, digits=1),
            "s, unresolved ", count(==(0), raymap.side))
    J.save_ppm(joinpath(root, "out", "fly_$k.ppm"), J.shade(raymap, sky))
end

# turn around (yaw 180 about up, frame-relative) and look back at the mouth
cam = J.maneuver(cam, J.SMatrix{3, 3, Float64}([-1 0 0; 0 1 0; 0 0 -1]))
raymap = J.keyframe_raymap(nothing, scene, budget, cam, thf, w, h)
println("frame 9 (", where(cam), ", looking back): unresolved ", count(==(0), raymap.side))
J.save_ppm(joinpath(root, "out", "fly_9.ppm"), J.shade(raymap, sky))
