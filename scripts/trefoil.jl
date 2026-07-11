using jupiter
const J = jupiter

m = J.load_obj(joinpath(@__DIR__, "..", "res", "models", "trefoil.obj"))
surf = J.Surface(m, J.fit_geometry(m))
# tube radius ~0.5, so the collar must stay well inside the cross-section curvature
th = J.Throat(surf, J.ThroatParams(0.4, 1.0, 0.2, 0.3))
mouths = (J.TessellatedMouth(nothing, J.HalfThroat(th, 1), 2),
          J.TessellatedMouth(nothing, J.HalfThroat(th, 2), 2))
scene = J.Scene(th, mouths, J.checker_sky)
budget = J.RayBudget(0.05, 400, 4)
cam = J.look_at_camera([0.5, -7.5, 3.5], zeros(3), [0.0, 0.0, 1.0], pi / 3.2)

w, h = parse.(Int, get(ARGS, 1, "192x144") |> s -> split(s, "x"))
println("tracing ", w, "x", h, " on ", Threads.nthreads(), " thread(s)...")
t0 = time()
raymap = J.render_raymap(nothing, scene, budget, cam, 1, w, h)
println("traced in ", round(time() - t0, digits=1), "s;  unresolved pixels: ",
        count(==(0), raymap.side))
mkpath(joinpath(@__DIR__, "..", "out"))
J.save_raymap(joinpath(@__DIR__, "..", "out", "trefoil.raymap"), raymap)
J.save_ppm(joinpath(@__DIR__, "..", "out", "trefoil.ppm"), J.shade(raymap, J.checker_sky))
