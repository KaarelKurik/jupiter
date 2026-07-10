using jupiter
const J = jupiter

m = J.cubemesh()
surf = J.Surface(m, J.fit_geometry(m))
th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75))
mouths = (J.TessellatedMouth(nothing, J.HalfThroat(th, 1), 10),
          J.TessellatedMouth(nothing, J.HalfThroat(th, 2), 10))
scene = J.Scene(th, mouths, J.checker_sky)
budget = J.RayBudget(0.05, 400, 4)
cam = J.look_at_camera([0.9, -2.8, 0.9], zeros(3), [0.0, 0.0, 1.0], pi / 3.2)

w, h = 192, 144
println("tracing ", w, "x", h, " on ", Threads.nthreads(), " thread(s)...")
t0 = time()
raymap = J.render_raymap(nothing, scene, budget, cam, 1, w, h)
println("traced in ", round(time() - t0, digits=1), "s;  unresolved pixels: ",
        count(==(0), raymap.side))
J.save_raymap(joinpath(@__DIR__, "..", "res", "first_light.raymap"), raymap)
J.save_ppm(joinpath(@__DIR__, "..", "res", "first_light.ppm"), J.shade(raymap, J.checker_sky))
