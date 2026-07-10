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
img = zeros(3, w, h)
cam_space = J.HalfThroat(th, 1)
aspect = w / h
t0 = time()
unresolved = 0
for j in 1:h
    for i in 1:w
        x = (2 * (i - 0.5) / w - 1) * aspect
        y = 1 - 2 * (j - 0.5) / h
        out = J.trace_ray(nothing, scene, budget, J.camera_ray(cam, cam_space, x, y))
        if out === nothing
            global unresolved += 1
            img[:, i, j] = J.unresolved_color()
        else
            img[:, i, j] = scene.sky(J.side(J.half_throat(out)), J.generic_normalize(out.vel))
        end
    end
    j % 12 == 0 && println("row ", j, "/", h, "  elapsed ", round(time() - t0, digits=1), "s")
end
println("render done in ", round(time() - t0, digits=1), "s;  unresolved pixels: ", unresolved)
J.save_ppm(joinpath(@__DIR__, "..", "res", "first_light.ppm"), img)
