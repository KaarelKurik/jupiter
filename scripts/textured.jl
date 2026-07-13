# re-shade a saved raymap under the res/textures skies — the cheap pass, no
# tracing. Usage: julia --project scripts/textured.jl [raymap] [out.ppm]
using jupiter
const J = jupiter

root = joinpath(@__DIR__, "..")
raymap_path = get(ARGS, 1, joinpath(root, "out", "trefoil.raymap"))
out_path = get(ARGS, 2, joinpath(root, "out", "textured.ppm"))

raymap = J.load_raymap(raymap_path)
sky = J.TexturedSky(J.load_ppm(joinpath(root, "res", "textures", "sky1.ppm")),
                    J.load_ppm(joinpath(root, "res", "textures", "sky2.ppm")))
t0 = time()
img = J.shade(raymap, sky)
println("shaded ", join(size(raymap.side), "x"), " raymap in ", round(time() - t0, digits=2), "s")
J.save_ppm(out_path, img)
println("wrote ", out_path)
