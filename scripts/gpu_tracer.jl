# device tracer kernel measurement (GPU step 5): sweep_ray in-kernel — one
# thread marches one ray through settled RK4 attempts, chart hops, half
# transitions, and the to_ambient exit — with the host running entry solves
# and compaction between launches (jupitergpu.DeviceSweep inside the same
# run_wavefront! driver the CPU uses). Parity is judged raymap-vs-raymap
# against the CPU wavefront at the same precision (transcendental intrinsics
# differ by ulps, so bit-identity is not expected — the 14c decomposition
# frame applies: side flips + conditional angle quantiles); perf is wall
# time against the recursive renderer and the CPU wavefront, plus
# kernel-only seconds, occupancy stats, and the bin lever.
#
#   julia --project=gpu --threads=auto scripts/gpu_tracer.jl [WxH]

using jupiter
using CUDA
using LinearAlgebra
include(joinpath(@__DIR__, "..", "gpu", "jupitergpu.jl"))
using .jupitergpu
const J = jupiter
const G = jupitergpu

f32_surface(s) = J.Surface(J.mesh(s), J.chart_polys(s), [Float32.(p) for p in J.packed_polys(s)])
f32_throat(th) = J.Throat(f32_surface(J.geometry(th)), J.params(th), th.placements)

function angle_between(v1, v2)
    n1, n2 = J.generic_normalize(v1), J.generic_normalize(v2)
    2 * asin(clamp(norm(n1 - n2) / 2, 0.0, 1.0))
end

quant(sorted, p) = sorted[clamp(ceil(Int, p * length(sorted)), 1, length(sorted))]
fmt(x) = x == 0 ? "0" : string(round(x, sigdigits = 3))

function compare(name, a::J.RayMap, b::J.RayMap)
    n = length(a.side)
    flips = count(@. a.side != b.side && a.side != 0 && b.side != 0)
    res_mismatch = count(@. (a.side == 0) != (b.side == 0))
    agree = [(i, j) for j in axes(a.side, 2), i in axes(a.side, 1) if a.side[i, j] == b.side[i, j] != 0]
    angles = sort([angle_between(a.vel[:, t[1], t[2]], b.vel[:, t[1], t[2]]) for t in agree])
    identical = count(t -> a.vel[:, t[1], t[2]] == b.vel[:, t[1], t[2]] &&
                           a.pos[:, t[1], t[2]] == b.pos[:, t[1], t[2]], agree)
    println("  ", rpad(name, 24), " flips ", flips, "/", n, "  res-mismatch ", res_mismatch,
            "  bitwise-equal ", identical, "/", length(agree),
            "  | angle q50/q99/max ",
            isempty(angles) ? "-" : join((fmt(quant(angles, 0.5)), fmt(quant(angles, 0.99)),
                                          fmt(angles[end])), " "))
end

function main()
    wh = length(ARGS) >= 1 ? parse.(Int, split(ARGS[1], "x")) : [192, 144]
    w, h = wh
    m = J.cubemesh()
    surf = J.Surface(m, J.fit_geometry(m))
    th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75)) # the physics_diff cube
    th32 = f32_throat(th)
    sky = J.checker_sky
    mouths(t) = (J.TessellatedMouth(nothing, J.HalfThroat(t, 1), 10),
                 J.TessellatedMouth(nothing, J.HalfThroat(t, 2), 10))
    scene = J.Scene(th, mouths(th), sky)
    scene32 = J.Scene(th32, mouths(th32), sky) # F64 entry solves against the F32-coefficient surface (the 14c seam)
    cam = J.look_at_camera([0.9, -2.8, 0.9], zeros(3), [0.0, 0.0, 1.0], pi / 3)
    budget = J.RayBudget(0.05, 400, 4)
    println("gpu_tracer: ", w, "x", h, " on ", CUDA.name(CUDA.device()),
            ", host ", Threads.nthreads(), " threads")

    # kernel occupancy stats per precision, on a token record vector
    for T in (Float64, Float32)
        r0 = J.WavefrontRay{T}(J.RAY_THROAT, Int8(1), Int32(m.vertex_he[1]), Int32(1), Int32(0),
                               T(0.05), J.SVector{3, T}(0.1, 0.05, 0.2), J.SVector{3, T}(0.05, 0.02, 0.4))
        dth = G.device_throat(T === Float32 ? th32 : th, T)
        kc = @cuda launch = false G.kernel_sweep(CuArray([r0]), dth, budget, 8)
        mem = CUDA.memory(kc)
        println(T, " kernel: ", CUDA.registers(kc), " registers, ", mem.local, " local (spill) bytes/thread")
    end

    # warm the device path small, then parity + perf
    G.device_raymap(nothing, scene, budget, cam, 1, 16, 12; T = Float64, sweep = 16)
    G.device_raymap(nothing, scene32, budget, cam, 1, 16, 12; T = Float32, sweep = 16)

    println("parity (device vs CPU wavefront, same precision + tables):")
    rm64 = J.wavefront_raymap(nothing, scene, budget, cam, 1, w, h)
    dev64, _ = G.device_raymap(nothing, scene, budget, cam, 1, w, h; T = Float64, sweep = 64)
    compare("F64 device vs cpu", rm64, dev64)
    rm32 = J.wavefront_raymap(nothing, scene32, budget, cam, 1, w, h; T = Float32)
    dev32, _ = G.device_raymap(nothing, scene32, budget, cam, 1, w, h; T = Float32, sweep = 64)
    compare("F32 device vs cpu", rm32, dev32)
    compare("F32 device vs F64 truth", rm64, dev32)

    println("perf (wall seconds, best of 2):")
    best2(f) = min(f(), f())
    tr = best2(() -> @elapsed J.render_raymap(nothing, scene, budget, cam, 1, w, h))
    println("  recursive cpu:            ", fmt(tr))
    tw = best2(() -> @elapsed J.wavefront_raymap(nothing, scene, budget, cam, 1, w, h))
    println("  wavefront cpu F64:        ", fmt(tw))
    tw32 = best2(() -> @elapsed J.wavefront_raymap(nothing, scene32, budget, cam, 1, w, h; T = Float32))
    println("  wavefront cpu F32:        ", fmt(tw32))
    for (T, sc) in ((Float64, scene), (Float32, scene32)), (sweep, bin) in
        ((64, true), (64, false), (256, true), (16, true))
        ks = 0.0
        td = best2() do
            t = @elapsed ((_, ks) = G.device_raymap(nothing, sc, budget, cam, 1, w, h; T, sweep, bin))
            t
        end
        println("  device ", T, " sweep=", rpad(sweep, 3), " bin=", rpad(bin, 5), ": ", fmt(td),
                "  (kernel-only ", fmt(ks), ", ", fmt(td / tr), "x recursive wall)")
    end
end

main()
