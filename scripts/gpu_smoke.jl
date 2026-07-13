# device smoke test (GPU groundwork step 3): compile the hot arithmetic on the
# GPU and compare against the CPU — a wall-finder, not a port. Kernel A:
# packed-poly Horner surface evaluation (eval_packed). Kernel B: the ad.jl
# dual pass (value_and_jacobian_columns) through the wedge transcendentals
# (square_coords_to_chart -> fake_complex_pow/safe_atan2) — ForwardDiff duals
# and StaticArrays on device. Float64 throughout; the Float32 policy is a
# later, deliberate step.
#
#   julia --project=gpu scripts/gpu_smoke.jl [npoints]

using jupiter
using CUDA
const J = jupiter

function kernel_surface(out, packed, uvs)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i > size(uvs, 2) && return nothing
    p = J.eval_packed(packed, J.SVector(uvs[1, i], uvs[2, i]))
    out[1, i] = p[1]; out[2, i] = p[2]; out[3, i] = p[3]
    nothing
end

function kernel_dualpass(vals, jacs, packed, n, wedge, sts)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i > size(sts, 2) && return nothing
    f = w -> J.eval_packed(packed, J.square_coords_to_chart(n, wedge, w))
    v, jac = J.value_and_jacobian_columns(f, J.SVector(sts[1, i], sts[2, i]), Val(2))
    for c in 1:3
        vals[c, i] = v[c]
        jacs[c, 1, i] = jac[c, 1]
        jacs[c, 2, i] = jac[c, 2]
    end
    nothing
end

function cpu_surface(packed, uvs)
    out = zeros(3, size(uvs, 2))
    for i in 1:size(uvs, 2)
        out[:, i] = J.eval_packed(packed, J.SVector(uvs[1, i], uvs[2, i]))
    end
    out
end

function cpu_dualpass(packed, n, wedge, sts)
    vals = zeros(3, size(sts, 2))
    jacs = zeros(3, 2, size(sts, 2))
    f = w -> J.eval_packed(packed, J.square_coords_to_chart(n, wedge, w))
    for i in 1:size(sts, 2)
        v, jac = J.value_and_jacobian_columns(f, J.SVector(sts[1, i], sts[2, i]), Val(2))
        vals[:, i] = v
        jacs[:, :, i] = jac
    end
    (vals, jacs)
end

function main()
    npoints = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 1_000_000
    m = J.cubemesh()
    surf = J.Surface(m, J.fit_geometry(m))
    packed = J.packed_polys(surf)[1] # vertex 1's chart
    n = J.valence(m, 1)
    wedge = 1

    golden = 2pi * 0.6180339887498949
    uvs = stack([0.5 * sqrt(k / npoints) .* [cos(golden * k), sin(golden * k)] for k in 1:npoints])
    sts = stack([[mod(k * 0.6180339887498949, 1.0), mod(k * 0.3819660112501051, 1.0)] for k in 1:npoints])

    d_packed = CuArray(packed)
    d_uvs = CuArray(uvs)
    d_sts = CuArray(sts)
    d_out = CUDA.zeros(Float64, 3, npoints)
    d_vals = CUDA.zeros(Float64, 3, npoints)
    d_jacs = CUDA.zeros(Float64, 3, 2, npoints)

    threads = 256
    blocks = cld(npoints, threads)

    ka = @cuda launch=false kernel_surface(d_out, d_packed, d_uvs)
    kb = @cuda launch=false kernel_dualpass(d_vals, d_jacs, d_packed, n, wedge, d_sts)
    println("kernel A registers: ", CUDA.registers(ka), "   kernel B registers: ", CUDA.registers(kb))

    ka(d_out, d_packed, d_uvs; threads, blocks); CUDA.synchronize() # warm
    ta = CUDA.@elapsed ka(d_out, d_packed, d_uvs; threads, blocks)
    kb(d_vals, d_jacs, d_packed, n, wedge, d_sts; threads, blocks); CUDA.synchronize()
    tb = CUDA.@elapsed kb(d_vals, d_jacs, d_packed, n, wedge, d_sts; threads, blocks)

    cpu_surface(packed, uvs[:, 1:10]) # warm
    ca = @elapsed cpu_out = cpu_surface(packed, uvs)
    cpu_dualpass(packed, n, wedge, sts[:, 1:10])
    cb = @elapsed (cpu_vals, cpu_jacs) = cpu_dualpass(packed, n, wedge, sts)

    da = maximum(abs, Array(d_out) - cpu_out)
    dv = maximum(abs, Array(d_vals) - cpu_vals)
    dj = maximum(abs, Array(d_jacs) - cpu_jacs)

    rate(t) = round(npoints / t / 1e6, sigdigits = 3)
    println("A (Horner surface eval), $npoints points:")
    println("  gpu ", round(ta * 1000, digits = 2), " ms (", rate(ta), " Meval/s)   cpu 1 thread ",
            round(ca * 1000, digits = 2), " ms (", rate(ca), " Meval/s)   speedup ", round(ca / ta, sigdigits = 3))
    println("  max |gpu - cpu| = ", da)
    println("B (dual pass through wedge map), $npoints points:")
    println("  gpu ", round(tb * 1000, digits = 2), " ms (", rate(tb), " Meval/s)   cpu 1 thread ",
            round(cb * 1000, digits = 2), " ms (", rate(cb), " Meval/s)   speedup ", round(cb / tb, sigdigits = 3))
    println("  max |gpu - cpu| value = ", dv, "   jacobian = ", dj)
end

main()
