# GPU groundwork, the simpler road (kaarel 2026-07-14): before tabulating Γ,
# see how the device does running the production AD christoffel *directly* —
# second-order nested duals through the full blended surface eval. The type
# layer is storage-parametric, so the real Throat is Adapt-ed onto the device
# (ragged construction-only Mesh fields dropped to nothing, packed_polys
# padded into one 4D array behind an indexable view — padding sits at the
# leading end of the descending Horner loop, so eval stays bit-identical) and
# the kernel calls jupiter.christoffel unmodified on the same SituatedPhase
# the CPU tracer uses. Points cover all charts, wedges, and all three depth
# branches (pure outer d<0, blend, pure inner past cylinder_depth).
#
#   julia --project=gpu scripts/gpu_christoffel.jl [npoints]

using jupiter
using CUDA
using Adapt
const J = jupiter

# device-facing packed_polys: every chart padded to a common size in one 4D
# array; getindex returns the [component, u-power, v-power] slice eval_packed
# expects. Adapt rules for jupiter types live here until a gpu/ module earns
# its keep.
struct PackedTable{A}
    a::A
end
Base.getindex(p::PackedTable, v::Int) = view(p.a, :, :, :, v)

function PackedTable(polys::Vector{Array{Float64, 3}})
    du = maximum(size.(polys, 2))
    dv = maximum(size.(polys, 3))
    a = zeros(3, du, dv, length(polys))
    for (v, p) in enumerate(polys)
        a[:, 1:size(p, 2), 1:size(p, 3), v] = p
    end
    PackedTable(a)
end

Adapt.@adapt_structure PackedTable
Adapt.@adapt_structure J.Mesh
Adapt.@adapt_structure J.Surface
Adapt.@adapt_structure J.Throat

# host-side device view: flat tables as CuArrays, host-only fields dropped;
# @cuda's cudaconvert recurses through the adapt rules above at launch.
# T picks the on-device scalar type; the eltype-honest christoffel path then
# computes wholly in T because the phase and the packed table carry it.
function device_throat(th::J.Throat, T)
    m = J.mesh(J.geometry(th))
    dmesh = J.Mesh(nothing, nothing, nothing,
                   (CuArray(getfield(m, f)) for f in
                    (:he_tail, :he_head, :he_next, :he_prev, :he_twin,
                     :he_offset, :vertex_he, :vertex_valence))...)
    padded = PackedTable(J.packed_polys(J.geometry(th)))
    dsurf = J.Surface(dmesh, nothing, PackedTable(CuArray(T.(padded.a))))
    J.Throat(dsurf, J.params(th), th.placements)
end

function kernel_christoffel(out, throat, hes, poss)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i > size(poss, 2) && return nothing
    T = eltype(poss)
    chart = J.Chart(J.HalfThroat(throat, 1), J.HalfEdgeHandle(J.mesh(J.HalfThroat(throat, 1)), hes[i]))
    v = J.SituatedPhase(chart, J.SVector(poss[1, i], poss[2, i], poss[3, i]), J.SVector(zero(T), zero(T), one(T)))
    Γ = J.christoffel(nothing, v)
    for k in 1:27
        out[k, i] = Γ[k]
    end
    nothing
end

function cpu_christoffels(th, hes, poss)
    m = J.mesh(J.geometry(th))
    out = zeros(27, size(poss, 2))
    for i in 1:size(poss, 2)
        chart = J.Chart(J.HalfThroat(th, 1), J.HalfEdgeHandle(m, hes[i]))
        v = J.SituatedPhase(chart, J.SVector(poss[1, i], poss[2, i], poss[3, i]), J.SVector(0.0, 0.0, 1.0))
        Γ = J.christoffel(nothing, v)
        out[:, i] = Γ[:]
    end
    out
end

frac(x) = x - floor(x)

function main()
    npoints = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 200_000
    m = J.cubemesh()
    surf = J.Surface(m, J.fit_geometry(m))
    th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75)) # the physics_diff cube

    φ = 0.6180339887498949
    hes = zeros(Int, npoints)
    poss = zeros(3, npoints)
    for k in 1:npoints
        vtx = mod1(k, length(m.vertex_he))
        n = J.valence(m, vtx)
        wedge = mod(k, n)
        st = J.SVector(0.02 + 0.96 * frac(k * φ), 0.02 + 0.96 * frac(k * φ^2))
        uv = J.square_coords_to_chart(n, wedge, st)
        t = frac(k * 0.7548776662466927)
        # 25% pure outer (d<0), the rest sweeps blend + pure inner (cylinder_depth 0.5)
        d = t < 0.25 ? (t - 0.25) * 0.8 : (t - 0.25) / 0.75 * 0.7
        hes[k] = m.vertex_he[vtx]
        poss[:, k] = [uv[1], uv[2], d]
    end

    d_hes = CuArray(hes)

    nsub = min(npoints, 20_000)
    cpu_christoffels(th, hes[1:10], poss[:, 1:10]) # warm
    tc = @elapsed cpu_out = cpu_christoffels(th, hes[1:nsub], poss[:, 1:nsub])
    println("cpu 1 thread: ", round(tc / nsub * 1e6, digits = 2), " µs/eval (truth on $nsub points)")

    for T in (Float64, Float32)
        dth = device_throat(th, T)
        d_poss = CuArray(T.(poss))
        d_out = CUDA.zeros(T, 27, npoints)

        kc = @cuda launch=false kernel_christoffel(d_out, dth, d_hes, d_poss)
        mem = CUDA.memory(kc)
        println(T, ":")
        println("  registers: ", CUDA.registers(kc), "   local (spill) bytes: ", mem.local)

        kc(d_out, dth, d_hes, d_poss; threads = 256, blocks = cld(npoints, 256)); CUDA.synchronize() # warm
        times = map((64, 128, 256)) do threads # 512 exceeds the 64K regs/block budget at 255 regs/thread
            CUDA.@elapsed kc(d_out, dth, d_hes, d_poss; threads, blocks = cld(npoints, threads))
        end
        tg = minimum(times)
        gpu_out = Float64.(Array(d_out)) # before the divergence launches overwrite it

        # branch-divergence datum: same lateral points, all three depth branches
        # uniform per launch (warps never split on the metric's depth branch)
        for (name, dfix) in (("pure outer (d=-0.1)", -0.1), ("blend (d=0.25)", 0.25), ("pure inner (d=0.6)", 0.6))
            d_pu = CuArray(T.(vcat(poss[1:2, :], fill(dfix, 1, npoints))))
            t = CUDA.@elapsed kc(d_out, dth, d_hes, d_pu; threads = 256, blocks = cld(npoints, 256))
            println("  ", name, ": ", round(t / npoints * 1e6, digits = 2), " µs/eval uniform-branch")
        end

        println("  mixed branches: ", round(tg * 1000, digits = 2), " ms  (",
                round(tg / npoints * 1e9, digits = 1), " ns/eval, ",
                round(npoints / tg / 1e6, sigdigits = 3), " Meval/s)")
        println("  speedup vs 1 cpu thread: ", round(tc / nsub / (tg / npoints), sigdigits = 3), "x")

        # accuracy vs the Float64 cpu truth: absolute, and per-point sup-relative
        adiff = vec(maximum(abs, gpu_out[:, 1:nsub] - cpu_out, dims = 1))
        rel = sort(adiff ./ max.(vec(maximum(abs, cpu_out, dims = 1)), 1e-12))
        println("  vs cpu truth: max abs ", round(maximum(adiff), sigdigits = 3),
                "   sup-relative median ", round(rel[(nsub + 1) ÷ 2], sigdigits = 3),
                "  p99 ", round(rel[ceil(Int, 0.99 * nsub)], sigdigits = 3),
                "  max ", round(rel[end], sigdigits = 3))
    end
end

main()
