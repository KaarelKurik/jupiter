# device side of the wavefront tracer (GPU step 5). Promoted from
# scripts/gpu_christoffel.jl when the tracer kernel landed, per plan: the
# PackedTable/adapt rules live here now, plus the kernel itself — sweep_ray
# unmodified per thread, so certification tiers hold (device runs the same
# production code the CPU wavefront runs, which is bit-identical to the
# recursive tracer, which is certified against reference). Loaded by
# `include` from scripts running under --project=gpu; src/ stays CUDA-free.

module jupitergpu

using jupiter
using CUDA
using Adapt
const J = jupiter

# device-facing packed_polys: every chart padded to a common size in one 4D
# array; getindex returns the [component, u-power, v-power] slice eval_packed
# expects. Padding sits at the leading end of the descending Horner loop, so
# eval stays bit-identical.
struct PackedTable{A}
    a::A
end
Base.getindex(p::PackedTable, v::Int) = view(p.a, :, :, :, v)

function PackedTable(polys::Vector{<:Array{<:AbstractFloat, 3}}) # F32-rounded tables pass through exactly (F32→F64→T)
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
# T picks the on-device scalar type; the eltype-honest trace path then
# computes wholly in T because the records and the packed table carry it.
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

# the tracer kernel: one thread marches one ray — sweep_ray's settled
# steps, chart hops, half transitions, and to_ambient exit, unmodified
function kernel_sweep(recs, env, throat, budget, sweep)
    i = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    i > length(recs) && return nothing
    recs[i] = J.sweep_ray(env, throat, budget, recs[i], sweep)
    nothing
end

"""
sweep stage for run_wavefront! that launches the active records through
kernel_sweep on a device throat (build it with the same T as the pool, or
mixed-precision promotion silently puts F64 math on the card). Gather →
launch → scatter; kernel wall time accumulates in `kernel_seconds`.

`handoff`: rounds with at most this many active records run on host threads
instead, to completion (sweep = the full step budget) — the straggler tail
is latency-bound on the device (a near-constant ~33 ms per round however few
rays remain), while a few hundred budget-capped rays are milliseconds on
host threads. The handed-off rays go through the same sweep_ray the kernel
runs, on the host tables of the same values, so they land in the standing
CPU-vs-device intrinsic-ulps band; host time is not counted in
`kernel_seconds`.
"""
struct DeviceSweep{DT}
    dthroat::DT
    threads::Int
    handoff::Int
    kernel_seconds::Base.RefValue{Float64}
end

DeviceSweep(dthroat; threads = 128, handoff = 1024) = DeviceSweep(dthroat, threads, handoff, Ref(0.0))

function (ds::DeviceSweep)(pool::Vector{J.WavefrontRay{T}}, active, env, throat, budget, sweep) where {T}
    if length(active) <= ds.handoff
        J.threaded_sweep!(pool, active, env, throat, budget, budget.max_steps)
        return nothing
    end
    recs = pool[active]
    drecs = CuArray(recs)
    t = CUDA.@elapsed @cuda threads = ds.threads blocks = cld(length(recs), ds.threads) kernel_sweep(
        drecs, env, ds.dthroat, budget, sweep)
    ds.kernel_seconds[] += t
    copyto!(recs, drecs)
    pool[active] = recs
    nothing
end

"""
wavefront_raymap with the sweep stage on the device: builds the T-typed
device throat, runs the host driver around kernel launches, returns
(raymap, accumulated kernel seconds). bin defaults on — depth grouping is
the whole point on a warp machine.
"""
function device_raymap(env, scene, budget, camera, cam_side, width, height;
                       T = Float32, sweep = 64, bin = true, threads = 128, handoff = 1024)
    ds = DeviceSweep(device_throat(scene.throat, T); threads, handoff)
    rm = J.wavefront_raymap(env, scene, budget, camera, cam_side, width, height;
                            T, sweep, bin, sweep_stage! = ds)
    (rm, ds.kernel_seconds[])
end

end
