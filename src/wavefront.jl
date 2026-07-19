# wavefront tracer: render.jl's per-ray passage recursion restructured as
# staged sweeps over a flat pool of isbits ray records (GPU step 5). The
# sum-type returns (Union{Nothing, AmbientRay}, Union{SituatedPhase,
# AmbientRay}) become a stage tag on one fixed-size record, so the in-throat
# stage is kernel-shaped: record in, record out, no boxes. Mouth entry stays
# a host-side Float64 solve — the precision seam precision_diff certified:
# the phase converts to the trace type T at entry, exits convert at
# to_ambient, and the initial camera ray reaches its first entry solve in
# full Float64 (it never lives in a T record). F64 wavefront raymaps are
# bit-identical to render_raymap's; the sweep granularity, the depth-branch
# binning of the active list, and thread scheduling never change any
# per-ray operation sequence.

const RAY_THROAT = Int8(0)     # in-throat, tracing
const RAY_ENTRY = Int8(1)      # ambient, awaiting an entry solve
const RAY_ESCAPED = Int8(2)    # ambient, missed every mouth: final
const RAY_UNRESOLVED = Int8(3) # budget exhausted: final

struct WavefrontRay{T}
    stage::Int8
    side::Int8
    he::Int32       # chart handle id while in-throat; 0 in ambient stages
    passages::Int32 # entry solves consumed
    steps::Int32    # integrator attempts consumed in the current passage
    h::T            # dopri's adaptive step, persisted across sweeps
    pos::SVector{3, T}
    vel::SVector{3, T}
end

retag(r::WavefrontRay{T}, stage) where {T} =
    WavefrontRay{T}(stage, r.side, r.he, r.passages, r.steps, r.h, r.pos, r.vel)

function throat_state(stage, v::SituatedPhase, passages, steps, h::T) where {T}
    WavefrontRay{T}(stage, Int8(side(half_throat(v.chart))),
                    Int32(half_edge_handle(v.chart).id), Int32(passages),
                    Int32(steps), h, v.pos, v.vel)
end

throat_record(::Type{T}, v::SituatedPhase, passages, budget) where {T} =
    throat_state(RAY_THROAT, SituatedPhase(v.chart, T.(v.pos), T.(v.vel)),
                 passages, Int32(0), T(budget.step_size))

function ambient_state(ray::AmbientRay, passages, h::T) where {T}
    WavefrontRay{T}(RAY_ENTRY, Int8(side(half_throat(ray))), Int32(0),
                    Int32(passages), Int32(0), h, ray.pos, ray.vel)
end

"""
the kernel body: advance one in-throat ray by up to `sweep` integrator
attempts in its own precision, replicating trace_geodesic's loops op-for-op
(settle, exit check, dopri accept/step-update order), so resuming across
sweeps is invisible to the ray. Flips to RAY_ENTRY on a mouth exit (pos/vel
become ambient) or RAY_UNRESOLVED when the passage's step budget runs out.
"""
function sweep_ray(env, throat, budget, r::WavefrontRay{T}, sweep) where {T}
    chart = Chart(HalfThroat(throat, Int(r.side)),
                  HalfEdgeHandle(mesh(geometry(throat)), Int(r.he)))
    v = SituatedPhase(chart, r.pos, r.vel)
    steps = r.steps
    h = r.h
    # a non-finite phase (F32 blowup, filament class) can never resolve: retire
    # it as RAY_UNRESOLVED instead of burning budget on NaN arithmetic (or, on
    # device, throwing in-kernel). Control flow only — never fires on finite
    # rays, so F64 bit-identity with trace_geodesic is untouched.
    finite(p) = isfinite(sum(p.pos) + sum(p.vel))
    if budget.tolerance > 0
        tol = T(budget.tolerance)
        for _ in 1:sweep
            finite(v) || return throat_state(RAY_UNRESOLVED, v, r.passages, steps, h)
            h = min(h, T(0.25) / (maximum(abs, v.vel) + T(1e-12)))
            u = vcat(SVector{3}(v.pos), SVector{3}(v.vel))
            u5, err = dopri_step(env, v.chart, u, h)
            scale = tol * (1 + maximum(abs, u))
            if err <= scale
                v = settle_phase(env, SituatedPhase(v.chart, SVector(u5[1], u5[2], u5[3]),
                                                    SVector(u5[4], u5[5], u5[6])))
                exits_mouth(v) && return ambient_state(to_ambient(env, v), r.passages, h)
            end
            h = h * clamp(T(0.9) * (scale / (err + floatmin(T)))^T(1/5), T(0.2), T(5.0))
            steps += Int32(1)
            steps >= budget.max_steps && return throat_state(RAY_UNRESOLVED, v, r.passages, steps, h)
        end
    else
        for _ in 1:sweep
            finite(v) || return throat_state(RAY_UNRESOLVED, v, r.passages, steps, h)
            v = settle_phase(env, geodesic_step(env, v, h))
            steps += Int32(1)
            exits_mouth(v) && return ambient_state(to_ambient(env, v), r.passages, h)
            steps >= budget.max_steps && return throat_state(RAY_UNRESOLVED, v, r.passages, steps, h)
        end
    end
    throat_state(RAY_THROAT, v, r.passages, steps, h)
end

"""
host stage: one entry attempt for a RAY_ENTRY record. The passage budget is
charged here (a ray that exits its last allowed passage never gets another
entry solve — trace_ray's loop shape); the solve runs in Float64 on the
record's values (exact: T embeds), and a hit converts the phase to T.
"""
function entry_stage(env, scene, budget, r::WavefrontRay{T}) where {T}
    r.passages >= budget.max_passages && return retag(r, RAY_UNRESOLVED)
    v = enter_mouth(env, scene.mouths[Int(r.side)], Float64.(r.pos), Float64.(r.vel))
    v === nothing && return retag(r, RAY_ESCAPED)
    throat_record(T, v, r.passages + Int32(1), budget)
end

pool_pixel(k, width) = ((k - 1) % width + 1, (k - 1) ÷ width + 1)

function finalize_exit!(raymap, i, j, s, pos, vel)
    raymap.side[i, j] = s
    raymap.pos[:, i, j] = pos
    raymap.vel[:, i, j] = vel
end

function depth_branch(throat, r::WavefrontRay) # the metric's three-way branch, the divergence axis
    d = r.pos[3]
    d < 0 ? Int8(0) : (d < params(throat).cylinder_depth ? Int8(1) : Int8(2))
end

"""
default sweep stage: sweep_ray over the active records on host threads. Any
callable with this signature can stand in (the gpu module's DeviceSweep
launches the same records through the kernel instead).
"""
function threaded_sweep!(pool, active, env, throat, budget, sweep)
    Threads.@threads for n in 1:length(active)
        k = active[n]
        pool[k] = sweep_ray(env, throat, budget, pool[k], sweep)
    end
end

"""
drive the pool to completion: sweep the active in-throat rays, route mouth
exits through the entry stage, drop finished rays from the lists. `bin`
groups the active list by depth branch before each sweep — the device
divergence lever; on the CPU it only permutes an order no result depends on.
"""
function run_wavefront!(raymap, pool, active, env, scene, budget, width, sweep, bin,
                        sweep_stage! = threaded_sweep!)
    pending = Int32[]
    while !isempty(active)
        bin && sort!(active, by = k -> depth_branch(scene.throat, pool[k]))
        sweep_stage!(pool, active, env, scene.throat, budget, sweep)
        empty!(pending)
        keep = 0
        for k in active
            if pool[k].stage == RAY_THROAT
                keep += 1
                active[keep] = k
            elseif pool[k].stage == RAY_ENTRY
                push!(pending, k)
            end
        end
        resize!(active, keep)
        isempty(pending) && continue
        Threads.@threads for n in 1:length(pending)
            k = pending[n]
            pool[k] = entry_stage(env, scene, budget, pool[k])
        end
        for k in pending
            r = pool[k]
            if r.stage == RAY_THROAT
                push!(active, k)
            elseif r.stage == RAY_ESCAPED
                i, j = pool_pixel(k, width)
                finalize_exit!(raymap, i, j, Int(r.side), Float64.(r.pos), Float64.(r.vel))
            end
        end
    end
end

"""
render_raymap through the wavefront driver: same contract, bit-identical at
T = Float64. T picks the in-throat trace precision (the entry solves stay
F64); `sweep` is the per-launch attempt granularity, `bin` the depth-branch
grouping — both scheduling-only.
"""
function wavefront_raymap(env, scene, budget, camera::Camera, cam_side::Int, width::Int, height::Int;
                          T::Type = Float64, sweep::Int = 64, bin::Bool = false,
                          sweep_stage! = threaded_sweep!)
    raymap = RayMap(zeros(Int, width, height), zeros(3, width, height), zeros(3, width, height))
    pool = Vector{WavefrontRay{T}}(undef, width * height)
    cam_space = HalfThroat(scene.throat, cam_side)
    Threads.@threads for k in 1:length(pool) # first entry straight from the F64 camera ray
        i, j = pool_pixel(k, width)
        ray = pixel_ray(camera, cam_space, i, j, width, height)
        v = budget.max_passages < 1 ? nothing :
            enter_mouth(env, scene.mouths[cam_side], ray.pos, ray.vel)
        if v === nothing
            budget.max_passages >= 1 &&
                finalize_exit!(raymap, i, j, cam_side, ray.pos, ray.vel)
            pool[k] = WavefrontRay{T}(RAY_UNRESOLVED, Int8(cam_side), Int32(0), Int32(0),
                                      Int32(0), zero(T), zero(SVector{3, T}), zero(SVector{3, T}))
        else
            pool[k] = throat_record(T, v, Int32(1), budget)
        end
    end
    active = Int32[k for k in 1:length(pool) if pool[k].stage == RAY_THROAT]
    run_wavefront!(raymap, pool, active, env, scene, budget, width, sweep, bin, sweep_stage!)
    raymap
end

"""
wavefront pass for a camera inside the throat: the initial segment is a
passage-0 in-throat record (emitted phase converted to T), so exits feed the
usual entry rounds with the full passage budget — trace_ray's situated shape.
"""
function wavefront_raymap(env, scene, budget, camera::SituatedCamera, width::Int, height::Int;
                          T::Type = Float64, sweep::Int = 64, bin::Bool = false,
                          sweep_stage! = threaded_sweep!)
    raymap = RayMap(zeros(Int, width, height), zeros(3, width, height), zeros(3, width, height))
    pool = Vector{WavefrontRay{T}}(undef, width * height)
    Threads.@threads for k in 1:length(pool)
        i, j = pool_pixel(k, width)
        x, y = pixel_ndc(i, j, width, height)
        pool[k] = throat_record(T, camera_ray(env, camera, x, y), Int32(0), budget)
    end
    active = Int32[Int32(k) for k in 1:length(pool)]
    run_wavefront!(raymap, pool, active, env, scene, budget, width, sweep, bin, sweep_stage!)
    raymap
end
