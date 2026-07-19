# dense flythrough with measured optical flow, two modes:
#  - adaptive (default): frame spacing chosen so peak apparent motion between
#    consecutive frames stays near a pixel-scale target, measured on cheap
#    sparse probe raymaps — accept/reject in the spirit of the DP5 stepper.
#    Perceptually uniform, time-warped. Emits out/frames/flyd_NNNN.ppm.
#  - uniform=Δτ: physically honest constant speed at fixed spacing; probes
#    only log the flow so the choice can be audited. Emits out/frames/flyu_NNNN.ppm.
#    Pick Δτ from an adaptive run: Δτ = θ_target / max(flow rate) — the whole
#    flight at the rate its most demanding moment requires.
# scene=trefoil flies the knotted throat (cube default); frames get the scene
# name prefixed (out/frames/trefoil_flyd_NNNN.ppm).
# skies: sky1=/sky2= point at equirect PPMs (default the res/textures test
# graticules). Real skies live as tracked JPGs in res/textures/: skybox[12].jpg
# (CC0 Poly Haven nature panoramas) and space[12].jpg (kaarel's manual_wormhole
# bg0/bg1 cubemaps, converted once via scripts/cube2equirect.jl); regenerate
# the render-ready PPMs per checkout with magick — see below
# unresolved pixels: each frame the raymap's unresolved pixels (limit-cycle
# boundary filaments) are retraced under a doubled budget until their fraction
# drops below unres= (default 0.005) or esc= doublings (default 3) are spent —
# only the unresolved pixels are retraced, so escalation costs what it chases.
# Assemble with e.g.
#   ffmpeg -framerate 24 -i out/frames/flyd_%04d.ppm -pix_fmt yuv420p out/flythrough.mp4
#   magick res/textures/space1.jpg out/space1.ppm   # once per checkout
#   magick res/textures/space2.jpg out/space2.ppm
#   julia --project --threads=auto scripts/flyvideo.jl [WxH] [uniform=Δτ] [scene=trefoil]
#       [sky1=out/space1.ppm sky2=out/space2.ppm] [unres=0.005] [esc=3]
# gpu=1 (run under --project=gpu): frame raymaps trace on the device in F32
# against F32-rounded tables (entry solves stay F64 against the F32-coefficient
# surface — the 14c seam). Flight dynamics, flow probes, and escalation retraces
# stay CPU F64, so the flight path is bit-identical to a cpu run and escalated
# pixels land at F64:
#   julia --project=gpu --threads=auto scripts/flyvideo.jl [args...] gpu=1
using jupiter
using LinearAlgebra
const J = jupiter

const root = joinpath(@__DIR__, "..")
argval(key, default) = (i = findfirst(a -> startswith(a, key * "="), ARGS);
                        i === nothing ? default : parse(Float64, split(ARGS[i], "=")[2]))
argstr(key, default) = (i = findfirst(a -> startswith(a, key * "="), ARGS);
                        i === nothing ? default : String(split(ARGS[i], "=", limit = 2)[2]))
const uniform = argval("uniform", 0.0)
const tau0 = argval("tau0", 0.0)     # resume: coast (no rendering) to this τ first
const frame0 = Int(argval("frame0", 0.0)) # resume: continue numbering from here
const dims = (d = findfirst(a -> occursin("x", a) && !occursin("=", a), ARGS);
              d === nothing ? "192x144" : ARGS[d])
const scene_name = argstr("scene", "cube")
const sky1_path = argstr("sky1", joinpath("res", "textures", "sky1.ppm"))
const sky2_path = argstr("sky2", joinpath("res", "textures", "sky2.ppm"))
const unres_max = argval("unres", 0.005) # per-frame unresolved-pixel fraction to escalate against
const esc_max = Int(argval("esc", 3.0))  # budget doublings before we accept the leftovers
const w, h = parse.(Int, split(dims, "x"))
const thf = tan(pi / 3.2 / 2)
const θ_pixel = (pi / 3.2) / w        # angular size of a render pixel
# pacing knobs, overridable per flight: a strong lens filling the view (the
# trefoil) holds median flow high through the whole approach, so the cube's
# floor would spend the frame budget before the crossing
const θ_target = argval("target", 3.0) * θ_pixel # aim for ~target px of median flow per frame
const flip_max = 0.15                 # probe fraction allowed to change exit side
const Δτ_floor = argval("floor", 0.004) # hard floor: below this, motion is shimmer, not flow
const frame_cap = 1200                # absolute runaway guard
const probe_w, probe_h = 16, 12
const coast_h = 0.02
const use_gpu = argval("gpu", 0.0) > 0
use_gpu && include(joinpath(root, "gpu", "jupitergpu.jl")) # brings CUDA; needs --project=gpu

# F32-rounded packed tables on the same mesh/params/placements: what the device
# traces against, and what the F64 entry solves see (the 14c precision seam)
f32_throat(th) = J.Throat(J.Surface(J.mesh(J.geometry(th)), J.chart_polys(J.geometry(th)),
                                    [Float32.(p) for p in J.packed_polys(J.geometry(th))]),
                          J.params(th), th.placements)

# angular feature motion between probes: median, q90, max, and the side-flip
# fraction. The controller steers on the MEDIAN: near the mouth the
# chaotically-sensitive pixel population grows past any fixed upper quantile
# (run 1: q90 pinned at ~2.5px independent of Δτ once >10% of the view was
# near-critical, and the controller ground to Δτ=0.0008 chasing it). Chaotic
# shimmer is image-space aliasing — no temporal density fixes it; only the
# coherent bulk motion is a temporal sampling signal.
function flow(rm1, rm2)
    angs = Float64[]
    flips = 0
    for j in 1:probe_h, i in 1:probe_w
        s1, s2 = rm1.side[i, j], rm2.side[i, j]
        if s1 != s2
            flips += 1
        elseif s1 != 0
            d1 = J.generic_normalize(rm1.vel[:, i, j])
            d2 = J.generic_normalize(rm2.vel[:, i, j])
            push!(angs, acos(clamp(d1' * d2, -1.0, 1.0)))
        end
    end
    sort!(angs)
    q(p) = isempty(angs) ? 0.0 : angs[clamp(ceil(Int, p * length(angs)), 1, length(angs))]
    (q(0.5), q(0.9), q(1.0), flips / (probe_w * probe_h))
end

# retrace just the raymap's unresolved pixels under a bigger budget, in place.
# Mirrors keyframe_raymap's camera dispatch and render_raymap's per-pixel body —
# the full-frame drivers stay untouched, escalation is this script's business.
function refine!(raymap, scene, budget, cam, tan_half_fov, w, h)
    idx = findall(==(0), raymap.side)
    if cam.state isa J.AmbientRay
        camera = J.Camera(cam.state.pos, cam.frame, tan_half_fov)
        cam_space = J.half_throat(cam.state)
        Threads.@threads for n in eachindex(idx)
            i, j = Tuple(idx[n])
            out = J.trace_ray(nothing, scene, budget, J.pixel_ray(camera, cam_space, i, j, w, h))
            out === nothing && continue
            raymap.side[i, j] = J.side(J.half_throat(out))
            raymap.pos[:, i, j] = out.pos
            raymap.vel[:, i, j] = out.vel
        end
    else
        camera = J.SituatedCamera(cam.state.chart, cam.state.pos, cam.frame, tan_half_fov)
        Threads.@threads for n in eachindex(idx)
            i, j = Tuple(idx[n])
            x, y = J.pixel_ndc(i, j, w, h)
            out = J.trace_ray(nothing, scene, budget, J.camera_ray(nothing, camera, x, y))
            out === nothing && continue
            raymap.side[i, j] = J.side(J.half_throat(out))
            raymap.pos[:, i, j] = out.pos
            raymap.vel[:, i, j] = out.vel
        end
    end
end

# ambient-2 placement for the trefoil flight: rotate the second mouth's
# embedding so the flight's outro flies toward the space2 sky's sun. Scene
# authoring, not physics: the sky asset stays as-shot, the chart-side flight
# is untouched, and mouth geometry + exit rays rotate rigidly together.
# Minimal rotation (Rodrigues) taking the exit camera fwd — measured at τ=14
# with the identity placement — to the sun's luminance centroid in
# res/textures/space2.jpg. Harmlessly reorients any other sky2 the flight is
# run under.
function exit_placement()
    f = J.generic_normalize(J.SVector(-0.7661, 0.0009, -0.6428)) # exit fwd, identity placement
    s = J.generic_normalize(J.SVector(-0.1144, -0.1365, 0.984))  # space2 sun direction
    a = cross(f, s)
    K = J.SMatrix{3, 3, Float64}(0, a[3], -a[2], -a[3], 0, a[1], a[2], -a[1], 0)
    J.Placement(J.SMatrix{3, 3, Float64}(I) + K + K * K / (1 + f' * s), zeros(3))
end

function build_scene() # (throat, mouth tessellation, campos, fwd, τ_total)
    if scene_name == "trefoil"
        m = J.load_obj(joinpath(root, "res", "models", "trefoil.obj"))
        th = J.Throat(J.Surface(m, J.fit_geometry(m)), J.ThroatParams(1.0, 1.0, 0.2, 0.3),
                      (J.identity_placement(), exit_placement()))
        campos = J.SVector(0.5, -7.5, 3.5) # the physics_diff approach vector: known to reach the mouth
        (th, 2, campos, J.generic_normalize(J.SVector(0.262, -1.0, 0.437) - campos), 14.0)
    else
        m = J.cubemesh()
        th = J.Throat(J.Surface(m, J.fit_geometry(m)), J.ThroatParams(1.0, 1.0, 0.5, 0.75))
        campos = J.SVector(0.9, -2.8, 0.9)
        (th, 10, campos, J.generic_normalize(-campos), 6.2)
    end
end

function main()
    th, tess, campos, fwd, τ_total = build_scene()
    mouths = (J.TessellatedMouth(nothing, J.HalfThroat(th, 1), tess),
              J.TessellatedMouth(nothing, J.HalfThroat(th, 2), tess))
    sky = J.TexturedSky(J.load_ppm(joinpath(root, sky1_path)),
                        J.load_ppm(joinpath(root, sky2_path)))
    scene = J.Scene(th, mouths, sky)
    budget = J.RayBudget(0.05, 400, 4)

    right = J.generic_normalize(cross(fwd, [0.0, 0.0, 1.0]))
    cam = J.FlyingCamera(J.AmbientRay(J.HalfThroat(th, 1), campos, fwd),
                         J.SMatrix{3, 3, Float64}([right cross(right, fwd) fwd]))

    # device render path (gpu=1): the full-res frame raymaps run the wavefront
    # driver with the sweep stage on the card, one DeviceSweep reused across
    # frames (device throat uploaded once). Everything else — coast, probes,
    # refine! escalation — stays on the CPU F64 scene.
    render_scene, sweep_stage = scene, nothing
    if use_gpu
        th32 = f32_throat(th)
        render_scene = J.Scene(th32, (J.TessellatedMouth(nothing, J.HalfThroat(th32, 1), tess),
                                      J.TessellatedMouth(nothing, J.HalfThroat(th32, 2), tess)), sky)
        sweep_stage = jupitergpu.DeviceSweep(jupitergpu.device_throat(th32, Float32))
    end
    function render_keyframe(c) # keyframe_raymap's camera dispatch, device edition
        use_gpu || return J.keyframe_raymap(nothing, scene, budget, c, thf, w, h)
        if c.state isa J.AmbientRay
            J.wavefront_raymap(nothing, render_scene, budget,
                               J.Camera(c.state.pos, c.frame, thf),
                               J.side(J.half_throat(c.state)), w, h;
                               T = Float32, sweep = 64, bin = true, sweep_stage! = sweep_stage)
        else # emission stays on the F64 chart (host); the record traces in F32 like the CPU F32 wavefront
            J.wavefront_raymap(nothing, render_scene, budget,
                               J.SituatedCamera(c.state.chart, c.state.pos, c.frame, thf),
                               w, h; T = Float32, sweep = 64, bin = true, sweep_stage! = sweep_stage)
        end
    end

    probe(c) = J.keyframe_raymap(nothing, scene, budget, c, thf, probe_w, probe_h)
    where(c) = c.state isa J.AmbientRay ?
        "ambient$(J.side(J.half_throat(c.state)))" : "chart d=$(round(c.state.pos[3], digits=2))"

    mkpath(joinpath(root, "out", "frames"))
    prefix = (scene_name == "cube" ? "" : scene_name * "_") * (uniform > 0 ? "flyu_" : "flyd_")
    if tau0 > 0 # resume: replay the flight to τ0 without rendering
        cam = J.coast(nothing, scene, cam, tau0, coast_h)
        println("resumed at τ=", tau0, " (", where(cam), "), frame numbering from ", frame0 + 1)
    end
    τ, Δτ, frame, t0 = tau0, uniform > 0 ? uniform : 0.05, frame0, time()
    while τ < τ_total && frame < frame_cap
        p1 = probe(cam)
        medflow, q90flow, maxflow, flipfrac = 0.0, 0.0, 0.0, 0.0
        if uniform > 0 # honest constant speed; the probe only audits the flow
            cam = J.coast(nothing, scene, cam, Δτ, coast_h)
            medflow, q90flow, maxflow, flipfrac = flow(p1, probe(cam))
        else
            while true # shrink Δτ until the probes say the step is smooth enough
                cam2 = J.coast(nothing, scene, cam, Δτ, coast_h)
                medflow, q90flow, maxflow, flipfrac = flow(p1, probe(cam2))
                if (medflow <= 2 * θ_target && flipfrac <= flip_max) || Δτ <= Δτ_floor
                    cam = cam2
                    break
                end
                Δτ = max(Δτ / 2, Δτ_floor)
            end
        end
        τ += Δτ
        frame += 1
        raymap = render_keyframe(cam)
        unres0, esc_budget, escs = count(==(0), raymap.side), budget, 0
        while count(==(0), raymap.side) > unres_max * w * h && escs < esc_max
            esc_budget = J.RayBudget(esc_budget.step_size, 2 * esc_budget.max_steps,
                                     2 * esc_budget.max_passages)
            refine!(raymap, scene, esc_budget, cam, thf, w, h)
            escs += 1
        end
        J.save_ppm(joinpath(root, "out", "frames", prefix * lpad(frame, 4, '0') * ".ppm"),
                   J.shade(raymap, sky))
        println("frame ", frame, "  τ=", round(τ, digits=3), "  Δτ=", round(Δτ, digits=4),
                "  medflow=", round(medflow / θ_pixel, digits=1), "px  q90flow=",
                round(q90flow / θ_pixel, digits=1), "px  maxflow=",
                round(maxflow / θ_pixel, digits=1), "px  flips=",
                round(flipfrac * 100, digits=1), "%  ", where(cam),
                "  unresolved=", unres0, "->", count(==(0), raymap.side),
                "  esc=", escs)
        flush(stdout)
        # steer the next leg on the coherent bulk motion, floored and capped
        uniform > 0 || (Δτ = clamp(Δτ * clamp(0.85 * θ_target / (medflow + 1e-9), 0.6, 1.8), Δτ_floor, 0.4))
    end
    println(frame, " frames over τ=", τ_total, " in ", round((time() - t0) / 60, digits=1), " min")
    use_gpu && println("device kernel total ", round(sweep_stage.kernel_seconds[], digits=1), " s")
end

main()
