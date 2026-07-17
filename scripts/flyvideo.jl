# dense flythrough with measured optical flow, two modes:
#  - adaptive (default): frame spacing chosen so peak apparent motion between
#    consecutive frames stays near a pixel-scale target, measured on cheap
#    sparse probe raymaps — accept/reject in the spirit of the DP5 stepper.
#    Perceptually uniform, time-warped. Emits out/flyd_NNNN.ppm.
#  - uniform=Δτ: physically honest constant speed at fixed spacing; probes
#    only log the flow so the choice can be audited. Emits out/flyu_NNNN.ppm.
#    Pick Δτ from an adaptive run: Δτ = θ_target / max(flow rate) — the whole
#    flight at the rate its most demanding moment requires.
# scene=trefoil flies the knotted throat (cube default); frames get the scene
# name prefixed (out/trefoil_flyd_NNNN.ppm).
# Assemble with e.g.
#   ffmpeg -framerate 24 -i out/flyd_%04d.ppm -pix_fmt yuv420p out/flythrough.mp4
#   julia --project --threads=auto scripts/flyvideo.jl [WxH] [uniform=Δτ] [scene=trefoil]
using jupiter
using LinearAlgebra
const J = jupiter

const root = joinpath(@__DIR__, "..")
argval(key, default) = (i = findfirst(a -> startswith(a, key * "="), ARGS);
                        i === nothing ? default : parse(Float64, split(ARGS[i], "=")[2]))
const uniform = argval("uniform", 0.0)
const tau0 = argval("tau0", 0.0)     # resume: coast (no rendering) to this τ first
const frame0 = Int(argval("frame0", 0.0)) # resume: continue numbering from here
const dims = (d = findfirst(a -> occursin("x", a) && !occursin("=", a), ARGS);
              d === nothing ? "192x144" : ARGS[d])
const scene_name = (i = findfirst(a -> startswith(a, "scene="), ARGS);
                    i === nothing ? "cube" : String(split(ARGS[i], "=")[2]))
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

function build_scene() # (throat, mouth tessellation, campos, fwd, τ_total)
    if scene_name == "trefoil"
        m = J.load_obj(joinpath(root, "res", "models", "trefoil.obj"))
        th = J.Throat(J.Surface(m, J.fit_geometry(m)), J.ThroatParams(1.0, 1.0, 0.2, 0.3))
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
    sky = J.TexturedSky(J.load_ppm(joinpath(root, "res", "textures", "sky1.ppm")),
                        J.load_ppm(joinpath(root, "res", "textures", "sky2.ppm")))
    scene = J.Scene(th, mouths, sky)
    budget = J.RayBudget(0.05, 400, 4)

    right = J.generic_normalize(cross(fwd, [0.0, 0.0, 1.0]))
    cam = J.FlyingCamera(J.AmbientRay(J.HalfThroat(th, 1), campos, fwd),
                         J.SMatrix{3, 3, Float64}([right cross(right, fwd) fwd]))

    probe(c) = J.keyframe_raymap(nothing, scene, budget, c, thf, probe_w, probe_h)
    where(c) = c.state isa J.AmbientRay ?
        "ambient$(J.side(J.half_throat(c.state)))" : "chart d=$(round(c.state.pos[3], digits=2))"

    mkpath(joinpath(root, "out"))
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
        raymap = J.keyframe_raymap(nothing, scene, budget, cam, thf, w, h)
        J.save_ppm(joinpath(root, "out", prefix * lpad(frame, 4, '0') * ".ppm"),
                   J.shade(raymap, sky))
        println("frame ", frame, "  τ=", round(τ, digits=3), "  Δτ=", round(Δτ, digits=4),
                "  medflow=", round(medflow / θ_pixel, digits=1), "px  q90flow=",
                round(q90flow / θ_pixel, digits=1), "px  maxflow=",
                round(maxflow / θ_pixel, digits=1), "px  flips=",
                round(flipfrac * 100, digits=1), "%  ", where(cam),
                "  unresolved=", count(==(0), raymap.side))
        flush(stdout)
        # steer the next leg on the coherent bulk motion, floored and capped
        uniform > 0 || (Δτ = clamp(Δτ * clamp(0.85 * θ_target / (medflow + 1e-9), 0.6, 1.8), Δτ_floor, 0.4))
    end
    println(frame, " frames over τ=", τ_total, " in ", round((time() - t0) / 60, digits=1), " min")
end

main()
