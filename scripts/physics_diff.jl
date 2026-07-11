# physics diff: trace a deterministic ray bundle through reference scenes and
# compare exit rays against committed baselines. The standing check that an
# optimization changed only the arithmetic, not the physics.
#
#   julia --project scripts/physics_diff.jl          compare against baselines
#   julia --project scripts/physics_diff.jl --save   (re)write the baselines
#
# Reading a failure: exit *sides* are discrete and should essentially never
# flip (a flip means a ray traversed differently — semantic change). pos/vel
# deviations up to ~1e-9 are consistent with pure arithmetic reordering
# amplified through the throat's chaotic stretches; well beyond that means
# investigate. Directions come from an explicit golden-angle spiral and
# baselines are plain text, so both survive Julia version changes.

using jupiter
using LinearAlgebra
const J = jupiter

const TOLERANCE = 1e-9
const BUNDLE_SIZE = 80
const BASELINE_DIR = joinpath(@__DIR__, "..", "res", "physics_diff")

function cube_config()
    m = J.cubemesh()
    surf = J.Surface(m, J.fit_geometry(m))
    th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75))
    mouths = (J.TessellatedMouth(nothing, J.HalfThroat(th, 1), 10),
              J.TessellatedMouth(nothing, J.HalfThroat(th, 2), 10))
    (scene = J.Scene(th, mouths, J.checker_sky), budget = J.RayBudget(0.05, 400, 4),
     origin = [0.9, -2.8, 0.9], aim = [0.0, 0.0, 0.0], spread = 0.5)
end

function trefoil_config()
    m = J.load_obj(joinpath(@__DIR__, "..", "res", "models", "trefoil.obj"))
    surf = J.Surface(m, J.fit_geometry(m))
    th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.2, 0.3))
    mouths = (J.TessellatedMouth(nothing, J.HalfThroat(th, 1), 2),
              J.TessellatedMouth(nothing, J.HalfThroat(th, 2), 2))
    (scene = J.Scene(th, mouths, J.checker_sky), budget = J.RayBudget(0.05, 400, 4),
     origin = [0.5, -7.5, 3.5], aim = [0.262, -1.0, 0.437], spread = 0.35)
end

const SCENES = [("cube", cube_config), ("trefoil", trefoil_config)]

"""
golden-angle spiral of directions in a cone around aim − origin; deterministic
by explicit construction, no RNG stream to drift between Julia versions
"""
function bundle_directions(config)
    forward = J.generic_normalize(config.aim - config.origin)
    e1 = J.generic_normalize(cross(forward, [0.0, 0.0, 1.0]))
    e2 = cross(forward, e1)
    golden = 2pi * 0.6180339887498949
    [J.generic_normalize(forward + config.spread * sqrt(k / BUNDLE_SIZE) *
                         (cos(golden * k) * e1 + sin(golden * k) * e2))
     for k in 1:BUNDLE_SIZE]
end

function bundle_exits(config)
    map(bundle_directions(config)) do dir
        ray = J.AmbientRay(J.HalfThroat(config.scene.throat, 1), config.origin, dir)
        out = J.trace_ray(nothing, config.scene, config.budget, ray)
        out === nothing ? (0, zeros(3), zeros(3)) : (J.side(J.half_throat(out)), out.pos, out.vel)
    end
end

baseline_path(name) = joinpath(BASELINE_DIR, name * ".baseline")

function save_baseline(name, exits)
    mkpath(BASELINE_DIR)
    open(baseline_path(name), "w") do io
        println(io, "# physics_diff baseline: ", name, "; per ray: side pos_xyz vel_xyz (side 0 = unresolved)")
        for (side, pos, vel) in exits
            println(io, side, " ", join(string.([pos; vel]), " "))
        end
    end
end

function load_baseline(name)
    map(filter(l -> !startswith(l, "#"), readlines(baseline_path(name)))) do line
        parts = split(line)
        nums = parse.(Float64, parts[2:7])
        (parse(Int, parts[1]), nums[1:3], nums[4:6])
    end
end

function compare(name, exits)
    ref = load_baseline(name)
    side_flips = count(a[1] != b[1] for (a, b) in zip(exits, ref))
    deviation = maximum((a[1] == b[1] != 0 ? max(maximum(abs, a[2] - b[2]), maximum(abs, a[3] - b[3])) : 0.0)
                        for (a, b) in zip(exits, ref))
    ok = side_flips == 0 && deviation < TOLERANCE
    println(rpad(name, 8), " side flips: ", side_flips, "/", length(ref),
            "  worst pos/vel deviation: ", round(deviation, sigdigits=3),
            "  ", ok ? "OK" : "DEVIATES (tolerance $(TOLERANCE))")
    ok
end

function main()
    save = "--save" in ARGS
    all_ok = true
    for (name, config_fn) in SCENES
        exits = bundle_exits(config_fn())
        if save
            save_baseline(name, exits)
            println(rpad(name, 8), " baseline written (", length(exits), " rays, ",
                    count(e -> e[1] == 0, exits), " unresolved)")
        else
            all_ok &= compare(name, exits)
        end
    end
    save || exit(all_ok ? 0 : 1)
end

main()
