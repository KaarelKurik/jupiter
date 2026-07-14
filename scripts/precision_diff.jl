# precision diff: the F32 acceptance instrument (step-5 Float32 policy).
# Traces deterministic ray bundles through the physics_diff scenes in four
# variants and reports the decomposition kaarel specified (2026-07-14b):
# side-flip rate (size of the undecidable band), deviation conditional on
# side agreement as quantiles / exceedance fractions (the conditional tail is
# power-law, so means need not exist), binned by passage count (the cheap
# proxy ordinal for distance to the limit-cycle boundary). Sup error over the
# bundle is deliberately NOT the headline: the exit map is discontinuous, so
# worst-case error is unbounded a priori for any finite-precision method.
#
# Variants, factoring method-ε into its sources:
#   truth  F64 arithmetic on F64 tables — the oracle
#   ulp    truth with the initial direction perturbed by ~1 ulp — the
#          "just rounding" baseline: quantile curves parallel to (offset
#          above) this baseline mean pure-precision scaling, bends mean
#          method pathology
#   tab32  F64 arithmetic on Float32-rounded polynomial tables — isolates
#          representation rounding (what the GPU's tables cost by themselves)
#   f32    Float32 phases on Float32 tables — the GPU candidate: mouth entry
#          stays an F64 solve (once per passage, well-conditioned), the
#          in-throat trace runs wholly in F32, exits convert at to_ambient
#
#   julia --project --threads=auto scripts/precision_diff.jl [nrays]

using jupiter
using LinearAlgebra
const J = jupiter

f32_surface(s) = J.Surface(J.mesh(s), J.chart_polys(s), [Float32.(p) for p in J.packed_polys(s)])
f32_throat(th) = J.Throat(f32_surface(J.geometry(th)), J.params(th), th.placements)

function cube_config()
    m = J.cubemesh()
    surf = J.Surface(m, J.fit_geometry(m))
    th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.5, 0.75))
    (throat = th, budget = J.RayBudget(0.05, 400, 4),
     origin = [0.9, -2.8, 0.9], aim = [0.0, 0.0, 0.0], spread = 0.5, mouth_res = 10)
end

function trefoil_config()
    m = J.load_obj(joinpath(@__DIR__, "..", "res", "models", "trefoil.obj"))
    surf = J.Surface(m, J.fit_geometry(m))
    th = J.Throat(surf, J.ThroatParams(1.0, 1.0, 0.2, 0.3))
    (throat = th, budget = J.RayBudget(0.05, 400, 4),
     origin = [0.5, -7.5, 3.5], aim = [0.262, -1.0, 0.437], spread = 0.35, mouth_res = 2)
end

const SCENES = [("cube", cube_config), ("trefoil", trefoil_config)]

mouths(env, th, res) = (J.TessellatedMouth(env, J.HalfThroat(th, 1), res),
                        J.TessellatedMouth(env, J.HalfThroat(th, 2), res))

function bundle_directions(config, n) # golden-angle spiral, no RNG stream to drift
    forward = J.generic_normalize(config.aim - config.origin)
    e1 = J.generic_normalize(cross(forward, [0.0, 0.0, 1.0]))
    e2 = cross(forward, e1)
    golden = 2pi * 0.6180339887498949
    [J.generic_normalize(forward + config.spread * sqrt(k / n) *
                         (cos(golden * k) * e1 + sin(golden * k) * e2))
     for k in 1:n]
end

"""
trace_ray with the phase converted to T at each mouth entry (the entry solve
itself stays F64) and the passage count reported; T = Float64 on the F64
throat reproduces render.jl's trace_ray exactly. Returns
(side, exit vel, exit pos, passages) with side 0 = unresolved within budget.
"""
function traced_exit(env, th, mts, budget, origin, dir, T)
    ray = J.AmbientRay(J.HalfThroat(th, 1), origin, dir)
    for passage in 1:budget.max_passages
        mouth = mts[J.side(J.half_throat(ray))]
        v = J.enter_mouth(env, mouth, ray.pos, ray.vel)
        v === nothing && return (J.side(J.half_throat(ray)), ray.vel, ray.pos, passage - 1)
        vT = J.SituatedPhase(v.chart, T.(v.pos), T.(v.vel))
        res = J.trace_geodesic(env, vT, budget.step_size, budget.max_steps)
        res isa J.SituatedPhase && return (0, zeros(3), zeros(3), passage) # out of steps inside
        ray = J.AmbientRay(J.half_throat(res), Float64.(res.pos), Float64.(res.vel))
    end
    (0, zeros(3), zeros(3), budget.max_passages) # out of passages: orbiting near a limit cycle
end

function bundle_exits(env, th, mts, budget, origin, dirs, T)
    out = Vector{Tuple{Int, Vector{Float64}, Vector{Float64}, Int}}(undef, length(dirs))
    Threads.@threads for i in eachindex(dirs)
        out[i] = traced_exit(env, th, mts, budget, origin, dirs[i], T)
    end
    out
end

# angular deviation of exit directions; robust for tiny angles unlike acos
function angle_between(v1, v2)
    n1, n2 = J.generic_normalize(v1), J.generic_normalize(v2)
    2 * asin(clamp(norm(n1 - n2) / 2, 0.0, 1.0))
end

quant(sorted, p) = sorted[clamp(ceil(Int, p * length(sorted)), 1, length(sorted))]

"""
side agreement + conditional deviation of a variant against truth; angles only
where both resolve to the same side (the decomposition's second factor)
"""
function against_truth(truth, variant)
    n = length(truth)
    flips = count(t[1] != v[1] && t[1] != 0 && v[1] != 0 for (t, v) in zip(truth, variant))
    res_mismatch = count((t[1] == 0) != (v[1] == 0) for (t, v) in zip(truth, variant))
    both_unresolved = count(t[1] == 0 && v[1] == 0 for (t, v) in zip(truth, variant))
    agree = [(t, v) for (t, v) in zip(truth, variant) if t[1] == v[1] != 0]
    angles = sort([angle_between(t[2], v[2]) for (t, v) in agree])
    (n = n, flips = flips, res_mismatch = res_mismatch, both_unresolved = both_unresolved,
     agree = length(agree), angles = angles)
end

fmt(x) = x == 0 ? "0" : string(round(x, sigdigits = 3))

function report(name, stats)
    a = stats.angles
    println("  ", rpad(name, 6),
            " flips ", stats.flips, "/", stats.n,
            "  res-mismatch ", stats.res_mismatch,
            "  unresolved(both) ", stats.both_unresolved,
            "  | angle q50/q90/q99/max ",
            join([fmt(quant(a, p)) for p in (0.5, 0.9, 0.99)], " "), " ", fmt(a[end]),
            "  exceed 1e-4/1e-3 rad ",
            fmt(count(>(1e-4), a) / length(a)), " ", fmt(count(>(1e-3), a) / length(a)))
end

function passage_binned(truth, variant)
    for b in 0:maximum(t[4] for t in truth)
        sel = [(t, v) for (t, v) in zip(truth, variant) if t[4] == b]
        isempty(sel) && continue
        flips = count(t[1] != v[1] for (t, v) in sel)
        agree = [(t, v) for (t, v) in sel if t[1] == v[1] != 0]
        angles = sort([angle_between(t[2], v[2]) for (t, v) in agree])
        println("    passages=", b, "  n=", length(sel), "  flips=", flips,
                isempty(angles) ? "" : "  angle q50 " * fmt(quant(angles, 0.5)) *
                                       "  q90 " * fmt(quant(angles, 0.9)) *
                                       "  max " * fmt(angles[end]))
    end
end

function main()
    nrays = length(ARGS) >= 1 ? parse(Int, ARGS[1]) : 2000
    println("precision_diff: ", nrays, " rays/scene on ", Threads.nthreads(), " threads")
    for (name, config_fn) in SCENES
        config = config_fn()
        th = config.throat
        th32 = f32_throat(th)
        mts64 = mouths(nothing, th, config.mouth_res)
        mts32 = mouths(nothing, th32, config.mouth_res)
        dirs = bundle_directions(config, nrays)
        dirs_ulp = [J.generic_normalize(nextfloat.(d)) for d in dirs]

        t = @elapsed begin
            truth = bundle_exits(nothing, th, mts64, config.budget, config.origin, dirs, Float64)
            ulp   = bundle_exits(nothing, th, mts64, config.budget, config.origin, dirs_ulp, Float64)
            tab32 = bundle_exits(nothing, th32, mts32, config.budget, config.origin, dirs, Float64)
            f32   = bundle_exits(nothing, th32, mts32, config.budget, config.origin, dirs, Float32)
        end
        println(name, " (", round(t, digits = 1), "s):")
        s_ulp, s_tab, s_f32 = against_truth(truth, ulp), against_truth(truth, tab32), against_truth(truth, f32)
        report("ulp", s_ulp)
        report("tab32", s_tab)
        report("f32", s_f32)
        # quantile-curve offset vs the ulp baseline: log2 ratio at matched
        # quantiles; flat ≈ 29 (= 52 − 23 mantissa bits) means "just rounding"
        qs = (0.25, 0.5, 0.75, 0.9, 0.95, 0.99)
        offs = [log2(quant(s_f32.angles, p) / max(quant(s_ulp.angles, p), 1e-300)) for p in qs]
        println("  f32/ulp quantile offset (log2, at ", join(qs, "/"), "): ",
                join([fmt(o) for o in offs], " "))
        println("  f32 by truth passage count:")
        passage_binned(truth, f32)
    end
end

main()
