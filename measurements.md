# Measurements log

Findings worth not re-deriving, but which don't change the working plan.
Probe scripts are deliberately disposable (they'd be rewritten against changed
code anyway); enough method detail lives here to reconstruct them.

## 2026-07-13b — half-edge flattening (GPU groundwork step 1): Dicts weren't a measurable CPU cost; immutable-Mesh inlining grows the passage boxes

**Motivation.** GPU port step 1: replace Mesh's name-keyed Dicts
(`half_edges`, `half_edge_offsets`) with flat id-indexed arrays so the
hot-loop connectivity is device-shaped. Certified as a pure representation
change (bit-identity), with before/after timing and allocation.

**Method.** Three deterministic knot-hitting trefoil rays (physics_diff aim +
two 0.01 tilts toward e1/e2), cold runs, min of 5; `@allocated` per ray;
`Profile.Allocs` at sample_rate 1 by type, before-side run via a jj workspace
at the parent commit.

**Findings.**

- **Timing unchanged**: 3.27/2.32/2.29 → 3.27/2.31/2.28 ms. The per-step Dict
  lookups (`next`/`twin`/`prev`/`half_edge_offset` in chart transitions) were
  not a measurable CPU cost — the step loop is dominated by collar/Γ
  arithmetic. The value of the flattening is representational (GPU-ready
  tables), not CPU speed.
- **Bit-identity trap found and dodged**: `catmullclark` numbers the refined
  mesh's faces by `collect(he_to_face)` — Dict hash-iteration order. The
  rebuild helpers (`half_edge_name_to_face_index`/`..._to_next_name`)
  reconstruct the transient Dict with the historical insertion sequence, so
  iteration order — hence CC numbering, hence fitted polynomials — is
  unchanged. Confirmed: physics_diff worst deviations bit-identical to the
  previous certified run (1.34e-14 cube, 1.84e-10 trefoil), 0 flips, 275/275
  cold.
- **Allocation +240–416 B/ray, zero new sites**: same 7 per-passage sum-type
  boxes (2-passage knot ray). Mechanism: Mesh is immutable and Julia inlines
  immutable structs into parents, so every box embedding a Mesh-carrying
  handle inlines Mesh's field block, now 11 pointer fields (88 B) vs 5
  (40 B). A `SituatedPhase` box holds a Chart = two embedded Mesh copies
  (HalfThroat chain + HalfEdgeHandle) → +96 B each, matching the profile
  (800 → 992 B for n=2). Bytes not boxes, per-passage not per-step; the
  device port carries only the flat tables. If it ever matters: a mutable
  Mesh (by-reference everywhere) or a slimmer hot handle.
- **Step 2 (mouth concretization): two BVH stack representations died before
  the threaded form.** Replacing the traversal recursion with an explicit
  stack for device parity: (1) `MVector{64,Int}` heap-allocates 544 B per
  entry attempt (~+1.6 kB/ray) — MArray `setindex!` goes through
  `pointer_from_objref`, which defeats escape analysis; (2) an
  `NTuple{64,Int}` stack via `Base.setindex` explodes to 330–500 kB/ray —
  the 66-argument `_setindex` splat exceeds vararg specialization and every
  push boxes the whole tuple (also ~+3% time). Resolution: **threaded BVH**
  (skip-pointer DFS, one extra Int per node computed in a reverse sweep at
  build) — no stack at all, allocation exactly back to the step-1 floor,
  timing unchanged, and it's the GPU-native shape (no local-memory stack).
  Right-child-first descent preserves the recursive visit order, so best-hit
  evolution is bit-identical: physics_diff again 1.34e-14/1.84e-10, 0 flips,
  275/275 cold.
- **Step 3 (device smoke), RTX 4070 SUPER, Float64, 1M points each**
  (scripts/gpu_smoke.jl, gpu/ subenv): eval_packed Horner kernel
  **bit-identical** to CPU (max |Δ| = 0.0), 0.78 ms = 1.29 Geval/s, 76x one
  CPU thread. The ad.jl dual pass (value_and_jacobian_columns through
  square_coords_to_chart → fake_complex_pow/safe_atan2) — ForwardDiff duals,
  StaticArrays, closures, SMatrix hcat — compiles clean on device: 5.8 ms =
  172 Meval/s, 31x one thread, max deviation 7.2e-16 value / 5.0e-15
  jacobian (fma/libm reordering). 56/78 registers — headroom for a full RK4
  step. **No compilability wall in the AD core.** Remaining port cost is
  what we already knew: device throat view (flattened packed_polys), and the
  per-passage sum-type boxes (kernels can't heap-allocate → wavefront
  restructure). eval_packed's signature loosened to AbstractArray{Float64,3}
  so device arrays dispatch to the same body.

## 2026-07-13 — temporal sampling for fly-through video: the flow-rate profile, and two ways a flow controller can die

**Motivation.** How densely to sample a fly-through in the flight parameter τ,
for both a perceptually-smooth adaptive video and an honest constant-speed one.

**Method.** Optical flow measured directly: sparse probe raymaps (16×12 rays)
before/after each candidate step; per-pixel angular change of exit directions
(same-side pixels), expressed in render pixels (192-wide, fov π/3.2 →
0.29°/px); side-flip fraction as a second signal. Adaptive controller in the
DP5 accept/reject shape, growth factor 0.85·target/flow. Cube scene, textured
skies, camera straight through the throat, τ 0→6.5. scripts/flyvideo.jl.

**Findings.**

- **Flow-rate profile** (px of bulk feature motion per unit τ, 3px/frame
  target → Δτ = 3/rate): far approach ~195; final exterior approach (τ
  1.5–2.0, mouth filling the view) ~610 median, 805 peak — the global worst;
  entry zone ~480; throat interior 115–160; departure ~0 (no parallax against
  the far sky — Δτ immediately maxes out; the whole exit took 7 frames).
  **Constant-speed answer: uniform Δτ = target/peak = 0.0037 for 3 px.**
- **Controller death #1 — max statistic.** Near-critical-ring pixels' exit
  directions change ~2.5px+ for ANY camera displacement (chaotic
  sensitivity, does not scale with Δτ). Steering on max flow pinned Δτ at
  ~0.002 from the very start, 100x oversampling.
- **Controller death #2 — any fixed upper quantile.** q90 worked until the
  mouth filled >10% of the probe with near-critical pixels (τ≈1.4+), then
  pinned exactly like max (2.5px at Δτ=0.0008; measured max there up to
  319px) and the un-floored growth formula ground to a halt: 630 frames to
  reach τ=1.77. Fix: steer on the **median** + a hard Δτ floor (0.004);
  q90/max stay in the log as an audit of the chaotic population.
- **Chaotic shimmer is image-space aliasing, not temporal undersampling** —
  no frame density smooths it (it needs ray bundles / supersampling, see
  roadmap). The sensitive population clusters on the **boundary between the
  two ambient spaces' images** in the camera's view; kaarel's correction of
  an earlier "grazing rays" reading: the truly *expensive* rays are the
  limit-cycle ones sitting on that boundary — they are the unresolved
  pixels. Whether that same set drives the flow-controller trouble is
  plausible but unconfirmed. This reading explains the entry behavior
  cleanly: the moment the camera crossed d=0 the chaos collapsed
  (q90≈max≈median, flips 0, unresolved 0) because the two-image boundary
  left the forward FOV — not because the interior is special. Consistently,
  the side-2 look-back frame (flythrough demo frame 9), which points at that
  boundary again, had 60 unresolved pixels.
- Run artifacts: 839 frames (652 from the q90 run whose tail oversamples,
  187 from the resumed median run — resume is exact because the pre-entry
  camera state is closed-form in flat space), resampled to 314 at uniform
  3px perceptual speed (scratch script; rate profile from the logs), videos
  in out/flythrough_adaptive{,_raw}.mp4, logs preserved as
  out/flyvideo_run{1,2}.log. Render cost ~19min for the median-run portion
  at 192×144.

## 2026-07-12 — jacobian_columns twin retired: the inner pass is honestly a different operation, and fusing it makes rays ~18% faster

**Motivation.** The `nested_jacobian_columns` twin (below) worked but smelled
like a workaround. Question: is there a way to work *with* inference instead?

**Findings.**

- No user-facing knob exists for the per-method recursion-widening heuristic
  (it fires when the same method appears on inference's own stack with a
  grown signature — nested `Dual{Tag2, Dual{Tag1}}` types guarantee growth).
  Alternatives weighed: a concrete typeassert at the inner call site stops
  the Any from propagating downstream but leaves a dynamic dispatch plus a
  boxed SMatrix-of-nested-Dual return per collar evaluation (~600 B × several
  per step ≈ right back at MB/ray); two `directional` passes in
  surface_normal_out breaks the cycle via a different method but pays an
  extra primal surface evaluation; OpaqueClosure barriers are exotic. The
  clean resolution: give the inner level *different semantics* — precedented
  by ForwardDiff itself, whose hessian is jacobian-over-gradient, distinct
  methods per nesting level, never jacobian twice.
- The different semantics were already wanted: collar evaluated the surface
  twice (plain value + surface_normal_out's dual pass over the same point).
  New primitive `value_and_jacobian_columns` returns both from one dual pass;
  the value lane performs exactly the plain evaluation's arithmetic, so the
  value is bit-identical and free. `jacobian_columns` and it keep separate
  bodies (any shared helper would re-form the method cycle).
- **Numbers** (three deterministic knot-hitting trefoil rays: physics_diff
  aim + two small tilts, cold runs, min of 5): 4.08/5.37/4.98 ms →
  3.32/4.40/4.06 ms (−18% each); allocation unchanged at 2208 B/ray (the
  sum-type boxes — confirms no inference regression). Certified: 219/219
  cold, physics_diff 0 flips with worst deviations bit-identical to the
  pre-change run (1.34e-14 cube, 1.84e-10 trefoil).
- Reference keeps the two-evaluation collar: that *is* its simplest correct
  form; the fused production collar is exactly the kind of elaboration the
  equivalence testset exists to hold in place.

## 2026-07-12 — per-ray allocation: ~1 MB was one inference failure, not many small ones

**Motivation.** After the StaticArrays pass, ~1 MB/ray remained (trefoil,
knot-hitting). A priori the tracer needs ~0 heap: a geodesic step is a
fixed-size state update, outputs land in the preallocated raymap, so the
algorithm's honest floor is a few boxed sum-type returns per passage
(`enter_mouth` → `SituatedPhase | nothing`, `trace_geodesic` →
`AmbientRay | SituatedPhase`), i.e. well under 1 kB/ray.

**Method.** `@allocated` on `trace_ray` for three camera rays (trefoil 
scene), `Profile.Allocs` at sample_rate 1 grouped by nearest jupiter frame
and by type; `Test.@inferred` + `code_warntype` bisection down the call tree.

**Findings.**

- Measured 0.7–1.4 MB/ray, ~7 kB and ~42 allocations per RK4 step — closures,
  `SMatrix`-of-`Dual` intermediates, and `SituatedPhase`s all heap-boxed.
- Root cause was singular: Julia's *method self-recursion widening*. The chain
  `jacobian_columns → collar → surface_normal_out → jacobian_columns` puts the
  same method on the inference stack twice with a grown signature; the
  termination heuristic widens the inner call to `Any`, which poisons
  `outer_metric → metric → christoffel → wvel_along_v → geodesic_step` — the
  entire step loop ran boxed. Every callee inferred fine in isolation; only
  the nested chain failed. Fix: a structurally identical twin
  (`nested_jacobian_columns`) for the inner level.
- Same trap, second appearance: rewriting the BVH traversal recursively with a
  `nothing`-then-tuple accumulator re-tripped the identical widening
  (signature grows Nothing → Tuple across the self-call) and re-poisoned
  `enter_mouth` to 54 kB/ray. Threading fixed-type state
  `(t, (u, v), tri_index; 0 = none)` fixes it. Moral: in hot recursive/nested
  code, keep self-call signatures *exactly* constant — Union accumulators and
  method reuse across nesting levels are both inference hazards that profile
  as "death by boxing" with no single hot allocation site.
- Cleanups riding along: `next`/`twin`/`prev` allocated a Vector per call via
  NTuple range-slicing (`vertices[3:-1:2]` hits a generic `map`); explicit
  index pairs now. `AmbientRay{P,V}` parametrized; `RayBudget`/`Camera` fields
  concretized; `Scene` parametric in mouths+sky; BVH stack Vector replaced by
  the recursion above (right-child-first preserves the old pop order exactly).
- **After:** 1.5 kB and 3 allocations per knot-hitting ray (the sum-type
  boxes; per-passage, step loop allocation-free). Knot-hitting ray 5.5 →
  2.8 ms; trefoil 384×288 render 14.8 → 10.1 s wall on 16 threads. Certified:
  219/219 cold, physics_diff 0 flips, worst deviation 1.8e-10 (FMA/inlining
  reordering, within the 1e-9 budget).

## 2026-07-11 — grazing-ray deflection vs entry angle

**Motivation.** Renders (cube first light, trefoil at cross_scale 1.0 and 0.4)
show lensing that reads as aggressively switched-on near the mouth silhouette.
Question: is deflection a continuous/smooth function of the grazing angle, as
the smooth-by-construction metric demands, or did the machinery break somewhere?

**Method.** Construct phases directly on the d=0 surface, bypassing
`enter_mouth` entirely: at a chart point (square coords [0.3, 0.3] of wedge 0),
vel = cos θ · t̂ + sin θ · ê₃ with t̂ a metric-normalized surface tangent and ê₃
the depth axis (exactly metric-orthonormal to the tangents at d=0, so θ is an
honest angle; checked g33 = 1, g13 = g23 = 0). Trace with the raw
`geodesic_step`/`settle_phase` loop; deflection δ = angle between ambient entry
and exit velocities via `to_ambient`. Two exit conventions: (a) `exits_mouth`
as in production; (b) for tangency studies, integrate through shallow d < 0
dips (metric there is exactly the flat collar pullback) and exit at d < −0.05 —
convention (a) truncates the θ < 0 family at concave points before the surface
curves back into the ray's path, which masquerades as a discontinuity at θ = 0.

**Findings.**

- Cube (params 1.0/1.0/0.5/0.75, mid-face, along-edge direction): δ < 0.001°
  until θ ≈ 0.3 rad, i.e. until the ray's deepest point reaches ~25% of
  cylinder_depth — the C∞-flat blend really does keep shallow dips in exactly
  flat geometry. Then a smooth rise (0.08° at θ=0.36, 2.7° at 0.47, 6° at 0.50)
  steepening logarithmically toward capture at θ_c ≈ 0.566; a chaotic
  escape-side band θ ∈ [0.566, 0.574]; then a smooth side-2 traversal branch.
  Step-halving (h = 0.02 → 0.01) moves smooth-branch deflections by ~0.001°
  and does not move band boundaries.
- Trefoil (1.0/1.0/0.2/0.3), convex directions: same shape, remarkably uniform
  across the mesh — δ > 0.5° from θ ≈ 0.36, capture at θ_c ≈ 0.48 (5 vertices).
- Trefoil, concave directions (inner tube side): tangential rays plunge to full
  transition depth — there is no small-deflection grazing regime. Under exit
  convention (b), δ(θ) is smooth straight across tangency (v=1000: δ = 130.69°/
  130.73°/130.77° at θ = −0.002/0/+0.002, dδ/dθ itself smooth). v=4000 sits in
  a genuine chaotic-scattering zone for |θ| ≲ 0.06 (all rays reach transition
  depth, ~1000-step windings, escape side alternating in bands whose locations
  are h-converged at h = 0.02/0.01/0.005), with no feature at θ = 0 and a
  smooth branch (incl. smooth max of δ) for θ > 0.016.

**Conclusions.** (See also the collar-caustic entry below: cylinder_depth is
not capped by collar self-intersection, so it is a genuinely free gentleness
knob, up to the in-principle limits in the last point.)

- The geodesic/chart-transition/integration machinery shows no discontinuity or
  step-size artifact anywhere: smooth branches are smooth and h-converged;
  the fine-scale structure is chaotic scattering near the trapped set (closed
  geodesics), which a smooth metric with a throat must have.
- The "aggressive" look is the metric design, not a bug: (1) exactly flat
  outside d=0, so near-miss rays bend not at all — lensing has compact support
  ending in a hard edge at convex silhouettes; (2) exp-flat blend, so entries
  shallower than ~10–25% of cylinder_depth exit visually straight — the visible
  bending is compressed into a thin annulus; (3) concave-side grazing rays
  plunge and scramble. Softening knobs are limited in principle: head-on rays
  must traverse, so all lensing has to happen between a grazing ray and the
  nearest head-on ray, and there is not a lot of daylight between them.

## Recorded 2026-07-11, measured in an earlier session — collar caustics are benign

**Question.** Past the limit surface's focal radius the collar map
self-intersects (a caustic; its Jacobian drops rank), which naively caps
cylinder_depth. Does pushing cylinder_depth past the caustic make the blended
geometry singular?

**Conclusion.** No — cylinder_depth can be made arbitrarily high. The final
space is a genuine Riemannian 3-manifold regardless: on the open band the
blend has strictly positive inner weight, and g_outer = JᵀJ stays smooth even
past the caustic (merely dropping rank on it), so the blend is smooth and
eigenvalue-floored by (1−w)·λ_min(g_inner). Empirically: with cylinder_depth
pushed to 0.6 (cube), det(g_outer) bottoms out at 3×10⁻⁸ at d = 0.5 (neatly
measuring the limit surface's focal radius) and re-grows, while the blend's
λ_min never dips below 0.017. The one true validity condition is immersedness
of the fitted YZ surface — unrelated to focal depth.
