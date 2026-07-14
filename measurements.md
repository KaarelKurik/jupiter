# Measurements log

Findings worth not re-deriving, but which don't change the working plan.
Probe scripts are deliberately disposable (they'd be rewritten against changed
code anyway); enough method detail lives here to reconstruct them.

## 2026-07-14c — F32 traced end-to-end: quantile curves parallel to the 1-ulp baseline (offset ≈ 2^28 = "just rounding"), the undecidable band is ~1e-4 of pixels and sits on the limit-cycle boundary, renders indistinguishable

**Motivation.** 2026-07-14b established F32 christoffel accuracy per *eval*
(median 9.5e-7 sup-rel) but flagged that Γ-eval error is not exit-map error:
F32 state accumulates over hundreds of steps and chaotic stretches amplify.
Kaarel's metric correction governs the design here: the exit map
image→(side, ray) is discontinuous at the limit-cycle boundary, so sup
exit-map error is vacuous for *any* finite-precision method — measure the
physics_diff decomposition instead (side-flip rate = size of the undecidable
band; deviation conditional on side agreement, in tail-safe aggregates,
since deviation ~ ε·κ(ray) with κ power-law near the boundary means even the
conditional *mean* need not exist), plus eyeballed render pairs as the
acceptance bar.

**Method.** Two parts.

- *Eltype-honesty completed past christoffel* (the 14b sweep's remainder):
  the dopri tableau ratios (previously `(1/5)k1`-style Float64 literals),
  the tracers' h/tol/`floatmin()`, RK4's h, `half_transition`'s
  transition_depth read, `to_ambient`'s Placement reads, `emit_ray`'s 1.0,
  `view_phase_at_target`'s rotation angle (kept exact in Float64, one
  rounding into the carrier — the wedge_square_coords idiom), and
  `directional`'s Float64 derivative seed in ad.jl all convert to
  `carrier(x)` instead of promoting. Float64 semantics certified unchanged —
  275 tests green, physics_diff deviation exactly 0.0 — since every
  conversion is F64→F64 identity there.
- *The F32 trace design* (scripts/precision_diff.jl — tracked, not
  disposable: it is the standing precision-acceptance instrument,
  physics_diff's sibling): mouth entry stays an F64 Newton solve against the
  F32-coefficient surface (once per passage, well-conditioned, and the
  plausible wavefront split anyway); the phase converts to F32 at entry; the
  in-throat trace runs wholly in F32 on Float32-rounded polynomial tables
  (mirroring device_throat — with F64 tables eval_packed would silently
  promote every step back to F64); exits convert at to_ambient and re-enter
  the next passage in F64. 8000-ray golden-angle bundles through the two
  physics_diff scenes (budget 0.05/400/4), four variants: **truth** (F64),
  **ulp** (F64, direction perturbed by ~1 ulp — the "just rounding"
  baseline), **tab32** (F64 arithmetic on F32 tables — isolates
  representation rounding from arithmetic rounding), **f32** (the GPU
  candidate). Angles via 2·asin(‖n̂₁−n̂₂‖/2); conditional stats over rays
  where both variants resolve to the same side (sky-miss rays are in the
  denominator and deviate by exactly 0 there — cube 67%, trefoil 35% of the
  bundle).

**Findings** (conditional angular deviation in radians; 1e-3 rad ≈ ⅓ pixel
at 384-wide/60°-fov, sub-pixel at anything plausible):

- cube: ulp q50/q90/q99/max 1.7e-16 / 3.3e-15 / 3.0e-14 / 7.5e-11;
  tab32 0 / 1.4e-7 / 1.4e-6 / 0.021; **f32 0 / 5.0e-7 / 2.9e-6 / 0.118**.
  Side flips 0/8000 every variant. Exceedance(1e-3) for f32: 0.04%
  (3 rays, the same grazing-rim rays that dominate tab32's tail).
- trefoil: ulp 9.1e-15 / 9.2e-14 / 3.8e-12 / 2.4e-9;
  tab32 8.7e-7 / 8.3e-6 / 3.2e-4 / 0.31; **f32 2.8e-6 / 3.1e-5 / 1.3e-3 /
  0.46**. Side flips: f32 1/8000 + 2 resolution mismatches (ulp and tab32:
  0). Exceedance for f32: 5.1% at 1e-4, **1.1% at 1e-3**.
- **The structural test passes**: the f32/ulp quantile offset on the
  trefoil is flat at 28.2–28.4 log2 across q50–q99 — parallel curves, no
  bends — i.e. F32 behaves as pure precision scaling of the same method,
  with offset ≈ the 2^29 mantissa-width ratio. (Cube offsets 26.6–27.7 on
  the entered subset; its bundle is mostly sky and single-passage.)
- **Error factorization** (tab32 vs f32, trefoil): representation rounding
  alone (F32 tables under F64 arithmetic) is ~⅓ of the f32 median
  (8.7e-7 of 2.8e-6) with the same tail shape; F32 step arithmetic
  multiplies ~3x on top. No arithmetic reordering can beat the tab32 curve
  while tables are F32 — that is the representation floor.
- **Passage count works as the boundary-distance ordinal**: trefoil f32
  medians by truth passage count 5.2e-6 (p=1, n=4201) → 4.0e-5 (p=2,
  n=728) → 2.0e-4 (p=3, n=207) — roughly ×8 per passage — and all three
  side flips live at p ≥ 3. The ulp baseline shows the same amplification
  from the other end: 1e-16 input noise reaches 2.4e-9 max — F64's own
  undecidable band is merely thinner, not absent.
- **Eyeball bar** (384×288 render pairs + diff images, both scenes):
  indistinguishable by eye. The diff structure is exactly the theory's
  picture — grayscale deviation on the knot silhouette ~1e-5, bright
  filaments along the boundary between the two ambient images, and the
  flipped pixels (cube 1/110592, trefoil 12/110592 = 0.011%) sitting *on*
  those filaments.

**Read.** F32 is accepted as the tracer arithmetic at the image standard
that already picked RK4: renders indistinguishable, wrong-for-precision
pixels ~1e-4 of the image concentrated where the F64-fallback flag was
already planned (the same filaments are where F64's own chaotic shimmer
lives — the ray-bundles item). The acceptance-aligned scalars to re-run as
the port proceeds: side-flip rate + exceedance-at-1-mrad, by passage count
(the cheap flag ordinal until Jacobi-field κ lands). The sweep also covers
emit_ray/transport, so in-throat cameras go F32 the same way; reference.jl
stays untouched as the F64 oracle. Remaining step-5 work is the wavefront
restructure (+ depth binning, ~2x from 14b's divergence datum).

## 2026-07-14b — production christoffel straight on the GPU: Float64 loses to the CPU, Float32 matches it naive and beats every tabulation variant on accuracy

**Motivation.** Kaarel's call for this session: before tabulating Γ (and
before the representation fork that tabulation accuracy opened), try the
simpler road — run the production AD christoffel *directly* in a kernel and
see what the card actually does with it. The step-3 smoke test only compiled
the first-order dual pass; christoffel is second-order (a 3-partial seed over
`metric`, whose collar/surface evals are themselves dual passes — innermost
scalars are (1+3)(1+3)(1+2) = 48 carrier floats), which was the untested
compilability wall.

**Method.** The type layer was made storage-parametric (pure type
generalization; 275 tests + physics_diff bit-identical), so the real cube
`Throat` Adapt-s onto the device: ragged construction-only Mesh fields drop
to `nothing`, packed_polys pad into one 4D array behind a view wrapper
(padding sits at the leading end of the descending Horner loop — evaluation
bit-identical). The kernel builds the tracer's own `SituatedPhase` from
half-edge ids and calls `jupiter.christoffel` unmodified
(scripts/gpu_christoffel.jl, kept). 200k points sweep all 8 charts, all
wedges, st ∈ [0.02, 0.98]², and all three depth branches (25% pure outer
d<0, the rest through blend and pure inner). CPU truth: the same production
christoffel, Float64, 1 thread (20–22 µs/eval, matching the step-4 26–32 µs
on fancier charts). For the Float32 leg, the christoffel path was first made
eltype-honest — Float64 literals (sqrt2, the 2π sincos wedge angles, wedge
powers n/4 and 4/n, 0.5) and ThroatParams field reads convert to
`carrier(x)` = the innermost primal type through the dual nesting — with
Float64 certified bit-identical (tests + physics_diff 0.0) since every
conversion is F64→F64 identity there.

**Findings** (RTX 4070 SUPER, 1:64 FP64; npoints=200k, threads 64/128/256
indistinguishable, 512 exceeds the 64K regs/block budget):

- **No compilability wall.** The full second-order nested-dual christoffel —
  wedge transcendentals, blend, collar, 3×3 inverse — compiles clean and
  matches CPU to 2.1e-13 abs / 2.1e-14 sup-relative (device libm, not bugs).
- **Float64 is spill-bound and loses to the multicore CPU**: 255 registers
  (the cap) + 52,648 B/thread local spill; 3.9–4.1 µs/eval mixed = 5.0–5.5x
  one CPU thread ≈ 1/3 of the 16-thread CPU. Not a usable port.
- **Float32 is 2.8x faster than F64-on-device and ≈ the whole CPU, naive**:
  1.43 µs/eval mixed = 14.3x one CPU thread (0.70 Meval/s). Still 255 regs +
  27,672 B spill — spill-bound, not flop-bound, so this is a floor, not the
  card's ceiling.
- **Branch divergence costs ~2.3x on both types**: uniform-depth-branch
  launches give 2.09–2.14 (outer) / 2.13–2.14 (blend) / 0.14–0.17 (inner)
  µs/eval at F64, 0.58 / 0.63 / 0.06 at F32, vs 3.9 / 1.43 mixed — warps
  near-fully serialize the three depth branches (blend ≈ outer + inner as
  expected: the blend evaluates both metrics; inner is 15x cheaper, no
  collar jacobian). Depth-binning rays in a wavefront structure recovers
  ~2x ⇒ F32 ~29x one CPU thread without touching physics.
- **Float32 Γ accuracy beats every tabulation variant probed**: vs Float64
  truth, sup-relative per point median 9.5e-7, p99 3.8e-6, max 1.1e-5 —
  better than the C⁴-blend table target (~1e-5), far better than exp-flat
  (~1e-3), and it is the *exact* blend semantics merely rounded, so nothing
  is re-fit, re-baselined, or re-rendered.

**Read.** The corner-blend accuracy problem that opened the representation
fork does not gate the GPU road: the card can afford to evaluate the truth.
The fork stays open as a CPU-side/representation question (tabulation still
wins 10–50x *on the CPU*, and affine atlases/level sets have independent
virtues), but the port no longer waits on it. What F32 does **not** yet
answer is end-to-end: Γ-eval error is not exit-map error — F32 positions
accumulate over hundreds of steps and chaotic bands amplify; the planned
policy (F32 with F64 fallback for flagged pixels) needs a traced comparison,
which wants the eltype-honesty sweep extended past christoffel to the
steppers/transitions (settle, half_transition, the integrators still carry
F64 literals). That, plus the wavefront restructure, is the remaining
step-5 work; register-pressure surgery (splitting the fused dual passes) is
the optimization deliberately not done today.

## 2026-07-14 — tabulated Γ (GPU step 4) probed: depth factors out exactly; lateral convergence is blend-limited, and the two blends fail in complementary places

**Motivation.** Step 4 planned "Chebyshev-tabulated Γ per face×depth box
(~8³ coeffs/component) … geometric convergence for a C^∞ metric, likely
10–50x per step". Before implementing: measure where integrator stages
actually evaluate Γ, and what Chebyshev convergence the blended metric really
delivers. Face-coordinate tables were rejected on paper first: face coords
are z^(4/n)-singular at irregular corners, while chart coords are smooth by
construction — so tables live in chart coords.

**Method.** Four probe rounds (scratchpad, disposable). (1) A recording env
(`christoffel(env::Recorder, v)` delegating to the AD path — the env slot is
a clean dispatch extension point) traced the physics_diff bundles and logged
per-eval wedge-square coords and depth. (2–4) Tensor-Chebyshev fits
(first-kind nodes, hand-rolled DCT + series-derivative recurrence) of Γ and
of exact-in-d metric pieces, on wedge-aligned domains; end-to-end assembled-Γ
error vs the AD christoffel on the needed region; a blend-parametrized local
copy of the surface→metric→Γ pipeline (reference-style ForwardDiff, verified
4.4e-16 against production) to compare corner blends; tiny-interior-box
anchor test (collapses to g~1e-14/Γ~1e-12, validating the pipeline — a
round-2 "anchor" that spanned the full wedge angle failed at 0.75 and was
meaningless, not a bug).

**Findings.**

- **Stage excursions are large; tables need margins + an AD fallback.**
  Settled phases live in st ∈ [0,1/2]², d ∈ [0, td], but RK/DP stages
  evaluate beyond: cube maxst q99/max 0.51/0.53, d ∈ [−0.036, 0.79]; trefoil
  maxst q99/max 0.80/1.12 with 19% of evals past 0.5 (fixed h=0.05 ≈ half a
  face on that mesh), d ∈ [−0.048, 0.35]. Rays genuinely dip to d<0 with
  inward velocity before `exits_mouth` triggers.
- **Depth is exactly solved, no approximation needed**: g_outer(u,v,d) is
  exactly quadratic in d (collar = s − d·n̂ ⇒ JᵀJ degree 2; residual 6e-16),
  g_inner is d-independent, and the depth blend w(d/cyl) with its flat ends
  reproduces the outer/inner branch structure of `metric` exactly (residual
  3e-16, including d<0 and d>cyl). So tabulate 2D lateral pieces A,B,C,im
  with g = w·(A+Bd+Cd²)+(1−w)·im and ∂_d g = w′·(om−im)+w·(B+2Cd); w, w′
  evaluated exactly at runtime. The depth blend never needs to change and
  never enters the approximation error.
- **Naive boxes fail; seam-aligned polar sectors are the right lateral
  domain.** Wedge-frame bounding boxes straddle the corner-blend seam rays
  and barely converge (cube Γ rel error 0.28→0.12 from N=6→20). Per-wedge
  (r,φ) sectors put seams and chart center on box edges; a small Cartesian
  core box avoids the 1/r series-derivative blowup at the center. Also the
  partition of unity dies outside the face ring (0/0), so domains must stay
  inside st ≲ 1 — box validity is a hard constraint, checked at fit nodes.
- **Convergence under the production exp-flat corner blend plateaus at
  ~1e-3 relative Γ.** Trefoil (n=4, wedge maps trivial): rel 2.8e-4 at N=24
  single sector. Whole-mesh at feasible N: q50 per-chart Γsup 3.3e-3 (N=8) /
  1.7e-3 (N=12), zero charts below 1e-4, 0.84/1.9 GB uniform. Cube (n=3):
  sectors stuck at ~0.013 sup even at N=24; the exp-flat *core* however
  converges cleanly when refined (1.2e-6 at Ncore=32, rc=0.05R — exp-flat is
  flat to all orders at the vertex).
- **A C⁴ polynomial corner blend (nonic smootherstep) flips the picture.**
  Cube sectors collapse 1000x to 1.2e-5 at N=24; trefoil reaches rel 8.4e-5
  (N=16) / 2.1e-5 (N=24) and still dropping steeply. But at irregular
  vertices its core stops converging (~6e-3 sup, insensitive to rc: the
  polynomial blend composed with the z^(n/4) wedge coords leaves a genuine
  fractional-power singularity at the vertex, where exp-flat was flat).
  Complementary failure modes: exp-flat is bad at seams / good at centers;
  polynomial is good at seams / bad at irregular centers. A C² blend
  (quintic) is strictly worse than both for Γ (needs ∂g ⇒ only C¹; rel
  errors 0.2–0.9, non-converging) — C⁴ is the minimum interesting order.
- **The AD christoffel costs 26–32 µs/eval warm** (cube/trefoil, single
  thread) — it dominates ray cost (4 evals/RK4 step, 7/DP5 attempt).
  Estimated table eval (values + derivatives via Clenshaw from one 24×N²
  coefficient set, N=8–12, plus ~300 flop Γ assembly) is ~1–4 µs: the
  roadmap's 10–50x per-step stands *if* the accuracy question is settled.
  Derivative tables need not be stored (evaluate d/du series from value
  coefficients), so memory is 24·N²·8B per wedge ≈ 12–28 kB.

**Conclusions.** Tabulation is worth it (26–32 µs → ~1–4 µs per Γ), the
exact-in-d factorization + seam-aligned sectors + tiny-core-AD-fallback is
the right architecture, and accuracy is decided by the corner blend, which
plan/claude.md records as "placeholder pending a deliberate choice":

- Keep exp-flat ⇒ tables certify at ~1e-3 relative Γ (visually sub-pixel-ish,
  but a real tolerance downgrade in the reference chain).
- Adopt a C⁴ polynomial corner blend ⇒ ~1e-5 relative Γ at N=12–24 with AD
  fallback inside r < ~0.1R of irregular vertices (regular-vertex meshes like
  the trefoil have no hard cores at all). This is a deliberate physics change:
  surface geometry shifts subtly, all baselines re-save, galleries re-render.
- Orthogonal lever either way: per-chart step-size discipline (h=0.05 is half
  a face on the trefoil) would shrink the needed margins and every error above.

## 2026-07-13c — de-pinning hash iteration order from mesh construction: the slack is ~1e-13/e-10, within budget

**Motivation.** Step 1 below froze Dict hash-iteration order as implicit
semantics (the CC rebuild helpers reconstructed the historical insertion
sequence so refined-mesh numbering stayed bit-identical). Kaarel's call:
incidental hash order must not be pinned — allow the slack downstream and
let construction be canonical instead.

**Method.** Pre/post pointwise comparison of fitted surfaces
(`sample_surface` grids are anchored at original-mesh face handles, whose
numbering never changes, so points pair up across the change regardless of
gauge); physics_diff against the *old* baselines before re-saving; knot-ray
timing/alloc probe.

**Findings.**

- What the order actually fed: refined-mesh vertex/face *numbering* (pure
  labels), FP summation order in the CC vertex means, and each vertex's
  canonical frame = first neighbor (a gauge choice — rotates chart
  coordinates, not the embedded surface). Nothing semantic. Construction is
  now deterministic without hash order: encounter-order neighbor lists,
  first-encounter edge enumeration, face-slot refined-face numbering; the
  HalfEdge struct, the half_edges Dict builder, and both order-freezing CC
  rebuild helpers are deleted. The one remaining construction Dict is
  lookup-only (order never observed).
- **Fitted surfaces moved 1.6e-14 (cube) / 5.8e-15 (trefoil) pointwise** —
  pure arithmetic reordering, geometry intact.
- **physics_diff passed against the old baselines**: 0 side flips, worst
  deviations 8.2e-13 cube / 2.7e-10 trefoil (vs 1.34e-14 / 1.84e-10
  previously certified; 1e-9 budget). Baselines re-saved to restore
  headroom (deliberate change; the 1 unresolved trefoil ray was unresolved
  before too). 275/275 cold; knot-ray timing and allocation flat.
- Side benefit: physics_diff's claim that baselines survive Julia version
  changes is now actually true of the construction path — hash order was
  the remaining version dependence.

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
