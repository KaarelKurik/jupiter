# Claude's working plan

Companion to kaarel's journal (`plan/kaarel.md`); holds current state, the
roadmap, and standing conventions. Dated session records live alongside as
`plan/YYYY-MM-DD.md`. Last updated 2026-07-14 (fourth session). **This file
is the resume point** — read it plus `jj log` before touching code.

## Where we are (2026-07-14, fourth session)

The POC pipeline is complete end-to-end and now includes cameras *in and
through* the wormhole: fitted YZ surface → blending → throat metric →
christoffels → geodesics with chart/half/mouth transitions → BVH mouth entry →
two-pass raymap renderer → textured equirect skies → parallel-transported
camera frames (reference + production, 0 B/step) → in-throat cameras
(ambient/chart continuity certified) → piecewise-geodesic flight with
mouth crossings → adaptive fly-through video (839 frames through the cube;
temporal-sampling profile in measurements.md 2026-07-13). 275 tests green
cold, physics_diff bit-stable all session. Perf state: knot-hitting ray
~3.3-4.4 ms (collar single-pass −18%, 2026-07-12b), tracer allocates ~3-4
boxed sum-type returns per ray, step loops heap-free. The honest
constant-speed video (uniform Δτ=0.0037 → ~1670 frames) is deliberately
deferred: CPU render cost is the binding constraint now. GPU groundwork
steps 1–3 landed 2026-07-13 (second session): half-edge connectivity is flat
id-indexed arrays, the mouth is isbits with a threaded (stackless) BVH, and
the device smoke test passed — eval_packed bit-identical on the RTX 4070
SUPER, the ad.jl dual pass compiles clean at 31–76x one CPU thread
(measurements.md 2026-07-13b). No compilability wall in the AD core.
Step 4 (tabulated Γ) probed 2026-07-14: architecture settled (exact-in-d
pieces, seam-aligned sectors), 10–50x confirmed available, but achievable
accuracy hinges on the corner-blend choice; that widened into the full
**Representation fork** section below (C⁴ blend vs manifold splines vs
level-set) — kaarel picks before step-4 implementation resumes.
2026-07-14 second session (kaarel's call: try the simple road before
tabulating): the type layer went storage-parametric (bit-identical) and the
**production christoffel now runs unmodified in a device kernel** via Adapt —
no second-order-AD compilability wall. Float64 is spill-bound and loses to
the multicore CPU (5.5x one thread ≈ 1/3 of 16); **Float32 is the road**:
14.3x one CPU thread naive (~29x with depth binning; still spill-bound, so a
floor), and Γ accuracy vs F64 truth (median 9.5e-7, max 1.1e-5 sup-relative)
**beats every tabulation variant at exact blend semantics — the
representation fork no longer gates the GPU port** (measurements.md
2026-07-14b). Third session: eltype-honesty extended to the whole trace
path (steppers, dopri tableau, settle/transitions, mouth exit, emit_ray —
F64 certified bit-identical again) and **F32 judged end-to-end and
accepted** by kaarel's decomposition metric (measurements.md 2026-07-14c):
quantile curves parallel to the F64 1-ulp baseline at offset ≈ 2^28 ("just
rounding", no pathology bends), side flips 0/8000 cube and 1/8000 trefoil
concentrated at passage ≥ 3, exceedance-at-1-mrad 0.04%/1.1%, F32/F64
render pairs indistinguishable with the diff living exactly on the
two-image boundary filaments. `scripts/precision_diff.jl` is the standing
acceptance instrument (physics_diff's sibling); mouth entry stays an F64
solve, exits convert at to_ambient. Fourth session: the **wavefront
restructure landed** (src/wavefront.jl, measurements.md 2026-07-14d) — the
passage loop's sum-type boxes became a stage tag on one isbits
WavefrontRay{T} record; sweep_ray is the kernel-shaped body (0 B, both
integrators, resumable across sweeps), mouth entry the host-side F64 stage,
the driver sweep→compact→entry rounds with optional depth-branch binning.
Certified bit-identical to render_raymap (289 tests cold, physics_diff
0.0), and the staged CPU driver came out ~1.6x *faster* than the recursive
renderer (self-load-balancing active list vs row threading). Step 5's
remainder is now the device tracer kernel itself: sweep_ray + settle/
transition/to_ambient in-kernel on the Adapt-ed throat (gpu_christoffel.jl
scaffolding), host entry stage, then the register-pressure lever.
Gallery: first_light, textured trefoil, cube first flythrough.

## Next candidates, kaarel picks the order

- **GPU port, steps 4–5** (steps 1–3 done 2026-07-13b; decision that session:
  CPU-first groundwork aiming at mechanical device translation; tiering is a
  chain — GPU certified against production, production against reference,
  reference the sole oracle). Step 4: tabulated Γ developed on CPU against
  reference — probed 2026-07-14 (measurements.md; depth factors out exactly,
  seam-aligned sectors, 10–50x per eval *on the CPU* holds, accuracy
  corner-blend-limited: exp-flat ~1e-3, C⁴ ~1e-5 with a deliberate physics
  change) — **demoted from the GPU critical path 2026-07-14b**: kaarel had
  the simpler road tried first, and direct in-kernel AD christoffel at
  Float32 out-accuracies the tables (median 9.5e-7 sup-rel) at exact blend
  semantics and ≈ full-CPU speed naive. Tabulation remains a CPU-side
  option gated on the representation fork. Step 5: the rest of the real
  port — **eltype-honesty + F32 judgment DONE 2026-07-14c** (whole trace
  path carrier-honest, F64 bit-identical; F32 accepted end-to-end by the
  decomposition metric kaarel specified 14b — quantile parallelism vs the
  1-ulp baseline at ≈2^28, undecidable band ~1e-4 of pixels at passage ≥ 3
  on the boundary filaments, render pairs indistinguishable; instrument:
  scripts/precision_diff.jl, findings measurements.md 2026-07-14c; the
  F64-fallback flag ordinal until Jacobi-field κ lands is passage count,
  and the endgame remains Jacobi normalization — same infrastructure as
  ray bundles). **Wavefront restructure DONE 2026-07-14d** (src/wavefront.jl:
  isbits WavefrontRay{T} stage-tagged records, sweep_ray the box-free kernel
  body, host F64 entry stage, depth-branch binning plumbed; bit-identical to
  render_raymap, ~1.6x faster on CPU for free). Remaining: **the device
  tracer kernel** — sweep_ray + settle/transitions/to_ambient in-kernel on
  the Adapt-ed throat, entry solves staying host-side, wavefront_raymap
  growing a device driver (promote gpu_christoffel.jl's PackedTable/adapt
  rules to a gpu/ module when it lands). Register pressure: the named
  lever is **hand-rolled derivatives in production, AD stays in reference
  as the oracle** (kaarel 2026-07-14b) — the 27KB/thread F32 spill is dual
  bloat (48 carrier floats per scalar, mostly structural zeros); Γ needs
  s..∂³s through polynomial∘wedge∘blend, all closed-form (poly derivatives
  are packable polys, wedge is powers of z, blend/depth are scalar 1D
  chains), certified against reference AD exactly like every production
  optimization. Doubles as the portability road: realtime-flying endpoints
  likely mean shader stacks with no Julia AD, and 255-reg kernels hurt on
  every card. Subsumes "split the fused dual passes" surgery. `gpu/`
  subenv holds CUDA.jl (+Adapt); `scripts/gpu_smoke.jl` is the compile
  check, `scripts/gpu_christoffel.jl` the christoffel bench (its
  PackedTable/adapt rules promote to a gpu/ module when the tracer kernel
  lands).
- **Silhouette tightening**: tessellation-only first-hit gives a polygonal
  wormhole outline (spurious rim misses). `Mouth` is an abstract interface
  ready for a second strategy (outward-offset tessellation, or a conservative
  bound + Newton). Related but distinct (kaarel, 2026-07-13): the *expensive*
  rays are limit-cycle rays on the boundary between the two ambient images —
  the unresolved pixels — not grazing hits; that boundary is also where the
  fly-through's chaotic shimmer lives (image-space aliasing → the ray-bundles
  item).
- ~~Textured sky~~ done 2026-07-12: `load_ppm` (P3/P6; `magick x.png x.ppm`
  for real images), `sample_equirect` (bilinear, azimuth-wrap), `TexturedSky`
  per-side callable; res/textures graticule test skies
  (scripts/make_test_skies.jl regenerates), scripts/textured.jl re-shades a
  raymap (~0.2s at 384×288).
- **David mesh path**: res/models/*.stl are triangles; fit_geometry assumes
  quads, so one CC pre-subdivision or a quad remesh comes first; also STL
  loading. (`scripts/david.jl` scaffolding is in kaarel's working copy.)
- **Camera transport / fly-through**: DONE 2026-07-12 end-to-end — reference
  impl + production port (transport_flow/step joint RK4, map_frame riding
  transitions, settle_transport, trace_transport, emit_ray; 0 B/step), camera
  inside the throat (enter_transport = exact inverse of to_ambient via
  mouth_entry; metric_normalize for geometry-pegged budgets, emit_ray stays
  raw for effort semantics; SituatedCamera render_raymap; ambient/chart
  continuity certified in tests), and piecewise-geodesic flight (flight.jl:
  FlyingCamera, coast crossing mouths both directions with continuous
  parameter, frame-relative steer/maneuver, keyframe_raymap;
  scripts/flythrough.jl demo — cube traversal under textured skies, incl. the
  look-back-from-side-2 frame). Camera = point + arbitrary frame (NOT
  necessarily orthonormal; loop holonomy may rescale/shear, and should).
  Longer-term vision (kaarel): smooth curve authoring should be real-time
  flying controls — either a GPU build or a lighter schematic view with
  realtime control; not a near-term design driver.
- **Perf**: see roadmap below. Ground truth for aggressive optimization
  (kaarel's requirement): deliberate semantic changes go to `src/reference.jl`
  first, then production follows; run `scripts/physics_diff.jl` around any
  optimization.
- Multi-mouth ambient spaces need nearest-entry selection across mouths
  (enter_mouth would need to expose hit distance).

## Representation fork (recorded 2026-07-14)

The step-4 probes (measurements.md 2026-07-14) reduced tabulated-Γ accuracy
to one question — the corner blend's smoothness class — and localized the
hard points at irregular (valence ≠ 4) vertices. That reframes the
representation choices; all plausible paths in one place, kaarel to pick.
**Scope narrowed 2026-07-14b**: the fork no longer gates the GPU port
(direct F32 AD christoffel beats the tables' accuracy in-kernel at
acceptable speed — measurements.md 2026-07-14b); it remains live as the
CPU-perf / representation-elegance question, on kaarel's clock:

1. **Stay YZ + C⁴ corner blend** (kaarel's lean *if* tabulation on the
   current representation proceeds). Nonic smootherstep corner blend: sector
   fits reach ~1e-5 relative Γ at N=12–24; AD fallback discs cover the
   irregular-vertex cores, where the polynomial blend composed with z^(n/4)
   leaves a genuine fractional-power singularity (exp-flat is the mirror
   image: flat-at-vertices, bad at seams, plateaus ~1e-3; C² worse than
   both). Deliberate physics change: re-fit, re-baseline, gallery re-render.
   The depth blend stays exp-flat forever (exact at runtime in the tables).
2. **Manifold splines / affine atlas** (steps away from YZ; discussed
   2026-07-14). Key mapping: their "extraordinary point" (affine structure
   obstructed, cone defect ≠ 0) = our valence ≠ 4 vertex = exactly the
   non-converging tabulation cores; the trefoil (knotted torus, χ=0,
   valence-4-only) is the zero-extraordinary-point case and is precisely
   where our probes showed everything converging. With affine transitions,
   polynomials are transition-closed ⇒ **corner blending dies entirely**
   (one global piecewise-polynomial spline, no partition of unity); Γ per
   patch is closed-form rational-with-sqrt — tabulation likely unnecessary.
   Design axis: EP count vs severity — Ricci-flow single-EP construction
   (Gu–He–Qin, CAD 2008) concentrates the whole 4π genus-0 defect in one
   savage point; mild π/2 cones at the 8 cube corners with rigid transitions
   elsewhere + spline caps (Reif TURBS / Peters–Karčiauskas guided splines /
   Prautzsch) is likely the numerically kinder variant. Migration cost:
   fit.jl + wedge machinery + transitions replaced; steppers/transport/mouth
   survive; reference migrates with it (certification across the change is
   geometric, not baseline-diff).
3. **Level-set throat** (kaarel's preferred endpoint if a good rep existed):
   fit smooth implicit φ, g = A(φ)(I − n̂n̂ᵀ) + B(φ)n̂n̂ᵀ in one global ambient
   chart; no atlas, no transitions, no collar caustics. Requirements from
   kaarel's past failed searches: cheap eval, cheap precompute, good
   locality, respects mesh topology (no metaball-style merging/pathology).
   Hard constraint from the physics: Γ needs the Hessian of φ, so φ must be
   ~C³ for well-behaved geodesics. Candidates: (a) **superfrusta +
   metaballs** (kaarel 2026-07-14): ResFit decomposition of the throat into
   analytic 8-parameter primitives (arXiv 2512.09201, CVPR 2026), log-sum-exp
   smooth union (analytic — the union is not the weak point); open question
   is the primitive's own creases ("differentiable a.e.") vs the C³ need —
   confine or smooth them. (b) polygon-soup IMLS over a dense CC-limit
   sampling: C^∞, local via BVH, no global solve, topology-faithful at small
   kernels; eval cost ~ current blend arithmetic. (c) compactly-supported
   HRBF: banded solve precompute, C^k, locality by kernel support.
4. **Orthogonal regardless of path**: per-chart step-size discipline
   (h=0.05 ≈ half a face on the trefoil drives the huge stage-excursion
   margins); wavefront/GPU port structure (steps 1–3) survives every option.

## Performance roadmap (2026-07-11 discussion)

Known remaining costs (allocation is done — see 2026-07-12 log): `metric`
evaluates the surface jacobian twice (once inside collar's dual pass, once in
inner_metric) — sharing them entangles the clean outer/inner split, kaarel's
call. (Collar's *internal* double evaluation is gone: 2026-07-12,
`value_and_jacobian_columns` fuses value + normal into one dual pass, −18%
per ray, and retired the `nested_jacobian_columns` twin along the way.)
Standing hazard worth remembering: nested/recursive hot code must keep
self-call signatures exactly constant or Julia's inference widens to Any and
boxes everything (measurements.md 2026-07-12); the dual-pass primitives in
ad.jl must keep separate bodies for the same reason. In rough order of
leverage:

- **Cylinder-region product structure**: for d ≥ cylinder_depth the metric is
  d-independent — integrate a 2D surface geodesic + linear depth motion; the
  expensive winding rays live there.
- **Chebyshev-tabulated Γ per face×depth box** (~8³ coeffs/component):
  geometric convergence for a C^∞ metric, likely 10–50x per step; cross-face
  continuity then only to approximation tolerance (~1e-12) — certify against
  Reference + physics_diff rather than assume. The single biggest
  single-machine lever.
- **Mouth scattering-map caching** (entry point × direction → exit ray, per
  mouth): amortizes fly-throughs; chaotic bands flagged for true integration.
- **Jacobi-field adaptive ray bundles** (DNGR-style): trace sparse rays +
  geodesic deviation, interpolate pixels where the exit map is smooth,
  subdivide near critical rings — 10–100x fewer rays and principled
  antialiasing.
- **GPU port**: the post-StaticArrays core is isbits/kernel-shaped; the real
  port cost is flattening the Dict-based half-edge/offset lookups to index
  arrays, concretizing Mouth, and the BVH stack. Float64 is 1/32–1/64 rate on
  consumer cards — decide Float32 policy (fine away from chaotic bands) with
  Float64 fallback for flagged pixels. Wavefront-style ray compaction for
  divergence. Port the winning algorithm (likely tabulated Γ), not the
  current loop.
- **Representation alternatives** (research fork): see the dedicated
  "Representation fork" section below (consolidated 2026-07-14 after the
  step-4 tabulation probes made the irregular-vertex/blend problem concrete).

## Current architecture (post-redesign)

- `Surface` = mesh + chart_polys (three polynomials per vertex, in vertex
  order, each fitted relative to the vertex's canonical handle).
- `Throat` = Surface + ThroatParams + two Placements; `HalfThroat` is a
  handle (throat, side); `Chart` is a thin handle: (half_throat,
  half_edge_handle). Thick data lives on `Surface`, avoiding Julia's
  circular-type-definition awkwardness.
- Canonical frame convention: a vertex's induced neighbor is the first entry
  in its neighbor list (recorded in meshy.jl); `half_edge_offsets` numbers
  outgoing half-edges ccw from that canonical handle.
- Chart transitions re-rotate through the shared-edge frame on both sides
  using `half_edge_offset`, landing in the neighbor's canonical frame.
  **Verified correct** (2026-07-08), as is the `half_edge_offsets`
  construction itself and the agreement between fit-time wedge numbering and
  runtime offsets.
- `wedge_map` reimplemented via `fake_complex_pow`/`safe_atan2` (principal
  complex power from AD-friendly real primitives); `reference_wedge_map` kept
  as a test oracle. Caveats: NaN at the exact chart origin; `safe_atan2`
  diverges near the negative real axis (points at/beyond the far vertex along
  the extended edge) — outside the intended transition domain, but worth a
  guard or comment.
- AD backend is ForwardDiff (migrated from TaylorDiff 2026-07-08;
  `taylordiff-bugs.md` preserves the harvest). The AD contact surface is
  `ad.jl` (directional, the two dual-pass primitives, situate) — swap the
  backend there if ForwardDiff ever disappoints.
- src layout since 2026-07-12: meshy.jl (half-edge mesh), cc.jl
  (Catmull-Clark), throat.jl (types + accessors), ad.jl (AD contact surface),
  chart.jl (wedge geometry, blending, surface eval, chart transitions),
  geodesic.jl (collar/metric/christoffel, steppers, tracers, transport),
  mouth.jl (Mouth/BVH/enter_mouth/enter_transport), fit.jl (YZ fitting),
  render.jl (cameras/raymaps/skies), flight.jl (piecewise-geodesic camera
  flight), reference.jl (formerly one wew.jl).

## Open items

- Chart polys under-shoot the CC limit position at the chart center by ~0.015
  on the unit-ish cube — expected, the YZ fit is least-squares, not
  interpolatory. If exact vertex interpolation turns out to matter, add an
  interpolation constraint at the center to the fit.
- jd daemon left to test: explicit `--restart`/`--stop` exercises, and a warm
  run of the ray-bundle compare against the cold baseline (warm≡cold
  semantics spot check).
- The 14-check AD battery (scratchpad battery.jl) is worth promoting to
  test/.

## Working conventions

jj (never bare git mutations); one described change per step, `jj new` at
seams, describe-as-intent up front; Claude does the jj bookkeeping at phase
boundaries. `res/` tracked inputs, `out/` ignored outputs, `gallery/` curated
milestone PNGs. Raymaps stay ignored (cache semantics). Findings go in
repo-root `measurements.md`; probe scripts are disposable. Invariant tests
accumulate in `test/` as pieces stabilize; certification (tests /
physics_diff / final renders) always runs cold, day-to-day iteration goes
through `scripts/jd` (see `.claude/skills/julia-workflow`).
