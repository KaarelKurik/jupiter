# Claude's working plan

Companion to kaarel's journal (`plan/kaarel.md`); holds current state, the
open decisions, and standing conventions. Dated session records live
alongside as `plan/YYYY-MM-DD.md`; findings with their why-chains live in
repo-root `measurements.md`; the code map is repo-root `atlas.md`. Last
updated 2026-07-17. **This file is the resume point** — read it plus `jj log`
before touching code. This file describes the present; the history of how we
got here lives in the session logs and jj descriptions, not here.

## Where we are

The POC pipeline is complete end-to-end, on CPU and GPU: fitted YZ surface →
C^∞ corner blending → collar/throat metric → closed-form christoffels
(order-3 surface jets, AD retained as in-tree oracle) → geodesics with
chart/half/mouth transitions → BVH mouth entry (always an F64 Newton solve)
→ raymap renderers (recursive and wavefront drivers, bit-identical) →
textured equirect skies → parallel-transported camera frames → in-throat
cameras → piecewise-geodesic flight with mouth crossings → adaptive/uniform
fly-through video.

Certification state: 451 tests cold green; `scripts/physics_diff.jl` passes
its baselines (0 side flips, worst deviation 1.26e-10, inside the 1e-9
arithmetic-reordering band); `scripts/precision_diff.jl` is the standing F32
acceptance instrument and reproduces its accepted shape. The oracle tiering
holds: reference.jl is the sole semantic ground truth; production is
certified against it; the jet christoffel against its AD twin
(`christoffel_ad`); the device kernel runs production `sweep_ray` unmodified
(85% of F64 rays bit-identical end-to-end, rest transcendental-intrinsic
ulps). F32 is accepted end-to-end by kaarel's decomposition metric
(measurements.md 2026-07-14c, re-confirmed post-jets 2026-07-16).

Perf state (measurements.md 2026-07-16): the hand-rolled derivative lever
paid out on both platforms — CPU christoffel 12.3x, whole renders 11.5–12.5x
(cube 192×144 raymap ~0.10–0.15 s at 16 threads), device kernel 21.5x F64 /
34x F32. Device F32 reaches 16-thread-CPU parity at 384×288 and is **2.1x
ahead at 768×576**: resolution alone feeds the card, so no GPU throughput
lever remains on the critical path. The production Γ path is AD-free
closed-form, which also settles shader portability.

Latest flight work (2026-07-16 second session, closed out post hoc): the
first trefoil fly-through — flyvideo.jl grew scene=/pacing knobs; 701
uniform frames at 192×144 in 6 minutes (pre-jets this was hour-class);
out/trefoil_flythrough.mp4, gallery promotion still kaarel's call. The
honest constant-speed cube video (uniform Δτ ≈ 0.0037, ~1670 frames)
stopped being deferred-for-cost; it just hasn't been wanted yet.

Gallery: first_light, textured trefoil, cube first flythrough.

## Next candidates, kaarel picks the order

- **Representation fork** (own section below): the live design decision.
  Not on any critical path since direct F32 christoffel beat tabulation
  accuracy in-kernel (2026-07-14b); remains the CPU-perf / elegance
  question, on kaarel's clock.
- **Silhouette tightening**: tessellation-only first-hit gives a polygonal
  wormhole outline (spurious rim misses). `Mouth` is an abstract interface
  ready for a second strategy (outward-offset tessellation, or a
  conservative bound + Newton). Related but distinct (kaarel, 2026-07-13):
  the *expensive* rays are limit-cycle rays on the boundary between the two
  ambient images — the unresolved pixels — not grazing hits; that boundary
  is also where fly-through chaotic shimmer lives (image-space aliasing →
  the ray-bundles item).
- **Jacobi-field adaptive ray bundles** (DNGR-style): trace sparse rays +
  geodesic deviation, interpolate where the exit map is smooth, subdivide
  near critical rings — 10–100x fewer rays and principled antialiasing.
  Same κ machinery is the endgame for the F32→F64 fallback flag (passage
  count is the cheap ordinal until then) and the shimmer answer.
- **David mesh path**: res/models/*.stl are triangles; fit_geometry assumes
  quads, so one CC pre-subdivision or a quad remesh comes first; also STL
  loading. (`scripts/david.jl` scaffolding is in kaarel's working copy.)
- **Camera/authoring vision** (kaarel, longer-term): smooth curve authoring
  should be real-time flying controls — either a GPU build or a lighter
  schematic view with realtime control; not a near-term design driver.
- **GPU optional refinements** (all off the critical path, wanted only if a
  realtime build happens): persistent threads pulling from the pool queue
  (fixes wave quantization + compaction tails in-frame), host/device
  overlap of entry solves (~20% of wall at 442k rays), frame batching for
  offline video only (flyvideo: cameras known ahead; latency-hostile in
  realtime). Provenance correction (kaarel, 2026-07-17): "mouth entry stays
  an F64 host solve" was a design convenience (2026-07-14c), not a kaarel
  requirement, and was never measured as necessary — an F32 or in-kernel
  entry solve is fair game if the host stage ever blocks full-device
  rendering. Budget-as-frame-deadline with progressive fallback on boundary
  filaments — semantics RayBudget already has.
- **Retire served bit-identity pins** (agreed 2026-07-17): bit-identity is
  a landing instrument, not a permanent semantics; once a restructure's
  equivalence is proven, either merge the twins (a shared body keeps
  identity by construction — sweep_ray ↔ trace_geodesic's duplicated
  integrator loops are the standing case) or downgrade the pinned equality
  test to the physics_diff tolerance band. Tied to it: decide whether the
  recursive renderer is a production path or a certification twin, and
  label it accordingly.
- **Remaining CPU levers**, if CPU cost ever binds again: cylinder-region
  product structure (for d ≥ cylinder_depth the metric is d-independent —
  integrate a 2D surface geodesic + linear depth motion; winding rays live
  there); mouth scattering-map caching (entry × direction → exit ray, per
  mouth — amortizes fly-throughs, chaotic bands flagged for true
  integration); tabulated Γ (probed 2026-07-14, demoted — accuracy is
  corner-blend-limited, gated on the representation fork).
- Multi-mouth ambient spaces need nearest-entry selection across mouths
  (enter_mouth would need to expose hit distance). This includes the
  both-mouths-of-one-throat-in-one-room case: conceptually the two ambient
  spaces need not be distinct, but the code currently identifies ambient
  space with side (Scene keys one mouth + the sky per side), so a co-placed
  second mouth is invisible to exiting rays (noted in atlas.md 2026-07-17).

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
   The depth blend stays exp-flat — no current reason to revisit (exact at
   runtime in tables, a cheap closed-form chain in jets), though "forever"
   was written when tabulation made it free; the same conditioning arguments
   would formally apply to it if the corner blend's class ever changes.
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
   Smoothness demand — downgraded to a probe question 2026-07-17: Γ needs
   the Hessian of φ, and the recorded requirement was φ ~C³; but the step-4
   data says conditioning at step scale, not smoothness class, predicted
   quality (exp-flat is C^∞ and plateaus at 1e-3; the C⁴ nonic reaches
   1e-5), and the standing acceptance instruments judge exit maps, not
   analytic class. So sub-C³ candidates get probed rather than gated out
   (kaarel skeptical it changes the answer, but probing is cheap).
   Candidates: (a) **superfrusta +
   metaballs** (kaarel 2026-07-14): ResFit decomposition of the throat into
   analytic 8-parameter primitives (arXiv 2512.09201, CVPR 2026), log-sum-exp
   smooth union (analytic — the union is not the weak point); open question
   is the primitive's own creases ("differentiable a.e.") vs the smoothness
   question above — confine or smooth them, or probe whether they matter. (b) polygon-soup IMLS over a dense CC-limit
   sampling: C^∞, local via BVH, no global solve, topology-faithful at small
   kernels; eval cost ~ current blend arithmetic. (c) compactly-supported
   HRBF: banded solve precompute, C^k, locality by kernel support.
4. **Orthogonal regardless of path**: per-chart step-size discipline
   (h=0.05 ≈ half a face on the trefoil drives the huge stage-excursion
   margins); wavefront/GPU structure survives every option.

## Current architecture

See `atlas.md` for the full map (math-to-code dictionary, a ray's life, the
twin/oracle structure, precision seams, per-file roles). Standing facts that
belong here:

- `Surface` = mesh + chart_polys (three polynomials per vertex, in vertex
  order, each fitted relative to the vertex's canonical handle) +
  packed_polys (dense Horner tables). `Throat` = Surface + ThroatParams +
  two Placements; `HalfThroat` is a handle (throat, side); `Chart` is a thin
  handle (half_throat, half_edge_handle). Thick data lives on `Surface`.
  The whole type layer is storage-parametric so an Adapt-ed device view is
  still the same types.
- Canonical frame convention: a vertex's induced neighbor is the first entry
  in its neighbor list (recorded in meshy.jl); `half_edge_offsets` numbers
  outgoing half-edges ccw from that canonical handle. Chart transitions
  re-rotate through the shared-edge frame on both sides using
  `half_edge_offset`. **Verified correct** (2026-07-08), as is the
  `half_edge_offsets` construction and the agreement between fit-time wedge
  numbering and runtime offsets.
- `wedge_map` via `fake_complex_pow`/`safe_atan2` (principal complex power
  from AD-friendly real primitives); `reference_wedge_map` is a test oracle.
  Caveats: NaN at the exact chart origin; `safe_atan2` diverges near the
  negative real axis (at/beyond the far vertex along the extended edge) —
  outside the intended transition domain, but worth a guard or comment.
- AD backend is ForwardDiff (migrated from TaylorDiff 2026-07-08;
  `taylordiff-bugs.md` preserves the harvest). The AD contact surface is
  `ad.jl`; production christoffel no longer uses AD (jets.jl), but AD stays
  the oracle layer (`christoffel_ad`, reference.jl) and the collar/metric
  path still runs dual passes.
- The `env` slot threaded through the geometry functions is deliberately
  unused: it's the dispatch extension point a recording/instrumenting env
  rides (confirmed useful 2026-07-14), and where a `christoffel(::ΓTables, v)`
  variant would hang.
- Standing inference hazard: nested/recursive hot code must keep self-call
  signatures exactly constant or Julia widens to Any and boxes everything;
  the dual-pass primitives in ad.jl keep separate bodies for the same
  reason (measurements.md 2026-07-12).

## Open items

- **No convergence-established ground truth exists** (kaarel, 2026-07-17):
  every instrument compares same-discretization twins — physics_diff runs
  reference and production at the same budgets, precision_diff's "truth"
  variant is F64 production at the production budget, and the F32 accuracy
  numbers are relative to F64 at the same h/tol. Right comparison for
  isolating precision (the ε-vs-κ factorization), but the discretization
  error itself (RK4 h=0.05 / dopri tol against converged geodesics) has
  never been quantified, so "truth"/"accuracy" in the record means "our
  best at the same discretization." Fix when wanted: a Richardson-style
  probe — tighten h/tol until exit-map quantiles stabilize (decomposition
  metric, not sup) — to bound discretization error at production budgets.
- Chart polys under-shoot the CC limit position at the chart center by ~0.015
  on the unit-ish cube — expected, the YZ fit is least-squares, not
  interpolatory. If exact vertex interpolation turns out to matter, add an
  interpolation constraint at the center to the fit.
- jd daemon left to test: explicit `--restart`/`--stop` exercises, and a warm
  run of the ray-bundle compare against the cold baseline (warm≡cold
  semantics spot check).
- The 14-check AD battery lived in a session scratchpad (2026-07-08) and is
  likely gone; if promoting it to test/ ever matters, reconstruct from the
  description in plan/2026-07-08.md rather than hunting for the file.
- physics_diff deviations are nonzero (≤1.26e-10) though passing; if a
  future optimization stacks more reordering on top, re-baselining then is
  the natural moment (2026-07-16 session note).

## Working conventions

jj (never bare git mutations); one described change per step, `jj new` at
seams, describe-as-intent up front; Claude does the jj bookkeeping at phase
boundaries. `res/` tracked inputs, `out/` ignored outputs, `gallery/` curated
milestone PNGs. Raymaps stay ignored (cache semantics). Findings go in
repo-root `measurements.md`; probe scripts are disposable (precision_diff
and physics_diff are standing instruments, not probes). Invariant tests
accumulate in `test/` as pieces stabilize; certification (tests /
physics_diff / final renders) always runs cold, day-to-day iteration goes
through `scripts/jd` (see `.claude/skills/julia-workflow`). Ground truth for
aggressive optimization (kaarel's requirement): deliberate semantic changes
go to `src/reference.jl` first, then production follows; run
`scripts/physics_diff.jl` around any optimization. `atlas.md` is the code
map — update it when file layout or the twin/oracle structure changes, not
for routine edits.
