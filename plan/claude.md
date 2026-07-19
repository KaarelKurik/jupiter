# Claude's working plan

Companion to kaarel's journal (`plan/kaarel.md`); holds current state, the
open decisions, and standing conventions. Dated session records live
alongside as `plan/YYYY-MM-DD.md`; findings with their why-chains live in
repo-root `measurements.md`; the code map is repo-root `atlas.md`. Last
updated 2026-07-19. **This file is the resume point** — read it plus `jj log`
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
its baselines (0 side flips; baselines re-saved 2026-07-19 at the
Γ-fusion re-baseline moment, deviations currently read 0);
`scripts/precision_diff.jl` is the standing F32
acceptance instrument and reproduces its accepted shape. The oracle tiering
holds: reference.jl is the sole semantic ground truth; production is
certified against it; the jet christoffel against its AD twin
(`christoffel_ad`); the device kernel runs production `sweep_ray` unmodified
(85% of F64 rays bit-identical end-to-end, rest transcendental-intrinsic
ulps). F32 is accepted end-to-end by kaarel's decomposition metric
(measurements.md 2026-07-14c, re-confirmed post-jets 2026-07-16).

Perf state (measurements.md 2026-07-16 → 2026-07-19b): hand-rolled jet
christoffels are 12x CPU / 21–34x device over the AD path they replaced;
the production Γ path is AD-free closed-form (settles shader portability);
device F32 ~2x 16-thread CPU at 768×576 (kernel ~0.36 s, local 14.0
KB/thread). Standing diagnosis: the kernel is bandwidth-bound on its own
local-memory traffic, NOT occupancy-limited (register caps always lose);
~90% of per-attempt local traffic is flow6's *interior* — the collapsed
jet tower's live width exceeds the 255-register cap and spills; call
structure is no longer a factor (07-18 @inline collapse + 07-19 leaf
completion removed it), and FLOP cuts alone don't move the kernel
(Γ-contraction fusion: CPU 2–4%, device neutral). Kernel time is two
regimes: a throughput-bound first round and a latency-bound straggler tail
(46% of kernel for <1% of rays, the budget-capped filament pixels).
Device benchmark protocol: kernel numbers don't reproduce across sessions
(clock state) — A/B same-day, interleaved, ≥2 full-size warm passes,
min-of-reps; workload-sensitive levers benchmark against flight segments,
not the static camera.

Flight/video state: flyvideo.jl is the authoring surface — adaptive/uniform
pacing, scene=/sky1=/sky2= knobs, unresolved-pixel budget escalation
(`unres=`/`esc=`), `exit_placement()` (author the embedding frame, never
the sky asset), and `gpu=1` (device F32 frame raymaps; coast/probes/
escalation stay CPU F64, so flight paths are bit-identical to CPU runs and
escalated pixels land F64; ~2.3x at 768×576, CPU wins small frames). Real
skies: res/textures/ skybox[12].jpg + space[12].jpg (render-ready PPMs
regenerated per checkout via magick — see flyvideo.jl header); frames land
in out/frames/. Flagship: gallery/trefoil_flythrough_space.mp4 (768×576,
701 frames). The honest constant-speed cube video (uniform Δτ ≈ 0.0037)
still just hasn't been wanted.

Gallery: first_light, textured trefoil, cube first flythrough, trefoil
flythrough (nature skies low-res, and the space flythrough at 768×576).
Oversized gallery files enter via `jj file track --include-ignored <path>`
(kaarel's preference: explicit track, never raise the auto-track guard).

## Next candidates

The perf/quality path is ordered (agreed with kaarel, 2026-07-17b; a priori
cost analysis in that session's log). Headline of the analysis: the physics
doesn't price out realtime — ~2 orders of magnitude of pure execution
headroom on the device path (spill, wave tails, host entry), another 1–2
algorithmic orders in smooth regions of the exit map, and a hard core
(first hit + chaotic filament integration) that fits interactive
budgets everywhere except mid-crossing frames, where budget-as-frame-
deadline degradation (semantics RayBudget already has) is available and
perceptually honest — the filament is where the image already aliases.
(The core is "irreducible" only up to Lyapunov-limited depth: the
section-cache idea under item 2 attacks it directly.)

Item 3 runs ahead of 1–2 (kaarel's reordering, 2026-07-18): bundles and
caches are set aside to see how far execution alone goes. Execution
engineering helps *both* the realtime-authoring endpoint and offline
renders (flyvideo gpu=1).

1. **Jacobi-field ray bundles (DNGR-style).** Trace sparse rays +
   geodesic deviation κ; interpolate where the exit map is smooth,
   subdivide near critical rings — 10–100x fewer rays and principled
   antialiasing. The keystone item: κ is simultaneously the shimmer answer
   (the filament is image-space aliasing), the trust criterion the caches
   below need, and the endgame F32→F64 fallback flag (passage count stays
   the cheap ordinal until then). Backend-agnostic (CPU video renders drop
   from hours toward minutes) and representation-agnostic (survives the
   fork) — highest information value, no dependencies. Design notes from
   the 2026-07-18 DNGR discussion (kaarel checked the paper: DNGR bundles
   are one-per-pixel for filtering/quality; the sparse-ray reading is
   ours, at a scale whose linearization validity must be *checked*): trust
   is pairwise Hermite consistency between neighboring bundles plus
   discrete agreement (side, passage count) — never pointwise κ alone,
   since folds can hide between samples; land the architecture with
   finite-difference bundles (2 offset rays, no new machinery) before
   extending jets to order-4 for the true Jacobi equation (κ needs
   accuracy only in smooth regions, where FD is fine; in the chaotic band
   it only needs to be large, to veto interpolation); the speedup claim is
   conditional — 10–100x on smooth-dominated frames, ~break-even-but-
   antialiased at peak crossing, frame budget as the floor.
2. **Scattering-map caches — second, pulled by demand** (they pay when
   geometry is static across many frames, i.e. the authoring loop; softly
   gated on the representation fork since tables bake the current charts).
   Mouth cache first: entry × direction → exit ray per mouth (enter_mouth
   would expose hit distance), chaotic bands flagged for true integration.
   In-throat adaptation (worked out 2026-07-17b): the cache *composes* — an
   in-throat camera ray needs honest integration only for its first leg (a
   2D family per frame; after first exit it's ambient flight + re-entries,
   which the mouth cache covers), and the first leg's winding tail lives in
   the cylinder region, where the product structure (metric d-free ⇒
   surface geodesic × linear depth, ḋ constant) reduces caching to the 2D
   surface-geodesic flow: a 4D table shared across rays, frames, and both
   halves, singular exactly on the ḋ≈0 trapped set (closed surface
   geodesics), flagged like the critical rings. Composition compounds
   interpolation error, so table cells carry the κ trust criterion from
   item 1. (Absorbs the former standalone cylinder-product CPU lever.)
   **Upgrade (kaarel's loop-cut idea, 2026-07-18): index the in-throat
   table by Poincaré sections instead of time steps.** Sections = (a few
   well-chosen loop cuts of the mesh) × (d-interval), chosen so every
   trapped geodesic crosses the family infinitely often (a complete
   system of sections, Birkhoff's construction; sufficient condition:
   each complement patch of the loop family contains no complete surface
   geodesic — patches inside convex normal balls guarantee it; empirical
   probe: seed long surface geodesics, measure max sojourn between
   crossings). This makes the ḋ≈0 trapped set *interior* to the table
   domain (the return map doesn't care that d isn't progressing — the
   singularity was an artifact of time parametrization), and the product
   structure collapses the table to a 2D-domain return map
   (s,θ) → (s′,θ′, next-loop-id, arc length σ), with d′ = d + ḋσ/√(1−ḋ²)
   exact and ḋ constant. Per-crossing error can be pushed below
   per-winding RK4 error at trivial memory cost, and the Lyapunov
   amplification per winding hits integration and table-iteration
   identically — so where the table wins the per-winding ε comparison,
   the cache is *more* accurate than integration, not a lossy shortcut.
   Caveats: composed expansion is still the depth limit (store the return
   map's Jacobian per cell — the item-1 bundle data — as the trust flag);
   tangency boundaries (grazing crossings) flag cells for integration
   fallback; valid only where the metric is d-free, which is exactly
   where the winding tail lives. Attacks the cost analysis's "irreducible
   core" directly — softens it to Lyapunov-limited depth.
3. **GPU execution engineering — ACTIVE.** Done so far (why-chains in
   measurements.md 2026-07-18 → 2026-07-19b): wall decomposition, the
   ABI-stack diagnosis and @inline tree collapse (−31% kernel / −15% CPU),
   occupancy ruled out, Γ-contraction fusion (reference geodesic_accel
   carries the meaning; physics_diff re-baselined), the byte-weighted
   traffic census, leaf-inline completion (−3.5% kernel / −3–5% CPU,
   bit-identical). Mapped constraints: structural-layer inlining OOMs at
   127–188 KB/thread; sweep_ray's 7.9K depot is footprint-not-traffic
   (retired as a lever); remaining flow6 spill is live-width-forced, so
   only narrower objects — not call restructuring — cut the dominant
   traffic. **Next moves agreed 2026-07-19: (b) first, then (d).**
   - (b) **CPU tail handoff**: the straggler tail is 46% of kernel for
     <1% of rays at ~0.5 ms serial latency per attempt. After round ~2,
     ship the few hundred stragglers to host threads via the sweep_stage!
     seam (backend-per-round is trivial there; a 400-attempt straggler is
     ~1 ms serial on host vs ~200 ms on device). Budget-as-deadline stays
     the realtime degradation semantics; threads=128 is a free ~6%
     (candidate DeviceSweep default).
   - (d) **Deep fusion** (kaarel's fuller reading of the opaque-block
     doctrine, promoted by the census): fuse chart logic + interpolation
     + contraction *algebraically*. With v constant per RHS evaluation,
     the acceleration is the Euler–Lagrange form a = g⁻¹(½∇(vᵀgv) −
     D_v(gv)) — all derivatives of contracted objects (with g = (∂c)ᵀ(∂c):
     vᵀgv = |D_v c|², gv = (∂c)ᵀD_v c) — so a directional-jet tower
     (v-contracted lanes from eval_packed_partials up, instead of all 10
     partials) shrinks exactly the live width that is the dominant local
     traffic, plus a real CPU FLOP cut. Cost: a second jet algebra beside
     Jet3 through the whole chain (packed eval → wedge/cpow → wirtinger →
     blend → metric); production goes genuinely opaque (doctrine covers
     it; reference states the EL form as the meaning-carrier).
   - (c) behind those: **pool-init overlap** — the wall-kernel gap is the
     initial 442k F64 entry solves + device throat build, pipelinable
     against round 1 (or an in-kernel entry solve; the F64 host solve was
     a design convenience, never a requirement).
   Frame batching for offline video unchanged (latency-hostile in
   realtime); persistent threads demoted (the tail is latency-bound, not
   wave-quantized). Tools ask pending: `pacman -S nsight-compute
   nsight-systems` (kaarel to install) — ncu counters and nsys timelines
   sharpen (b) and (c).

Unordered, as wanted:

- **Representation fork** (own section below): the live design decision.
  Not on any critical path since direct F32 christoffel beat tabulation
  accuracy in-kernel (2026-07-14b); remains the CPU-perf / elegance
  question, on kaarel's clock. Item 2 above leans on it only softly.
- **Silhouette tightening**: tessellation-only first-hit gives a polygonal
  wormhole outline (spurious rim misses). `Mouth` is an abstract interface
  ready for a second strategy (outward-offset tessellation, or a
  conservative bound + Newton). Related but distinct (kaarel, 2026-07-13):
  the *expensive* rays are limit-cycle rays on the boundary between the two
  ambient images — the unresolved pixels — not grazing hits; that boundary
  is also where fly-through chaotic shimmer lives (image-space aliasing →
  the ray-bundles item).
- **David mesh path**: res/models/*.stl are triangles; fit_geometry assumes
  quads, so one CC pre-subdivision or a quad remesh comes first; also STL
  loading. (`scripts/david.jl` scaffolding is in kaarel's working copy.)
- **Camera/authoring vision** (kaarel, longer-term): smooth curve authoring
  should be real-time flying controls — either a GPU build or a lighter
  schematic view with realtime control; not a near-term design driver, but
  see the reordering trigger in item 3.
- **Retire served bit-identity pins** (agreed 2026-07-17): bit-identity is
  a landing instrument, not a permanent semantics; once a restructure's
  equivalence is proven, either merge the twins (a shared body keeps
  identity by construction — sweep_ray ↔ trace_geodesic's duplicated
  integrator loops are the standing case) or downgrade the pinned equality
  test to the physics_diff tolerance band. Tied to it: decide whether the
  recursive renderer is a production path or a certification twin, and
  label it accordingly.
- **Tabulated Γ** (probed 2026-07-14, demoted — accuracy is
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
