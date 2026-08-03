# Claude's working plan

Companion to kaarel's journal (`plan/kaarel.md`); holds current state, the
open decisions, and standing conventions. Dated session records live
alongside as `plan/YYYY-MM-DD.md`; findings with their why-chains live in
repo-root `measurements.md`; the code map is repo-root `atlas.md`. Last
updated 2026-08-01. **This file is the resume point** — read it plus `jj log`
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

Certification state: 782 tests cold green (451 + the 2026-08-01 pullback
and directional-jet batteries); `scripts/physics_diff.jl` passes its
baselines (0 side flips; baselines re-saved 2026-08-01 at the fused-accel
moment with kaarel's blessing, deviations read 0);
`scripts/precision_diff.jl` is the standing F32
acceptance instrument and reproduces its accepted shape. The oracle tiering
holds: reference.jl is the sole semantic ground truth; production is
certified against it; the jet christoffel against its AD twin
(`christoffel_ad`); the device kernel runs production `sweep_ray` unmodified
(85% of F64 rays bit-identical end-to-end, rest transcendental-intrinsic
ulps). F32 is accepted end-to-end by kaarel's decomposition metric
(measurements.md 2026-07-14c, re-confirmed post-jets 2026-07-16).

Perf state (measurements.md 2026-07-16 → 2026-07-19c): hand-rolled jet
christoffels are 12x CPU / 21–34x device over the AD path they replaced;
the production Γ path is AD-free closed-form (settles shader portability).
**Production integrates with fixed-step RK4**: every production RayBudget
is the 3-arg constructor (tolerance = 0), so the dopri/error-control
branch executes only in tests — discovered 2026-07-19c; it inverts item
3(e)'s framing and corrects the 19b census (per-attempt traffic is 4 RK4
stages of the same jet tower, not 7 dopri stages; the shared-tower
conclusions and all landed wins survive). Device F32 at 768×576: wall
0.40 s = kernel 0.21 s + host straggler tail + ~0.19 s pool-init gap
(lever c, now the probe workload's biggest single item); local 14.1
KB/thread. The straggler tail (a ~33 ms/round latency floor, 37% of
kernel for <1% of rays) is handed to host threads since 2026-07-19c
(DeviceSweep defaults handoff=1024, threads=128), which also flips the
small-frame conclusion: device F32 beats CPU at 192×144 (0.051 vs 0.078
wall, raymap level). Standing diagnosis (ncu, measurements.md 2026-08-01c
— supersedes 07-19b's "bandwidth-bound on local traffic"): the kernel is
**latency-bound at register-pinned occupancy** — 255 regs/thread admits
only 8 of the SM's 48 warp slots (12% achieved occupancy), 61% of warp
wait is spill-load scoreboard, 20.5/32 lanes active under divergence, and
no memory unit is near saturation (DRAM 35%, L2 49% at production size).
This mechanically explains the two kernel-NEUTRAL fusions: FLOP/lane cuts
feed pipes that idle anyway; the byte-traffic census priced a
bandwidth-bound kernel we don't have. Kernel-side FLOP/lane micro-fusion
is CLOSED as a lever. The corner-loop transient set (pd tuples, wirtinger
intermediates) remains the culprit, recast: it is the *live width* that
pins registers at the cap and forces the spill chains. The fused
Gauss-form accel stays for the CPU win, the simpler meaning-carrier, and
as item-1 substrate. Device
benchmark protocol: kernel numbers don't reproduce across sessions (clock
state) — A/B same-day, interleaved, ≥2 full-size warm passes, min-of-reps;
wall is the honest metric (host tail time is deliberately outside
kernel_seconds); workload-sensitive levers benchmark against flight
segments, not just the static camera (the mid-crossing regime is a single
all-throat round — no tail, latency levers are no-ops there).

Flight/video state: flyvideo.jl is the authoring surface — adaptive/uniform
pacing, scene=/sky1=/sky2= knobs, unresolved-pixel budget escalation
(`unres=`/`esc=`), `exit_placement()` (author the embedding frame, never
the sky asset), and `gpu=1` (device F32 frame raymaps; coast/probes/
escalation stay CPU F64, so flight paths are bit-identical to CPU runs and
escalated pixels land F64; ~2.3x at 768×576, and since the 2026-07-19c tail
handoff the device also wins small frames at the raymap level — 0.051 vs
0.078 wall at 192×144; flyvideo-level confirmation rides the next render). Real
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
   the cheap ordinal until then). This is the item that owns the expensive
   rays (kaarel, 2026-07-13): they are limit-cycle rays on the boundary
   between the two ambient images — the unresolved pixels — not grazing
   hits, and that same boundary is where fly-through chaotic shimmer lives. Backend-agnostic (CPU video renders drop
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
   measurements.md 2026-07-18 → 2026-07-19c): wall decomposition, the
   ABI-stack diagnosis and @inline tree collapse (−31% kernel / −15% CPU),
   occupancy ruled out, Γ-contraction fusion (reference geodesic_accel
   carries the meaning; physics_diff re-baselined), the byte-weighted
   traffic census, leaf-inline completion (−3.5% kernel / −3–5% CPU,
   bit-identical), the CPU tail handoff (item b, served 2026-07-19c:
   static-probe wall −20%, small frames flip to device, mid-crossing
   correctly neutral; DeviceSweep defaults handoff=1024/threads=128), and
   the fixed-RK4 discovery (production never runs dopri — see perf state).
   Mapped constraints: structural-layer inlining OOMs at 127–188
   KB/thread; sweep_ray's 7.9K depot is footprint-not-traffic (retired);
   remaining stage-tower spill is live-width-forced, so only narrower
   objects — not call restructuring — cut the dominant traffic. **(d) deep
   fusion SERVED 2026-08-01** (measurements.md 2026-08-01 + 01b): the EL
   contraction collapses to the Gauss formula w = Jᵀ(D_vD_v c) (reference
   pullback_accel carries it), the 8-lane DJet tower + fused accel landed
   certified (782 cold, physics_diff in-band unsaved, precision/parity in
   character) — kernel NEUTRAL, CPU −3%: the spill gate turned out to be
   the corner-loop transients, not output lanes; kept for the CPU win,
   the simpler meaning-carrier, and as item-1 substrate. kernel_sweep now
   threads env (the seam reaches the device). A/B protocol refinement:
   in-process wall comparisons need ABBA ordering (position artifact
   ±5–10%; kernel numbers are position-robust). The 2026-08-01c ncu
   profile (first direct look; permission setup is boot-scoped, see
   measurements method note) re-ranked the remaining levers for the
   offline-throughput endpoint (kaarel's priority, 2026-08-01):
   **Next: (c), in its cross-frame form.**
   - (c) **cross-frame pipeline** (generalizes the old pool-init-overlap
     item; single-frame overlap is its degenerate case). The ncu profile
     shows the device idle or near-idle across round-2-class launches
     (29 blocks, DRAM 6%), inter-round host gaps, and the ~0.19 s
     pool-init gap — all fillable with *other frames'* independent work
     in offline rendering (per-task streams; host entry solves for frame
     N+1 under device rounds of frame N). Honest expectation ~1.5–2x at
     the flyvideo level; co-resident blocks share the register file, so
     concurrency fills tails/gaps but cannot raise warps/SM during
     round-1 saturation. nsys timelines are the instrument. (An
     in-kernel entry solve remains an alternative for the single-frame/
     realtime case; the F64 host solve was a design convenience, never a
     requirement.)
   - (f) **divergence / sweep granularity probe** (new, 2026-08-01c):
     20.5/32 lanes active with bin=true — early-exiting rays idle until
     their sweep window closes. Up to ~1.5x of issued work recoverable;
     knobs are sweep size vs round overhead and compaction cadence.
     Cheap A/B probes before any restructuring.
   - (g) **rematerialization probe** (post-wrap addendum, 2026-08-01c
     discussion): trade spilled *holds* for *recomputation* in the
     corner loop — the pebbling space-time tradeoff, distinct from the
     ruled-out register caps (those convert holds to spills, not
     recomputes; the allocator never leaves that corner and can't — FP
     algebraic rewrites are off-limits to it under IEEE semantics, so
     narrowing liveness is a source-level move). The profile makes the
     exchange rate unusually favorable: ALUs at 20% while spill-load
     scoreboard is 61% of warp wait — recompute is paid in an idle
     currency. Open question is whether the corner-loop DAG has
     profitable pebbling moves at all; probe at the source level,
     certify bit-identical (recomputation of the same expressions is
     bit-stable).
   - (e) **Integrator economy — reframed 2026-07-19c**: production is
     *already* the cheap fixed-step RK4 integrator (h = 0.05 constant, no
     h cap, no error controller — those live in the test-only dopri
     branch). The morning's slack analysis (h-cap binding, wasted order-5
     stages) described the non-executing branch and is void. The live
     question inverts: would dopri-with-tolerance buy *accuracy* worth
     its cost anywhere (chaotic band? tightened budgets?), or is fixed
     RK4 at chart-scale h simply correct here? The Richardson
     h-refinement probe — still the instrument that closes the standing
     no-converged-ground-truth open item — now judges both at once:
     quantify fixed-RK4's discretization error at production h, and
     dopri's marginal value over it, by the decomposition metric. The
     dopri arm is fair now: k1 reuse (reject cache + FSAL, 2026-07-19c)
     is in both twin loops, bit-identical, footprint-only cost.
   Frame batching for offline video unchanged (latency-hostile in
   realtime); persistent threads demoted (the tail is latency-bound and
   now handled at the seam). Tools: nsight-compute + nsight-systems
   installed (2026-07-19) — ncu stall-reason counters to confirm the
   local-latency story directly, nsys timelines for (c).

Unordered, as wanted:

- **Representation fork** (own section below): the live design decision.
  Not on any critical path since direct F32 christoffel beat tabulation
  accuracy in-kernel (2026-07-14b); remains the CPU-perf / elegance
  question, on kaarel's clock — but re-coupled to GPU perf 2026-08-01c
  as the sole lever on kernel occupancy (see the fork section). Item 2
  above leans on it only softly.
- **Silhouette tightening**: tessellation-only first-hit gives a polygonal
  wormhole outline (spurious rim misses). `Mouth` is an abstract interface
  ready for a second strategy (outward-offset tessellation, or a
  conservative bound + Newton). Purely an outline-fidelity item, with no
  bearing on unresolved pixels or shimmer (kaarel, 2026-08-03): what the
  silhouette resolves is grazing rays, and on convex silhouettes those are
  barely deflected however they are resolved — measured, δ < 0.001° until
  the ray's deepest point reaches ~25% of cylinder_depth (measurements.md
  2026-07-11). The expensive rays live elsewhere; see item 1.
- **David mesh path**: res/models/*.stl are triangles; fit_geometry assumes
  quads, so one CC pre-subdivision or a quad remesh comes first; also STL
  loading. (`scripts/david.jl` scaffolding is in kaarel's working copy.)
- **Camera/authoring vision** (kaarel, longer-term): smooth curve authoring
  should be real-time flying controls — either a GPU build or a lighter
  schematic view with realtime control; not a near-term design driver, but
  see the reordering trigger in item 3.
- **Retire served bit-identity pins / prune the certification surface**
  (agreed 2026-07-17; broadened 2026-08-01 by kaarel's growing unease at
  the code and test volume): bit-identity is a landing instrument, not a
  permanent semantics; once a restructure's equivalence is proven, either
  merge the twins (a shared body keeps identity by construction —
  sweep_ray ↔ trace_geodesic's duplicated integrator loops are the
  standing case) or downgrade the pinned equality test to the physics_diff
  tolerance band. Tied to it: decide whether the recursive renderer is a
  production path or a certification twin, and label it accordingly. The
  2026-08-01 additions are candidates once the fused path has settled: the
  328 DJet op-level contraction tests were landing scaffolding (a thin
  surface_djet-vs-surface_jet + accel-vs-oracle set may suffice), and the
  accel now has three bodies (fused / gradient oracle / AD twin) where
  two tiers might do.
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
CPU-perf / representation-elegance question, on kaarel's clock.
**Re-coupled to GPU perf 2026-08-01c**: the ncu profile pins kernel
occupancy at 8 of 48 warps/SM via 255 regs/thread, and the register
demand is the corner-blend live width — option 2's affine transitions
eliminate corner blending entirely, the only identified path to
≤170 regs (3 blocks/SM, +50% warps) or ≤128 (4 blocks, 2x). The fork is
now the sole lever on kernel occupancy, though not on any near-term
critical path:

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
