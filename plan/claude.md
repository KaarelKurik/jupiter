# Claude's working plan

Companion to kaarel's journal (`plan/kaarel.md`); holds current state, the
roadmap, and standing conventions. Dated session records live alongside as
`plan/YYYY-MM-DD.md`. Last updated 2026-07-13 (second session). **This file
is the resume point** — read it plus `jj log` before touching code.

## Where we are (2026-07-13)

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
Gallery: first_light, textured trefoil, cube first flythrough.

## Next candidates, kaarel picks the order

- **GPU port, steps 4–5** (steps 1–3 done 2026-07-13b; decision that session:
  CPU-first groundwork aiming at mechanical device translation; tiering is a
  chain — GPU certified against production, production against reference,
  reference the sole oracle). Step 4: tabulated Γ developed on CPU against
  reference (port the winning algorithm, not the current loop — see the
  roadmap bullet below). Step 5: the real port — device throat view
  (packed_polys still Vector{Array}, needs one padded/offset array), wavefront
  restructure of the passage loop (kernels can't heap-allocate the per-passage
  sum-type boxes), Float32 policy measured against Float64-on-device.
  `gpu/` subenv holds CUDA.jl; `scripts/gpu_smoke.jl` is the compile check.
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
- **Representation alternatives** (research fork): Prautzsch freeform splines /
  Reif TURBS / Peters–Karčiauskas guided splines give C^k with finite
  polynomial patches (no transcendental blends — cheaper derivatives,
  GPU-native, fiddlier construction); YZ is the canonical manifold-based
  choice, and valence-4-only meshes nearly degenerate its transcendental part
  anyway. More radical: **level-set throat** — fit a smooth implicit φ and
  define g = A(φ)(I − n̂n̂ᵀ) + B(φ)n̂n̂ᵀ in one global ambient chart; no atlas,
  no transitions, no collar caustics; a different but arguably more natural
  wormhole; current machinery would remain as the physics reference.

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
