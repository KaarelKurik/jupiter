# Claude's working plan

Companion to kaarel's journal (`plan/kaarel.md`); holds current state, the
roadmap, and standing conventions. Dated session records live alongside as
`plan/YYYY-MM-DD.md`. Last updated 2026-07-12. **This file is the resume
point** — read it plus `jj log` before touching code.

## Where we are (2026-07-11)

The POC pipeline is complete end-to-end: fitted YZ surface → blending → throat
metric → christoffels → geodesics with chart/half/mouth transitions → BVH mouth
entry → renderer with two-pass raymap output. First image lives in
`gallery/first_light.png`; the trefoil knot renders end-to-end via
`scripts/trefoil.jl` (heavily scrambled side-2 sky, as a knotted throat
should). 209 tests green via `julia --project -e 'using Pkg; Pkg.test()'`,
including reference-equivalence against the `Reference` submodule. The hot
loop has been through two perf passes (7.5x, then 3.1x from StaticArrays);
adaptive DP5(4) tracing exists but is opt-in (`RayBudget.tolerance` — see
2026-07-11 log for why fixed-step stays the render default).

## Next candidates, kaarel picks the order

- **Silhouette tightening**: tessellation-only first-hit gives a polygonal
  wormhole outline (spurious rim misses). `Mouth` is an abstract interface
  ready for a second strategy (outward-offset tessellation, or a conservative
  bound + Newton).
- **Textured sky** from an image file (raymap makes iterating on this ~0.6s).
- **David mesh path**: res/models/*.stl are triangles; fit_geometry assumes
  quads, so one CC pre-subdivision or a quad remesh comes first; also STL
  loading. (`scripts/david.jl` scaffolding is in kaarel's working copy.)
- **Camera transport / fly-through**: camera = point + arbitrary frame (NOT
  necessarily orthonormal — no global metric with non-isometric placements;
  loop holonomy may rescale/shear, and should). Needs a parallel-transport
  companion to geodesic_step; ray emission from inside via SituatedPhase.
- **Perf**: see roadmap below. Ground truth for aggressive optimization
  (kaarel's requirement): deliberate semantic changes go to `src/reference.jl`
  first, then production follows; run `scripts/physics_diff.jl` around any
  optimization.
- Multi-mouth ambient spaces need nearest-entry selection across mouths
  (enter_mouth would need to expose hit distance).

## Performance roadmap (2026-07-11 discussion)

Known remaining costs: ~1 MB allocated per ray (Placement in to_ambient and
mouth tessellation, BVH traversal stack, enter_mouth Newton on plain arrays);
`metric` evaluates the surface jacobian twice (once inside collar's dual pass,
once in inner_metric) — sharing them entangles the clean outer/inner split,
kaarel's call; `surface_normal_out`'s nested pass inside collar is ~45% of
trace time. In rough order of leverage:

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
  `directional(f, x, v)` plus `jacobian_columns(f, x, Val(N))` in wew.jl —
  swap the backend there if ForwardDiff ever disappoints.

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
