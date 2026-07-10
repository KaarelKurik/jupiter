# Claude's working plan

Companion to kaarel's `plan`; tracks review findings and the proposed implementation
sequence. Last updated 2026-07-10. **This file is the resume point** — read it plus
`jj log` before touching code.

## Where we are (2026-07-10)

The POC pipeline is complete end-to-end: fitted YZ surface → blending → throat
metric → christoffels → geodesics with chart/half/mouth transitions → BVH mouth
entry → renderer with two-pass raymap output. First image lives in
`gallery/first_light.png` (blue side-1 sky, orange side-2 sky visible through the
throat with concentric lensing rings; 18 unresolved pixels sit on the critical
ring, as limit-cycle theory predicts). 144 tests green via
`julia --project -e 'using Pkg; Pkg.test()'`.

Next candidates, kaarel picks the order:
- **Silhouette tightening**: tessellation-only first-hit gives a polygonal
  wormhole outline (spurious rim misses). `Mouth` is an abstract interface
  ready for a second strategy (outward-offset tessellation, or a conservative
  bound + Newton).
- **Textured sky** from an image file (raymap makes iterating on this ~0.6s).
- **David mesh path**: res/models/*.stl are triangles; fit_geometry assumes
  quads, so one CC pre-subdivision or a quad remesh comes first; also STL
  loading and probably per-valence caching of `yz_fitting_matrix`.
- **Camera transport / fly-through**: camera = point + arbitrary frame (NOT
  necessarily orthonormal — no global metric with non-isometric placements;
  loop holonomy may rescale/shear, and should). Needs a parallel-transport
  companion to geodesic_step; ray emission from inside via SituatedPhase.
- **Perf**: render is threaded (63s for 192×144 on 16 threads, byte-identical
  to serial); big wins left in christoffel (ForwardDiff.jacobian chunking) and
  precomputing limit positions / fitting matrices.

Working conventions: jj (never bare git mutations); one described change per
step, `jj new` at seams, describe-as-intent up front; Claude does the jj
bookkeeping at phase boundaries. `res/` tracked inputs, `out/` ignored outputs,
`gallery/` curated milestone PNGs. Raymaps stay ignored (cache semantics).

## Current architecture (post-redesign)

- `Surface` = mesh + chart_polys (three polynomials per vertex, in vertex order,
  each fitted relative to the vertex's canonical handle).
- `HalfThroat` wraps a `Surface`; `Chart` is a thin handle: (half_throat,
  half_edge_handle). Thick data lives on `Surface`, avoiding Julia's
  circular-type-definition awkwardness.
- Canonical frame convention: a vertex's induced neighbor is the first entry in its
  neighbor list (recorded in meshy.jl); `half_edge_offsets` numbers outgoing
  half-edges ccw from that canonical handle.
- Chart transitions re-rotate through the shared-edge frame on both sides using
  `half_edge_offset`, landing in the neighbor's canonical frame. **Verified correct**
  (2026-07-08), as is the `half_edge_offsets` construction itself and the agreement
  between fit-time wedge numbering and runtime offsets. This resolves the
  reference-frame design gap from the 2026-07-04 review.
- `wedge_map` reimplemented via `fake_complex_pow`/`safe_atan2` (principal complex
  power from AD-friendly real primitives) after TaylorDiff experiments;
  `reference_wedge_map` kept as a test oracle. Caveats: NaN at the exact chart
  origin; `safe_atan2` diverges near the negative real axis (points at/beyond the
  far vertex along the extended edge) — outside the intended transition domain, but
  worth a guard or comment.

## Open items

Observation from smoke testing: chart polys under-shoot the CC limit position at the
chart center by ~0.015 on the unit-ish cube — expected, the YZ fit is least-squares,
not interpolatory. If exact vertex interpolation turns out to matter, add an
interpolation constraint at the center to the fit.

## AD backend (migrated 2026-07-08)

Migrated TaylorDiff → ForwardDiff after the nested-use bug harvest (details and fix
ideas preserved in `taylordiff-bugs.md`; bugs 1 and 2 there are upstream-worthy).
The AD contact surface is now a single function, `directional(f, x, v)` in wew.jl —
swap the backend there if ForwardDiff ever disappoints. Nested duals (christoffel)
work out of the box via ForwardDiff's tag system; no patches needed. Full 14-check
battery (scratchpad battery.jl — worth promoting to test/) passes identically to
the TaylorDiff figures: christoffel-vs-FD 4.7e-9, tensoriality 2.2e-15, round
trips ~3e-16. Perf note for the renderer hot loop: `directional` is 3 calls per
Jacobian; ForwardDiff.jacobian with chunking would do it in one pass.

For rendering speed later: cache `yz_fitting_matrix` per valence and precompute
limit positions over the mesh (existing comments in wew.jl note both).

## Implementation sequence

1. ~~**Wiring fix pass**~~ Done 2026-07-08, smoke-tested: `fit_geometry(cubemesh())`
   end-to-end; poly centers vs limit positions; `wedge_map` ≡ `reference_wedge_map`
   to 6 digits across valences; chart-transition round trips at machine precision on
   pos *and* vel (so TaylorDiff differentiates `fake_complex_pow` fine).
2. ~~**YZ blending**~~ Done 2026-07-08. `surface(env, chart, uv)` blends the 4
   corner charts of the containing face, with per-wedge weights
   `blend_scalar(u)·blend_scalar(v)` in conformal unit-square coordinates
   (`wedge_square_coords` / `square_coords_to_chart`; corner frames related by
   (s,t)→(t,1−s)). `blend_scalar` is a swappable placeholder
   (`H(1−x)/(H(x)+H(1−x))`, H = e^{−1/x}): contract is f(0)=1, f(1)=0, C^∞-flat at
   BOTH ends (0: wedge seams + centers; 1: chart support boundary); its
   f(x)+f(1−x)=1 property makes weights sum to exactly 1 (blend divides by the sum
   anyway, so any conforming f works). Smoke-tested: same face point via all 4
   corner charts agrees to 5e-16; weights sum to 1; no C0/C1 kink crossing a mesh
   edge; blend ≈ chart poly at center; TaylorDiff differentiates through the blend
   (`surface_normal_out` works). Wedge selection branches on `primal(x)`
   (TaylorScalar-aware) since `floor` fails on TaylorScalar.
   Caveats: NaN at the exact chart center (safe_atan2(0,0)) — guard when it bites;
   2-arg `atan` on TaylorScalar silently RETURNS ONLY THE PRIMAL (derivative info
   dropped) — never use it in a differentiated path; cubemesh winds ccw-from-inside,
   so `surface_normal_out` currently points INTO the cube — decide the orientation
   convention before wiring up `collar`'s depth sign.
3. ~~**Sampling + export**~~ Done 2026-07-08 except the Blender eyeball (kaarel's).
   `sample_surface(half_throat, samples_per_edge)` samples each face on a grid
   (boundary vertices duplicated — raw (verts, faces), not a Mesh, since the
   duplicates make it non-manifold); `save_obj` in meshy.jl. `res/yz_cube.obj`
   exported (12×12 per face, NaN-free, bbox ±0.85). Numeric smoothness check:
   finite differences up to 4th order along an arc crossing a mesh edge are ≤ a
   no-seam baseline arc at every order — no discontinuity signature, C3 survives.
   Fixes en route: `fake_complex_pow` guarded at the origin (grid corners are chart
   centers); `safe_atan2` got the complementary half-angle form for x<0 (roundoff
   at face-corner samples put corner coords at (−6.6e-17, 0), i.e. 0/0 on the old
   form's branch cut, with blend weight 1).
4. ~~**Metric layer**~~ Done 2026-07-08. Vocabulary settled with kaarel:
   - `Throat` = Surface + ThroatParams + two Placements; `HalfThroat` is now a
     handle (throat, side). Params are per-Throat: cross_scale, depth_scale,
     cylinder_depth (metric fully cylindrical beyond it), transition_depth
     (half-to-half handover; constructor asserts ≥ cylinder_depth).
   - Orientation bookkeeping lives in Placements; the internal half-to-half map is
     natural (identity in (u,v), reversal in d). Same-orientation placements ⇒
     traversal reverses orientation, by design. Differing placement scales in a
     shared ambient space ⇒ no global metric (intended; test path is same-scale).
   - `depth_interpolate` reuses `blend_scalar(d/cylinder_depth)`: C^∞, exactly
     outer at d ≤ 0, exactly cylindrical at d ≥ cylinder_depth (verified exact).
   - `christoffel` works via NESTED TaylorDiff (TaylorScalar of TaylorScalar); the
     earlier symbolic-differentiate idea was wrong (blend weights are
     position-dependent). Required: seeds built with `basis_direction`
     (eltype-matched — TaylorDiff's make_seed demands it), recursive `primal`, and
     three small method patches on TaylorDiff.TaylorScalar constructors for nested
     use (one genuine TaylorDiff bug: the generic constructor routes nested
     constants through the truncate/extend path — candidate for upstreaming).
   - Verified: metric SPD + symmetric across the band; Γ symmetric in ij (exact);
     Γ matches central finite differences of the metric to 5e-9; metric is C0
     across chart transitions in the tensorial sense (|Jᵀ g_far J − g_near| ≈ 3e-15
     against entries of size 2.5).
5. ~~**Geodesic integration + transitions**~~ Done 2026-07-10.
   - `geodesic_step` (RK4 on SituatedPhase, fixed chart; acceleration =
     `wvel_along_v(v, v)`), `settle_phase` (re-homes after each step: lateral
     `chart_transition` hop when the face square-coords leave [0,1/2]², natural
     half-to-half `half_transition` past transition_depth — identity in (u,v),
     d ↦ 2·td − d, vel_d flipped; closed on phases by design), `exits_mouth`
     predicate + `to_ambient` (placement ∘ collar; exact, since the metric is
     the flat collar pullback for d ≤ 0), `trace_geodesic` driver returning
     AmbientRay or a still-inside SituatedPhase — the one place the sum type
     appears, as the honest codomain (kaarel raised type-stability/reasoning
     concerns; resolved by confining the union here).
   - Verified (in test/runtests.jl "geodesics"): energy g(v,v) conserved to 1e-6
     across lateral hops and a full side-1 → side-2 traversal; settling changes
     description but not the physical point; half_transition is an involution
     with exact energy; ambient exit speed² = conserved energy (isometry);
     time-reversal over a hop-crossing arc closes to ~1e-11.
   - Note: "settled" ≠ arbitrary — a phase at square-coords > 1/2 belongs to a
     neighbor chart; compare trajectories only between settled phases.
   - Observed physics: laterally-launched geodesics near the mouth get pulled
     into the throat — plausible for this geometry, and energy/reversibility
     checks pass, so not treated as a bug.
   - Deliberately deferred to the renderer: the *inverse* mouth problem
     (ambient ray hits the mouth surface → initial SituatedPhase), i.e.
     ambient-side ray/collar intersection and chart lookup.
6. **Renderer**: camera in one flat patch, integrate geodesics until they exit to
   flat space on either side, hit-test against a simple environment (textured sky
   sphere per side).
   - *Inverse mouth problem* done 2026-07-10: `Mouth(env, half_throat,
     samples_per_edge)` tessellates the d=0 surface in ambient coords with chart
     provenance per triangle (face handle + square-coord corners); `enter_mouth`
     takes the nearest Möller–Trumbore hit as a guess, Newton-refines in
     (s, t, ray-parameter) against the exact blended surface (~1e-13), and pulls
     the ray direction back through the inverse collar Jacobian (exact — collar
     is an isometry onto the flat outside), returning a settled SituatedPhase or
     nothing on a miss. Verified: entry inverts `to_ambient` to machine
     precision on pos and vel; flat speed² = chart-metric energy; full pipeline
     ambient-in side 1 → ambient-out side 2 passes in the suite.
   - *Generic entry interface + BVH* done 2026-07-10 (kaarel's requirements:
     David mesh is O(1e5) tris, and other entry strategies may be needed for
     visual fidelity — tessellation-only first-hit gives spurious silhouette
     misses, i.e. a polygonal wormhole outline). `Mouth` is now an abstract type
     with contract enter_mouth(env, mouth, origin, dir) → SituatedPhase |
     nothing, plus half_throat(mouth); `TessellatedMouth <: Mouth` carries a
     median-split BVH (`build_bvh`/`nearest_mouth_hit`, verified identical to
     brute force). A tighter entry strategy (e.g. outward-offset tessellation +
     Newton, or sphere-tracing a conservative bound) can slot in later.
   - *Renderer assembly* done 2026-07-10 in src/render.jl: `Scene` (throat, one
     mouth per side — each side its own universe for now), `RayBudget`
     (step_size, max_steps per passage, max_passages) as the user-facing knob
     for limit-cycle rays (all wormholes admit them), `trace_ray` (enter → trace
     → re-enter loop; escaping AmbientRay or nothing when budget runs out),
     pinhole `look_at_camera`/`camera_ray`, `checker_sky` (30° checkers, blue
     side 1 / orange side 2), `render` image loop, `save_ppm`.
   - *Raymap two-pass rendering* done 2026-07-10 (kaarel's request): `RayMap`
     stores per-pixel exit (side, pos, vel; side 0 = unresolved) so skies swap
     without retracing. `render_raymap` (expensive pass, Threads.@threads over
     rows) → `shade(raymap, sky)` (cheap pass); `render` = the composition;
     `save_raymap`/`load_raymap` via Serialization. Exit *positions* stored too,
     so finite-distance environments stay possible later.
   - *Camera stance settled with kaarel*: a camera is a point + an arbitrary
     frame (basis), NOT necessarily orthonormal — with non-isometric placements
     in one ambient space there is no global metric, only a connection, so loop
     transport may rescale/shear a returning camera; images should honestly
     show that. Orthonormality is just look_at_camera's initial condition; the
     pinhole formula is indifferent (and homogeneous in overall frame scale).
     Documented on the Camera struct; transport code itself is future work
     (parallel-transport companion to geodesic_step for the fly-through shot).
   - Remaining/deferred: multi-mouth ambient spaces need nearest-entry selection
     across mouths (enter_mouth would need to expose hit distance); textured sky
     from an image; David mesh needs a quad remesh or one CC pre-subdivision
     (fit_geometry assumes quads). Perf notes: cache `yz_fitting_matrix` per
     valence; precompute limit positions; ForwardDiff.jacobian chunking in
     christoffel.

## Testing infrastructure

Add `test/` with invariant checks as pieces stabilize: CC subdivision against known
cube values (done manually vs Blender), limit stencil weights, fitting-matrix shape
consistency (1 + 12·valence points vs monomial count), `half_edge_offsets` vs
repeated `ccw`, wedge-map equivalence with the complex reference, transition
round-trips, flat-mesh Christoffels.

## Resolved (for the record)

- 2026-07-04 review: mechanical bugs in `half_edge_ccw`, `chart_transition`,
  `christoffel` tensor assembly (index permutation re-verified), `inner_metric`
  fundamental form + eltype promotion, normal normalization — all fixed.
- Reference-frame design gap: resolved by `half_edge_offsets` + two-sided rotation
  in `chart_transition` (verified 2026-07-08).
