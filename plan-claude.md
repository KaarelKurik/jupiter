# Claude's working plan

Companion to kaarel's `plan`; tracks review findings and the proposed implementation
sequence. Last updated 2026-07-08 after the Chart-as-handle redesign.

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

**Critical path**: blending done, so the metric → christoffel → geodesic stack is
unblocked; next up is sampling/export (step 3) for visual + C3 verification, or
straight to the metric layer (step 4). For rendering speed later: cache
`yz_fitting_matrix` per valence and precompute limit positions over the mesh
(existing comments in wew.jl note both).

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
   - Remaining: camera model, per-ambient-space mouth set (a ray must be tested
     against every mouth placed in its ambient space), sky-sphere hit-test +
     texture, image loop. Perf notes for then: `enter_mouth` is brute force over
     triangles (fine at 768 tris; add a bbox/BVH if it hurts); cache
     `yz_fitting_matrix` per valence; precompute limit positions; consider
     ForwardDiff.jacobian chunking in christoffel.

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
