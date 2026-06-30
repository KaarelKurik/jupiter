# Claude's working plan

Companion to kaarel's `plan`; this tracks the review findings and the proposed
implementation sequence from the 2026-07-04 review session.

## Review findings

### Fixed (by kaarel, 2026-07-05)
- `half_edge_ccw` loop advanced `half_edge` instead of `cur`
- `chart_transition`: undefined `chart_valence`; `newpos` dropped the depth coordinate
- `neighbor_chart`: passed the function `half_throat` instead of `half_throat(chart)`
- `christoffel`: derivatives now `stack`ed into a 3×3×3 tensor; index permutation in the
  `@tensor` expression re-derived and verified correct for derivative-index-last layout
- `inner_metric`: bracketing of the Jacobian; first fundamental form `sf_jac' * sf_jac`
- `surface_normal_out` now normalized

### Still open
- `surface` (wew.jl): `u`, `v` not in scope (fix: `@polyvar u v` at use site — see notes
  below), and `chart.poly` is a `Vector` of 3 polynomials, not callable — evaluate
  componentwise
- `vertex_index(chart::Chart)` (wew.jl): `HalfEdge` has no `vertex_index` field
  (fields are `vertices`, `face_index`); for a `HalfEdgeHandle` it's `name[1]`.
  Depends on deciding whether `Chart.half_edge` is a `HalfEdge` or a handle.
- `chart_by_half_edge` constructs `Chart` with 2 args; `Chart` has 3 fields (`poly`
  missing). Tied to the chart-identity design decision below.
- `inner_metric`: `out = zeros(3,3)` is Float64-rigid; MethodErrors once `christoffel`
  pushes TaylorDiff duals through `pos`. Fix: promote the eltype from the values, e.g.
  `zeros(promote_type(eltype(g), typeof(params.depth_scale)), 3, 3)`.
- `inner_metric_params` and `depth_interpolate` don't exist yet. `depth_interpolate`
  should be at least C1 in depth so the Christoffels come out C0 (smoothstep-family
  interpolant gives margin).
- `wedge_map` + TaylorDiff: complex powers of dual-carrying complex numbers are
  unlikely to differentiate cleanly. The map is holomorphic, so the pushforward is
  available in closed form: d(a3)/d(a0) = (source_n/target_n) · a1/a0 · a3/a2 · (−1),
  applied to velocity as complex multiplication. Sidesteps AD entirely.
- Naming note (not a bug): `wedge_index_and_angle`'s "index" is a float and can be
  negative for α < 0; the remainder is still correct. Rename so nobody indexes a
  vector with it, or normalize with `mod(floor(Int, …), n_wedges)` if it will ever
  be used as an index.

### Design gap: chart identity and reference frames
`fit_geometry` returns polys in vertex order but discards the reference handle each
chart was fitted against. Chart coordinates are only meaningful relative to that edge.
`chart_transition` needs the neighbor's poly expressed relative to the shared edge
`twin(ce)`, which generally differs from the neighbor's stored reference edge.

Decision: store per-vertex `(poly, reference_handle)`; `neighbor_chart` composes the
wedge map with a rotation by 2πk/n in chart coordinates to land in the neighbor's
stored frame (k = wedge offset between shared edge and reference edge). Cheap, avoids
storing n rotated copies of every poly.

Related: `Chart` currently bundles `half_throat + half_edge + poly` — decide whether a
Chart *is* the per-vertex fitted object (reference handle fixed at fit time) or a view
relative to an arbitrary edge. The rotation-composition design suggests the former.

## Implementation sequence

1. **Finish the fix pass** (items above), then smoke test in the REPL:
   `fit_geometry(cubemesh())`, evaluate a chart on a grid, check the center value
   against `limit_position`.
2. **YZ blending**: per-vertex bump weights in chart coordinates; evaluate a point
   inside a face as the weighted blend of its 4 corner charts. Needs face-local
   (s,t) → corner-chart-coordinate maps (the z^{4/n} maps already exist) plus
   partition-of-unity normalization.
3. **Sampling + export**: sample the blended surface over each base face, write OBJ
   (MeshIO already a dep), eyeball in Blender. Numeric continuity check: sample along
   a curve crossing an edge/vertex, finite-difference up to 3rd derivatives to verify
   the C3 claim survives the implementation.
4. **Metric layer**: `inner_metric_params`, `depth_interpolate`, dual-number
   genericity throughout. Consider symbolic `TypedPolynomials.differentiate` for the
   surface Jacobian instead of nested AD (charts are polynomials — exact derivatives
   are free, and it avoids TaylorDiff-inside-TaylorDiff in `christoffel`).
   Unit test: flat quad grid ⇒ Christoffels ≈ 0 everywhere.
5. **Geodesic integration + transitions**: RK4 on (pos, vel); chart-switch trigger via
   wedge angle/depth thresholds. Round-trip test: transition to neighbor and back ≈
   identity on pos and vel; geodesic crossing a chart boundary stays continuous.
6. **Renderer**: camera in one flat patch, integrate geodesics until they exit to flat
   space on either side, hit-test against a simple environment (textured sky sphere
   per side).

## Testing infrastructure
Add `test/` with invariant checks as pieces stabilize: CC subdivision against known
cube values (done manually vs Blender), limit stencil weights, fitting-matrix shape
consistency (1 + 12·valence points vs monomial count), transition round-trips,
flat-mesh Christoffels.
