# Measurements log

Findings worth not re-deriving, but which don't change the working plan.
Probe scripts are deliberately disposable (they'd be rewritten against changed
code anyway); enough method detail lives here to reconstruct them.

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

**Conclusions.**

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
