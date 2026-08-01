# Measurements log

Findings worth not re-deriving, but which don't change the working plan.
Probe scripts are deliberately disposable (they'd be rewritten against changed
code anyway); enough method detail lives here to reconstruct them.

## 2026-08-01 — deep-fusion opener: the EL contraction collapses further to the Gauss formula; the fused tower needs 8 lanes, not 9–10

**The identity.** Working item 3(d)'s EL form through the collar structure
before building: with g = JᵀJ (J = ∂c, the collar jacobian — outer metric),
the contraction w_u = v^i v^j (∂_i g_uj − ½ ∂_u g_ij) collapses past the
plan's sketch. Using ∂_i e_u = ∂_u e_i (mixed partials of c),
v^iv^j ∂_i g_uj = (∂_u D_v c)·(D_v c) + e_u·(D_v D_v c) while
½ ∂_u(vᵀgv) = (∂_u D_v c)·(D_v c) — the first-derivative-of-J terms cancel
exactly, leaving

    w = Jᵀ (D_v D_v c),    a = −(JᵀJ)⁻¹ w

— the Gauss formula: geodesic acceleration of a flat pullback is the
tangential projection of the ambient second derivative. Pointwise algebra on
g = JᵀJ, no immersion assumption, so it holds at any d where outer_metric is
in play (including as the g_o arm inside the blend).

**Probe** (disposable, method): cube throat with non-unit ThroatParams
(1.3, 0.9, 0.5, 0.75), charts 1/3/6 × 8 random (uv, vel) × depths spanning
all three branches; reference-AD implementations of (1) the EL restatement
a = g⁻¹(½∇(vᵀgv) − D_v(gv)) everywhere, (2) the Gauss form at d ≤ 0,
(3) the blend-region decomposition w = ω·w_o + (1−ω)·w_i +
(D_vω)(g_o v − g_i v) − ½(ω′/cd)(vᵀg_o v − vᵀg_i v)·ê_d with w_o the Gauss
form and w_i the inner contraction. Worst relative deviation vs
Reference.geodesic_accel: EL 1.4e-13, Gauss 6.2e-15, blend 7.3e-15.

**What it does to the lane count.** The plan's "v-contracted lanes instead
of all 10 partials" over-promised and its 9-lane refinement under-promised;
the true requirement per surface component is 8 lanes: full Jet2
{f, fu, fv, fuu, fuv, fvv} (J needs n̂, nu, nv fully, hence full seconds —
that floor is real: g must invert) + two *doubly*-contracted third lanes
{D_wD_w su-lane, D_wD_w sv-lane}, since thirds enter only through
D_v D_v c = D_wD_w s − d·D_wD_w n̂ − 2 v_d·D_w n̂ (w = uv-part of vel; d
enters exactly as always). The (g, dgu, dgv, dgd) sym-matrix layer of
outer_metric_gradient disappears from the geodesic path entirely; the inner
arm needs no third derivatives at all (g_i is d-free and built from firsts,
so its contracted gradient closes on seconds); the blend arm is the linear
combination above, whose extra ingredients (g_o v, g_i v, vᵀg_o v, vᵀg_i v,
ω, ω′) all come from objects already in hand. Landed in reference.jl as
`pullback_accel` (the meaning-carrier for the fusion) + 3 outer-region
equivalence tests; 454 cold green.

**Calibrated expectation for the build** (recorded before measuring, per
protocol): tower width 10 → 8 lanes (−20%) at the fat layers (pd, wirtinger
outputs, accumulators, jdiv), with the census's spill leverage — traffic
rides the excess over the 255-register cap, not total width — amplifying
that by roughly ×(width/excess); plus the FLOP cut is concentrated in the
most expensive lanes (order-3 compose1/leibniz/jdiv terms, wirtinger
D30/D21). The census's "~90% of traffic is the tower" bounds the win from
above; this is not a 90% lever. Honest guess: kernel −15–30%, CPU −10–20%;
wall moves half the kernel's move on the static probe (kernel is ~0.21 of
0.40 s post-handoff).

## 2026-07-19c — production has always run fixed-RK4 (dopri is test-only); the tail handoff lands: static-probe wall −20%, small frames flip to the device; k1 reuse is free but optimizes a dead branch

**The discovery.** The k1-reuse rider went in first (dopri_step gains a k1
argument and returns the FSAL k7; both twin loops carry k1 across rejects
and reuse k7 on no-hop accepts, the no-hop test being componentwise scalar
egal — whole-struct `===` lowers to a memcmp runtime call the device
compiler rejects, today's portability lesson). Certified bit-identical:
451 cold green, physics_diff exactly 0/0, parity counts in character. Then
every A/B read dead neutral — device 0.366/0.375/0.368/0.366, CPU
recursive 0.1177/0.1175/0.1174 — and the why is the finding: **every
production RayBudget is the 3-arg constructor, tolerance = 0.0, so
trace_ray and sweep_ray take the fixed-step RK4 branch everywhere;
dopri executes only in tests** (the 1e-8 cases). Corrections that follow:
- The 19b census's "7 × flow6 per dopri attempt" counted the dopri call
  graph; the executing loop is geodesic_step's 4 RK4 stages. The
  *conclusions* survive — both branches walk the same jet-tower working
  set, which is why the @inline collapse, Γ-fusion, and leaf-inline wins
  were real — but per-attempt arithmetic based on 7 stages should read 4.
- Item 3(e) inverts: the "cheaper integrator" is the status quo, dopri-
  with-tolerance is the untaken *expensive* path, and the h-cap/error-
  controller slack analysis described a branch production never enters.
  The live (e) question is now whether error control would buy accuracy
  worth its cost — the Richardson probe (still the ground-truth item's
  instrument) judges it.
- k1 reuse kept: +104 B/thread of *footprint* (occupancy doesn't bind),
  zero runtime cost, and it makes the dopri arm fair when (e) runs it as
  the accuracy reference.

**The round census** (768×576 static probe, F32, sweep=64, logging
sweep_stage! wrapper): rounds go 65114 → 3649 → 616 → 449 → 308 → 239 →
186 active; kernel 0.198 / 0.033 / 0.033 / 0.033 / 0.033 / 0.033 / 0.008 s.
Round 1 is 54% of kernel for 98% of the work; **rounds 3–7 are 37% for
≤616 rays — a ~33 ms latency floor per round that doesn't care how few
rays remain**. This is the 07-18 "46% tail" localized per-round, and it
prices the handoff exactly: a few hundred budget-capped rays are
milliseconds on 16 host threads.

**The handoff** (plan item 3b, agreed 2026-07-19 morning): DeviceSweep
gains `handoff` — a round with ≤ handoff active records runs
threaded_sweep! to completion (sweep = budget.max_steps) instead of
launching the kernel. Records already live in the host pool: pure control
flow at the sweep_stage! seam, no record change, flyvideo gpu=1 inherits
through the constructor defaults. Correctness: with handoff > pool size
every round runs on host and the raymap is **bit-identical to the CPU
wavefront (diff 0, F64 and F32)** — the "sweep granularity is
scheduling-only" claim survives run-to-completion sweeps; with the default
threshold, hybrid-vs-pure moves 219/27648 pixels (0.8%, the straggler
class), gpu_tracer parity stays 0 flips with quantiles in the 14c band,
and 384×288 shaded renders are eyeball-identical (855/110592 moved).

**Numbers** (same-day interleaved, min-of-reps; wall is the honest metric
now — host tail time is deliberately outside kernel_seconds):
- static probe 768×576 F32: wall 0.500 → 0.401 (**−20%**, threads=128 +
  handoff=1024), kernel 0.368 → 0.212. handoff=4096 within noise of 1024;
  threads=128 alone −7% kernel.
- mid-crossing frame (cube flight coasted to d=0.38, SituatedCamera):
  **neutral, correctly** — the census there is a single 442k-ray round;
  no tail exists, the handoff never fires. The lever is regime-matched:
  it pays on approach/ambient frames, vanishes at crossing.
- 192×144: device F32 wall 0.18 → **0.051, now beating CPU F32 (0.078)**
  — the per-round latency floor was exactly what starved small frames.
  07-18b's "CPU wins small frames" flips at the raymap level (flyvideo-
  level confirmation rides the next real render).
- remaining static-probe wall is ~0.19 s of wall–kernel gap (initial-pool
  entry solves + upload): lever (c) is now the biggest single item on the
  probe workload, ahead of deep fusion's share of the remaining 0.21 s
  kernel.

## 2026-07-19b — the traffic census overturns the lever ranking: flow6's interior is ~90% of per-attempt local traffic, sweep_ray's 7.9K frame is footprint-not-traffic; leaf-inline completion buys 3.5% kernel / 3–5% CPU, bit-identical

**Method.** Lever (a) opened with the promised census, byte-weighted this
time: PTX dump of kernel_sweep (F32), awk over per-function `__local_depot`
sizes and ld/st.local *bytes* (width-decoded), plus the call graph from
`call.uni` targets. The 07-18 census ranked by static depot size; this one
weights by dynamic execution (calls per attempt × bytes per body).

**The correction.** Per dopri attempt, local traffic is: 7 × flow6 body
(1752 ld + 2128 st B each) ≈ 27 KB, vs sweep_ray's own ≤ 1 KB, dopri_step
336 B, settle ~200 B (chart_transition ~1 KB but only on hops). The two
monster closures (`__56_*`, 2.5 KB/2 KB traffic each) are the ForwardDiff
`situate` closures — exit-path only, per passage not per attempt.
**sweep_ray's 7856 B depot is stack *footprint*, not traffic** — it gates
streamed bytes/wave only, and occupancy was already ruled out as the
binding constraint. So lever (a) as written (slim what crosses the
dopri_step boundary) attacks 4% of the traffic; the gate is flow6's
*interior* — the spilled working set of the collapsed jet tower.

**What still crossed ABI inside flow6.** eval_packed_partials' ntuple
do-block lowers to a closure too big for Julia's inliner despite the
enclosing @inline (6 calls/body, 10-lane sret round-trips through local);
same for eval_packed's; compose1 and vjet3 simply lacked annotations; plus
a small [2×f32] sincos-class math leaf (left alone). Fix: per-component
bodies extracted into named @inline functions (throat.jl, jets.jl),
annotations added.

**Result.** flow6 depot 1680 → 1360 B, body traffic 3880 → 3560 B (−8%);
whole-kernel local 13,976 → 13,952 B/thread. Interleaved A/B (protocol per
2026-07-19): kernel min 0.378/0.375 baseline vs 0.363/0.365 leaf-inlined —
**−3.5% kernel**; CPU recursive −4–5% (1.95 → 1.87 at 768×576, 0.123 →
0.116 at 192×144), wavefront F32 −3% at 192×144, neutral-noisy at 768.
physics_diff deviations **exactly 0** against the morning's baselines (the
extraction moved no arithmetic), 451 cold green, device parity counts
character-identical to the pre-leaf run (F64 0 flips, 377268/442187
bit-identical).

**Why the win is 3.5% and not 8%.** The sret slots vanished but most of
flow6's spill stayed: with registers capped at 255 the compiler spills
whatever live state exceeds them, call structure or not. Local traffic
tracks *live aggregate width*, not call count — which is exactly the case
for kaarel's deep-fusion vision (plan item 3d): shrinking the jet tower's
live lanes algebraically (directional jets, the Euler–Lagrange form) is
the only remaining way to cut the dominant traffic. Ranking after this
session: (d) deep fusion attacks the ~90% (throughput regime), (b) CPU
tail handoff attacks the straggler tail (latency regime, 46% of kernel),
(c) pool-init overlap attacks the wall-kernel gap; "slim sweep_ray's
frame" is retired as a traffic lever (it remains relevant only if
footprint ever binds).

## 2026-07-19 — Γ-contraction fusion lands: CPU 2–4%, device kernel NEUTRAL — post-collapse the kernel doesn't pay for FLOPs, and cross-session kernel numbers carry GPU clock state

**Context.** The agreed opener (2026-07-18 doctrine discussion, first
go-wild-with-reference move): never materialize the 27-component Γ on the
geodesic hot path. Contracted against the symmetric v^i v^j, the two
symmetric-derivative terms of Γ are the same number, so

    w_u = v^i v^j (∂_i g_uj − ½ ∂_u g_ij),    a = −g⁻¹ w

reference.jl carries the meaning (`geodesic_accel`, next to `wvel_along_v`
which stays as the general transport law — w ≢ v does not symmetrize);
production fuses it via t_c = (∂_c g)·v (three matvecs feed both
contractions), `metric_gradient` extracted from christoffel's branch body so
`christoffel = christoffel_from ∘ metric_gradient` survives for camera-frame
transport. geodesic_flow and flow6 switch on both twins; production
wvel_along_v deleted (dead). Contraction-stage FLOPs roughly 2.5x down
(~250 → ~95 mults+adds), but the stage is small next to the jet assembly.

**Certification.** Warm probe: fused vs materialized-Γ contraction ≤ 8.7e-17
relative over every cube chart × depth branch; production vs reference
≤ 4.4e-16; F32 eltype honest. 451 cold green. physics_diff 0 flips, worst
deviations moved in-band (cube 6.61e-15 → 8.55e-15, trefoil 1.61e-10 →
1.65e-10) — **re-baselined with --save** (the pre-agreed moment; baselines
now read 0). gpu_tracer 768×576 parity: F64 0 flips, 377k/442k bit-identical,
max 8.6e-8 rad; F32 1 flip + 3 res-mismatch in the 14c band.

**Perf: the benchmark trap first.** First gpu_tracer read 0.421 s kernel
(F32 sweep=64 bin) against the recorded 0.386 — an apparent 9% regression.
A/B against a pre-fusion jj workspace killed it: the *baseline* also
measures 0.42–0.46 today when run cold-ish, and both checkouts settle near
0.38 after two full-size warm passes. Last session's 0.386 was a
warm-clocks number. **Protocol from here: device A/B comparisons must be
same-day, interleaved, after ≥2 full-size warm passes; min-of-reps is the
statistic** (noise is one-sided). Alternated 15-rep runs
(baseline/fused/baseline/fused), kernel min: 0.382 / 0.375 / 0.382 / 0.380.

**Result: device kernel neutral (±1%), CPU wins 2–4%** (768×576 wavefront
F32 1.29 → 1.26, recursive 2.00 → 1.93; 192×144 wavefront 0.081 → 0.078,
recursive 0.129 → 0.123 — CPU numbers, unlike GPU ones, reproduce across
sessions). Registers 255 unchanged; local 13,968 → 13,976 B/thread (+8,
noise). Why neutral: post-collapse Γ already lived entirely inside flow6's
inlined frame — the fusion removes arithmetic but not one byte of the local
ABI traffic the kernel is bandwidth-bound on (2026-07-18 diagnosis,
now confirmed from the FLOP side). Corollary that re-ranks the levers:
**FLOP reduction does not pay on device until the frame-traffic gate is
broken** — sweep_ray's 7.9K frame (lever a) is the confirmed next target,
and the fusion's FLOP savings should be re-credited on device *after* the
kernel stops being traffic-bound. No flight-segment run: a lever neutral in
the static probe (which already contains the 46% tail-round regime) has
nothing for workload composition to amplify.

## 2026-07-18b — flyvideo grows gpu=1; the first F32 flight finds a latent crash in the certified CPU wavefront; GPU 2.3x CPU on the production-resolution flythrough, CPU wins at 192×144

**Context.** Kaarel's catch: the 07-17 space flythrough was assumed
GPU-rendered; in fact keyframe_raymap → render_raymap = the *recursive*
renderer, in *F64*, on *CPU* — the production video path had never used the
wavefront driver, F32, or the card. flyvideo now takes `gpu=1` (under
--project=gpu): frame raymaps run the wavefront driver with DeviceSweep
(one device throat upload, reused across frames), F32 against F32-rounded
tables, entry solves F64 against the F32-coefficient surface (the 14c
seam). Deliberately CPU F64: coast (flight path bit-identical to CPU runs),
flow probes, and refine! escalation — so escalated pixels land at F64.

**The surprise (kaarel predicted one).** Frame 4 of the first GPU cube
flight: KernelException. With -g2: InexactError at wedge_index —
`Int(...)` of NaN — via christoffel ← geodesic_flow. An F32 ray blows up
mid-step (filament class): vel goes non-finite, the NaN reaches a *stage
position* within one RK4 step, and wedge_index throws. On device one
throwing thread kills the whole launch. Reproduced identically on the
**CPU F32 wavefront** (cube flight camera at τ=1.6): the certified F32
path has carried this latent crash since it landed — 14c acceptance ran
static cameras; a flight's worth of viewpoints reaches the measure-zero
cases. Fix (control-flow only, never fires on finite rays — F64
bit-identity with trace_geodesic untouched): wedge_index returns wedge 0
on non-finite angle instead of throwing (any wedge is as good as any other
for a ray that is already NaN), and sweep_ray retires non-finite phases as
RAY_UNRESOLVED at the top of both loop branches — NaN rays die into the
standing escalation machinery, which retraces them at F64. Certified: 451
cold green; physics_diff worst deviations **exactly** unchanged
(6.61e-15 / 1.61e-10); cube-flight frames GPU vs CPU differ on ~4–6% of
pixels by ≤1 8-bit step (PAE 257/65535), eyeball-identical.

**Numbers.** Cube flight 192×144 uniform=0.4 (16 frames, crossing
included): GPU kernel total 14.6 s vs CPU total wall ~6 s — **CPU wins
small flights** (small frames starve the card; crossing frames are
straggler-latency-bound). Trefoil 768×576 uniform=0.02 from τ=6 (the
lens-heavy approach), equal 420 s wall windows including startup: **GPU 86
frames vs CPU 37** (~2.3x; the GPU window carries heavier startup — kernel
compile — so the steady-state per-frame advantage is larger). ~4 s/frame
GPU at production resolution vs the 07-17 render's 10.7 s/frame average.

**Workload shift, confirmed.** The flythrough regime is not the static
probe regime: mid-crossing frames put essentially all pixels in-throat
(the static 768×576 camera enters only 15%), situated cameras start
passage-0 legs in-throat, and escalation adds CPU F64 tail work per frame.
Execution-engineering measurements from here on should benchmark against a
flight segment, not the static camera. Per-frame device profiling of a
real flight is the next instrument to build when needed.

## 2026-07-18 — execution engineering opens: the 19.2K local was ABI call-frame stack, not live state; @inline tree collapse buys 31% kernel / 15% CPU; occupancy is NOT the gate

**Context.** Kaarel's call: set bundles/caches aside, see how far pure
execution engineering goes. Attack order agreed: profile first, rank levers
with data. No ncu/nsys on the box, so instrumentation was host-side (an
instrumented copy of run_wavefront! timing every stage per round) plus PTX
census. Workload throughout: the physics_diff cube, 768×576, F32, sweep=64
bin, RayBudget(0.05, 400, 4), RTX 4070 SUPER (56 SMs), 16 host threads.

**Wall decomposition (pre-change baseline).** Wall 0.649 s excluding pool
init: kernel 90%, per-round entry solves 0.5%, transfers+gather+scatter+sort
+compact ≈ 1.5%. Two corrections to the 2026-07-16 addendum's reading:

- The "~20% host" gap between wall and kernel-only is the **initial** pool
  fill (442k F64 first-entry solves + device throat build + transfers), not
  the per-round entry stage. Still pipelinable, now correctly attributed.
- Only 65,114 of 442,368 pixels enter the throat at all (this camera); 94%
  of entrants finish their passage within the first 64-attempt sweep.

The kernel time splits into two regimes with different physics:

    round  active  entries  kernel_ms      round  active  entries  kernel_ms
        1   65114    61465     313.5           5     309       69      50.4
        2    3649     3033      54.7           6     240       54      51.5
        3     616      167      49.6           7     186        5      13.2
        4     449      140      50.5

  Rounds 2–7 = 46% of kernel time for <6% of rays, ~50 ms/round *regardless
  of active count*: a 64-attempt round costs the same for 616 rays as for
  3649. That is per-attempt serial latency ~0.8 ms — local-memory traffic at
  raw L2/DRAM latency, unhidden. The stragglers are almost exactly the 181
  budget-capped unresolved pixels (side==0) burning 400 sequential attempts.
  Ray lifetimes: q50 1 round, q99 2, max 7 (the budget ceiling).

**Occupancy is not the gate** (twice measured). maxregs ∈ {128, 96, 64} ×
threads ∈ {128, 256, 512}: every register cap makes the kernel *slower*
(freed registers just move to local: 255→128 regs adds +496 B spill), both
at the 19.2K baseline and at the 14.0K post-change footprint. Doubling
theoretical occupancy does nothing because the kernel is bandwidth-bound on
its own local traffic (19.2 KB × 14.3k resident threads ≈ 275 MB streamed
per wave). One real scheduling nugget: threads=128 beats 256 by ~6% kernel
(0.364 vs 0.388 post-change) — finer blocks spread the tiny tail rounds
across more SMs. Candidate DeviceSweep default.

**Where the 19.2 KB actually lived.** PTX census (`CUDA.@device_code`; the
.asm dump is full PTX): Julia-codegen `__local_depot`s total just 656 B —
everything else is **ABI call-frame stack** along the deepest chain:
sweep_ray 7856 B → dopri_step 6776 → settle_phase 4064 → geodesic_step 3024
→ surface_jet 2128 → surface 1968 → chart_transition 1320 → christoffel 900
+ jet leaves (corner_contribution_jet 616, normal_jet 336, …). 224 call
sites, ~1300 st.local + ~1100 ld.local static: every call passes/returns
aggregates (Jet3 of SVector{3} = 30 floats, metric triples, the ~500 B
by-value Chart handle) through local memory. The 2026-07-16 attribution
("remaining local is integrator state + record") was wrong — it was call
ABI, mostly eliminable at the language level. Also seen: 49 call sites to
throw_boundserror; `--check-bounds=no` was measured and is a dud (−176 B,
−3% kernel, raymaps bit-identical) — bounds checks are not the cost.

**The lever that paid: @inline the christoffel tree into flow6.** Annotated
(src/chart.jl, jets.jl, throat.jl, geodesic.jl): the wedge chain
(safe_atan2, fake_complex_pow, wedge_index/square coords/to_chart + cjet
variants), corner/surface evaluators (polynomial_surface,
corner_contribution[_jet], surface, surface_jet), jet algebra (jadd, jdu/v,
leibniz, jdiv, cjet ops, cpow_jet, wirtinger_compose, flat_bump[_jet],
blend_scalar/blend_jet, eval_packed[_partials]), metric assembly (collar,
normal_jet, sym3, outer/inner_metric_gradient, christoffel_from,
christoffel, christoffel_pull, wvel_along_v). The structural frames
(sweep_ray, dopri_step, settle_phase, flow6, chart_transition, to_ambient)
stay ABI calls; flow6 remains the single cheap crossing (6 floats in/out).
Result: local 19,232 → **13,968 B/thread**, kernel 0.560 → **0.386 s**
(−31%), both regimes uniformly (round 1: 313→209 ms; tail rounds 50→33 ms).
Wall 0.712 → 0.535–0.57. The tree now sits in flow6's 1808 B frame (was
~4.5K spread over five layers). CPU wins too (192×144: recursive
0.154→0.129, wavefront F32 0.095→0.081 — same call-frame arithmetic, x86
edition).

**The blowup boundary, mapped.** `always_inline=true`: 187,936 B/thread,
launch OOMs — LLVM/ptxas cannot share stack slots across branches at that
scale. Inlining the next layer up (dopri_step/settle_phase/to_ambient/
chart_transition into sweep_ray): 126,672 B, same OOM — to_ambient drags the
AD collar tree (the dual-pass exit path, still ForwardDiff by design) into
the union. Inlining just flow6 into dopri_step: 21,480 B, kernel unchanged —
footprint up, no traffic win. The sweet spot is exactly "leaves collapsed
into flow6, structural frames separate," and it is now occupied.

**Certification.** 451 tests cold green. physics_diff: 0 side flips, worst
deviation 6.61e-15 cube / 1.61e-10 trefoil — *moved* from 8.1e-15 / 1.26e-10
(CPU arithmetic shifted at the muladd/fusion level once inlined; StaticArrays
uses muladd internally), still inside the 1e-9 band; this is the "stacking
reordering" case the 2026-07-16 note anticipated — re-baseline when the next
optimization lands on top. gpu_tracer parity re-run: F64 device vs CPU 0
flips, 23,581/27,639 bit-identical (same count as 2026-07-16), max 3.6e-11
rad; F32 in the 14c band. Device pre-vs-post (same precision, same tables):
0 flips, 94% bit-identical, q99 4e-7, max 1.4e-3 rad on one
boundary-filament ray — the NVPTX backend contracts/schedules differently
across former call boundaries, so device-side bit-stability across kernel
versions is not a thing @inline preserves (it was never a certified
invariant; parity is judged by the 14c decomposition).

**Lever ranking after this session** (next targets, in order):

1. **sweep_ray's own 7.9K frame** — unchanged by the collapse and now the
   biggest; its ld/st traffic runs per attempt iteration (record staging,
   chart handle, h/steps state around the dopri call). Census its body,
   then either slim what crosses the dopri_step boundary or restructure the
   loop so the hot state stays in registers.
2. **Tail rounds: 46% of kernel serves <1% of rays.** Per-attempt serial
   latency improved 50→33 ms/round but the structure stands. Candidates:
   CPU tail handoff — after round ~2, ship the few hundred stragglers to
   host threads (the sweep_stage! seam makes backend-per-round trivial; a
   400-attempt straggler is ~1 ms serial on CPU vs ~200 ms on device);
   budget-as-deadline for realtime; threads=128 as a free 6%.
3. **Pool-init overlap** — the ~0.15 s of initial F64 entry solves +
   device-throat build, pipelinable against round 1 (or an in-kernel entry
   solve; the F64 host solve was a design convenience, not a requirement).

## 2026-07-17 — unresolved-pixel budget escalation: retracing only the side==0 pixels under doubled budgets clears the approach for ~30% extra wall, and the mid-crossing filament band shrinks ~2x per doubling

**What was measured.** flyvideo.jl now escalates per frame: after the normal
raymap, the unresolved (side==0) pixels — and only those — are retraced from
scratch under a budget with max_steps and max_passages both doubled, repeating
until the unresolved fraction drops below `unres=` (default 0.005) or `esc=`
doublings are spent. Base budget RayBudget(0.05, 400, 4), trefoil flight,
uniform Δτ=0.02, 701 frames.

**Numbers.** 192×144, esc=3: whole-flight avg unresolved 302→37 px/frame
(1.09%→0.13%), worst frame 1950→262 (7.1%→0.95%); wall 8.3 min vs 6.0 min
baseline (+30%). 768×576, esc=4: avg 4852→550 (1.10%→0.12%), worst frame
31110→2210 (7.0%→0.50%); wall 125.5 min (10.7 s/frame — 16x the pixels of
192×144 at ~15x the time, so escalation keeps the scaling linear). Escalation
histogram at 768×576: esc=0 on 315 frames, 1 on 343, 2 on 18, 3 on 2, 4 on 23.

**Why the shape.** The unresolved population is kaarel's limit-cycle boundary
filament (2026-07-13): rays out of passages, not steps. One passage doubling
(4→8) resolves the great majority everywhere in the approach — esc=1 frames
go from ~1-2% to ~0.05% in one round. The stubborn band is exactly the
mid-crossing frames (chart d≈0.2-0.3, τ≈9.3-9.7, ~25 frames): there ~7% of
pixels start unresolved and each doubling only halves-ish the survivors
(frame 467: 27060→1623 after four doublings ≈ 2.0x/doubling) — consistent
with a measure-zero critical set being approached geometrically from both
sides. Since retracing costs only the pixels chased, the escalated frames pay
proportionally to their filament content; brute-force full-frame re-renders
would have paid the resolved 93% again each round.

**Escalation interacts with nothing.** The base raymap is untouched (frames
with esc=0 are bit-identical to pre-escalation flyvideo), refine! lives in the
script and mirrors keyframe_raymap's dispatch + render_raymap's per-pixel body,
and identical unresolved counts under nature and space skies re-confirm shading
and tracing are independent.

**Also recorded (same session, scene authoring, not physics):** first
non-identity `Placement` use — the trefoil flight's side-2 embedding is
rotated (Rodrigues, exit-fwd → destination-sky hero feature) so the outro
looks at something worth seeing; probe flow/flip stats are identical to
identity-placement runs (rigid rotation preserves both angles and sides), and
export/entry both route through pl.linear, so cross-and-return stays exact.
manual_wormhole's cubemap skyboxes convert seam-free via
scripts/cube2equirect.jl (face conventions replicated from deadsimple.rs
Skybox::sample; bilinear where the Rust is nearest).

## 2026-07-16 — hand-rolled christoffel derivatives land: 12x CPU, 21-34x device, device F32 reaches CPU parity at 384×288; F32 accuracy *improves*; all certification instruments green without re-baselining

**What changed** (src/jets.jl + surface_jet in chart.jl + christoffel assembly
in geodesic.jl). Production christoffel no longer runs three nested ForwardDiff
dual passes (3×3×2 seeding = 48 carrier floats per scalar, mostly structural
zeros — the spill source named 14b/14e). It now assembles g and ∂g closed-form
from one order-3 (u,v) jet of the blended surface. The old dual-pass body
survives verbatim as `christoffel_ad` (the in-tree oracle twin, pinned in
test/); reference.jl is untouched. Structure of the hand derivation:

- **The wedge chain is holomorphic** in z = u+iv: wedge_square_coords is
  rotation ∘ z^(n/4) ∘ scale, next_corner_coords is w ↦ i + (−i)·w, and
  square_coords_to_chart is scale ∘ z^(4/m) ∘ rotation. So its whole order-3
  jet is 4 complex numbers (CJet) with 1D chain rules; z^p derivatives are
  the recurrence φ⁽ᵏ⁺¹⁾ = (p−k)φ⁽ᵏ⁾/z reusing fake_complex_pow for the value.
  Real (u,v) partials fall out by Cauchy–Riemann at the seams to real code.
- **Corner polynomial ∘ holomorphic via Wirtinger calculus**: ∂z acts only
  through w, ∂z̄ only through w̄, so the bivariate Faà di Bruno factorizes
  into two 1D chains; h real ⇒ D_ab = conj(D_ba) leaves five D terms. The
  ten poly partials come from one table read (derivative-carrying double
  Horner: A′ₖ₊₁ = A′ₖx + Aₖ, A″ₖ₊₁ = A″ₖx + 2A′ₖ, …).
- **Blends are scalar exp(−1/x) chains** with closed-form b′..b‴, gated on
  b == 0 so F32 near-seam points can't make 0·Inf (the 1/x powers overflow
  F32 long before exp(−1/x) stops underflowing to zero — checked: past the
  gate w⁶ ≤ ~1.2e12 in F32, no overflow possible).
- **Depth enters exactly**: collar jacobian columns are linear in d
  (e_a = ∂_a s − d ∂_a n̂, ∂d-column = −n̂), the inner metric is d-free, the
  depth blend is a 1D chain — no third derivative direction ever needed, the
  step-4 "depth factors out exactly" observation cashed in. n̂ needs only an
  order-2 jet (cross + √ chains); Γ contraction is the same final block.
- **Value lanes reproduce the plain evaluators op-for-op** (same divisions,
  same addition order): surface_jet's value is bit-identical to surface(),
  which made oracle debugging trivial (any mismatch is in a derivative lane).

**Correctness.** surface_jet's 10 coefficients vs third-order *nested*
ForwardDiff of production surface(): worst 5.3e-13 rel (third-order lanes,
trefoil), first-order lanes ~1e-15; value lane exactly bit-identical, both
scenes, all charts (probe swept every vertex × 12 spiral points). christoffel
vs the AD twin and vs Reference.christoffel over 258k stratified points (both
scenes, both halves, depths spanning outer / t=0 / blend / t=1 / inner):
median ~2e-14, max 5.9e-13 rel. Cold: **451 tests green** (new testset pins
jet-vs-AD-twin at 1e-11 across all cube charts × depth branches + F32 eltype
honesty), and **physics_diff passes against the existing baselines** — 0 side
flips, worst deviation 8.1e-15 cube / 1.26e-10 trefoil, inside the 1e-9
arithmetic-reordering band. No re-baseline needed; the derivative-lane
arithmetic change lands within the tolerance physics_diff was built to absorb.

**F32 accuracy improves** (fewer, better-conditioned ops): jet christoffel on
F32 tables/phases vs F64 truth median 4.5e-7 / q99 2.3e-6 / max 4.3e-6
sup-rel — vs AD-F32's 9.5e-7 / 3.8e-6 / 1.1e-5 (14b). End-to-end,
precision_diff reproduces the 14c acceptance shape: f32/ulp quantile offset
flat at ~2^27–28.6, **0 side flips both scenes**, passage count still the
deviation ordinal. One 110k-ray render shows a single boundary-filament ray
at large F32 deviation (max angle 1.69 rad, q99 1.4e-6) — the 14c
undecidable-band class, not a regression.

**Throughput** (scripts/gpu_tracer.jl, cube, RayBudget(0.05,400,4), RTX 4070
SUPER, 16 host threads; wall seconds best-of-2, device best config):

    192×144 (27.6k rays)       14e (AD)    now (jets)   speedup
    recursive cpu               1.93        0.154        12.5x
    wavefront cpu F64           1.21        0.102        11.9x
    wavefront cpu F32           1.09        0.095        11.5x
    device F64 best            16.2         0.752        21.5x
    device F32 best             8.35        0.245        34x

    384×288 (110.6k rays):  wavefront cpu F32 0.365,  device F32 best 0.366
    96×72   (6.9k rays):    wavefront cpu F32 0.027,  device F32 best 0.222

CPU 1-thread christoffel microbench (cold): jet 1.49 µs/eval vs AD 18.3
µs/eval = **12.3x** — the register lever is equally a CPU lever, and the
end-to-end 11.5-12.5x render speedups confirm christoffel was ~90%+ of trace
cost. Kernel footprint: still 255 registers, but local (spill) bytes/thread
drop 42.5K → **19.2K** F32 and 68.9K → **23.6K** F64 (the remaining local is
integrator state + record, not duals). Device F64 is now 2.9x behind
16-thread CPU (1:64-FP64 card), device F32 **reaches CPU parity at 384×288**
(0.61x the recursive-renderer wall) and its per-ray cost is still falling
with ray count (32 → 8.9 → 3.3 µs/ray at 6.9k/27.6k/110.6k): with per-thread
traffic cut, the remaining gap to the card's arithmetic capacity is
**parallelism supply** — exactly 14e's second gate. Frame batches (flat pool
+ frame id) are the named next lever; persistent threads behind that.

**Why-chain.** 14e diagnosed spill-bound (F64 only 2x F32 on a 1:64 card;
kernel-only ≈ wall). Cutting per-scalar dual freight 48 → 10 real
coefficients (and n̂'s to 6) cut local traffic ~2.2-2.9x and bought 21-34x
kernel time — super-linear because spilled bytes were being re-read per
RK4/dopri stage, per passage. The F64-vs-F32 device ratio moving from 2.0x
to 2.9x says the kernel is now partly compute-bound: the lever did its job;
what's left is feeding the card.

**Addendum — kaarel's pushback measured: resolution alone feeds the card;
"rays in flight" is not a lever to build.** Ladder extended to 768×576
(442k rays), device F32 best wall µs/ray: 32.1 / 8.86 / 3.31 / 1.53 at
6.9k / 27.6k / 110.6k / 442k rays (kernel-only 1.18 at 442k) — still
falling, and the device F32 wall is now **2.1x faster than the 16-thread
CPU wavefront** (0.679 vs 1.44 s; 0.28x the recursive wall). F64 device
2.05 s vs CPU F64 1.66 s — within 1.25x of the CPU on a 1:64-FP64 card.
So production-scale resolutions supply the parallelism by themselves;
frame batching demotes to the offline-video corner (flyvideo: cameras all
known ahead — legitimate batch regime) and small-res test runs. Residual
starvation at fixed resolution decomposes into (a) wave quantization —
255 regs ⇒ one 256-thread block/SM ⇒ 14.3k resident threads/wave — and
(b) compaction-tail rounds running below wave size (the limit-cycle
stragglers); both shrink as fractions of the work with resolution. The
realtime-idiomatic fix for both, if ever needed, is **persistent threads**
pulling rays from the pool queue (same driver semantics, different launch
shape — the standard wavefront-path-tracing idiom), not cross-frame
batching, which costs input latency realtime can't pay. Realtime frame
deadlines map to budget caps + progressive fallback on unresolved boundary
filaments — semantics RayBudget/unresolved already has. Two scale notes:
wall > kernel-only at 442k (0.679 vs 0.52 — host entry solves/transfers
~20%, pipelinable by overlapping round k+1 entry solves with kernel round
k); and 1-2/442k F32 side flips vs truth — the boundary-filament class
appearing at scale, consistent with 14c density.

## 2026-07-14e — the device tracer kernel runs and is correct (0 flips, 85% of F64 rays bit-identical); naive throughput loses to the CPU ~5-8x, and the gate is spill, exactly where 14b pointed

**What ran** (scripts/gpu_tracer.jl; gpu/jupitergpu.jl is the promoted
module — PackedTable/adapt/device_throat moved out of gpu_christoffel.jl,
which was re-verified after the dedup). `sweep_ray` executes in-kernel
*unmodified* — settled RK4/dopri attempts, chart hops, half transitions,
to_ambient exits — one thread per ray over `WavefrontRay` records, inside
the same run_wavefront! driver the CPU uses with only the sweep stage
swapped (jupitergpu.DeviceSweep: gather active records → launch → scatter).
Cube scene, 192×144, RayBudget(0.05, 400, 4), RTX 4070 SUPER.

**Correctness — the tiering holds on device.** Zero side flips and zero
resolution mismatches in 27,648 rays for all three comparisons. F64 device
vs F64 CPU wavefront: 23,580/27,639 rays **bit-identical end-to-end**
through hundreds of steps; the rest deviate by transcendental-intrinsic
ulps only (angle q99 3.1e-15, max 2.9e-12 rad — CUDA's sincos/atan/pow are
not bit-equal to glibc's). F32 device vs F32 CPU: q99 1.5e-6, max 3.1e-3
rad (one near-boundary ray amplified, the 14c picture). F32 device vs F64
truth: q99 1.4e-6, max 1.4e-2 — mirrors 14c's acceptance data. Since the
device runs the same production sweep_ray that is bit-identical to the
recursive tracer, which physics_diff certifies against reference, the
oracle chain reaches the card.

**Throughput — the honest number.** Wall seconds (kernel-only ≈ wall, so
host entry solves + transfers + compaction are *negligible* — the wavefront
split is structurally right):

    recursive cpu (16 thr)      1.93
    wavefront cpu F64           1.21
    wavefront cpu F32           1.09
    device F64  best (sweep=256, bin)   16.2
    device F32  best (sweep=16, bin)     8.35
    device F32  sweep=64 bin/nobin      12.4 / 14.4

Naive device F32 is ~7.7x *slower* than the 16-thread CPU wavefront.
Binning is worth ~15% at sweep=64; sweep granularity moves things ~30%
non-monotonically (scheduling noise, not read into). F64 is only ~2x F32
on a 1:64-FP64 card — confirmation the kernel is **spill-bound, not
FLOP-bound**: 255 registers with 68.9KB (F64) / 42.5KB (F32) local
bytes/thread, up from 52K/27K for bare christoffel (integrator state on
top of nested-dual bloat).

**Why-chain, and why this was predictable.** 14b already said bare F32
christoffel naive ≈ break-even with the full CPU *at 200k-point
parallelism*. The tracer kernel is worse on both factors: ~1.6x more spill
per thread, and a 192×144 image supplies only 27.6k rays — shrinking as
compaction retires them — which cannot hide local-memory latency. The gate
is per-thread local traffic, so the ranked levers are: (1) **hand-rolled
closed-form derivatives in production** (the named register lever — each
nested-dual scalar hauls 48 carrier floats, mostly structural zeros;
s..∂³s through polynomial∘wedge∘blend is closed-form), (2) more rays in
flight (frame batches / bigger tiles — the fly-through use case supplies
them naturally), (3) persistent-thread scheduling. The port is
functionally complete; the remaining work is arithmetic density, not
structure.

**Addendum — resolution scaling (kaarel's question: doesn't more rays =
higher resolution too? Yes, measured).** Device F32 (sweep=16, bin) per-ray
wall time vs CPU wavefront F32, same scene/budget:

    96x72     6,912 rays   848 µs/ray   19.4x behind cpu
    192x144  27,648 rays   285 µs/ray    7.3x behind cpu
    384x288 110,592 rays   109 µs/ray    2.8x behind cpu (still improving)

The 255-regs-→-one-256-thread-block-per-SM arithmetic suggested saturation
near 14k rays (2 waves of 108 blocks over 56 SMs); the data refutes that —
wave quantization plus compaction-shrunk tail rounds leave SMs idle far
past nominal fill, and per-ray cost is still dropping at 110k rays. So
parallelism starvation is a real, resolution-curable chunk of the gap: at
production sizes the naive kernel is already only ~3x behind the CPU, and
frame batching (pool records are flat — batching frames is just a bigger
pool with a frame id) buys the same without touching resolution. The
CPU's per-ray cost is flat (~39-44 µs/ray) across the sweep, as expected.

## 2026-07-14d — wavefront restructure certified bit-identical, and the staged CPU driver is ~1.6x the recursive renderer for free

**Motivation.** Step 5's remaining blocker: the passage loop's sum-type
returns (`Union{Nothing, AmbientRay}` from trace_ray, `Union{SituatedPhase,
AmbientRay}` from trace_geodesic) heap-box on every passage and cannot exist
in a kernel; divergence wants rays grouped by depth branch (~2.3x,
2026-07-14). The restructure had to be certifiable as pure control-flow
surgery: same per-ray operation sequence, F64 bit-identical.

**Shape** (src/wavefront.jl). One isbits `WavefrontRay{T}` record — stage
tag (in-throat / at-entry / escaped / unresolved) + side, chart handle id,
passage and step counters, dopri's adaptive h, pos, vel. `sweep_ray` is the
kernel body: up to K settled integrator attempts in the ray's own precision,
replicating both trace_geodesic loops op-for-op (settle → exit-check order,
dopri's cap/accept/h-update order, exact carrier conversions), resumable
across sweeps because h and the attempt count live in the record. Mouth
entry stays a host-side F64 stage (the seam precision_diff certified;
passage budget charged there, reproducing trace_ray's exits-on-last-passage
shape). Driver: pool parallel to pixels, active-index list, sweep → compact
→ entry rounds; `bin` sorts the active list by the metric's three-way depth
branch — scheduling-only on CPU, the divergence lever on device. The
initial camera ray never lives in a T record, so its first entry solve sees
full F64; all later T→F64→T round trips are exact embeddings.

**Certification.** wavefront_raymap ≡ render_raymap bit-identically (side,
pos, vel) on the cube scene: RK4, starved budget (exercises out-of-steps
and out-of-passages), and dopri; sweep granularities 1/7/64; bin on/off;
ambient and in-throat cameras. `sweep_ray` allocates 0 bytes in both
integrator modes behind a function barrier (the naive `@allocated` on
non-const globals reports 80 phantom bytes of caller dispatch — barrier it).
Now pinned in test/ (289 cold); physics_diff exactly 0.0.

**Surprise datum.** The staged driver beats the recursive renderer on CPU:
96×72 cube, RayBudget(0.05, 400, 4), 16 threads — recursive 0.538 s,
wavefront 0.32–0.36 s (0.59–0.68x) across sweep/bin settings. Plausible
why-chain, not isolated: render_raymap threads over image *rows*, and ray
cost is wildly row-uneven (mouth-crossing rays cluster), so threads idle at
row barriers; the wavefront pool re-@threads over the compacted active list
each round, which is self-load-balancing, and the per-passage boxes are
gone. Not worth optimizing further on CPU — the point was the structure —
but it means the restructure costs nothing to adopt everywhere.

## 2026-07-14c — F32 traced end-to-end: quantile curves parallel to the 1-ulp baseline (offset ≈ 2^28 = "just rounding"), the undecidable band is ~1e-4 of pixels and sits on the limit-cycle boundary, renders indistinguishable

**Motivation.** 2026-07-14b established F32 christoffel accuracy per *eval*
(median 9.5e-7 sup-rel) but flagged that Γ-eval error is not exit-map error:
F32 state accumulates over hundreds of steps and chaotic stretches amplify.
Kaarel's metric correction governs the design here: the exit map
image→(side, ray) is discontinuous at the limit-cycle boundary, so sup
exit-map error is vacuous for *any* finite-precision method — measure the
physics_diff decomposition instead (side-flip rate = size of the undecidable
band; deviation conditional on side agreement, in tail-safe aggregates,
since deviation ~ ε·κ(ray) with κ power-law near the boundary means even the
conditional *mean* need not exist), plus eyeballed render pairs as the
acceptance bar.

**Method.** Two parts.

- *Eltype-honesty completed past christoffel* (the 14b sweep's remainder):
  the dopri tableau ratios (previously `(1/5)k1`-style Float64 literals),
  the tracers' h/tol/`floatmin()`, RK4's h, `half_transition`'s
  transition_depth read, `to_ambient`'s Placement reads, `emit_ray`'s 1.0,
  `view_phase_at_target`'s rotation angle (kept exact in Float64, one
  rounding into the carrier — the wedge_square_coords idiom), and
  `directional`'s Float64 derivative seed in ad.jl all convert to
  `carrier(x)` instead of promoting. Float64 semantics certified unchanged —
  275 tests green, physics_diff deviation exactly 0.0 — since every
  conversion is F64→F64 identity there.
- *The F32 trace design* (scripts/precision_diff.jl — tracked, not
  disposable: it is the standing precision-acceptance instrument,
  physics_diff's sibling): mouth entry stays an F64 Newton solve against the
  F32-coefficient surface (once per passage, well-conditioned, and the
  plausible wavefront split anyway); the phase converts to F32 at entry; the
  in-throat trace runs wholly in F32 on Float32-rounded polynomial tables
  (mirroring device_throat — with F64 tables eval_packed would silently
  promote every step back to F64); exits convert at to_ambient and re-enter
  the next passage in F64. 8000-ray golden-angle bundles through the two
  physics_diff scenes (budget 0.05/400/4), four variants: **truth** (F64),
  **ulp** (F64, direction perturbed by ~1 ulp — the "just rounding"
  baseline), **tab32** (F64 arithmetic on F32 tables — isolates
  representation rounding from arithmetic rounding), **f32** (the GPU
  candidate). Angles via 2·asin(‖n̂₁−n̂₂‖/2); conditional stats over rays
  where both variants resolve to the same side (sky-miss rays are in the
  denominator and deviate by exactly 0 there — cube 67%, trefoil 35% of the
  bundle).

**Findings** (conditional angular deviation in radians; 1e-3 rad ≈ ⅓ pixel
at 384-wide/60°-fov, sub-pixel at anything plausible):

- cube: ulp q50/q90/q99/max 1.7e-16 / 3.3e-15 / 3.0e-14 / 7.5e-11;
  tab32 0 / 1.4e-7 / 1.4e-6 / 0.021; **f32 0 / 5.0e-7 / 2.9e-6 / 0.118**.
  Side flips 0/8000 every variant. Exceedance(1e-3) for f32: 0.04%
  (3 rays, the same grazing-rim rays that dominate tab32's tail).
- trefoil: ulp 9.1e-15 / 9.2e-14 / 3.8e-12 / 2.4e-9;
  tab32 8.7e-7 / 8.3e-6 / 3.2e-4 / 0.31; **f32 2.8e-6 / 3.1e-5 / 1.3e-3 /
  0.46**. Side flips: f32 1/8000 + 2 resolution mismatches (ulp and tab32:
  0). Exceedance for f32: 5.1% at 1e-4, **1.1% at 1e-3**.
- **The structural test passes**: the f32/ulp quantile offset on the
  trefoil is flat at 28.2–28.4 log2 across q50–q99 — parallel curves, no
  bends — i.e. F32 behaves as pure precision scaling of the same method,
  with offset ≈ the 2^29 mantissa-width ratio. (Cube offsets 26.6–27.7 on
  the entered subset; its bundle is mostly sky and single-passage.)
- **Error factorization** (tab32 vs f32, trefoil): representation rounding
  alone (F32 tables under F64 arithmetic) is ~⅓ of the f32 median
  (8.7e-7 of 2.8e-6) with the same tail shape; F32 step arithmetic
  multiplies ~3x on top. No arithmetic reordering can beat the tab32 curve
  while tables are F32 — that is the representation floor.
- **Passage count works as the boundary-distance ordinal**: trefoil f32
  medians by truth passage count 5.2e-6 (p=1, n=4201) → 4.0e-5 (p=2,
  n=728) → 2.0e-4 (p=3, n=207) — roughly ×8 per passage — and all three
  side flips live at p ≥ 3. The ulp baseline shows the same amplification
  from the other end: 1e-16 input noise reaches 2.4e-9 max — F64's own
  undecidable band is merely thinner, not absent.
- **Eyeball bar** (384×288 render pairs + diff images, both scenes):
  indistinguishable by eye. The diff structure is exactly the theory's
  picture — grayscale deviation on the knot silhouette ~1e-5, bright
  filaments along the boundary between the two ambient images, and the
  flipped pixels (cube 1/110592, trefoil 12/110592 = 0.011%) sitting *on*
  those filaments.

**Read.** F32 is accepted as the tracer arithmetic at the image standard
that already picked RK4: renders indistinguishable, wrong-for-precision
pixels ~1e-4 of the image concentrated where the F64-fallback flag was
already planned (the same filaments are where F64's own chaotic shimmer
lives — the ray-bundles item). The acceptance-aligned scalars to re-run as
the port proceeds: side-flip rate + exceedance-at-1-mrad, by passage count
(the cheap flag ordinal until Jacobi-field κ lands). The sweep also covers
emit_ray/transport, so in-throat cameras go F32 the same way; reference.jl
stays untouched as the F64 oracle. Remaining step-5 work is the wavefront
restructure (+ depth binning, ~2x from 14b's divergence datum).

## 2026-07-14b — production christoffel straight on the GPU: Float64 loses to the CPU, Float32 matches it naive and beats every tabulation variant on accuracy

**Motivation.** Kaarel's call for this session: before tabulating Γ (and
before the representation fork that tabulation accuracy opened), try the
simpler road — run the production AD christoffel *directly* in a kernel and
see what the card actually does with it. The step-3 smoke test only compiled
the first-order dual pass; christoffel is second-order (a 3-partial seed over
`metric`, whose collar/surface evals are themselves dual passes — innermost
scalars are (1+3)(1+3)(1+2) = 48 carrier floats), which was the untested
compilability wall.

**Method.** The type layer was made storage-parametric (pure type
generalization; 275 tests + physics_diff bit-identical), so the real cube
`Throat` Adapt-s onto the device: ragged construction-only Mesh fields drop
to `nothing`, packed_polys pad into one 4D array behind a view wrapper
(padding sits at the leading end of the descending Horner loop — evaluation
bit-identical). The kernel builds the tracer's own `SituatedPhase` from
half-edge ids and calls `jupiter.christoffel` unmodified
(scripts/gpu_christoffel.jl, kept). 200k points sweep all 8 charts, all
wedges, st ∈ [0.02, 0.98]², and all three depth branches (25% pure outer
d<0, the rest through blend and pure inner). CPU truth: the same production
christoffel, Float64, 1 thread (20–22 µs/eval, matching the step-4 26–32 µs
on fancier charts). For the Float32 leg, the christoffel path was first made
eltype-honest — Float64 literals (sqrt2, the 2π sincos wedge angles, wedge
powers n/4 and 4/n, 0.5) and ThroatParams field reads convert to
`carrier(x)` = the innermost primal type through the dual nesting — with
Float64 certified bit-identical (tests + physics_diff 0.0) since every
conversion is F64→F64 identity there.

**Findings** (RTX 4070 SUPER, 1:64 FP64; npoints=200k, threads 64/128/256
indistinguishable, 512 exceeds the 64K regs/block budget):

- **No compilability wall.** The full second-order nested-dual christoffel —
  wedge transcendentals, blend, collar, 3×3 inverse — compiles clean and
  matches CPU to 2.1e-13 abs / 2.1e-14 sup-relative (device libm, not bugs).
- **Float64 is spill-bound and loses to the multicore CPU**: 255 registers
  (the cap) + 52,648 B/thread local spill; 3.9–4.1 µs/eval mixed = 5.0–5.5x
  one CPU thread ≈ 1/3 of the 16-thread CPU. Not a usable port.
- **Float32 is 2.8x faster than F64-on-device and ≈ the whole CPU, naive**:
  1.43 µs/eval mixed = 14.3x one CPU thread (0.70 Meval/s). Still 255 regs +
  27,672 B spill — spill-bound, not flop-bound, so this is a floor, not the
  card's ceiling.
- **Branch divergence costs ~2.3x on both types**: uniform-depth-branch
  launches give 2.09–2.14 (outer) / 2.13–2.14 (blend) / 0.14–0.17 (inner)
  µs/eval at F64, 0.58 / 0.63 / 0.06 at F32, vs 3.9 / 1.43 mixed — warps
  near-fully serialize the three depth branches (blend ≈ outer + inner as
  expected: the blend evaluates both metrics; inner is 15x cheaper, no
  collar jacobian). Depth-binning rays in a wavefront structure recovers
  ~2x ⇒ F32 ~29x one CPU thread without touching physics.
- **Float32 Γ accuracy beats every tabulation variant probed**: vs Float64
  truth, sup-relative per point median 9.5e-7, p99 3.8e-6, max 1.1e-5 —
  better than the C⁴-blend table target (~1e-5), far better than exp-flat
  (~1e-3), and it is the *exact* blend semantics merely rounded, so nothing
  is re-fit, re-baselined, or re-rendered.

**Read.** The corner-blend accuracy problem that opened the representation
fork does not gate the GPU road: the card can afford to evaluate the truth.
The fork stays open as a CPU-side/representation question (tabulation still
wins 10–50x *on the CPU*, and affine atlases/level sets have independent
virtues), but the port no longer waits on it. What F32 does **not** yet
answer is end-to-end: Γ-eval error is not exit-map error — F32 positions
accumulate over hundreds of steps and chaotic bands amplify; the planned
policy (F32 with F64 fallback for flagged pixels) needs a traced comparison,
which wants the eltype-honesty sweep extended past christoffel to the
steppers/transitions (settle, half_transition, the integrators still carry
F64 literals). That, plus the wavefront restructure, is the remaining
step-5 work; register-pressure surgery (splitting the fused dual passes) is
the optimization deliberately not done today.

## 2026-07-14 — tabulated Γ (GPU step 4) probed: depth factors out exactly; lateral convergence is blend-limited, and the two blends fail in complementary places

**Motivation.** Step 4 planned "Chebyshev-tabulated Γ per face×depth box
(~8³ coeffs/component) … geometric convergence for a C^∞ metric, likely
10–50x per step". Before implementing: measure where integrator stages
actually evaluate Γ, and what Chebyshev convergence the blended metric really
delivers. Face-coordinate tables were rejected on paper first: face coords
are z^(4/n)-singular at irregular corners, while chart coords are smooth by
construction — so tables live in chart coords.

**Method.** Four probe rounds (scratchpad, disposable). (1) A recording env
(`christoffel(env::Recorder, v)` delegating to the AD path — the env slot is
a clean dispatch extension point) traced the physics_diff bundles and logged
per-eval wedge-square coords and depth. (2–4) Tensor-Chebyshev fits
(first-kind nodes, hand-rolled DCT + series-derivative recurrence) of Γ and
of exact-in-d metric pieces, on wedge-aligned domains; end-to-end assembled-Γ
error vs the AD christoffel on the needed region; a blend-parametrized local
copy of the surface→metric→Γ pipeline (reference-style ForwardDiff, verified
4.4e-16 against production) to compare corner blends; tiny-interior-box
anchor test (collapses to g~1e-14/Γ~1e-12, validating the pipeline — a
round-2 "anchor" that spanned the full wedge angle failed at 0.75 and was
meaningless, not a bug).

**Findings.**

- **Stage excursions are large; tables need margins + an AD fallback.**
  Settled phases live in st ∈ [0,1/2]², d ∈ [0, td], but RK/DP stages
  evaluate beyond: cube maxst q99/max 0.51/0.53, d ∈ [−0.036, 0.79]; trefoil
  maxst q99/max 0.80/1.12 with 19% of evals past 0.5 (fixed h=0.05 ≈ half a
  face on that mesh), d ∈ [−0.048, 0.35]. Rays genuinely dip to d<0 with
  inward velocity before `exits_mouth` triggers.
- **Depth is exactly solved, no approximation needed**: g_outer(u,v,d) is
  exactly quadratic in d (collar = s − d·n̂ ⇒ JᵀJ degree 2; residual 6e-16),
  g_inner is d-independent, and the depth blend w(d/cyl) with its flat ends
  reproduces the outer/inner branch structure of `metric` exactly (residual
  3e-16, including d<0 and d>cyl). So tabulate 2D lateral pieces A,B,C,im
  with g = w·(A+Bd+Cd²)+(1−w)·im and ∂_d g = w′·(om−im)+w·(B+2Cd); w, w′
  evaluated exactly at runtime. The depth blend never needs to change and
  never enters the approximation error.
- **Naive boxes fail; seam-aligned polar sectors are the right lateral
  domain.** Wedge-frame bounding boxes straddle the corner-blend seam rays
  and barely converge (cube Γ rel error 0.28→0.12 from N=6→20). Per-wedge
  (r,φ) sectors put seams and chart center on box edges; a small Cartesian
  core box avoids the 1/r series-derivative blowup at the center. Also the
  partition of unity dies outside the face ring (0/0), so domains must stay
  inside st ≲ 1 — box validity is a hard constraint, checked at fit nodes.
- **Convergence under the production exp-flat corner blend plateaus at
  ~1e-3 relative Γ.** Trefoil (n=4, wedge maps trivial): rel 2.8e-4 at N=24
  single sector. Whole-mesh at feasible N: q50 per-chart Γsup 3.3e-3 (N=8) /
  1.7e-3 (N=12), zero charts below 1e-4, 0.84/1.9 GB uniform. Cube (n=3):
  sectors stuck at ~0.013 sup even at N=24; the exp-flat *core* however
  converges cleanly when refined (1.2e-6 at Ncore=32, rc=0.05R — exp-flat is
  flat to all orders at the vertex).
- **A C⁴ polynomial corner blend (nonic smootherstep) flips the picture.**
  Cube sectors collapse 1000x to 1.2e-5 at N=24; trefoil reaches rel 8.4e-5
  (N=16) / 2.1e-5 (N=24) and still dropping steeply. But at irregular
  vertices its core stops converging (~6e-3 sup, insensitive to rc: the
  polynomial blend composed with the z^(n/4) wedge coords leaves a genuine
  fractional-power singularity at the vertex, where exp-flat was flat).
  Complementary failure modes: exp-flat is bad at seams / good at centers;
  polynomial is good at seams / bad at irregular centers. A C² blend
  (quintic) is strictly worse than both for Γ (needs ∂g ⇒ only C¹; rel
  errors 0.2–0.9, non-converging) — C⁴ is the minimum interesting order.
- **The AD christoffel costs 26–32 µs/eval warm** (cube/trefoil, single
  thread) — it dominates ray cost (4 evals/RK4 step, 7/DP5 attempt).
  Estimated table eval (values + derivatives via Clenshaw from one 24×N²
  coefficient set, N=8–12, plus ~300 flop Γ assembly) is ~1–4 µs: the
  roadmap's 10–50x per-step stands *if* the accuracy question is settled.
  Derivative tables need not be stored (evaluate d/du series from value
  coefficients), so memory is 24·N²·8B per wedge ≈ 12–28 kB.

**Conclusions.** Tabulation is worth it (26–32 µs → ~1–4 µs per Γ), the
exact-in-d factorization + seam-aligned sectors + tiny-core-AD-fallback is
the right architecture, and accuracy is decided by the corner blend, which
plan/claude.md records as "placeholder pending a deliberate choice":

- Keep exp-flat ⇒ tables certify at ~1e-3 relative Γ (visually sub-pixel-ish,
  but a real tolerance downgrade in the reference chain).
- Adopt a C⁴ polynomial corner blend ⇒ ~1e-5 relative Γ at N=12–24 with AD
  fallback inside r < ~0.1R of irregular vertices (regular-vertex meshes like
  the trefoil have no hard cores at all). This is a deliberate physics change:
  surface geometry shifts subtly, all baselines re-save, galleries re-render.
- Orthogonal lever either way: per-chart step-size discipline (h=0.05 is half
  a face on the trefoil) would shrink the needed margins and every error above.

## 2026-07-13c — de-pinning hash iteration order from mesh construction: the slack is ~1e-13/e-10, within budget

**Motivation.** Step 1 below froze Dict hash-iteration order as implicit
semantics (the CC rebuild helpers reconstructed the historical insertion
sequence so refined-mesh numbering stayed bit-identical). Kaarel's call:
incidental hash order must not be pinned — allow the slack downstream and
let construction be canonical instead.

**Method.** Pre/post pointwise comparison of fitted surfaces
(`sample_surface` grids are anchored at original-mesh face handles, whose
numbering never changes, so points pair up across the change regardless of
gauge); physics_diff against the *old* baselines before re-saving; knot-ray
timing/alloc probe.

**Findings.**

- What the order actually fed: refined-mesh vertex/face *numbering* (pure
  labels), FP summation order in the CC vertex means, and each vertex's
  canonical frame = first neighbor (a gauge choice — rotates chart
  coordinates, not the embedded surface). Nothing semantic. Construction is
  now deterministic without hash order: encounter-order neighbor lists,
  first-encounter edge enumeration, face-slot refined-face numbering; the
  HalfEdge struct, the half_edges Dict builder, and both order-freezing CC
  rebuild helpers are deleted. The one remaining construction Dict is
  lookup-only (order never observed).
- **Fitted surfaces moved 1.6e-14 (cube) / 5.8e-15 (trefoil) pointwise** —
  pure arithmetic reordering, geometry intact.
- **physics_diff passed against the old baselines**: 0 side flips, worst
  deviations 8.2e-13 cube / 2.7e-10 trefoil (vs 1.34e-14 / 1.84e-10
  previously certified; 1e-9 budget). Baselines re-saved to restore
  headroom (deliberate change; the 1 unresolved trefoil ray was unresolved
  before too). 275/275 cold; knot-ray timing and allocation flat.
- Side benefit: physics_diff's claim that baselines survive Julia version
  changes is now actually true of the construction path — hash order was
  the remaining version dependence.

## 2026-07-13b — half-edge flattening (GPU groundwork step 1): Dicts weren't a measurable CPU cost; immutable-Mesh inlining grows the passage boxes

**Motivation.** GPU port step 1: replace Mesh's name-keyed Dicts
(`half_edges`, `half_edge_offsets`) with flat id-indexed arrays so the
hot-loop connectivity is device-shaped. Certified as a pure representation
change (bit-identity), with before/after timing and allocation.

**Method.** Three deterministic knot-hitting trefoil rays (physics_diff aim +
two 0.01 tilts toward e1/e2), cold runs, min of 5; `@allocated` per ray;
`Profile.Allocs` at sample_rate 1 by type, before-side run via a jj workspace
at the parent commit.

**Findings.**

- **Timing unchanged**: 3.27/2.32/2.29 → 3.27/2.31/2.28 ms. The per-step Dict
  lookups (`next`/`twin`/`prev`/`half_edge_offset` in chart transitions) were
  not a measurable CPU cost — the step loop is dominated by collar/Γ
  arithmetic. The value of the flattening is representational (GPU-ready
  tables), not CPU speed.
- **Bit-identity trap found and dodged**: `catmullclark` numbers the refined
  mesh's faces by `collect(he_to_face)` — Dict hash-iteration order. The
  rebuild helpers (`half_edge_name_to_face_index`/`..._to_next_name`)
  reconstruct the transient Dict with the historical insertion sequence, so
  iteration order — hence CC numbering, hence fitted polynomials — is
  unchanged. Confirmed: physics_diff worst deviations bit-identical to the
  previous certified run (1.34e-14 cube, 1.84e-10 trefoil), 0 flips, 275/275
  cold.
- **Allocation +240–416 B/ray, zero new sites**: same 7 per-passage sum-type
  boxes (2-passage knot ray). Mechanism: Mesh is immutable and Julia inlines
  immutable structs into parents, so every box embedding a Mesh-carrying
  handle inlines Mesh's field block, now 11 pointer fields (88 B) vs 5
  (40 B). A `SituatedPhase` box holds a Chart = two embedded Mesh copies
  (HalfThroat chain + HalfEdgeHandle) → +96 B each, matching the profile
  (800 → 992 B for n=2). Bytes not boxes, per-passage not per-step; the
  device port carries only the flat tables. If it ever matters: a mutable
  Mesh (by-reference everywhere) or a slimmer hot handle.
- **Step 2 (mouth concretization): two BVH stack representations died before
  the threaded form.** Replacing the traversal recursion with an explicit
  stack for device parity: (1) `MVector{64,Int}` heap-allocates 544 B per
  entry attempt (~+1.6 kB/ray) — MArray `setindex!` goes through
  `pointer_from_objref`, which defeats escape analysis; (2) an
  `NTuple{64,Int}` stack via `Base.setindex` explodes to 330–500 kB/ray —
  the 66-argument `_setindex` splat exceeds vararg specialization and every
  push boxes the whole tuple (also ~+3% time). Resolution: **threaded BVH**
  (skip-pointer DFS, one extra Int per node computed in a reverse sweep at
  build) — no stack at all, allocation exactly back to the step-1 floor,
  timing unchanged, and it's the GPU-native shape (no local-memory stack).
  Right-child-first descent preserves the recursive visit order, so best-hit
  evolution is bit-identical: physics_diff again 1.34e-14/1.84e-10, 0 flips,
  275/275 cold.
- **Step 3 (device smoke), RTX 4070 SUPER, Float64, 1M points each**
  (scripts/gpu_smoke.jl, gpu/ subenv): eval_packed Horner kernel
  **bit-identical** to CPU (max |Δ| = 0.0), 0.78 ms = 1.29 Geval/s, 76x one
  CPU thread. The ad.jl dual pass (value_and_jacobian_columns through
  square_coords_to_chart → fake_complex_pow/safe_atan2) — ForwardDiff duals,
  StaticArrays, closures, SMatrix hcat — compiles clean on device: 5.8 ms =
  172 Meval/s, 31x one thread, max deviation 7.2e-16 value / 5.0e-15
  jacobian (fma/libm reordering). 56/78 registers — headroom for a full RK4
  step. **No compilability wall in the AD core.** Remaining port cost is
  what we already knew: device throat view (flattened packed_polys), and the
  per-passage sum-type boxes (kernels can't heap-allocate → wavefront
  restructure). eval_packed's signature loosened to AbstractArray{Float64,3}
  so device arrays dispatch to the same body.

## 2026-07-13 — temporal sampling for fly-through video: the flow-rate profile, and two ways a flow controller can die

**Motivation.** How densely to sample a fly-through in the flight parameter τ,
for both a perceptually-smooth adaptive video and an honest constant-speed one.

**Method.** Optical flow measured directly: sparse probe raymaps (16×12 rays)
before/after each candidate step; per-pixel angular change of exit directions
(same-side pixels), expressed in render pixels (192-wide, fov π/3.2 →
0.29°/px); side-flip fraction as a second signal. Adaptive controller in the
DP5 accept/reject shape, growth factor 0.85·target/flow. Cube scene, textured
skies, camera straight through the throat, τ 0→6.5. scripts/flyvideo.jl.

**Findings.**

- **Flow-rate profile** (px of bulk feature motion per unit τ, 3px/frame
  target → Δτ = 3/rate): far approach ~195; final exterior approach (τ
  1.5–2.0, mouth filling the view) ~610 median, 805 peak — the global worst;
  entry zone ~480; throat interior 115–160; departure ~0 (no parallax against
  the far sky — Δτ immediately maxes out; the whole exit took 7 frames).
  **Constant-speed answer: uniform Δτ = target/peak = 0.0037 for 3 px.**
- **Controller death #1 — max statistic.** Near-critical-ring pixels' exit
  directions change ~2.5px+ for ANY camera displacement (chaotic
  sensitivity, does not scale with Δτ). Steering on max flow pinned Δτ at
  ~0.002 from the very start, 100x oversampling.
- **Controller death #2 — any fixed upper quantile.** q90 worked until the
  mouth filled >10% of the probe with near-critical pixels (τ≈1.4+), then
  pinned exactly like max (2.5px at Δτ=0.0008; measured max there up to
  319px) and the un-floored growth formula ground to a halt: 630 frames to
  reach τ=1.77. Fix: steer on the **median** + a hard Δτ floor (0.004);
  q90/max stay in the log as an audit of the chaotic population.
- **Chaotic shimmer is image-space aliasing, not temporal undersampling** —
  no frame density smooths it (it needs ray bundles / supersampling, see
  roadmap). The sensitive population clusters on the **boundary between the
  two ambient spaces' images** in the camera's view; kaarel's correction of
  an earlier "grazing rays" reading: the truly *expensive* rays are the
  limit-cycle ones sitting on that boundary — they are the unresolved
  pixels. Whether that same set drives the flow-controller trouble is
  plausible but unconfirmed. This reading explains the entry behavior
  cleanly: the moment the camera crossed d=0 the chaos collapsed
  (q90≈max≈median, flips 0, unresolved 0) because the two-image boundary
  left the forward FOV — not because the interior is special. Consistently,
  the side-2 look-back frame (flythrough demo frame 9), which points at that
  boundary again, had 60 unresolved pixels.
- Run artifacts: 839 frames (652 from the q90 run whose tail oversamples,
  187 from the resumed median run — resume is exact because the pre-entry
  camera state is closed-form in flat space), resampled to 314 at uniform
  3px perceptual speed (scratch script; rate profile from the logs), videos
  in out/flythrough_adaptive{,_raw}.mp4, logs preserved as
  out/flyvideo_run{1,2}.log. Render cost ~19min for the median-run portion
  at 192×144.

## 2026-07-12 — jacobian_columns twin retired: the inner pass is honestly a different operation, and fusing it makes rays ~18% faster

**Motivation.** The `nested_jacobian_columns` twin (below) worked but smelled
like a workaround. Question: is there a way to work *with* inference instead?

**Findings.**

- No user-facing knob exists for the per-method recursion-widening heuristic
  (it fires when the same method appears on inference's own stack with a
  grown signature — nested `Dual{Tag2, Dual{Tag1}}` types guarantee growth).
  Alternatives weighed: a concrete typeassert at the inner call site stops
  the Any from propagating downstream but leaves a dynamic dispatch plus a
  boxed SMatrix-of-nested-Dual return per collar evaluation (~600 B × several
  per step ≈ right back at MB/ray); two `directional` passes in
  surface_normal_out breaks the cycle via a different method but pays an
  extra primal surface evaluation; OpaqueClosure barriers are exotic. The
  clean resolution: give the inner level *different semantics* — precedented
  by ForwardDiff itself, whose hessian is jacobian-over-gradient, distinct
  methods per nesting level, never jacobian twice.
- The different semantics were already wanted: collar evaluated the surface
  twice (plain value + surface_normal_out's dual pass over the same point).
  New primitive `value_and_jacobian_columns` returns both from one dual pass;
  the value lane performs exactly the plain evaluation's arithmetic, so the
  value is bit-identical and free. `jacobian_columns` and it keep separate
  bodies (any shared helper would re-form the method cycle).
- **Numbers** (three deterministic knot-hitting trefoil rays: physics_diff
  aim + two small tilts, cold runs, min of 5): 4.08/5.37/4.98 ms →
  3.32/4.40/4.06 ms (−18% each); allocation unchanged at 2208 B/ray (the
  sum-type boxes — confirms no inference regression). Certified: 219/219
  cold, physics_diff 0 flips with worst deviations bit-identical to the
  pre-change run (1.34e-14 cube, 1.84e-10 trefoil).
- Reference keeps the two-evaluation collar: that *is* its simplest correct
  form; the fused production collar is exactly the kind of elaboration the
  equivalence testset exists to hold in place.

## 2026-07-12 — per-ray allocation: ~1 MB was one inference failure, not many small ones

**Motivation.** After the StaticArrays pass, ~1 MB/ray remained (trefoil,
knot-hitting). A priori the tracer needs ~0 heap: a geodesic step is a
fixed-size state update, outputs land in the preallocated raymap, so the
algorithm's honest floor is a few boxed sum-type returns per passage
(`enter_mouth` → `SituatedPhase | nothing`, `trace_geodesic` →
`AmbientRay | SituatedPhase`), i.e. well under 1 kB/ray.

**Method.** `@allocated` on `trace_ray` for three camera rays (trefoil 
scene), `Profile.Allocs` at sample_rate 1 grouped by nearest jupiter frame
and by type; `Test.@inferred` + `code_warntype` bisection down the call tree.

**Findings.**

- Measured 0.7–1.4 MB/ray, ~7 kB and ~42 allocations per RK4 step — closures,
  `SMatrix`-of-`Dual` intermediates, and `SituatedPhase`s all heap-boxed.
- Root cause was singular: Julia's *method self-recursion widening*. The chain
  `jacobian_columns → collar → surface_normal_out → jacobian_columns` puts the
  same method on the inference stack twice with a grown signature; the
  termination heuristic widens the inner call to `Any`, which poisons
  `outer_metric → metric → christoffel → wvel_along_v → geodesic_step` — the
  entire step loop ran boxed. Every callee inferred fine in isolation; only
  the nested chain failed. Fix: a structurally identical twin
  (`nested_jacobian_columns`) for the inner level.
- Same trap, second appearance: rewriting the BVH traversal recursively with a
  `nothing`-then-tuple accumulator re-tripped the identical widening
  (signature grows Nothing → Tuple across the self-call) and re-poisoned
  `enter_mouth` to 54 kB/ray. Threading fixed-type state
  `(t, (u, v), tri_index; 0 = none)` fixes it. Moral: in hot recursive/nested
  code, keep self-call signatures *exactly* constant — Union accumulators and
  method reuse across nesting levels are both inference hazards that profile
  as "death by boxing" with no single hot allocation site.
- Cleanups riding along: `next`/`twin`/`prev` allocated a Vector per call via
  NTuple range-slicing (`vertices[3:-1:2]` hits a generic `map`); explicit
  index pairs now. `AmbientRay{P,V}` parametrized; `RayBudget`/`Camera` fields
  concretized; `Scene` parametric in mouths+sky; BVH stack Vector replaced by
  the recursion above (right-child-first preserves the old pop order exactly).
- **After:** 1.5 kB and 3 allocations per knot-hitting ray (the sum-type
  boxes; per-passage, step loop allocation-free). Knot-hitting ray 5.5 →
  2.8 ms; trefoil 384×288 render 14.8 → 10.1 s wall on 16 threads. Certified:
  219/219 cold, physics_diff 0 flips, worst deviation 1.8e-10 (FMA/inlining
  reordering, within the 1e-9 budget).

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

**Conclusions.** (See also the collar-caustic entry below: cylinder_depth is
not capped by collar self-intersection, so it is a genuinely free gentleness
knob, up to the in-principle limits in the last point.)

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

## Recorded 2026-07-11, measured in an earlier session — collar caustics are benign

**Question.** Past the limit surface's focal radius the collar map
self-intersects (a caustic; its Jacobian drops rank), which naively caps
cylinder_depth. Does pushing cylinder_depth past the caustic make the blended
geometry singular?

**Conclusion.** No — cylinder_depth can be made arbitrarily high. The final
space is a genuine Riemannian 3-manifold regardless: on the open band the
blend has strictly positive inner weight, and g_outer = JᵀJ stays smooth even
past the caustic (merely dropping rank on it), so the blend is smooth and
eigenvalue-floored by (1−w)·λ_min(g_inner). Empirically: with cylinder_depth
pushed to 0.6 (cube), det(g_outer) bottoms out at 3×10⁻⁸ at d = 0.5 (neatly
measuring the limit surface's focal radius) and re-grows, while the blend's
λ_min never dips below 0.017. The one true validity condition is immersedness
of the fitted YZ surface — unrelated to focal depth.
