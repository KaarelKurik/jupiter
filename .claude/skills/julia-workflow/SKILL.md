---
name: julia-workflow
description: How to run, test, render, and benchmark Julia code in this repo — warm daemon (scripts/jd) vs cold runs, quoting pitfalls, verification rules. Read before running any Julia here.
---

# Running Julia in this repo

## Golden rules

- Write `.jl` script files and run them; avoid `julia -e '...'` one-liners.
  Julia's `'` transpose operator fights shell quoting, and top-level loops
  assigning globals hit soft-scope errors outside the REPL. Scratch scripts go
  in the session scratchpad, tracked ones in `scripts/`.
- Inner-loop work (probes, quick benchmarks, iteration) → warm daemon via
  `scripts/jd` (~1s overhead instead of ~20s package-load + JIT per call).
- Anything that **certifies** a change runs cold, never through the daemon:
  the test suite, `scripts/physics_diff.jl`, and final renders.

## Warm daemon: scripts/jd

- `scripts/jd script.jl [args...]` — run in the warm process; first call
  cold-starts the daemon (~30s). Args must not contain spaces (DaemonMode
  joins them with spaces). Relative paths resolve against the caller's cwd.
- `scripts/jd --restart | --stop | --status` — the escape hatch. Restart on
  any suspicion of staleness; if weirdness survives a cold run, it's real.
- **Benchmarks and any number you'll act on: `--restart` first** (or run
  cold). 2026-07-11: Revise twice silently failed to apply an edit (a new
  function, then a redefined convenience constructor) without the die-loudly
  shim firing — the runs used stale code and produced wrong timings. Cheap
  insurance: a fresh daemon is ~30s; a wrong conclusion is not.
- Function edits in `src/` are picked up automatically (Revise). Struct or
  const redefinitions **cannot** be hot-applied: the daemon detects this,
  reports it, and kills itself — the next `jd` call cold-starts. That message
  after editing a struct is normal, not a failure.
- `Project.toml`/`Manifest.toml` changes auto-restart the daemon. jupiter's
  memoization caches are cleared before every run.

## Standard cold commands

- Tests: `julia --project -e 'using Pkg; Pkg.test()'` (~20s)
- Physics regression: `julia --project scripts/physics_diff.jl` — deterministic
  ray bundles vs committed baselines in `res/physics_diff/`; side flips are the
  hard signal. `--save` only after a *deliberate* physics change.
- Renders: `julia --project --threads=auto scripts/trefoil.jl 192x144` (or
  `scripts/first_light.jl`); PPM output in `out/`; view via
  `magick out/x.ppm out/x.png` then Read the PNG.

## Repo conventions (pointers)

- `plan-claude.md` in the repo root is the session resume point; read it first.
- VCS is jj, never bare git mutations; one described change per step.
- `res/` tracked inputs, `out/` ignored outputs, `gallery/` curated renders.
