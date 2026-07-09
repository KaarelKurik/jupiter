# TaylorDiff v0.3.5 — bugs and limitations found in jupiter

Found 2026-07-08 while computing Christoffel symbols via nested TaylorDiff
(TaylorScalar-of-TaylorScalar: outer AD over a metric whose entries are themselves
computed with inner AD). Recorded before migrating the codebase to ForwardDiff, so
the findings survive the migration. All line references are to
`~/.julia/packages/TaylorDiff/F2JtC/` (v0.3.5).

## Bug 1 (silent!): 2-arg `atan` drops derivatives

`atan(y::TaylorScalar, x::TaylorScalar)` returns a plain `Float64` — the primal
value with all derivative information silently discarded (no error). Presumably it
falls through `isinf/isnan`-style "operate on value, drop partials" definitions or
a Real-fallback path.

Danger class: worst possible — differentiated code containing 2-arg atan produces
*wrong derivatives*, not exceptions.

Workaround used here: `safe_atan2` built from 1-arg `atan` via half-angle forms,
branching on the primal sign of x:

```julia
safe_atan2(y, x) = primal(x) >= 0 ? 2 * atan(y / (sqrt(x^2 + y^2) + x)) :
                                    2 * atan((sqrt(x^2 + y^2) - x) / y)
```

(first form is 0/0 near the negative real axis, second near the positive one;
the two together cover everything but the origin).

Proper fix idea: define the binary primitive `atan(y, x)` with the standard
∂atan2 rules in TaylorDiff's primitive table.

## Bug 2 (upstream-worthy): nested constant construction returns the wrong type

`scalar.jl:23`:

```julia
TaylorScalar{T, P}(x) where {T, P} = TaylorScalar{P}(T(x))
```

For nested `T <: TaylorScalar` (e.g. `TaylorScalar{TaylorScalar{Float64,1},1}(2)`),
`T(x)` produces an inner-type scalar, and `TaylorScalar{P}(::TaylorScalar)` then
dispatches to the **truncate/extend** method (`scalar.jl:44`) instead of the
"constant with zero partials" method — so the constructor returns the *inner* type
instead of the requested nested type. Manifests as
`TypeError: in typeassert, expected TaylorScalar{TaylorScalar{Float64,1},1}, got
TaylorScalar{Float64,1}` inside `convert` during arithmetic promotion.

Fix (what we patched locally):

```julia
TaylorScalar{T, P}(x::Number) where {T <: TaylorScalar, P} =
    TaylorScalar(convert(T, x), ntuple(i -> zero(T), Val(P)))
```

## Bug 3: constructor ambiguity with Core's identity constructor

`TaylorScalar{T,P}(x::TaylorScalar{T,P})` (identity call, arises inside primitive
rules, e.g. `atan` at the nesting boundary) is ambiguous between the generic
constructor above and Core's `(::Type{T})(x::T) where T<:Number`. Fix — the
identity method Julia's error message suggests, plus a doubly-constrained variant
so it stays unambiguous against the Bug-2 fix:

```julia
TaylorScalar{T, P}(x::TaylorScalar{T, P}) where {T, P} = x
TaylorScalar{T, P}(x::TaylorScalar{T, P}) where {T <: TaylorScalar, P} = x
```

## Limitation: `make_seed` requires seed eltype == point eltype

`derivative(f, x::Vector{TaylorScalar}, l::Vector{Float64}, Val(1))` throws
`MethodError: no method matching make_seed(...)` — `make_seed(::A, ::A, ...)` wants
both arguments of identical array type. Nested use therefore requires building
direction vectors in `eltype(x)` (our `basis_direction` helper did
`convert(eltype(x), ...)` per component). Fix idea: promote seed eltype in
`make_seed`.

## Limitation: `floor` (and friends) missing

`floor(::TaylorScalar)` is a MethodError. Locally-constant functions should
arguably return the floor of the primal (derivative zero a.e.). We routed branch
selection through a `primal` extractor instead — which is the safer pattern anyway.

## Limitation: no complex support

Complex powers of TaylorScalar-carrying complex numbers don't differentiate
(motivated `fake_complex_pow`, a principal-branch complex power from real
primitives — kept regardless of AD library, since ForwardDiff duals are also
`<: Real`).

## General lesson

Single-level TaylorDiff was solid; every problem above except the atan bug lives on
the nested path, where we found two genuine bugs within an hour. Robustness tracks
traffic on the specific code path. Bugs 1 and 2 are worth reporting upstream
(https://github.com/JuliaDiff/TaylorDiff.jl).
