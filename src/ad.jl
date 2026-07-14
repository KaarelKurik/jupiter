# the AD contact surface — the only file that talks to the AD library; see
# taylordiff-bugs.md for why we left TaylorDiff. Swap the backend here if
# ForwardDiff ever disappoints.

primal(x) = x
primal(d::ForwardDiff.Dual) = primal(ForwardDiff.value(d)) # recurse: christoffel nests duals

# the innermost scalar type through any dual nesting: what Float64 constants
# convert to so they scale with the phase's precision instead of promoting it
carrier(x) = typeof(primal(x))

directional(f, x, v) = ForwardDiff.derivative(t -> f(x .+ t .* v), zero(eltype(x))) # seed in x's eltype or the dual pass promotes everything to Float64

# The two dual-pass primitives below must keep separate bodies: when one pass
# runs inside another (jacobian_columns → collar → the surface pass), any
# shared method would meet itself on inference's stack with a grown dual-type
# signature, and the termination heuristic widens that call to Any, boxing the
# entire geodesic step loop (measurements.md 2026-07-12). Nesting levels must
# be distinct methods — the same reason ForwardDiff's hessian is jacobian over
# gradient rather than jacobian twice.

"""
all coordinate directional derivatives of f at x in one dual pass; returns an
SMatrix of columns [df/dx_1 ... df/dx_n] like the hcat-of-directionals it
replaces. N must equal length(x) (passed as Val so everything stack-allocates).
"""
function jacobian_columns(f, x, ::Val{N}) where {N}
    tag = typeof(ForwardDiff.Tag(f, eltype(x)))
    seeded = SVector{N}(ntuple(i -> ForwardDiff.Dual{tag}(x[i], ntuple(j -> eltype(x)(i == j), Val(N))), Val(N)))
    fd = f(seeded)
    hcat(ntuple(j -> ForwardDiff.partials.(fd, j), Val(N))...)
end

"""
f's value together with its jacobian columns, from the same dual pass — the
value lane performs exactly the plain evaluation's arithmetic, so the value is
bit-identical to f(x) and free. This is the pass to use under jacobian_columns
(see the inference note above); that it is a different operation, not a copy,
is what keeps the nesting honest.
"""
function value_and_jacobian_columns(f, x, ::Val{N}) where {N}
    tag = typeof(ForwardDiff.Tag(f, eltype(x)))
    seeded = SVector{N}(ntuple(i -> ForwardDiff.Dual{tag}(x[i], ntuple(j -> eltype(x)(i == j), Val(N))), Val(N)))
    fd = f(seeded)
    (ForwardDiff.value.(fd), hcat(ntuple(j -> ForwardDiff.partials.(fd, j), Val(N))...))
end

# f(pos)-shaped view of f(env, chart, pos); Fix1 (not a closure) so the callable's
# type stays concrete for AD tags and inference in the hot paths
situate(f, env, chart) = Base.Fix1(Base.Fix1(f, env), chart)

generic_normalize(x) = x ./ sqrt(sum(abs2, x)) # LinearAlgebra.normalize has scaling branches unfriendly to dual numbers
