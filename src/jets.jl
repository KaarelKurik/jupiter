# closed-form derivative jets — the hand-rolled register lever (2026-07-16).
# christoffel needs s..∂³s of the blended surface in (u,v); nested ForwardDiff
# duals haul 48 carrier floats per scalar, mostly structural zeros, which is
# what spills on device (measurements.md 2026-07-14e). These jets carry exactly
# the 10 distinct partials. ad.jl stays the AD contact surface and reference.jl
# the oracle; chart.jl/geodesic.jl assemble surface_jet/christoffel from these.
#
# Value lanes throughout reproduce the plain evaluation's arithmetic op-for-op
# (same divisions, same addition order), so values stay bit-identical to the
# production evaluators; only derivative lanes are new arithmetic.

# order-k truncated Taylor jets in the two chart variables (u,v). V is the
# coefficient: a scalar or an SVector{3} (a vector jet is three scalar jets).
struct Jet1{V}
    f::V
    fu::V
    fv::V
end

struct Jet2{V}
    f::V
    fu::V
    fv::V
    fuu::V
    fuv::V
    fvv::V
end

struct Jet3{V}
    f::V
    fu::V
    fv::V
    fuu::V
    fuv::V
    fvv::V
    fuuu::V
    fuuv::V
    fuvv::V
    fvvv::V
end

@inline jadd(a::Jet3, b::Jet3) = Jet3(a.f + b.f, a.fu + b.fu, a.fv + b.fv,
                              a.fuu + b.fuu, a.fuv + b.fuv, a.fvv + b.fvv,
                              a.fuuu + b.fuuu, a.fuuv + b.fuuv, a.fuvv + b.fuvv, a.fvvv + b.fvvv)

# ∂u/∂v as jets one order down: how a jet's own derivative enters chain rules
@inline jdu(a::Jet3) = Jet2(a.fu, a.fuu, a.fuv, a.fuuu, a.fuuv, a.fuvv)
@inline jdv(a::Jet3) = Jet2(a.fv, a.fuv, a.fvv, a.fuuv, a.fuvv, a.fvvv)
@inline jdu(a::Jet2) = Jet1(a.fu, a.fuu, a.fuv)
@inline jdv(a::Jet2) = Jet1(a.fv, a.fuv, a.fvv)

# Leibniz rule for any bilinear op (*, dot, cross): the one product rule,
# parameterized by the pairing
@inline leibniz(op, a::Jet1, b::Jet1) = Jet1(
    op(a.f, b.f),
    op(a.fu, b.f) + op(a.f, b.fu),
    op(a.fv, b.f) + op(a.f, b.fv))

@inline leibniz(op, a::Jet2, b::Jet2) = Jet2(
    op(a.f, b.f),
    op(a.fu, b.f) + op(a.f, b.fu),
    op(a.fv, b.f) + op(a.f, b.fv),
    op(a.fuu, b.f) + 2 * op(a.fu, b.fu) + op(a.f, b.fuu),
    op(a.fuv, b.f) + op(a.fu, b.fv) + op(a.fv, b.fu) + op(a.f, b.fuv),
    op(a.fvv, b.f) + 2 * op(a.fv, b.fv) + op(a.f, b.fvv))

@inline leibniz(op, a::Jet3, b::Jet3) = Jet3(
    op(a.f, b.f),
    op(a.fu, b.f) + op(a.f, b.fu),
    op(a.fv, b.f) + op(a.f, b.fv),
    op(a.fuu, b.f) + 2 * op(a.fu, b.fu) + op(a.f, b.fuu),
    op(a.fuv, b.f) + op(a.fu, b.fv) + op(a.fv, b.fu) + op(a.f, b.fuv),
    op(a.fvv, b.f) + 2 * op(a.fv, b.fv) + op(a.f, b.fvv),
    op(a.fuuu, b.f) + 3 * op(a.fuu, b.fu) + 3 * op(a.fu, b.fuu) + op(a.f, b.fuuu),
    op(a.fuuv, b.f) + op(a.fuu, b.fv) + 2 * op(a.fuv, b.fu) + 2 * op(a.fu, b.fuv) + op(a.fv, b.fuu) + op(a.f, b.fuuv),
    op(a.fuvv, b.f) + op(a.fvv, b.fu) + 2 * op(a.fuv, b.fv) + 2 * op(a.fv, b.fuv) + op(a.fu, b.fvv) + op(a.f, b.fuvv),
    op(a.fvvv, b.f) + 3 * op(a.fvv, b.fv) + 3 * op(a.fv, b.fvv) + op(a.f, b.fvvv))

# quotient a/b (b scalar-valued), solved order by order from a = q·b; the value
# lane is the same division the plain evaluator performs
@inline function jdiv(a::Jet2, b::Jet2)
    q = a.f / b.f
    qu = (a.fu - q * b.fu) / b.f
    qv = (a.fv - q * b.fv) / b.f
    Jet2(q, qu, qv,
         (a.fuu - 2 * (qu * b.fu) - q * b.fuu) / b.f,
         (a.fuv - qu * b.fv - qv * b.fu - q * b.fuv) / b.f,
         (a.fvv - 2 * (qv * b.fv) - q * b.fvv) / b.f)
end

@inline function jdiv(a::Jet3, b::Jet3)
    q = a.f / b.f
    qu = (a.fu - q * b.fu) / b.f
    qv = (a.fv - q * b.fv) / b.f
    quu = (a.fuu - 2 * (qu * b.fu) - q * b.fuu) / b.f
    quv = (a.fuv - qu * b.fv - qv * b.fu - q * b.fuv) / b.f
    qvv = (a.fvv - 2 * (qv * b.fv) - q * b.fvv) / b.f
    Jet3(q, qu, qv, quu, quv, qvv,
         (a.fuuu - 3 * (quu * b.fu) - 3 * (qu * b.fuu) - q * b.fuuu) / b.f,
         (a.fuuv - quu * b.fv - 2 * (quv * b.fu) - 2 * (qu * b.fuv) - qv * b.fuu - q * b.fuuv) / b.f,
         (a.fuvv - qvv * b.fu - 2 * (quv * b.fv) - 2 * (qv * b.fuv) - qu * b.fvv - q * b.fuvv) / b.f,
         (a.fvvv - 3 * (qvv * b.fv) - 3 * (qv * b.fvv) - q * b.fvvv) / b.f)
end

# scalar chain rule (Faà di Bruno to order 3): outer 1D f given by its
# derivatives f1..f3 at g.f, inner a (u,v) jet
@inline compose1(f0, f1, f2, g::Jet2) = Jet2(
    f0, f1 * g.fu, f1 * g.fv,
    f2 * g.fu^2 + f1 * g.fuu,
    f2 * (g.fu * g.fv) + f1 * g.fuv,
    f2 * g.fv^2 + f1 * g.fvv)

@inline compose1(f0, f1, f2, f3, g::Jet3) = Jet3(
    f0, f1 * g.fu, f1 * g.fv,
    f2 * g.fu^2 + f1 * g.fuu,
    f2 * (g.fu * g.fv) + f1 * g.fuv,
    f2 * g.fv^2 + f1 * g.fvv,
    f3 * g.fu^3 + 3 * (f2 * (g.fu * g.fuu)) + f1 * g.fuuu,
    f3 * (g.fu^2 * g.fv) + f2 * (g.fuu * g.fv + 2 * (g.fu * g.fuv)) + f1 * g.fuuv,
    f3 * (g.fu * g.fv^2) + f2 * (g.fvv * g.fu + 2 * (g.fv * g.fuv)) + f1 * g.fuvv,
    f3 * g.fv^3 + 3 * (f2 * (g.fv * g.fvv)) + f1 * g.fvvv)

# ---- complex jets: the wedge chain is holomorphic in z = u + iv, so its
# order-3 jet is one value and three complex z-derivatives, and every chain
# rule is 1D. Real (u,v) partials fall out by Cauchy–Riemann: ∂u ↦ f⁽ᵏ⁾,
# ∂v ↦ iᵏ-weighted f⁽ᵏ⁾.

struct CJet{T}
    w0::Complex{T}
    w1::Complex{T}
    w2::Complex{T}
    w3::Complex{T}
end

@inline cjet_identity(uv) = CJet(complex(uv[1], uv[2]), complex(one(uv[1]), zero(uv[1])),
                         complex(zero(uv[1]), zero(uv[1])), complex(zero(uv[1]), zero(uv[1])))

@inline cjet_scale(a, j::CJet) = CJet(a * j.w0, a * j.w1, a * j.w2, a * j.w3)
@inline cjet_rdiv(j::CJet, a) = CJet(j.w0 / a, j.w1 / a, j.w2 / a, j.w3 / a)

@inline re_jet(j::CJet) = Jet3(real(j.w0),
                       real(j.w1), -imag(j.w1),
                       real(j.w2), -imag(j.w2), -real(j.w2),
                       real(j.w3), -imag(j.w3), -real(j.w3), imag(j.w3))

@inline im_jet(j::CJet) = Jet3(imag(j.w0),
                       imag(j.w1), real(j.w1),
                       imag(j.w2), real(j.w2), -imag(j.w2),
                       imag(j.w3), real(j.w3), -imag(j.w3), -real(j.w3))

"""
principal power ζ^p as a complex jet: the value is fake_complex_pow's exact
arithmetic (chart.jl); derivatives are the recurrence φ⁽ᵏ⁺¹⁾ = (p−k)·φ⁽ᵏ⁾/ζ,
chained with the incoming jet by 1D Faà di Bruno. Zero at the origin like
fake_complex_pow (same guard, derivative lanes zero there too).
"""
@inline function cpow_jet(j::CJet, power)
    z0 = j.w0
    if iszero(primal(real(z0))) && iszero(primal(imag(z0)))
        zc = zero(z0)
        return CJet(zc, zc, zc, zc)
    end
    val = fake_complex_pow(SVector(real(z0), imag(z0)), power)
    φ0 = complex(val[1], val[2])
    p = carrier(real(z0))(power)
    φ1 = p * φ0 / z0
    φ2 = (p - 1) * φ1 / z0
    φ3 = (p - 2) * φ2 / z0
    CJet(φ0,
         φ1 * j.w1,
         φ2 * j.w1^2 + φ1 * j.w2,
         φ3 * j.w1^3 + 3 * (φ2 * (j.w1 * j.w2)) + φ1 * j.w3)
end

"""
order-3 jet of real-valued p∘ζ with ζ holomorphic: p's ten (x,y) partials at
ζ's value plus ζ's three complex derivatives. In Wirtinger form ∂z acts only
through w and ∂z̄ only through w̄, so the bivariate chain rule factorizes into
two 1D chains; p real makes D_ab = conj(D_ba), leaving five D terms.
"""
@inline function wirtinger_compose(pd, z1, z2, z3)
    f, fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy = pd
    P10 = complex(fx, -fy) / 2
    P20 = complex(fxx - fyy, -2 * fxy) / 4
    P11 = (fxx + fyy) / 4
    P30 = complex(fxxx - 3 * fxyy, fyyy - 3 * fxxy) / 8
    P21 = complex(fxxx + fxyy, -(fxxy + fyyy)) / 8
    z1sq = z1 * z1
    D10 = P10 * z1
    D20 = P20 * z1sq + P10 * z2
    D30 = P30 * (z1sq * z1) + 3 * (P20 * (z1 * z2)) + P10 * z3
    D11 = P11 * abs2(z1)
    D21 = P21 * (z1sq * conj(z1)) + P11 * (z2 * conj(z1))
    Jet3(f,
         2 * real(D10), -2 * imag(D10),
         2 * real(D20) + 2 * D11, -2 * imag(D20), -2 * real(D20) + 2 * D11,
         2 * real(D30) + 6 * real(D21),
         -2 * (imag(D30) + imag(D21)),
         -2 * real(D30) + 2 * real(D21),
         2 * imag(D30) - 6 * imag(D21))
end

# ---- scalar 1D chains (blend, flat bump) as (f, f', f'', f''') tuples

"""
flat_bump and its first three derivatives: exp(-1/x)·(rational in 1/x),
closed form. Gated on the value: once exp underflows (or x ≤ 0) every
derivative is zero too — and the gate keeps 0·Inf out of Float32 near seams,
where 1/x powers overflow long before the product stops being zero.
"""
@inline function flat_bump_jet(x)
    z = zero(x)
    primal(x) > 0 || return (z, z, z, z)
    b = exp(-1 / x)
    iszero(b) && return (b, z, z, z)
    w = 1 / x
    (b,
     b * w^2,
     b * (w^3 * (w - 2)),
     b * (w^4 * ((w - 6) * w + 6)))
end

"""
blend_scalar and its first three derivatives: fb(1−x)/(fb(x)+fb(1−x)) with the
(1−x) sign flips folded in, then the 1D quotient solved order by order. The
value lane is blend_scalar's exact arithmetic.
"""
@inline function blend_jet(x)
    a = flat_bump_jet(1 - x)
    b = flat_bump_jet(x)
    n1 = a[1]; n2 = -a[2]; n3 = a[3]; n4 = -a[4]
    d1 = b[1] + a[1]; d2 = b[2] - a[2]; d3 = b[3] + a[3]; d4 = b[4] - a[4]
    q0 = n1 / d1
    q1 = (n2 - q0 * d2) / d1
    q2 = (n3 - 2 * (q1 * d2) - q0 * d3) / d1
    q3 = (n4 - 3 * (q2 * d2) - 3 * (q1 * d3) - q0 * d4) / d1
    (q0, q1, q2, q3)
end

"""
value and all (x,y) partials to order 3 of one packed chart polynomial, per
component, in a single table read: double Horner where each accumulator drags
its derivatives along (Aₖ₊₁ = Aₖ·x + c ⇒ A′ₖ₊₁ = A′ₖ·x + Aₖ, A″ₖ₊₁ = A″ₖ·x + 2A′ₖ, …).
Returns three (f, fx, fy, fxx, fxy, fyy, fxxx, fxxy, fxyy, fyyy) tuples; the
f lane is eval_packed's exact arithmetic.
"""
@inline function eval_packed_partials_component(packed, component, x, y)
    # own @inline function, not an ntuple do-block: the lowered closure is too
    # big for the inliner's own judgment, and as an ABI call its 10-lane return
    # round-trips through device local memory per corner (census 2026-07-19)
    z = zero(x) * zero(y)
    m00 = z; m10 = z; m20 = z; m30 = z
    m01 = z; m11 = z; m21 = z
    m02 = z; m12 = z; m03 = z
    for i in size(packed, 2):-1:1
        r0 = z; r1 = z; r2 = z; r3 = z
        for j in size(packed, 3):-1:1
            r3 = r3 * y + 3 * r2
            r2 = r2 * y + 2 * r1
            r1 = r1 * y + r0
            r0 = r0 * y + packed[component, i, j]
        end
        m30 = m30 * x + 3 * m20
        m21 = m21 * x + 2 * m11
        m20 = m20 * x + 2 * m10
        m12 = m12 * x + m02
        m11 = m11 * x + m01
        m10 = m10 * x + m00
        m03 = m03 * x + r3
        m02 = m02 * x + r2
        m01 = m01 * x + r1
        m00 = m00 * x + r0
    end
    (m00, m10, m01, m20, m11, m02, m30, m21, m12, m03)
end

@inline function eval_packed_partials(packed::AbstractArray{<:AbstractFloat, 3}, x, y)
    ntuple(component -> eval_packed_partials_component(packed, component, x, y), Val(3))
end

# three scalar jets to one vector jet (surface components back into SVectors)
@inline vjet3(a::Jet3, b::Jet3, c::Jet3) = Jet3(
    SVector(a.f, b.f, c.f), SVector(a.fu, b.fu, c.fu), SVector(a.fv, b.fv, c.fv),
    SVector(a.fuu, b.fuu, c.fuu), SVector(a.fuv, b.fuv, c.fuv), SVector(a.fvv, b.fvv, c.fvv),
    SVector(a.fuuu, b.fuuu, c.fuuu), SVector(a.fuuv, b.fuuv, c.fuuv),
    SVector(a.fuvv, b.fuvv, c.fuvv), SVector(a.fvvv, b.fvvv, c.fvvv))
