//
//  DoubleDouble.swift
//  MandelbrotMetal
//
//  Created by Michael Stebel on 9/8/25.
//
//  Double-double (≈106-bit) arithmetic for CPU-side Swift.
//  Representation: x = hi + lo, with |lo| <= 0.5 ulp(hi)
//
//  Notes / Safety against regressions:
//  - We *must* use fma() to get error-free product residuals in TwoProd.
//  - Do not rewrite `fma(a,b,c)` as `(a*b)+c`; that breaks error accounting.
//  - The renormalization step ensures |lo| is small relative to hi.
//  - Subnormals: keep denormals enabled on CPU; flushing may lose a bit of robustness.
//  - Build with at least -O; inlining hints included.
//
//  If you port similar logic to Metal, keep fma(...) and renorm patterns identical.
//

import Foundation
import simd

@frozen
public struct DD: Sendable, Equatable {
    @usableFromInline var v: simd_double2   // v.x = hi, v.y = lo

    @inlinable public var hi: Double { v.x }
    @inlinable public var lo: Double { v.y }

    // MARK: - Inits
    @inlinable public init(_ hi: Double, _ lo: Double) {
        self.v = simd_double2(hi, lo)
        self = DD.renorm(self) // ensure canonical form
    }

    @inlinable public init(_ x: Double) {
        self.v = simd_double2(x, 0.0)
    }

    // MARK: - Constants
    public static let zero = DD(0.0)
    public static let one  = DD(1.0)
    public static let two  = DD(2.0)

    // MARK: - Core error-free transforms

    /// Knuth/Dekker TwoSum: returns s,e with a+b = s+e exactly (assuming RN)
    @inlinable @inline(__always)
    static func twoSum(_ a: Double, _ b: Double) -> (s: Double, e: Double) {
        let s = a + b
        // Recover roundoff: (s - a) cancels high bits of s leaving contribution of b
        let bb = s - a
        let e = (a - (s - bb)) + (b - bb)
        return (s, e)
    }

    /// FastTwoSum when |a| >= |b| (slightly faster); caller must ensure precondition.
    @inlinable @inline(__always)
    static func fastTwoSum(_ a: Double, _ b: Double) -> (s: Double, e: Double) {
        let s = a + b
        let e = b - (s - a)
        return (s, e)
    }

    /// Error-free product: p = a*b (rounded), e = exact residual so that a*b = p + e
    @inlinable @inline(__always)
    static func twoProd(_ a: Double, _ b: Double) -> (p: Double, e: Double) {
        let p = a * b
        // FMA gives exact (a*b - p) in one rounded step
        let e = fma(a, b, -p)
        return (p, e)
    }

    // MARK: - Renormalization

    /// Normalize (hi, lo) so that hi is the rounded sum and |lo| <= 0.5 ulp(hi)
    @inlinable @inline(__always)
    static func renorm(_ x: DD) -> DD {
        // Prefer fastTwoSum if |hi| >= |lo|; otherwise fall back to general twoSum.
        if abs(x.hi) >= abs(x.lo) {
            let (h, l) = fastTwoSum(x.hi, x.lo)
            return DD.unchecked(h, l)
        } else {
            let (h, l) = twoSum(x.hi, x.lo)
            return DD.unchecked(h, l)
        }
    }

    /// Construct without renormalization (internal use).
    @inlinable @inline(__always)
    static func unchecked(_ hi: Double, _ lo: Double) -> DD {
        var r = DD(0.0)
        r.v = simd_double2(hi, lo)
        return r
    }

    // MARK: - Basic ops

    /// x + y (double-double)
    @inlinable @inline(__always)
    public static func + (x: DD, y: DD) -> DD {
        // Add high parts, capture residual
        let (s, e) = twoSum(x.hi, y.hi)
        // Add low parts, fold in carefully
        let t = x.lo + y.lo
        let (s2, e2) = twoSum(s, t)
        // Combine all errors
        let z = unchecked(s2, e + e2)
        return renorm(z)
    }

    /// x - y
    @inlinable @inline(__always)
    public static func - (x: DD, y: DD) -> DD {
        return x + (-y)
    }

    /// Unary minus
    @inlinable @inline(__always)
    public static prefix func - (x: DD) -> DD {
        return unchecked(-x.hi, -x.lo)
    }

    /// x * y (double-double)
    @inlinable @inline(__always)
    public static func * (x: DD, y: DD) -> DD {
        // Main product and its rounding error
        let (p, pe) = twoProd(x.hi, y.hi)

        // Cross terms (smaller in magnitude)
        let c = x.hi * y.lo + x.lo * y.hi

        // Accumulate main + cross with error tracking
        let (s, e2) = twoSum(p, c)

        // Add all small pieces: product residual + cross residual + lo*lo
        let loAll = pe + e2 + (x.lo * y.lo)

        // Final renormalization
        let (hi2, lo2) = twoSum(s, loAll)
        return renorm(unchecked(hi2, lo2))
    }

    /// Square (slightly cheaper than x*x)
    @inlinable @inline(__always)
    public func squared() -> DD {
        // (hi + lo)^2 = hi^2 + 2*hi*lo + lo^2
        let (p, pe) = DD.twoProd(self.hi, self.hi)
        let cross = 2.0 * self.hi * self.lo
        let (s, e2) = DD.twoSum(p, cross)
        let loAll = pe + e2 + (self.lo * self.lo)
        let (hi2, lo2) = DD.twoSum(s, loAll)
        return DD.renorm(DD.unchecked(hi2, lo2))
    }

    /// x / y using one Newton refinement (good to ~106 bits)
    @inlinable @inline(__always)
    public static func / (x: DD, y: DD) -> DD {
        // Initial double quotient
        let q0 = x.hi / y.hi

        // r = x - y*q0 (DD-accurate)
        let r = x - (y * DD(q0))

        // Correction term ≈ (r.hi + r.lo)/y.hi (double is fine; then fold in)
        let corr = (r.hi + r.lo) / y.hi

        // q = q0 + corr, normalize
        let q = DD(q0, 0.0) + DD(corr, 0.0)
        return renorm(q)
    }

    // MARK: - Mixed ops

    @inlinable @inline(__always)
    public static func + (x: DD, a: Double) -> DD {
        let (s, e) = twoSum(x.hi, a)
        let (hi2, lo2) = twoSum(s, x.lo + e)
        return renorm(unchecked(hi2, lo2))
    }

    @inlinable @inline(__always)
    public static func - (x: DD, a: Double) -> DD { x + (-a) }

    @inlinable @inline(__always)
    public static func * (x: DD, a: Double) -> DD {
        let (p, pe) = twoProd(x.hi, a)
        let loAll = pe + x.lo * a
        let (hi2, lo2) = twoSum(p, loAll)
        return renorm(unchecked(hi2, lo2))
    }

    @inlinable @inline(__always)
    public static func / (x: DD, a: Double) -> DD {
        // Note: this is safe since ‘a’ is a double; we normalize result
        let q0 = x.hi / a
        // r = x - a*q0
        let r = x - DD(q0) * a
        let corr = (r.hi + r.lo) / a
        return renorm(DD(q0) + DD(corr))
    }

    // MARK: - Utilities

    /// ulp of the high part (for diagnostics)
    @inlinable public var ulpHi: Double { Double.ulpOfOne * abs(hi) }

    /// Convert back to Double (rounded)
    @inlinable public var asDouble: Double { hi + lo }

    /// Exact-ish string (for debugging / logging)
    public var description: String {
        // Compact, avoids scientific unless large/small
        return String(format: "DD(hi=%0.17g, lo=%0.17g)", hi, lo)
    }
}

// MARK: - Quick self-checks (leave in Debug builds)
#if DEBUG
@usableFromInline
func __dd_self_test() {
    // Addition vs. plain double
    let a = DD(1.0e16) + DD(1.0)   // standard double would lose the +1 here
    assert(a.asDouble - 1.0e16 == 1.0, "DD add lost low bit")

    // Multiplication sanity
    let x = DD(1.234567890123456)  // ~53-bit payload
    let y = DD(9.876543210987654)
    let z = x * y
    let ref = x.asDouble * y.asDouble
    // DD should be within a handful of double ulps of the *true* 106-bit product projection
    precondition(abs(z.asDouble - ref) <= max(1.0, abs(ref)) * 1e-30, "DD mul suspicious")

    // Division round-trip
    let q = z / y
    let err = abs((q.asDouble - x.asDouble) / x.asDouble)
    precondition(err < 1e-30, "DD div too imprecise")
}
#endif
