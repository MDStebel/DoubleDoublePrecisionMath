//  DoubleDoublePrecisionMath.swift
//  Mandelbrot Metal
//
//  Created by Michael Stebel on 8/8/25.
//  Updated on 5/15/25.
//

import Foundation

// MARK: - Double-Double Arithmetic for High-Precision Computing

/// A double-double precision floating-point number representation.
///
/// Represents a high-precision real value as the unevaluated sum of two
/// `Double` values: `x ≈ hi + lo`, where `hi` holds the leading bits and
/// `lo` holds the trailing error term.
///
/// This gives ~106 bits of precision (~31–32 decimal digits), which is
/// useful for extreme zoom levels in Mandelbrot Metal where FP64 is not
/// sufficient.
internal struct DD: Equatable, CustomStringConvertible {
    /// Leading (high) part of the value.
    var hi: Double
    /// Trailing (low) part, storing the rounding error of `hi`.
    var lo: Double
    
    /// Initializes a `DD` from a regular `Double`, treating it as exact
    /// (all bits go into `hi`, and `lo` is set to zero).
    @inline(__always)
    init(_ x: Double = 0.0) {
        self.hi = x
        self.lo = 0.0
    }
    
    /// Initializes a `DD` from explicit high and low components.
    /// Use this when you already have a compensated representation.
    @inline(__always)
    init(hi: Double, lo: Double) {
        self.hi = hi
        self.lo = lo
    }

    /// Debug-friendly representation that keeps both components visible.
    var description: String {
        "DD(hi: \(hi), lo: \(lo))"
    }
}

// MARK: - Low-Level Error-Free Transformations

/// Computes the exact sum of two finite `Double` values when the rounded sum
/// is finite, returning both the rounded sum and the rounding error.
///
/// This is Knuth's *two-sum* algorithm, an error-free transformation:
///   - `s` is the correctly rounded sum `a + b` in `Double`.
///   - `e` is the residual such that `a + b = s + e` exactly.
///
/// Non-finite or overflowing sums return a zero residual instead of attempting
/// an error transform on infinities.
///
/// This forms the foundation of double-double addition and normalization.
@inline(__always)
internal func twoSum(_ a: Double, _ b: Double) -> (Double, Double) {
    let s  = a + b
    guard s.isFinite else { return (s, 0.0) }

    let bb = s - a
    let e  = (a - (s - bb)) + (b - bb)
    return (s, e)
}

/// A faster variant of `twoSum` used when we know that `|a| >= |b|`.
///
/// This is the standard *quick-two-sum* algorithm:
///   - `s` is the correctly rounded sum `a + b`.
///   - `e` is the small residual such that `a + b = s + e` exactly.
///
/// Non-finite or overflowing sums return a zero residual instead of attempting
/// an error transform on infinities.
///
/// Precondition: `abs(a) >= abs(b)` must hold for the error formula to
/// be stable. In this file we only use `quickTwoSum` in contexts where
/// `a` is the large leading term and `b` is a small correction.
@inline(__always)
internal func quickTwoSum(_ a: Double, _ b: Double) -> (Double, Double) {
    let s = a + b
    guard s.isFinite else { return (s, 0.0) }

    let e = b - (s - a)
    return (s, e)
}

/// Computes the exact product of two finite `Double` values when the rounded
/// product is finite, returning both the rounded product and the rounding error.
///
/// When the rounded product is finite, `fma(a, b, -p)` gives the residual
/// of `a * b - p` in one rounded step, so:
///   - `p` is the correctly rounded product `a * b`.
///   - `e` is the residual such that `a * b = p + e` exactly.
///
/// Non-finite or overflowing products return a zero residual instead of
/// producing a NaN from `fma(a, b, -p)`.
@inline(__always)
internal func twoProd(_ a: Double, _ b: Double) -> (Double, Double) {
    let p = a * b
    guard p.isFinite else { return (p, 0.0) }

    let e = fma(a, b, -p)
    return (p, e)
}

/// Normalizes a high/low pair into a canonical double-double value.
///
/// Unlike `quickTwoSum`, this does not assume the high component already
/// dominates the low component, which keeps cancellation-heavy results safe.
@inline(__always)
internal func ddRenormalize(_ hi: Double, _ lo: Double) -> DD {
    let (normalizedHi, normalizedLo) = twoSum(hi, lo)
    return DD(hi: normalizedHi, lo: normalizedLo)
}

// MARK: - Double-Double Operations

/// Adds two double-double numbers with extended precision.
///
/// We:
///  1. Use `twoSum` on the high parts to obtain a rounded sum `s` and
///     an error term `e1`.
///  2. Use `twoSum` on the low parts so their own residual is not discarded.
///  3. Fold the high- and low-part residuals together.
///  4. Renormalize without relying on `quickTwoSum`'s magnitude precondition.
///
/// This is a slightly more conservative DD addition than the common fast
/// "sloppy add" pattern, which helps when the high parts nearly cancel.
@inline(__always)
internal func ddAdd(_ x: DD, _ y: DD) -> DD {
    // Sum leading parts and track the rounding error.
    let (s, e1) = twoSum(x.hi, y.hi)
    
    // Sum lower parts separately so their residual can still contribute.
    let (t, e2) = twoSum(x.lo, y.lo)
    
    // Fold the high-part residual and low-part sum together before the final
    // two-step normalization.
    let (correction, e3) = twoSum(e1, t)
    let partial = ddRenormalize(s, correction)
    return ddRenormalize(partial.hi, partial.lo + e2 + e3)
}

/// Multiplies two double-double numbers with extended precision.
///
/// We:
///  1. Use `twoProd` on the high parts to obtain `p ≈ x.hi * y.hi` and
///     its error `e1`.
///  2. Add the cross terms `x.hi * y.lo + x.lo * y.hi` into `e2`.
///  3. Renormalize `p + e1 + e2` into `(hi, lo)`.
///
/// This yields a full double-double product with ~106 bits of precision,
/// suitable for extreme Mandelbrot iterations.
@inline(__always)
internal func ddMul(_ x: DD, _ y: DD) -> DD {
    // Product of high parts plus its error term.
    let (p, e1)  = twoProd(x.hi, y.hi)
    
    // Cross terms (the low*low term is small enough to ignore here).
    let e2       = x.hi * y.lo + x.lo * y.hi
    
    // Normalize the accumulated product.
    return ddRenormalize(p, e1 + e2)
}

/// Subtracts two double-double numbers with extended precision.
///
/// Implemented in terms of `ddAdd` by negating both components of `y`.
/// The optimizer typically eliminates the temporary `DD` in hot paths.
@inline(__always)
internal func ddSub(_ x: DD, _ y: DD) -> DD {
    ddAdd(x, DD(hi: -y.hi, lo: -y.lo))
}

// MARK: - DD Convenience Operators & Helpers

extension DD {
    // MARK: Basic arithmetic operators (DD × DD)
    /// Adds two `DD` values using `ddAdd`.
    @inline(__always)
    static func + (lhs: DD, rhs: DD) -> DD {
        ddAdd(lhs, rhs)
    }

    /// Subtracts two `DD` values using `ddSub`.
    @inline(__always)
    static func - (lhs: DD, rhs: DD) -> DD {
        ddSub(lhs, rhs)
    }

    /// Multiplies two `DD` values using `ddMul`.
    @inline(__always)
    static func * (lhs: DD, rhs: DD) -> DD {
        ddMul(lhs, rhs)
    }

    // MARK: Mixed arithmetic with `Double`
    /// Adds a `Double` to a `DD` by promoting the `Double` to `DD`.
    @inline(__always)
    static func + (lhs: DD, rhs: Double) -> DD {
        ddAdd(lhs, DD(rhs))
    }

    /// Adds a `DD` to a `Double` by promoting the `Double` to `DD`.
    @inline(__always)
    static func + (lhs: Double, rhs: DD) -> DD {
        ddAdd(DD(lhs), rhs)
    }

    /// Subtracts a `Double` from a `DD` by promoting the `Double` to `DD`.
    @inline(__always)
    static func - (lhs: DD, rhs: Double) -> DD {
        ddSub(lhs, DD(rhs))
    }

    /// Subtracts a `DD` from a `Double` by promoting the `Double` to `DD`.
    @inline(__always)
    static func - (lhs: Double, rhs: DD) -> DD {
        ddSub(DD(lhs), rhs)
    }

    /// Multiplies a `DD` by a `Double` by promoting the `Double` to `DD`.
    @inline(__always)
    static func * (lhs: DD, rhs: Double) -> DD {
        ddMul(lhs, DD(rhs))
    }

    /// Multiplies a `Double` by a `DD` by promoting the `Double` to `DD`.
    @inline(__always)
    static func * (lhs: Double, rhs: DD) -> DD {
        ddMul(DD(lhs), rhs)
    }

    // MARK: Compound assignment operators
    /// Compound-adds a `DD` to this value.
    @inline(__always)
    static func += (lhs: inout DD, rhs: DD) {
        lhs = ddAdd(lhs, rhs)
    }

    /// Compound-subtracts a `DD` from this value.
    @inline(__always)
    static func -= (lhs: inout DD, rhs: DD) {
        lhs = ddSub(lhs, rhs)
    }

    /// Compound-multiplies this value by a `DD`.
    @inline(__always)
    static func *= (lhs: inout DD, rhs: DD) {
        lhs = ddMul(lhs, rhs)
    }
}

/// Squares a double-double number (`x * x`) with extended precision.
///
/// This is a small convenience wrapper around `ddMul(x, x)` used
/// frequently in inner loops (e.g., Mandelbrot iteration).
@inline(__always)
internal func ddSquare(_ x: DD) -> DD {
    ddMul(x, x)
}

// MARK: - Complex Double-Double Type for Mandelbrot

/// A complex number with double-double real and imaginary parts.
///
/// This is tailored for Mandelbrot Metal's deep-zoom kernel, where
/// both the real and imaginary components need ~106 bits of precision.
internal struct DDComplex: Equatable, CustomStringConvertible {
    /// Real part of the complex value.
    var re: DD
    /// Imaginary part of the complex value.
    var im: DD

    /// Initializes a complex value from `DD` components.
    @inline(__always)
    init(re: DD = DD(), im: DD = DD()) {
        self.re = re
        self.im = im
    }

    /// Initializes a complex value from `Double` components, promoted to `DD`.
    @inline(__always)
    init(re: Double, im: Double) {
        self.re = DD(re)
        self.im = DD(im)
    }

    /// Debug-friendly representation that keeps both DD components visible.
    var description: String {
        "DDComplex(re: \(re), im: \(im))"
    }
}

// MARK: Complex Arithmetic Helpers

/// Adds two `DDComplex` numbers component-wise.
@inline(__always)
internal func ddComplexAdd(_ z1: DDComplex, _ z2: DDComplex) -> DDComplex {
    DDComplex(re: ddAdd(z1.re, z2.re),
              im: ddAdd(z1.im, z2.im))
}

/// Subtracts two `DDComplex` numbers component-wise.
@inline(__always)
internal func ddComplexSub(_ z1: DDComplex, _ z2: DDComplex) -> DDComplex {
    DDComplex(re: ddSub(z1.re, z2.re),
              im: ddSub(z1.im, z2.im))
}

/// Multiplies two `DDComplex` numbers with full double-double precision.
///
/// (a + bi) * (c + di) = (ac - bd) + (ad + bc)i
@inline(__always)
internal func ddComplexMul(_ z1: DDComplex, _ z2: DDComplex) -> DDComplex {
    let ac = ddMul(z1.re, z2.re)
    let bd = ddMul(z1.im, z2.im)
    let ad = ddMul(z1.re, z2.im)
    let bc = ddMul(z1.im, z2.re)

    let real = ddSub(ac, bd)
    let imag = ddAdd(ad, bc)

    return DDComplex(re: real, im: imag)
}

/// Returns |z|^2 for a complex double-double number.
///
/// This is often used in Mandelbrot iteration to test bailout
/// without needing a true square root.
@inline(__always)
internal func ddComplexAbsSquared(_ z: DDComplex) -> DD {
    let re2 = ddMul(z.re, z.re)
    let im2 = ddMul(z.im, z.im)
    return ddAdd(re2, im2)
}

// MARK: DDComplex Operator overloads

extension DDComplex {
    /// Adds two complex double-double values using `ddComplexAdd`.
    @inline(__always)
    static func + (lhs: DDComplex, rhs: DDComplex) -> DDComplex {
        ddComplexAdd(lhs, rhs)
    }

    /// Subtracts two complex double-double values using `ddComplexSub`.
    @inline(__always)
    static func - (lhs: DDComplex, rhs: DDComplex) -> DDComplex {
        ddComplexSub(lhs, rhs)
    }

    /// Multiplies two complex double-double values using `ddComplexMul`.
    @inline(__always)
    static func * (lhs: DDComplex, rhs: DDComplex) -> DDComplex {
        ddComplexMul(lhs, rhs)
    }

    /// Compound-adds another complex value.
    @inline(__always)
    static func += (lhs: inout DDComplex, rhs: DDComplex) {
        lhs = ddComplexAdd(lhs, rhs)
    }

    /// Compound-subtracts another complex value.
    @inline(__always)
    static func -= (lhs: inout DDComplex, rhs: DDComplex) {
        lhs = ddComplexSub(lhs, rhs)
    }

    /// Compound-multiplies by another complex value.
    @inline(__always)
    static func *= (lhs: inout DDComplex, rhs: DDComplex) {
        lhs = ddComplexMul(lhs, rhs)
    }
}
