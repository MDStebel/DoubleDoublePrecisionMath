# Double-Double Precision Math

High-precision floating-point and complex arithmetic for **Mandelbrot Metal**, implemented using **double-double (DD)** arithmetic in pure Swift.

This module provides:

- `DD`: a high-precision real type (~106 bits, ~31–32 decimal digits)
- `DDComplex`: a complex number type using `DD` for real + imaginary parts
- Exact error-free transforms (`twoSum`, `quickTwoSum`, `twoProd`)
- DD renormalization (`ddRenormalize`)
- High-precision operations (`ddAdd`, `ddSub`, `ddMul`, `ddSquare`)
- Complex arithmetic (`ddComplexAdd`, `ddComplexMul`, etc.)
- Operator overloads for natural mathematical syntax

These operations power Mandelbrot Metal's **CPU Deep Mode** when zooming past the precision limits of standard 64-bit floating-point. CPU Deep Mode is currently focused on Mandelbrot Set exploration; Julia Set rendering uses the GPU path while sharing the app's palette, SSAA, lighting, capture, and bookmark controls.

---

## 1. Double-Double Real Representation: `DD`

A double-double value stores a real number as:

**x ≈ hi + lo**

Where:

- `hi` — the leading significand bits (normal Double)
- `lo` — the trailing error term from previous operations

Together they yield ~106 bits of precision — roughly **twice** the precision of native FP64.

### Type Definition

```swift
internal struct DD: Equatable, CustomStringConvertible {
    var hi: Double
    var lo: Double
}
```

### Initializers

```swift
DD(1.0)
DD(hi: x, lo: y)
```

---

## 2. Error-Free Transformations

### 2.1 `twoSum(a, b)`

Computes an **exact sum** for finite, non-overflowing inputs using Knuth’s
algorithm:

```swift
let (s, e) = twoSum(a, b)
```

---

### 2.2 `quickTwoSum(a, b)`

A faster version of `twoSum` when `|a| >= |b|`.

---

### 2.3 `twoProd(a, b)`

Computes an **exact product** using fused-multiply-add when the rounded
product is finite. Overflow and other non-finite products return a zero tail
instead of producing a NaN residual.

---

### 2.4 `ddRenormalize(hi, lo)`

Normalizes a high/low pair using `twoSum`, without requiring the high
component to dominate the low component. This is used after cancellation-heavy
operations where `quickTwoSum`'s precondition may not hold.

---

## 3. Double-Double Arithmetic

### 3.1 `ddAdd(x, y)`

High-precision addition using separate `twoSum` passes for the high and low
components, followed by conservative renormalization. This is a little more
defensive than the common fast "sloppy add" pattern and behaves better when
the high components nearly cancel.

### 3.2 `ddMul(x, y)`

High-precision multiplication using:

- `twoProd`
- Cross terms
- Final `ddRenormalize`

### 3.3 `ddSub(x, y)`

Subtracts by negating `y` and using `ddAdd`.

### 3.4 `ddSquare(x)`

Convenience wrapper around `ddMul(x, x)`.

---

## 4. Operator Overloads for `DD`

Supports:

```swift
+  -  *
+= -= *=
Double mixed arithmetic
```

---

## 5. Complex Double-Double: `DDComplex`

```swift
internal struct DDComplex: Equatable, CustomStringConvertible {
    var re: DD
    var im: DD
}
```

---

## 6. Complex Arithmetic

Supports:

- Addition
- Subtraction
- Multiplication
- Squared magnitude

---

## 7. Operator Overloads for `DDComplex`

Allows clean expressions like:

```swift
z = z * z + c
```

---

## 8. Mandelbrot Example

```swift
func mandelbrotIterationsDD(...) -> Int {
    var z = DDComplex(re: 0.0, im: 0.0)
    let c = DDComplex(re: cr, im: ci)

    for i in 0..<maxIter {
        z = z * z + c
        if ddComplexAbsSquared(z).hi > 4.0 { return i }
    }
    return maxIter
}
```

### Julia Mode Note

The same `z = z * z + c` recurrence is used for Julia sets, but the roles differ:

- Mandelbrot Set: `z` starts at `0`, and each pixel supplies `c`.
- Julia Set: each pixel supplies the starting `z`, and the renderer uses a fixed Julia constant.

This README documents the double-double primitives used by the Mandelbrot deep path. If CPU deep rendering is later extended to Julia mode, the same `DD` and `DDComplex` types can represent the Julia starting point and fixed constant.

---

## 9. Performance Notes

- All math primitives use `@inline(__always)`
- No heap allocation
- `twoProd` uses hardware FMA
- Optimized for deep fractal zooming
- Current CPU deep-render integration is Mandelbrot-focused; Julia mode remains GPU-rendered.

---

## 10. Future Extensions

- Division
- Reciprocal
- Trigonometric functions
- Additional utilities
- Julia-specific CPU deep-render integration

---

## 11. License

Part of **Mandelbrot Metal**  
© 2025-2026 Michael Stebel. All rights reserved.
