# Double-Double Precision Math

High-precision floating-point and complex arithmetic for **Mandelbrot Metal**, implemented using **double-double (DD)** arithmetic in pure Swift.

This module provides:

- `DD`: a high-precision real type (~106 bits, ~31–32 decimal digits)
- `DDComplex`: a complex number type using `DD` for real + imaginary parts
- Exact error-free transforms (`twoSum`, `quickTwoSum`, `twoProd`)
- High-precision operations (`ddAdd`, `ddSub`, `ddMul`, `ddSquare`)
- Complex arithmetic (`ddComplexAdd`, `ddComplexMul`, etc.)
- Operator overloads for natural mathematical syntax

These operations power Mandelbrot Metal’s **CPU Deep Mode** when zooming past the precision limits of standard 64-bit floating-point.

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
internal struct DD {
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

Computes an **exact sum** using Knuth’s algorithm:

```swift
let (s, e) = twoSum(a, b)
```

---

### 2.2 `quickTwoSum(a, b)`

A faster version of `twoSum` when `|a| >= |b|`.

---

### 2.3 `twoProd(a, b)`

Computes an **exact product** using fused-multiply-add.

---

## 3. Double-Double Arithmetic

### 3.1 `ddAdd(x, y)`

High-precision addition using `twoSum` + `quickTwoSum`.

### 3.2 `ddMul(x, y)`

High-precision multiplication using:

- `twoProd`
- Cross terms
- Final normalization

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
internal struct DDComplex {
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

---

## 9. Performance Notes

- All math primitives use `@inline(__always)`
- No heap allocation
- `twoProd` uses hardware FMA
- Optimized for deep fractal zooming

---

## 10. Future Extensions

- Division
- Reciprocal
- Trigonometric functions
- Additional utilities

---

## 11. License

Part of **Mandelbrot Metal**  
© 2025 Michael Stebel. All rights reserved.

