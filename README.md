# Double-Double Precision Arithmetic

## 1. What Is Double-Double Precision?

A **double-double number** is a software-emulated floating-point type that extends the precision of a standard IEEE-754 double (64-bit, ~53 bits of mantissa ≈ 16 decimal digits) by representing one logical number as the unevaluated sum of two doubles:

$$
x = x_{hi} + x_{lo}, \quad |x_{lo}| \le \tfrac{1}{2} \, \mathrm{ulp}(x_{hi})
$$

- **`hi`** holds the main value (nearest double).  
- **`lo`** carries the round-off error or residual.  
- Together, they give roughly **106 bits of mantissa (~31–32 decimal digits)**.

This is *not* a new IEEE format — it’s a software trick relying on carefully ordered standard double arithmetic.

---

## 2. Core Tools: Error-Free Transforms

The entire scheme depends on algorithms that decompose floating-point operations into exact parts.

### TwoSum
Computes the exact sum of two doubles:

$$
(a + b) = s + e \quad \text{exactly.}
$$

### TwoProd
Computes the exact product of two doubles (using FMA if available):

$$
(a \cdot b) = p + e \quad \text{exactly.}
$$

These **error-free transforms** guarantee correctness within IEEE-754 double precision — essential for stable double-double arithmetic.

---

## 3. Double-Double Addition

Given two double-double numbers  
\( x = (x_h, x_l) \), \( y = (y_h, y_l) \):

1. \( (s, e) = \text{TwoSum}(x_h, y_h) \)  
   Add high parts, capture residual.  
2. \( t = x_l + y_l \)  
   Add low parts.  
3. \( (s_2, e_2) = \text{TwoSum}(s, t) \)  
   Fold in the low parts.  
4. Result:

$$
z_h = s_2, \quad z_l = e + e_2
$$

After a brief **renormalization** step (ensuring \( |z_l| < \mathrm{ulp}(z_h) \)), this yields ~106-bit precision.

---

## 4. Double-Double Multiplication

Given \( x = (x_h, x_l) \), \( y = (y_h, y_l) \):

1. \( (p, e) = \text{TwoProd}(x_h, y_h) \)  
   Main product and rounding error.  
2. \( p_1 = p \)  
   Start with the main product.  
3. \( p_2 = (x_h \cdot y_l) + (x_l \cdot y_h) \)  
   Cross terms.  
4. \( (s, e_2) = \text{TwoSum}(p_1, p_2) \)  
   Combine partial results.  
5. Accumulate small terms:

$$
\text{lo} = e + e_2 + (x_l \cdot y_l)
$$

6. Renormalize:

$$
(\text{hi}, \text{lo}_2) = \text{TwoSum}(s, \text{lo})
$$

Final result:

$$
z_h = \text{hi}, \quad z_l = \text{lo}_2
$$

This achieves full **double-double precision multiplication**, with error bounded to about one ulp of the 106-bit result.

---

## 5. Why It Works

- IEEE-754 doubles include guard digits and FMA, enabling recovery of round-off errors exactly.  
- Keeping the small remainder separately yields roughly \( 53 + 53 = 106 \) bits of precision.  
- Performance: about **10–20× slower** than native double operations, but still **much faster** than arbitrary-precision libraries.

---

## 6. Other Operations

- **Division:** Compute \( \text{hi} = x_h / y_h \), then refine via Newton iteration using DD multiply.  
- **Square root:** Use a double-precision initial guess, refine with Newton’s method.  
- **Normalization:** After each operation, ensure `hi` is the rounded double and `lo` is the small correction.

---

## 7. Applications

- Deep-zoom fractal rendering (Mandelbrot / Julia).  
- High-precision geometry kernels.  
- Numerical analysis and compensated algorithms.  
- Validation of IEEE quad-precision implementations.

---

### References

- Dekker, T. J. (1971). *A Floating-Point Technique for Extending the Available Precision.*  
- Knuth, D. E. (1997). *The Art of Computer Programming, Vol. 2: Seminumerical Algorithms.*  
- Bailey, D. H. (2001). *High-Precision Floating-Point Arithmetic in Scientific Computation.*
