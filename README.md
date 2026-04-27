# Picard Curve Enumeration

This project enumerates Picard curves with small discriminant, both over **Q** (the rationals)
and over **Q(ζ₃)** (the Eisenstein integers).

A **Picard curve** is a genus-3 curve of the form

```
y³ = f(x)
```

where `f(x)` is a degree-4 polynomial.

---

## Original Version — Enumeration over Q (`curves/`)

### What it does

Searches for all Picard curves `y³ = f(x)` with **integer coefficients**

```
f(x) = f4·x⁴ + f3·x³ + f2·x² + f1·x + f0
```

in a coefficient box `|fi| ≤ c`, whose **discriminant** is small.

The discriminant of the Picard curve `y³ = f(x)` is

```
Δ = 3⁹ · f4³ · disc(f)²
```

where `disc(f)` is the classical polynomial discriminant of `f`.
The program outputs every non-singular curve (disc ≠ 0) with `Δ ≤ disc_bound`.

### Key algorithms

| Component | File(s) | Purpose |
|---|---|---|
| Monomial tree | `mpoly128.h`, `mpoly128.c` | Stores `disc(f)` as a tree of 128-bit nodes; evaluating at a new coefficient only recomputes affected branches |
| Discriminant polynomial | `mpolydisc128.c` | Encodes the hundreds of monomials of `disc(f0,f1,f2,f3,f4)` into the tree |
| Finite differences | `polyenum128.h` | In the innermost `f0` loop, advances `disc` by **only additions** (method of finite differences for polynomials) |
| Exact arithmetic | `mpzpolyutil.c` + **GMP** | Once a candidate passes the 128-bit filter, recomputes `disc(f)` in arbitrary precision to verify and format the output |
| Output | `pdisc_box.c` | Main driver; multi-threaded via **OpenMP** |

### Libraries required

| Library | Why |
|---|---|
| **GMP** (`libgmp`) | Arbitrary-precision integers for the final discriminant computation |
| **OpenMP** (`-fopenmp`) | Multi-threaded parallelism across the coefficient box |

### Compile

```bash
cd curves
gcc -O2 -fopenmp pdisc_box.c mpolydisc128.c mpoly128.c \
    mpzpolyutil.c polyparse.c cstd.c -lgmp -lm -o pdisc_box
```

### Usage

```bash
cd curves
./pdisc_box coeff-bound disc-bound [threads [instances instance-id]]
```

Example:
```bash
./pdisc_box 3 1000000 > results.txt 2> results_log.txt
```

### Output format (`results.txt`)

```
# Picard curves y^3 = f(x),  f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0
# Output format:  Delta:[f(x)]  where  Delta = 3^9 * f4^3 * disc(f)^2
Delta:f(x)
...
```

Example line:
```
729:[x^4+x^3-x]
```

- **Left of `:`** — the curve discriminant `Δ = 3⁹ · f4³ · disc(f)²`
- **Right of `:`** — the polynomial `f(x)` defining the Picard curve `y³ = f(x)`

A smaller `Δ` means the curve has more special arithmetic structure (analogous to
small-discriminant elliptic curves, which tend to have extra symmetry or CM).

### Log format (`results_log.txt`)

Stderr is redirected here. Example:
```
Using 8 threads
Beginning scan of 4116 of 4116 (2^12.01, 10^3.61) curve equations with |f_i| <= 3...
Found 3995 of 4116 non-singular curves in 0.006 secs using 8 threads (12548.4 ns/curve)
```

---

## New Version — Enumeration over Q(ζ₃) (`curves_zeta3/`)

### What is Q(ζ₃)?

`Q(ζ₃)` is the number field obtained by adjoining a primitive cube root of unity
`ω = e^(2πi/3) = (-1 + √-3) / 2` to the rationals. Its ring of integers is the
**Eisenstein integers** `Z[ω]`, where every element is written as `a + b·ω` with
`a, b ∈ Z`.

Multiplication rule: `(a₁ + b₁ω)(a₂ + b₂ω) = (a₁a₂ − b₁b₂) + (a₁b₂ + a₂b₁ − b₁b₂)·ω`

Norm: `N(a + bω) = a² − ab + b²` (always a non-negative integer)

### What it does

Searches for monic depressed Picard curves over `Q(ζ₃)`:

```
y³ = f(x),   f(x) = x⁴ + c2·x² + c1·x + c0
```

where each coefficient is an **Eisenstein integer**: `ci = ai + bi·ω ∈ Z[ω]`.

The search box is `|ai|, |bi| ≤ c`, and the program outputs every curve with

```
0 < N(disc(f)) ≤ norm_bound
```

where `N(disc(f))` is the **norm** of the polynomial discriminant from `Q(ζ₃)` to `Q`.

### Why monic and depressed (no x³ term)?

- **Monic** (`f4 = 1`): removes one free parameter without loss of generality
  (any Picard curve can be scaled to be monic over Q(ζ₃)).
- **No x³ term** (`f3 = 0`): by a substitution `x ↦ x + α` for a suitable
  `α ∈ Z[ω]`, any quartic can be put in this "depressed" form. This reduces the
  search from 5 Eisenstein-integer coefficients to 3, making enumeration tractable.

### Discriminant formula

The discriminant of `x⁴ + c2·x² + c1·x + c0` over `Z[ω]`, with `ci = ai + bi·ω`, is

```
disc(f) = disc_re + disc_im · ω
```

where `disc_re` and `disc_im` are explicit degree-5 polynomials in
`(a0, b0, a1, b1, a2, b2)` with integer coefficients, derived in
`norm_of_disc.ipynb` (Sage computation).

When all `bi = 0` (rational coefficients), `disc_re` reduces to the classical formula

```
256·a0³ − 128·a0²·a2² + 144·a0·a1²·a2 − 27·a1⁴ + 16·a0·a2⁴ − 4·a1²·a2³
```

and `disc_im = 0`, as expected.

The norm is then:
```
N(disc(f)) = disc_re² − disc_re·disc_im + disc_im²
```

This is always a non-negative integer, and equals zero only when `f` is singular.

### Compile

```bash
cd curves_zeta3
gcc -O2 -fopenmp pdisc_box_zeta3.c -o pdisc_box_zeta3 -lm
```

### Usage

```bash
cd curves_zeta3
./pdisc_box_zeta3 coeff-bound norm-bound [threads]
```

Example:
```bash
./pdisc_box_zeta3 2 1000000 > results_zeta3.txt 2> log_zeta3.txt
```

### Output format (`results_zeta3.txt`)

```
# Picard curves over Q(zeta_3):  y^3 = f(x),  f(x) = x^4 + (a2+b2*w)*x^2 + (a1+b1*w)*x + (a0+b0*w),  w = e^(2*pi*i/3)
# Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2]  where  N(disc(f)) = disc_re^2 - disc_re*disc_im + disc_im^2
N(disc(f)):[a0,b0,a1,b1,a2,b2]
...
```

Example line:
```
729:[0,0,0,-1,0,0]
```

- **Left of `:`** — the norm `N(disc(f)) ∈ Z`, a positive integer
- **Right of `:`** — the six Eisenstein-integer parameters; the curve is
  ```
  y³ = x⁴ + (a2+b2·ω)x² + (a1+b1·ω)x + (a0+b0·ω)
  ```

### Log format (`log_zeta3.txt`)

Stderr is redirected here. Example:
```
Using 8 threads
Scanning 15625 curves with |ai|,|bi| <= 2, N(disc) <= 1000000...
Found 2987 curves in 0.016 secs using 8 threads
```

---

## General Version — Enumeration over Q(ζ₃) with arbitrary coefficients (`curves_zeta3/`)

### What it does

Searches for **all** Picard curves over `Q(ζ₃)` with general Eisenstein integer coefficients:

```
y³ = f(x),   f(x) = f4·x⁴ + f3·x³ + f2·x² + f1·x + f0,   fi = ai + bi·ω ∈ Z[ω]
```

This is the direct Q(ζ₃) analogue of the original `pdisc_box` — no restrictions on f4 or f3.

### Symmetry reductions

To avoid enumerating isomorphic curves twice:

| Symmetry | Isomorphism | Normalization applied |
|---|---|---|
| `y → −y` | `f(x) ~ −f(x)`, i.e. `f4 ~ −f4` | Require f4 in positive half: `a4 > 0`, or `a4 = 0` and `b4 > 0` |
| `x → −x` | `(f4,f3,f2,f1,f0) ~ (f4,−f3,f2,−f1,f0)` | Require f3 in positive half: `a3 > 0`, or `a3 = 0` and `b3 ≥ 0` |

### Discriminant formula

The classical 16-term discriminant formula for `f4·x⁴ + f3·x³ + f2·x²+ f1·x + f0` is evaluated entirely in `Z[ω]` using `__int128_t` arithmetic:

```
disc(f) = 256·f4³·f0³ − 192·f4²·f3·f1·f0² − 128·f4²·f2²·f0²
        + 144·f4²·f2·f1²·f0 − 27·f4²·f1⁴ + 144·f4·f3²·f2·f0²
        − 6·f4·f3²·f1²·f0 − 80·f4·f3·f2²·f1·f0 + 18·f4·f3·f2·f1³
        + 16·f4·f2⁴·f0 − 4·f4·f2³·f1² − 27·f3⁴·f0²
        + 18·f3³·f2·f1·f0 − 4·f3³·f1³ − 4·f3²·f2³·f0 + f3²·f2²·f1²
```

The output quantity is `N(disc(f)) = disc_re² − disc_re·disc_im + disc_im²`.

### Compile

```bash
cd curves_zeta3
gcc -O2 -fopenmp pdisc_box_zeta3_general.c -o pdisc_box_zeta3_general -lm
```

### Usage

```bash
cd curves_zeta3
./pdisc_box_zeta3_general coeff-bound norm-bound [threads]
```

Example:
```bash
./pdisc_box_zeta3_general 2 1000000
```

Results are written **automatically** to files in the same directory as the executable — no shell redirection needed:

| File | Contents |
|---|---|
| `results_zeta3_general.txt` | One curve per line |
| `log_zeta3_general.txt` | Thread count, scan size, timing |

### Output format (`results_zeta3_general.txt`)

```
# Picard curves over Q(zeta_3):  y^3 = f(x)
# f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,  fi = ai+bi*w,  w = e^(2*pi*i/3)
# Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
...
```

Example line:
```
729:[0,0,0,-1,0,0,0,0,1,0]
```

- **Left of `:`** — `N(disc(f))`, a positive integer
- **Right of `:`** — the ten Eisenstein-integer parameters defining `f(x)`

---

## Fast General Version — Monomial Tree + Finite Differences (`curves_zeta3/`)

A faster version of `pdisc_box_zeta3_general` for large coefficient searches (`c ≥ 5`).
Uses the same two techniques as the rational `pdisc_box`. Design notes: `FAST_GENERAL_DESIGN.md`.

### Key idea: disc is cubic in f0

Regroup the 16-term quartic discriminant by power of f0:

```
disc(f) = P3·f0³ + P2·f0² + P1·f0 + P0
```

where P0…P3 ∈ Z[ω] depend only on (f1, f2, f3, f4):

```
P3 = 256·f4³

P2 = -27·f3⁴ + (-128·f4²·f2² + 144·f4·f3²·f2) + (-192·f4²·f3)·f1

P1 = (16·f4·f2⁴ - 4·f3²·f2³) + (144·f4²·f2 - 6·f4·f3²)·f1²
   + (-80·f4·f3·f2² + 18·f3³·f2)·f1

P0 = -27·f4²·f1⁴ + (18·f4·f3·f2 - 4·f3³)·f1³ + (-4·f4·f2³ + f3²·f2²)·f1²
```

For fixed b0, disc is cubic in the integer a0, enabling finite differences.

### Speedups

| Technique | Where applied | Saving |
|---|---|---|
| Monomial tree | Outer loops (f4→f3→f2→f1) | Cache P3, partial P2/P1/P0 at each level |
| Finite differences | Innermost a0 loop | 3 `zw_t` additions per step vs ~16 multiplications |

Measured speedup: **~16× faster** than `pdisc_box_zeta3_general` at c=4 (8 threads).

### Compile

```bash
cd curves_zeta3
gcc -O2 -fopenmp pdisc_box_zeta3_general_fast.c -o pdisc_box_zeta3_general_fast -lm
```

### Usage

```bash
cd curves_zeta3
./pdisc_box_zeta3_general_fast coeff-bound norm-bound [threads]
```

Example:
```bash
./pdisc_box_zeta3_general_fast 4 1000000
```

Results written automatically to the same directory:

| File | Contents |
|---|---|
| `results_zeta3_general_fast.txt` | One curve per line |
| `log_zeta3_general_fast.txt` | Thread count, scan size, timing |

### Output format

Same as `pdisc_box_zeta3_general`:

```
# Picard curves over Q(zeta_3):  y^3 = f(x)
# f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,  fi = ai+bi*w,  w = e^(2*pi*i/3)
# Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
...
```

Same symmetry reductions as `pdisc_box_zeta3_general` (f4 and f3 positive-half normalisation).

---

## Prerequisites (Ubuntu / WSL)

```bash
sudo apt-get install libgmp-dev   # for the original version only
# GCC with OpenMP is included by default
```

---

## References

- Paper: *Picard curves with small conductor* — see `1806.06289v5.pdf`
- Notebook: `curves_zeta3/norm_of_disc.ipynb` — Sage derivation of the discriminant formulas over Q(ζ₃)
