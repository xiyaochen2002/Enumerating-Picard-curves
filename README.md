# Picard Curve Enumeration

This project enumerates Picard curves with small discriminant, both over **Q** (the rationals)
and over **Q(ζ₃)** (the Eisenstein integers).

A **Picard curve** is a genus-2 curve of the form

```
y³ = f(x)
```

where `f(x)` is a degree-4 polynomial.

---

## Original Version — Enumeration over Q (`curves/`)

### What it does

Searches for all Picard curves `y³ = f(x)` with **integer coefficients**

```
f(x) = f₄·x⁴ + f₃·x³ + f₂·x² + f₁·x + f₀
```

in a coefficient box `|fᵢ| ≤ c`, whose **discriminant** is small.

The discriminant of the Picard curve `y³ = f(x)` is

```
Δ = 3⁹ · f₄³ · disc(f)²
```

where `disc(f)` is the classical polynomial discriminant of `f`.  
The program outputs every non-singular curve (disc ≠ 0) with `Δ ≤ norm_bound`.

### Key algorithms

| Component | File(s) | Purpose |
|---|---|---|
| Monomial tree | `mpoly128.h`, `mpoly128.c` | Stores `disc(f)` as a tree of 128-bit nodes; evaluating at a new coefficient only recomputes affected branches |
| Discriminant polynomial | `mpolydisc128.c` | Encodes the hundreds of monomials of `disc(f₀,f₁,f₂,f₃,f₄)` into the tree |
| Finite differences | `polyenum128.h` | In the innermost `f₀` loop, advances `disc` by **only additions** (method of finite differences for polynomials) |
| Exact arithmetic | `mpzpolyutil.c` + **GMP** | Once a candidate passes the 128-bit filter, recomputes `disc(f)` in arbitrary precision to verify and format the output |
| Output | `pdisc_box.c` | Main driver; multi-threaded via **OpenMP** |

### Libraries required

| Library | Why |
|---|---|
| **GMP** (`libgmp`) | Arbitrary-precision integers for the final discriminant computation |
| **OpenMP** (`-fopenmp`) | Multi-threaded parallelism across the coefficient box |

### Executables

| Program | Description |
|---|---|
| `pdisc_box` | Enumerates all non-singular curves (no smoothness filter) |
| `psmoothdisc_box` | Same, but additionally filters for 7-smooth discriminants (requires `smooth128.c`) |

### Usage

```bash
./pdisc_box coeff-bound norm-bound [threads [instances instance-id]]
```

Example:
```bash
./pdisc_box 3 1000000 > results.txt
```

### Output format (`results.txt`)

Each line:
```
Δ:f(x)
```

Example:
```
18915363:[x⁴ - 3x² - 3x - 1]
```

- **Left of `:`** — the curve discriminant `Δ = 3⁹ · f₄³ · disc(f)²`
- **Right of `:`** — the polynomial `f(x)` defining the Picard curve `y³ = f(x)`

A smaller `Δ` means the curve has more special arithmetic structure (analogous to
small-discriminant elliptic curves, which tend to have extra symmetry or CM).

---

## New Version — Enumeration over Q(ζ₃) (`curves_zeta3/`)

### What is Q(ζ₃)?

`Q(ζ₃)` is the number field obtained by adjoining a primitive cube root of unity
`ω = e^(2πi/3) = (-1 + √-3) / 2` to the rationals.  Its ring of integers is the
**Eisenstein integers** `Z[ω]`, where every element is written as `a + b·ω` with
`a, b ∈ Z`.

Multiplication rule: `(a₁ + b₁ω)(a₂ + b₂ω) = (a₁a₂ − b₁b₂) + (a₁b₂ + a₂b₁ − b₁b₂)·ω`

Norm: `N(a + bω) = a² − ab + b²` (always a non-negative integer)

### What it does

Searches for monic depressed Picard curves over `Q(ζ₃)`:

```
y³ = f(x),   f(x) = x⁴ + c₂x² + c₁x + c₀
```

where each coefficient is an **Eisenstein integer**: `cᵢ = aᵢ + bᵢ·ω ∈ Z[ω]`.

The search box is `|aᵢ|, |bᵢ| ≤ c`, and the program outputs every curve with

```
0 < N(disc(f)) ≤ norm_bound
```

where `N(disc(f))` is the **norm** of the discriminant from `Q(ζ₃)` to `Q`.

### Why monic and depressed (no x³ term)?

- **Monic** (`f₄ = 1`): removes one free parameter without loss of generality
  (any Picard curve can be scaled to be monic over Q(ζ₃)).
- **No x³ term** (`f₃ = 0`): by a substitution `x ↦ x + α` for a suitable
  `α ∈ Z[ω]`, any quartic can be put in this "depressed" form. This reduces the
  search from 5 Eisenstein-integer coefficients to 3, making enumeration tractable.

### Discriminant formula

The discriminant of `x⁴ + c₂x² + c₁x + c₀` over `Z[ω]`, with `cᵢ = aᵢ + bᵢω`, is

```
disc(f) = disc_re + disc_im · ω
```

where `disc_re` and `disc_im` are explicit degree-5 polynomials in
`(a₀, b₀, a₁, b₁, a₂, b₂)` with integer coefficients, derived in
`norm_of_disc.ipynb` (Sage computation). 

When all `bᵢ = 0` (rational coefficients), `disc_re` reduces to the classical formula

```
256·a₀³ − 128·a₀²·a₂² + 144·a₀·a₁²·a₂ − 27·a₁⁴ + 16·a₀·a₂⁴ − 4·a₁²·a₂³
```

and `disc_im = 0`, as expected.

The norm is then:
```
N(disc) = disc_re² − disc_re·disc_im + disc_im²
```

This is always a non-negative integer, and equals zero only when `f` is singular.

### Implementation

Everything is computed in a single file using **`__int128_t`** (128-bit integers
built into GCC), which is sufficient for coefficient bounds up to roughly `c ≤ 10⁶`
before overflow. No GMP is needed.

| Aspect | Original (Q) | New (Q(ζ₃)) |
|---|---|---|
| Coefficient space | 5 integers `(f₀…f₄)` | 6 integers `(a₀,b₀,a₁,b₁,a₂,b₂)` |
| Discriminant evaluation | mpoly128 tree + finite differences | Direct formula (35 terms each) |
| Arithmetic | 128-bit integers + GMP | 128-bit integers only |
| Libraries | GMP, OpenMP | OpenMP only |
| Parallelism | OpenMP over `(f₂, f₃)` slice | OpenMP over `(a₂, b₂)` slice |

### Usage

```bash
cd curves_zeta3
./pdisc_box_zeta3 coeff-bound norm-bound [threads]
```

Example:
```bash
./pdisc_box_zeta3 2 1000000 > results_zeta3.txt
```

### Output format (`results_zeta3.txt`)

Each line:
```
N(disc):[a0,b0,a1,b1,a2,b2]
```

Example:
```
18915363:[0,0,-1,0,0,0]
```

- **Left of `:`** — the norm `N(disc(f)) ∈ Z`, a positive integer
- **Right of `:`** — the six Eisenstein-integer parameters; the curve is
  ```
  y³ = x⁴ + (a₂+b₂ω)x² + (a₁+b₁ω)x + (a₀+b₀ω)
  ```

Curves with small `N(disc)` are the Eisenstein-integer analogues of
small-discriminant curves over Q, and are expected to have exceptional
arithmetic properties (e.g. complex multiplication by `Z[ω]`).

---

## Building

### Prerequisites (Ubuntu / WSL)

```bash
sudo apt-get install libgmp-dev   # for the original version only
# GCC with OpenMP is included by default
```

### Compile

```bash
# Original version (Q)
cd curves
gcc -O2 -fopenmp pdisc_box.c mpolydisc128.c mpoly128.c \
    mpzpolyutil.c polyparse.c cstd.c -lgmp -lm -o pdisc_box

# New version (Q(zeta_3))
cd curves_zeta3
gcc -O2 -fopenmp pdisc_box_zeta3.c -o pdisc_box_zeta3 -lm
```

---

## References

- Paper: *Picard curves with small conductor* — see `1806.06289v5.pdf`
- Notebook: `curves_zeta3/norm_of_disc.ipynb` — Sage derivation of the discriminant formulas over Q(ζ₃)
