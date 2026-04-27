# Fast General Version Design Notes
## `pdisc_box_zeta3_general_fast.c`

Goal: monomial tree + finite differences for the general Z[ω] Picard curve enumerator,
analogous to what `pdisc_box.c` does for the rational case.

---

## Key Observation: disc is cubic in f0

Regroup the 16-term quartic discriminant by power of f0:

```
disc(f) = P3·f0³ + P2·f0² + P1·f0 + P0
```

where each Pi ∈ Z[ω] depends only on (f1, f2, f3, f4):

```
P3 = 256·f4³

P2 = -192·f4²·f3·f1
   - 128·f4²·f2²
   + 144·f4·f3²·f2
   -  27·f3⁴

P1 = 144·f4²·f2·f1²
   -   6·f4·f3²·f1²
   -  80·f4·f3·f2²·f1
   +  16·f4·f2⁴
   +  18·f3³·f2·f1
   -   4·f3²·f2³

P0 = -27·f4²·f1⁴
   +  18·f4·f3·f2·f1³
   -   4·f4·f2³·f1²
   -   4·f3³·f1³
   +      f3²·f2²·f1²
```

Since f0 = a0 + b0·ω and each fi = ai + bi·ω, disc is a cubic in a0 for fixed b0 and
the outer variables. This enables finite differences in the a0 inner loop.

---

## Finite Differences in a0 (Inner Loop)

For fixed (b0, f1, f2, f3, f4), disc(a0) is a cubic in the integer a0. Set up 4
state variables d[0..3] ∈ Z[ω]:

```
Evaluate disc at a0 = -c, -c+1, -c+2, -c+3:
  v[k] = disc({-c+k, b0})   for k = 0,1,2,3

d[0] = v[0]
d[1] = v[1] - v[0]
d[2] = v[2] - 2*v[1] + v[0]
d[3] = v[3] - 3*v[2] + 3*v[1] - v[0]   (constant; = 6*P3 as a zw_t)
```

Then for each a0 = -c+4, -c+5, ..., c (advance one step):
```
d[0] += d[1]
d[1] += d[2]
d[2] += d[3]
```
Each step: 3 zw_t additions = 6 i128 additions. No multiplications.

After each step: N = zw_norm(d[0]), check bound and output.

Cost of setup: 4 cubic evaluations (each ~12 zw_t muls) + 3 differences = ~48 zw_t muls.
Cost per a0 step: 6 i128 additions + 1 norm (3 i128 muls).

Compare to current brute force: ~16 zw_t muls = 32 i128 muls per step.
Speedup on a0 loop: ~5x for large c (saves multiplications).

---

## Monomial Tree for Outer Variables

The monomial tree amortizes recomputation of P0, P1, P2, P3 as we iterate over
the outer variables (a4,b4,a3,b3,a2,b2,a1,b1).

### Decomposition by outer variable dependency

```
P3        depends on: f4 only
P2_f4f3   depends on: f4, f3        (terms: -27*f3^4)
P2_f4f3f2 depends on: f4, f3, f2   (terms: -128*f4^2*f2^2 + 144*f4*f3^2*f2)
P2_f4f3f1 depends on: f4, f3, f1   (terms: -192*f4^2*f3*f1)
P2 = P2_f4f3 + P2_f4f3f2 + P2_f4f3f1    (combine after f2 and f1 loops)

P1 similarly splits into:
  P1_f4f3f2   = -4*f3^2*f2^3 + 16*f4*f2^4   (f2-only terms, given f4,f3)
  P1_f4f3f2f1 = 144*f4^2*f2*f1^2 - 6*f4*f3^2*f1^2 - 80*f4*f3*f2^2*f1 + 18*f3^3*f2*f1
                                                                  (needs f1)

P0 similarly.
```

### Loop structure with tree

```
for (a4,b4):  f4 = a4+b4*w                    [symmetry: f4 in positive half]
  P3 = 256*zw_cube(f4)                         [once per f4]

  for (a3,b3):  f3 = a3+b3*w                  [symmetry: f3 in positive half]
    f3sq = zw_sq(f3); f3cu = ...; f3qu = ...
    P2_f4f3 = -27*f3qu                         [once per (f4,f3)]
    P1_f4f3 = 18*f3cu*f1... wait, depends on f1 too

    for (a2,b2):  f2 = a2+b2*w
      f2sq = ...; f2cu = ...; f2qu = ...
      P2_f4f2 = -128*f4sq*f2sq + 144*f4*f3sq*f2   [update per f2 change]
      P1_f4f2 = 16*f4*f2qu - 4*f3sq*f2cu          [update per f2 change]

      for (a1,b1):  f1 = a1+b1*w
        f1sq = ...; f1cu = ...; f1qu = ...
        P2 = P2_f4f3 + P2_f4f2 - 192*f4sq*f3*f1
        P1 = P1_f4f2 + (144*f4sq*f2 - 6*f4*f3sq)*f1sq
                      + (-80*f4*f3*f2sq + 18*f3cu*f2)*f1
        P0 = -27*f4sq*f1qu + (18*f4*f3*f2 - 4*f3cu)*f1cu
             + (-4*f4*f2cu + f3sq*f2sq)*f1sq

        for b0 = -c to c:
          // evaluate disc at a0 = -c, -c+1, -c+2, -c+3 using P0..P3
          // and f0 = {a0, b0}
          setup_finite_diff(d, P0,P1,P2,P3, b0, -c)

          for a0 = -c to c:
            check N(d[0]) against bound, output if within
            // advance
            d[0] += d[1]; d[1] += d[2]; d[2] += d[3]
```

### Key savings from tree

- P3: computed once per (a4,b4), shared across all (a3,b3,a2,b2,a1,b1,b0,a0)
- P2_f4f3: computed once per (a4,b4,a3,b3)
- P2_f4f2, P1_f4f2: computed once per (a4,b4,a3,b3,a2,b2)
- Powers of f2 (f2sq, f2cu, f2qu): computed once per (a2,b2)
- Powers of f3 (f3sq, f3cu, f3qu): computed once per (a3,b3)
- Powers of f4 (f4sq, f4cu): computed once per (a4,b4)

This is the monomial tree: each node cached at the loop level of its deepest variable.

---

## setup_finite_diff

```c
static void setup_finite_diff(zw_t d[4],
                               zw_t P0, zw_t P1, zw_t P2, zw_t P3,
                               long b0, long a0_start)
{
    // evaluate cubic at a0 = a0_start, ..., a0_start+3
    zw_t v[4];
    for (int k = 0; k < 4; k++) {
        zw_t f0 = {a0_start + k, b0};
        zw_t f0sq = zw_mul(f0,f0), f0cu = zw_mul(f0sq,f0);
        v[k] = zw_add(zw_add(zw_add(zw_smul_zw(P3,f0cu),
                                     zw_smul_zw(P2,f0sq)),
                              zw_smul_zw(P1,f0)),
                       P0);
    }
    d[0] = v[0];
    d[1] = zw_sub(v[1], v[0]);
    d[2] = zw_sub(zw_sub(v[2], zw_add(v[1],v[1])), zw_neg(v[0]));
    // i.e. v[2] - 2*v[1] + v[0]
    d[3] = zw_sub(zw_add(zw_sub(v[3], zw_smul(3,v[2])),
                          zw_smul(3,v[1])), v[0]);
    // Note: d[3] = 6*P3 (constant for cubic) — can also set directly
}
```

Actually d[3] = 6*P3 exactly (3rd finite difference of a monic cubic is 6 times leading coeff).
So: d[3] = zw_smul(6, P3).  This avoids 3 of the 4 evaluations.

---

## File name and output

New file: `curves_zeta3/pdisc_box_zeta3_general_fast.c`
Executable: `pdisc_box_zeta3_general_fast`
Output files (auto in same directory):
- `results_zeta3_general_fast.txt`
- `log_zeta3_general_fast.txt`

Same output format as the slow version:
```
N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
```

---

## Helper needed: zw_smul_zw (multiply two zw_t)

In current code, `zw_smul(i128 s, zw_t x)` multiplies a scalar by a zw_t.
We also need `zw_mul(zw_t, zw_t)` for P*f0 type products — already exists.

For the cubic evaluation:
```c
P3*f0^3 = zw_mul(P3, f0cu)    // zw_t * zw_t
P2*f0^2 = zw_mul(P2, f0sq)
P1*f0   = zw_mul(P1, f0)
```

---

## Status

- [ ] Write `pdisc_box_zeta3_general_fast.c` following the loop structure above
- [ ] Verify output matches slow version for c=2, norm_bound=1000000
- [ ] Compile and benchmark vs slow version

---

## Compile command (planned)

```bash
cd ~/mxm/curves_zeta3
gcc -O2 -fopenmp pdisc_box_zeta3_general_fast.c -o pdisc_box_zeta3_general_fast -lm
```
