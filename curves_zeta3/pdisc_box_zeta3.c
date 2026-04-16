/*
 * pdisc_box_zeta3.c
 *
 * Enumerate monic depressed Picard curves  y^3 = f(x)  over Q(omega),
 * where omega = e^(2*pi*i/3)  (primitive cube root of unity).
 *
 * f(x) = x^4 + c2*x^2 + c1*x + c0,   ci = ai + bi*omega  in Z[omega]
 *
 * We search  |ai|, |bi| <= c  and output every curve whose discriminant
 * norm  N(disc(f)) = disc_re^2 - disc_re*disc_im + disc_im^2  satisfies
 * 0 < N(disc(f)) <= norm_bound.
 *
 * Discriminant formulas come from norm_of_disc.ipynb:
 *   disc(x^4 + x2*x^2 + x1*x + x0) over Z[omega],  xi = ai + bi*omega
 *   gives disc_re and disc_im as explicit polynomials in (a0,b0,a1,b1,a2,b2).
 *
 * When all bi = 0 these reduce to the classical rational discriminant
 *   256*a0^3 - 128*a0^2*a2^2 + 144*a0*a1^2*a2 - 27*a1^4 + 16*a0*a2^4 - 4*a1^2*a2^3
 *
 * Build (from curves_zeta3/ directory):
 *   gcc -O2 -fopenmp pdisc_box_zeta3.c -o pdisc_box_zeta3
 *
 * Usage:
 *   ./pdisc_box_zeta3 coeff-bound norm-bound [threads]
 *
 * Output lines:  N(disc):[a0,b0,a1,b1,a2,b2]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <omp.h>
#include "../curves/cstd.h"


/* -----------------------------------------------------------------------
 * Real part of disc(x^4 + c2*x^2 + c1*x + c0)  over Z[omega],
 * where ci = ai + bi*omega.
 * All arithmetic in __int128_t to avoid overflow (coeffs fit in long,
 * but degree-5 products can reach ~c^5, and we need exact values).
 * ----------------------------------------------------------------------- */
static inline __int128_t
disc_re(__int128_t a0, __int128_t b0,
        __int128_t a1, __int128_t b1,
        __int128_t a2, __int128_t b2)
{
    return
        /* a2/b2 degree-5 terms */
        -4*a1*a1*a2*a2*a2   + 4*b1*b1*a2*a2*a2
        + 16*a0*a2*a2*a2*a2
        + 24*a1*b1*a2*a2*b2 - 12*b1*b1*a2*a2*b2
        - 64*b0*a2*a2*a2*b2
        + 12*a1*a1*a2*b2*b2 - 24*a1*b1*a2*b2*b2
        - 96*a0*a2*a2*b2*b2 + 96*b0*a2*a2*b2*b2
        - 4*a1*a1*b2*b2*b2  + 4*b1*b1*b2*b2*b2
        + 64*a0*a2*b2*b2*b2
        - 16*b0*b2*b2*b2*b2
        /* a1/b1 degree-4 terms */
        - 27*a1*a1*a1*a1 + 162*a1*a1*b1*b1 - 108*a1*b1*b1*b1
        /* mixed degree-4 terms */
        + 144*a0*a1*a1*a2   - 144*a0*b1*b1*a2
        - 128*a0*a0*a2*a2   + 128*b0*b0*a2*a2
        - 144*b0*a1*a1*b2   + 144*b0*b1*b1*b2
        + 128*a0*a0*b2*b2   - 128*b0*b0*b2*b2
        /* a0/b0 degree-3 terms (from 256*x0^3) */
        + 256*a0*a0*a0 - 768*a0*b0*b0 + 256*b0*b0*b0
        /* degree-3 cross terms */
        + 256*a0*b0*b2 - 128*b0*b0*b2
        - 288*a1*b1*b2 + 144*b1*b1*b2;
}


/* -----------------------------------------------------------------------
 * Imaginary part (omega-coefficient) of disc(x^4 + c2*x^2 + c1*x + c0).
 * ----------------------------------------------------------------------- */
static inline __int128_t
disc_im(__int128_t a0, __int128_t b0,
        __int128_t a1, __int128_t b1,
        __int128_t a2, __int128_t b2)
{
    return
        /* a2/b2 degree-5 terms */
        -8*a1*b1*a2*a2*a2   + 4*b1*b1*a2*a2*a2
        + 16*b0*a2*a2*a2*a2
        - 12*a1*a1*a2*a2*b2 + 24*a1*b1*a2*a2*b2
        + 64*a0*a2*a2*a2*b2 - 64*b0*a2*a2*a2*b2
        + 12*a1*a1*a2*b2*b2 - 12*b1*b1*a2*b2*b2
        - 96*a0*a2*a2*b2*b2
        - 8*a1*b1*b2*b2*b2  + 4*b1*b1*b2*b2*b2
        + 64*b0*a2*b2*b2*b2
        + 16*a0*b2*b2*b2*b2 - 16*b0*b2*b2*b2*b2
        /* a1/b1 degree-4 terms */
        - 108*a1*a1*a1*b1 + 162*a1*a1*b1*b1 - 27*b1*b1*b1*b1
        /* mixed degree-4 terms */
        + 288*a0*a1*b1*a2   - 144*a0*b1*b1*a2
        - 256*a0*b0*a2*a2   + 128*b0*b0*a2*a2
        - 288*b0*a1*b1*b2   + 144*b0*b1*b1*b2
        + 256*a0*b0*b2*b2   - 128*b0*b0*b2*b2
        /* a0/b0 degree-3 terms (from 256*x0^3 via im part) */
        + 768*a0*a0*b0 - 768*a0*b0*b0
        /* degree-3 cross terms */
        - 128*a0*a0*b2 + 256*a0*b0*b2
        + 144*a1*a1*b2 - 288*a1*b1*b2;
}


/* -----------------------------------------------------------------------
 * Format a curve as  x^4+(a2+b2*w)*x^2+(a1+b1*w)*x+(a0+b0*w)
 * Short form: [a0,b0,a1,b1,a2,b2]
 * ----------------------------------------------------------------------- */
static inline void
curve_string(char *buf, long a0, long b0, long a1, long b1, long a2, long b2)
{
    sprintf(buf, "[%ld,%ld,%ld,%ld,%ld,%ld]", a0, b0, a1, b1, a2, b2);
}


/* -----------------------------------------------------------------------
 * main
 * ----------------------------------------------------------------------- */
int main(int argc, char *argv[])
{
    if (argc < 3) {
        puts("Usage: pdisc_box_zeta3 coeff-bound norm-bound [threads]");
        puts("");
        puts("  Enumerates monic depressed Picard curves over Q(zeta_3):");
        puts("    y^3 = x^4 + c2*x^2 + c1*x + c0,  ci = ai+bi*omega in Z[omega]");
        puts("  with |ai|,|bi| <= coeff-bound and 0 < N(disc(f)) <= norm-bound.");
        puts("");
        puts("  Output format per line:  N(disc):[a0,b0,a1,b1,a2,b2]");
        return 0;
    }

    long c    = atol(argv[1]);  assert(c > 0);
    long Nbnd = atol(argv[2]);  assert(Nbnd >= 0);

    int threads = 0;
    if (argc > 3) threads = atoi(argv[3]);
    if (!threads) threads = omp_get_max_threads();
    fprintf(stderr, "Using %d threads\n", threads);

    /* total iterations = (2c+1)^6 */
    {
        long side = 2*c + 1;
        double tot = 1;
        for (int i = 0; i < 6; i++) tot *= side;
        fprintf(stderr,
            "Scanning %.0f curves with |ai|,|bi| <= %ld, N(disc) <= %ld...\n",
            tot, c, Nbnd);
    }

    long cnt   = 0;
    double t0  = get_time();

    /*
     * Parallelise over the two outermost loops (a2, b2).
     * Inner four loops run sequentially per thread.
     */
    #pragma omp parallel for num_threads(threads) schedule(dynamic) \
        reduction(+:cnt) collapse(2)
    for (long a2 = -c; a2 <= c; a2++) {
    for (long b2 = -c; b2 <= c; b2++) {

        for (long a1 = -c; a1 <= c; a1++) {
        for (long b1 = -c; b1 <= c; b1++) {
        for (long a0 = -c; a0 <= c; a0++) {
        for (long b0 = -c; b0 <= c; b0++) {

            __int128_t dr = disc_re(a0, b0, a1, b1, a2, b2);
            __int128_t di = disc_im(a0, b0, a1, b1, a2, b2);
            if (!dr && !di) continue;   /* degenerate (singular) curve */

            /* N(u + v*omega) = u^2 - u*v + v^2  (always >= 0) */
            __int128_t N = dr*dr - dr*di + di*di;
            /* N is always >= 0, but guard just in case */
            if (N <= 0 || N > (__int128_t)Nbnd) continue;

            char nbuf[64], cbuf[128];
            curve_string(cbuf, a0, b0, a1, b1, a2, b2);

            #pragma omp critical(output)
            {
                fprintf(stdout, "%s:%s\n", itoa128(nbuf, N), cbuf);
                fflush(stdout);
                cnt++;
            }

        }}}}
    }}

    double elapsed = get_time() - t0;
    fprintf(stderr,
        "Found %ld curves in %.3f secs using %d threads\n",
        cnt, elapsed, threads);
    return 0;
}
