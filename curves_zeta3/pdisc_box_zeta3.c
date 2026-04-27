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
 * disc(x^4 + c2*x^2 + c1*x + c0)  is computed directly in Z[omega] using
 * the 6-term monic depressed quartic formula:
 *   256*c0^3 - 128*c2^2*c0^2 + 144*c2*c1^2*c0 - 27*c1^4 + 16*c2^4*c0 - 4*c2^3*c1^2
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

typedef __int128_t i128;

/* Z[omega] element: re + im*omega,  omega^2 = -omega - 1 */
typedef struct { i128 re, im; } zw_t;

static inline zw_t zw_mul(zw_t x, zw_t y) {
    return (zw_t){ x.re*y.re - x.im*y.im,
                   x.re*y.im + x.im*y.re - x.im*y.im };
}
static inline zw_t zw_add(zw_t x, zw_t y) { return (zw_t){ x.re+y.re, x.im+y.im }; }
static inline zw_t zw_sub(zw_t x, zw_t y) { return (zw_t){ x.re-y.re, x.im-y.im }; }
static inline zw_t zw_smul(i128 s, zw_t x) { return (zw_t){ s*x.re, s*x.im }; }

/* N(a + b*omega) = a^2 - a*b + b^2 */
static inline i128 zw_norm(zw_t x) { return x.re*x.re - x.re*x.im + x.im*x.im; }


/* -----------------------------------------------------------------------
 * disc(x^4 + c2*x^2 + c1*x + c0) over Z[omega]
 * = 256*c0^3 - 128*c2^2*c0^2 + 144*c2*c1^2*c0 - 27*c1^4 + 16*c2^4*c0 - 4*c2^3*c1^2
 * ----------------------------------------------------------------------- */
static inline zw_t disc_monic(zw_t c0, zw_t c1, zw_t c2)
{
    zw_t c0sq = zw_mul(c0,c0), c0cu = zw_mul(c0sq,c0);
    zw_t c1sq = zw_mul(c1,c1), c1qu = zw_mul(c1sq,c1sq);
    zw_t c2sq = zw_mul(c2,c2), c2cu = zw_mul(c2sq,c2), c2qu = zw_mul(c2cu,c2);
    zw_t d = {0,0};
    d = zw_add(d, zw_smul( 256, c0cu));
    d = zw_sub(d, zw_smul( 128, zw_mul(c2sq, c0sq)));
    d = zw_add(d, zw_smul( 144, zw_mul(c2, zw_mul(c1sq, c0))));
    d = zw_sub(d, zw_smul(  27, c1qu));
    d = zw_add(d, zw_smul(  16, zw_mul(c2qu, c0)));
    d = zw_sub(d, zw_smul(   4, zw_mul(c2cu, c1sq)));
    return d;
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

    fprintf(stdout, "# Picard curves over Q(zeta_3):  y^3 = f(x),  f(x) = x^4 + (a2+b2*w)*x^2 + (a1+b1*w)*x + (a0+b0*w),  w = e^(2*pi*i/3)\n");
    fprintf(stdout, "# Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2]  where  N(disc(f)) = disc_re^2 - disc_re*disc_im + disc_im^2\n");
    fflush(stdout);

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

            zw_t c0={a0,b0}, c1={a1,b1}, c2={a2,b2};
            zw_t disc = disc_monic(c0, c1, c2);
            if (!disc.re && !disc.im) continue;   /* degenerate (singular) curve */

            /* N(u + v*omega) = u^2 - u*v + v^2  (always >= 0) */
            i128 N = zw_norm(disc);
            if (N <= 0 || N > (i128)Nbnd) continue;

            char nbuf[64], cbuf[128];
            sprintf(cbuf, "[%ld,%ld,%ld,%ld,%ld,%ld]", a0, b0, a1, b1, a2, b2);

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
