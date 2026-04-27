/*
 * pdisc_box_zeta3_general_fast.c
 *
 * Fast enumeration of Picard curves  y^3 = f(x)  over Q(omega),
 * omega = e^(2*pi*i/3), with general Eisenstein integer coefficients:
 *
 *   f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,   fi = ai+bi*omega in Z[omega]
 *
 * Uses two speedups over pdisc_box_zeta3_general:
 *
 *   1. Monomial tree: partial sums of the discriminant are cached at the
 *      loop level of their deepest variable, avoiding redundant recomputation.
 *
 *   2. Finite differences in the innermost a0 loop: disc(a0) is cubic in
 *      a0 (for fixed b0 and outer variables), so each a0 step costs only
 *      3 zw_t additions instead of ~16 zw_t multiplications.
 *
 * Key decomposition (regroup 16-term quartic disc by power of f0):
 *
 *   disc(f) = P3*f0^3 + P2*f0^2 + P1*f0 + P0
 *
 *   P3 = 256*f4^3
 *
 *   P2 = -27*f3^4
 *      + (-128*f4^2*f2^2 + 144*f4*f3^2*f2)
 *      + (-192*f4^2*f3) * f1
 *
 *   P1 = (16*f4*f2^4 - 4*f3^2*f2^3)
 *      + (144*f4^2*f2 - 6*f4*f3^2) * f1^2
 *      + (-80*f4*f3*f2^2 + 18*f3^3*f2) * f1
 *
 *   P0 = -27*f4^2 * f1^4
 *      + (18*f4*f3*f2 - 4*f3^3) * f1^3
 *      + (-4*f4*f2^3 + f3^2*f2^2) * f1^2
 *
 * Build:
 *   gcc -O2 -fopenmp pdisc_box_zeta3_general_fast.c -o pdisc_box_zeta3_general_fast -lm
 *
 * Usage:
 *   ./pdisc_box_zeta3_general_fast coeff-bound norm-bound [threads]
 *
 * Output (same format as pdisc_box_zeta3_general):
 *   N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <libgen.h>
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

/* Evaluate P3*f0^3 + P2*f0^2 + P1*f0 + P0 at a given f0 */
static inline zw_t eval_cubic_at(zw_t P3, zw_t P2, zw_t P1, zw_t P0, zw_t f0)
{
    zw_t f0sq = zw_mul(f0, f0);
    zw_t f0cu = zw_mul(f0sq, f0);
    return zw_add(zw_add(zw_add(zw_mul(P3, f0cu),
                                 zw_mul(P2, f0sq)),
                          zw_mul(P1, f0)),
                   P0);
}


int main(int argc, char *argv[])
{
    if (argc < 3) {
        puts("Usage: pdisc_box_zeta3_general_fast coeff-bound norm-bound [threads]");
        puts("");
        puts("  Enumerates Picard curves over Q(zeta_3):  y^3 = f(x)");
        puts("  f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,  fi = ai+bi*w in Z[w]");
        puts("  with |ai|,|bi| <= coeff-bound and 0 < N(disc(f)) <= norm-bound.");
        puts("  Uses monomial tree + finite differences for speed.");
        puts("");
        puts("  Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]");
        return 0;
    }

    long c    = atol(argv[1]); assert(c > 0);
    long Nbnd = atol(argv[2]); assert(Nbnd >= 0);

    int threads = argc > 3 ? atoi(argv[3]) : 0;
    if (!threads) threads = omp_get_max_threads();

    /* Open output files in the same directory as the executable */
    char exe_path[4096], dir_buf[4096];
    ssize_t exe_len = readlink("/proc/self/exe", exe_path, sizeof(exe_path)-1);
    if (exe_len < 0) { perror("readlink"); return 1; }
    exe_path[exe_len] = '\0';
    strncpy(dir_buf, exe_path, sizeof(dir_buf));
    char *exe_dir = dirname(dir_buf);

    char results_path[4096], log_path[4096];
    snprintf(results_path, sizeof(results_path),
             "%s/results_zeta3_general_fast.txt", exe_dir);
    snprintf(log_path, sizeof(log_path),
             "%s/log_zeta3_general_fast.txt", exe_dir);

    FILE *out_fp = fopen(results_path, "w");
    if (!out_fp) { perror(results_path); return 1; }
    FILE *log_fp = fopen(log_path, "w");
    if (!log_fp) { perror(log_path); return 1; }

    fprintf(stderr, "Writing results to %s\n", results_path);
    fprintf(stderr, "Writing log     to %s\n", log_path);
    fprintf(log_fp, "Using %d threads\n", threads);
    fprintf(stderr, "Using %d threads\n", threads);

    {
        double s = 2*c+1;
        fprintf(log_fp,
            "Scanning ~%.0f curves with |ai|,|bi| <= %ld, N(disc) <= %ld...\n",
            0.5*(c+1)*s*s*s*s*s*s*s*s*s, c, Nbnd);
        fflush(log_fp);
    }

    fprintf(out_fp, "# Picard curves over Q(zeta_3):  y^3 = f(x)\n");
    fprintf(out_fp, "# f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,"
                    "  fi = ai+bi*w,  w = e^(2*pi*i/3)\n");
    fprintf(out_fp, "# Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]\n");
    fflush(out_fp);

    long cnt = 0;
    double t0 = get_time();

    /*
     * Symmetry reductions (same as pdisc_box_zeta3_general):
     *   f4 positive half (y->-y): a4 > 0, or a4==0 and b4 > 0
     *   f3 positive half (x->-x): a3 > 0, or a3==0 and b3 >= 0
     *
     * Loop order: (a4,b4) parallelised; inner loops sequential per thread.
     * Monomial tree: cache partial sums at the level of the deepest variable.
     */
    #pragma omp parallel for num_threads(threads) schedule(dynamic) \
        reduction(+:cnt) collapse(2)
    for (long a4 = 0; a4 <= c; a4++) {
    for (long b4 = -c; b4 <= c; b4++) {

        if (a4 == 0 && b4 <= 0) continue;   /* f4 positive-half normalisation */

        /* ---- cached per (a4, b4) ---- */
        zw_t f4    = {a4, b4};
        zw_t f4sq  = zw_mul(f4, f4);
        zw_t f4cu  = zw_mul(f4sq, f4);
        zw_t P3    = zw_smul(256, f4cu);
        zw_t d3    = zw_smul(6, P3);       /* constant 3rd finite difference */
        zw_t P0_c4 = zw_smul(-27, f4sq);   /* coeff of f1^4 in P0 */

        for (long a3 = 0; a3 <= c; a3++) {
        for (long b3 = (a3 == 0 ? 0 : -c); b3 <= c; b3++) {

            /* ---- cached per (a3, b3) ---- */
            zw_t f3    = {a3, b3};
            zw_t f3sq  = zw_mul(f3, f3);
            zw_t f3cu  = zw_mul(f3sq, f3);
            zw_t f3qu  = zw_mul(f3cu, f3);

            zw_t P2_A      = zw_smul(-27, f3qu);                   /* -27*f3^4 */
            zw_t tP2_f1    = zw_smul(-192, zw_mul(f4sq, f3));      /* mult by f1 -> -192*f4^2*f3*f1 */
            zw_t tP2_f2    = zw_smul(144,  zw_mul(f4, f3sq));      /* mult by f2 -> 144*f4*f3^2*f2 */
            zw_t tP1_c2_b  = zw_smul(-6,   zw_mul(f4, f3sq));      /* -6*f4*f3^2: part of P1 f1^2 coeff */
            zw_t tP1_c1_a  = zw_smul(-80,  zw_mul(f4, f3));        /* mult by f2^2 -> -80*f4*f3*f2^2 */
            zw_t tP1_c1_b  = zw_smul(18,   f3cu);                  /* mult by f2 -> 18*f3^3*f2 */
            zw_t tP0_c3_a  = zw_smul(18,   zw_mul(f4, f3));        /* mult by f2 -> 18*f4*f3*f2 */
            zw_t P0_c3_b   = zw_smul(-4,   f3cu);                  /* -4*f3^3: part of P0 f1^3 coeff */

        for (long a2 = -c; a2 <= c; a2++) {
        for (long b2 = -c; b2 <= c; b2++) {

            /* ---- cached per (a2, b2) ---- */
            zw_t f2   = {a2, b2};
            zw_t f2sq = zw_mul(f2, f2);
            zw_t f2cu = zw_mul(f2sq, f2);
            zw_t f2qu = zw_mul(f2cu, f2);

            /* P2 partial (all terms except -192*f4^2*f3*f1):
             *   -27*f3^4 - 128*f4^2*f2^2 + 144*f4*f3^2*f2  */
            zw_t P2_AB = zw_add(P2_A,
                          zw_add(zw_smul(-128, zw_mul(f4sq, f2sq)),
                                 zw_mul(tP2_f2, f2)));

            /* P1 constant part (no f1): 16*f4*f2^4 - 4*f3^2*f2^3 */
            zw_t P1_0  = zw_add(zw_smul(16, zw_mul(f4, f2qu)),
                                 zw_smul(-4, zw_mul(f3sq, f2cu)));

            /* P1 coefficient of f1^2: 144*f4^2*f2 - 6*f4*f3^2 */
            zw_t P1_c2 = zw_add(zw_smul(144, zw_mul(f4sq, f2)), tP1_c2_b);

            /* P1 coefficient of f1: -80*f4*f3*f2^2 + 18*f3^3*f2 */
            zw_t P1_c1 = zw_add(zw_mul(tP1_c1_a, f2sq), zw_mul(tP1_c1_b, f2));

            /* P0 coefficient of f1^3: 18*f4*f3*f2 - 4*f3^3 */
            zw_t P0_c3 = zw_add(zw_mul(tP0_c3_a, f2), P0_c3_b);

            /* P0 coefficient of f1^2: -4*f4*f2^3 + f3^2*f2^2 */
            zw_t P0_c2 = zw_add(zw_smul(-4, zw_mul(f4, f2cu)), zw_mul(f3sq, f2sq));

        for (long a1 = -c; a1 <= c; a1++) {
        for (long b1 = -c; b1 <= c; b1++) {

            /* ---- assembled per (a1, b1) ---- */
            zw_t f1   = {a1, b1};
            zw_t f1sq = zw_mul(f1, f1);
            zw_t f1cu = zw_mul(f1sq, f1);
            zw_t f1qu = zw_mul(f1cu, f1);

            /* P2 = P2_AB + (-192*f4^2*f3) * f1 */
            zw_t P2 = zw_add(P2_AB, zw_mul(tP2_f1, f1));

            /* P1 = P1_0 + P1_c2*f1^2 + P1_c1*f1 */
            zw_t P1 = zw_add(P1_0, zw_add(zw_mul(P1_c2, f1sq),
                                           zw_mul(P1_c1, f1)));

            /* P0 = -27*f4^2*f1^4 + P0_c3*f1^3 + P0_c2*f1^2 */
            zw_t P0 = zw_add(zw_mul(P0_c4, f1qu),
                      zw_add(zw_mul(P0_c3, f1cu),
                             zw_mul(P0_c2, f1sq)));

        for (long b0 = -c; b0 <= c; b0++) {

            /* Set up finite differences for the a0 loop.
             * disc(a0) for fixed b0 is cubic in a0.
             * Evaluate at a0 = -c, -c+1, -c+2 to get d0, d1, d2.
             * d3 = 6*P3 is the constant 3rd difference (pre-computed per f4).
             * Each advance step: d0+=d1; d1+=d2; d2+=d3. */
            zw_t v0 = eval_cubic_at(P3, P2, P1, P0, (zw_t){-c,   b0});
            zw_t v1 = eval_cubic_at(P3, P2, P1, P0, (zw_t){-c+1, b0});
            zw_t v2 = eval_cubic_at(P3, P2, P1, P0, (zw_t){-c+2, b0});

            zw_t d0 = v0;
            zw_t d1 = zw_sub(v1, v0);
            zw_t d2 = zw_add(zw_sub(v2, zw_smul(2, v1)), v0);  /* v2-2*v1+v0 */
            /* d3 already set per (a4,b4) */

            for (long a0 = -c; a0 <= c; a0++) {

                /* d0 = disc at (a0, b0); check if non-singular and within bound */
                if (d0.re || d0.im) {
                    i128 N = zw_norm(d0);
                    if (N > 0 && N <= (i128)Nbnd) {
                        char nbuf[64], cbuf[256];
                        sprintf(cbuf, "[%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld]",
                                a0, b0, a1, b1, a2, b2, a3, b3, a4, b4);
                        #pragma omp critical(output)
                        {
                            fprintf(out_fp, "%s:%s\n", itoa128(nbuf, N), cbuf);
                            fflush(out_fp);
                            cnt++;
                        }
                    }
                }

                /* Advance finite differences to next a0 */
                d0 = zw_add(d0, d1);
                d1 = zw_add(d1, d2);
                d2 = zw_add(d2, d3);

            }   /* a0 */
        }   /* b0 */
        }}  /* a1, b1 */
        }}  /* a2, b2 */
        }}  /* a3, b3 */
    }}  /* a4, b4 */

    fprintf(log_fp, "Found %ld curves in %.3f secs using %d threads\n",
            cnt, get_time()-t0, threads);
    fclose(out_fp);
    fclose(log_fp);
    return 0;
}
