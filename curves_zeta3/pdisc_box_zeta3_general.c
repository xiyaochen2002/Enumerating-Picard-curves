/*
 * pdisc_box_zeta3_general.c
 *
 * Enumerate Picard curves  y^3 = f(x)  over Q(omega),  omega = e^(2*pi*i/3),
 * with GENERAL Eisenstein integer coefficients:
 *
 *   f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,   fi = ai + bi*omega in Z[omega]
 *
 * Search box: |ai|, |bi| <= c.
 * Output every non-singular curve with 0 < N(disc(f)) <= norm_bound.
 *
 * Symmetry reductions used (same as the rational pdisc_box):
 *   (1) y -> -y  maps  f(x) -> -f(x),  so f4 and -f4 give isomorphic curves.
 *       We fix f4 in the "positive half" of Z[omega]:
 *         a4 > 0,  or  a4 == 0 and b4 > 0.
 *   (2) x -> -x  maps  (f4,f3,f2,f1,f0) -> (f4,-f3,f2,-f1,f0),  so f3 and -f3
 *       give isomorphic curves (for fixed f4,f2,f0 with f1 -> -f1).
 *       We fix f3 in the "positive half":
 *         a3 > 0,  or  a3 == 0 and b3 >= 0.
 *       (When f3=0 there is residual f1 <-> -f1 symmetry; not exploited here.)
 *
 * Discriminant formula for f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0  (16 terms):
 *   disc(f) =
 *      256*f4^3*f0^3
 *    - 192*f4^2*f3*f1*f0^2
 *    - 128*f4^2*f2^2*f0^2
 *    + 144*f4^2*f2*f1^2*f0
 *    -  27*f4^2*f1^4
 *    + 144*f4*f3^2*f2*f0^2
 *    -   6*f4*f3^2*f1^2*f0
 *    -  80*f4*f3*f2^2*f1*f0
 *    +  18*f4*f3*f2*f1^3
 *    +  16*f4*f2^4*f0
 *    -   4*f4*f2^3*f1^2
 *    -  27*f3^4*f0^2
 *    +  18*f3^3*f2*f1*f0
 *    -   4*f3^3*f1^3
 *    -   4*f3^2*f2^3*f0
 *    +      f3^2*f2^2*f1^2
 * All arithmetic performed in Z[omega] using __int128_t pairs.
 *
 * Build (from curves_zeta3/ directory):
 *   gcc -O2 -fopenmp pdisc_box_zeta3_general.c -o pdisc_box_zeta3_general -lm
 *
 * Usage:
 *   ./pdisc_box_zeta3_general coeff-bound norm-bound [threads]
 *
 * Output lines:  N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]
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

/* Z[omega] element: re + im*omega */
typedef struct { i128 re, im; } zw_t;

/* omega^2 = -omega - 1, so:
 * (a + b*w)(c + d*w) = (ac - bd) + (ad + bc - bd)*w  */
static inline zw_t zw_mul(zw_t x, zw_t y)
{
    return (zw_t){ x.re*y.re - x.im*y.im,
                   x.re*y.im + x.im*y.re - x.im*y.im };
}
static inline zw_t zw_add(zw_t x, zw_t y) { return (zw_t){ x.re+y.re, x.im+y.im }; }
static inline zw_t zw_sub(zw_t x, zw_t y) { return (zw_t){ x.re-y.re, x.im-y.im }; }
static inline zw_t zw_smul(i128 s, zw_t x) { return (zw_t){ s*x.re, s*x.im }; }

/* N(a + b*omega) = a^2 - a*b + b^2 */
static inline i128 zw_norm(zw_t x) { return x.re*x.re - x.re*x.im + x.im*x.im; }


/* Discriminant of f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0 in Z[omega] */
static inline zw_t quartic_disc(zw_t f0, zw_t f1, zw_t f2, zw_t f3, zw_t f4)
{
    zw_t f0sq = zw_mul(f0,f0), f0cu = zw_mul(f0sq,f0);
    zw_t f1sq = zw_mul(f1,f1), f1cu = zw_mul(f1sq,f1), f1qu = zw_mul(f1cu,f1);
    zw_t f2sq = zw_mul(f2,f2), f2cu = zw_mul(f2sq,f2), f2qu = zw_mul(f2cu,f2);
    zw_t f3sq = zw_mul(f3,f3), f3cu = zw_mul(f3sq,f3), f3qu = zw_mul(f3cu,f3);
    zw_t f4sq = zw_mul(f4,f4), f4cu = zw_mul(f4sq,f4);

    zw_t d = {0,0};
    d = zw_add(d, zw_smul( 256, zw_mul(f4cu, f0cu)));
    d = zw_sub(d, zw_smul( 192, zw_mul(zw_mul(f4sq, zw_mul(f3,f1)), f0sq)));
    d = zw_sub(d, zw_smul( 128, zw_mul(zw_mul(f4sq, f2sq), f0sq)));
    d = zw_add(d, zw_smul( 144, zw_mul(zw_mul(f4sq, zw_mul(f2,f1sq)), f0)));
    d = zw_sub(d, zw_smul(  27, zw_mul(f4sq, f1qu)));
    d = zw_add(d, zw_smul( 144, zw_mul(zw_mul(f4, zw_mul(f3sq,f2)), f0sq)));
    d = zw_sub(d, zw_smul(   6, zw_mul(zw_mul(f4, zw_mul(f3sq,f1sq)), f0)));
    d = zw_sub(d, zw_smul(  80, zw_mul(zw_mul(f4, zw_mul(f3, zw_mul(f2sq,f1))), f0)));
    d = zw_add(d, zw_smul(  18, zw_mul(f4, zw_mul(f3, zw_mul(f2,f1cu)))));
    d = zw_add(d, zw_smul(  16, zw_mul(zw_mul(f4, f2qu), f0)));
    d = zw_sub(d, zw_smul(   4, zw_mul(f4, zw_mul(f2cu, f1sq))));
    d = zw_sub(d, zw_smul(  27, zw_mul(f3qu, f0sq)));
    d = zw_add(d, zw_smul(  18, zw_mul(zw_mul(f3cu, zw_mul(f2,f1)), f0)));
    d = zw_sub(d, zw_smul(   4, zw_mul(f3cu, f1cu)));
    d = zw_sub(d, zw_smul(   4, zw_mul(zw_mul(f3sq, f2cu), f0)));
    d = zw_add(d,              zw_mul(zw_mul(f3sq, f2sq), f1sq));
    return d;
}


int main(int argc, char *argv[])
{
    if (argc < 3) {
        puts("Usage: pdisc_box_zeta3_general coeff-bound norm-bound [threads]");
        puts("");
        puts("  Enumerates Picard curves over Q(zeta_3):  y^3 = f(x)");
        puts("  f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,  fi = ai+bi*w in Z[w]");
        puts("  with |ai|,|bi| <= coeff-bound and 0 < N(disc(f)) <= norm-bound.");
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
    snprintf(results_path, sizeof(results_path), "%s/results_zeta3_general.txt", exe_dir);
    snprintf(log_path,     sizeof(log_path),     "%s/log_zeta3_general.txt",     exe_dir);

    FILE *out_fp = fopen(results_path, "w");
    if (!out_fp) { perror(results_path); return 1; }
    FILE *log_fp = fopen(log_path, "w");
    if (!log_fp) { perror(log_path); return 1; }

    fprintf(stderr, "Writing results to %s\n", results_path);
    fprintf(stderr, "Writing log     to %s\n", log_path);

    fprintf(log_fp, "Using %d threads\n", threads);

    /* Rough curve count after symmetry reductions: ~(2c+1)^8 * (c+1)*(2c+1) / 2 */
    {
        double s = 2*c+1;
        fprintf(log_fp,
            "Scanning ~%.0f curves with |ai|,|bi| <= %ld, N(disc) <= %ld...\n",
            0.5*(c+1)*s*s*s*s*s*s*s*s*s, c, Nbnd);
        fflush(log_fp);
    }

    fprintf(out_fp, "# Picard curves over Q(zeta_3):  y^3 = f(x)\n");
    fprintf(out_fp, "# f(x) = f4*x^4 + f3*x^3 + f2*x^2 + f1*x + f0,  fi = ai+bi*w,  w = e^(2*pi*i/3)\n");
    fprintf(out_fp, "# Output format:  N(disc(f)):[a0,b0,a1,b1,a2,b2,a3,b3,a4,b4]\n");
    fflush(out_fp);

    long cnt = 0;
    double t0 = get_time();

    /*
     * Outer two loops (a4, b4) are parallelised.
     * Symmetry reductions:
     *   f4 positive half  (y->-y):  a4 > 0,  or  a4==0 and b4 > 0
     *   f3 positive half  (x->-x):  a3 > 0,  or  a3==0 and b3 >= 0
     */
    #pragma omp parallel for num_threads(threads) schedule(dynamic) \
        reduction(+:cnt) collapse(2)
    for (long a4 = 0; a4 <= c; a4++) {
    for (long b4 = -c; b4 <= c; b4++) {

        /* f4 positive-half normalisation */
        if (a4 == 0 && b4 <= 0) continue;

        zw_t f4 = { a4, b4 };

        /* f3 positive-half normalisation (x -> -x symmetry) */
        for (long a3 = 0; a3 <= c; a3++) {
        for (long b3 = (a3 == 0 ? 0 : -c); b3 <= c; b3++) {

            zw_t f3 = { a3, b3 };

        for (long a2 = -c; a2 <= c; a2++) {
        for (long b2 = -c; b2 <= c; b2++) {
            zw_t f2 = { a2, b2 };

        for (long a1 = -c; a1 <= c; a1++) {
        for (long b1 = -c; b1 <= c; b1++) {
            zw_t f1 = { a1, b1 };

        for (long a0 = -c; a0 <= c; a0++) {
        for (long b0 = -c; b0 <= c; b0++) {
            zw_t f0 = { a0, b0 };

            zw_t disc = quartic_disc(f0, f1, f2, f3, f4);
            if (!disc.re && !disc.im) continue;    /* singular curve */

            i128 N = zw_norm(disc);
            if (N <= 0 || N > (i128)Nbnd) continue;

            char nbuf[64], cbuf[256];
            sprintf(cbuf, "[%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld,%ld]",
                    a0, b0, a1, b1, a2, b2, a3, b3, a4, b4);

            #pragma omp critical(output)
            {
                fprintf(out_fp, "%s:%s\n", itoa128(nbuf, N), cbuf);
                fflush(out_fp);
                cnt++;
            }
        }}
        }}
        }}
        }}
    }}

    fprintf(log_fp, "Found %ld curves in %.3f secs using %d threads\n",
            cnt, get_time()-t0, threads);
    fclose(out_fp);
    fclose(log_fp);
    return 0;
}
