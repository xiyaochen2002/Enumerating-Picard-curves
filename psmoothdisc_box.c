#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <gmp.h>
#include <omp.h>
#include <ctype.h>
#include "mpzpolyutil.h"
#include "polyparse.h"
#include "mpoly128.h"
#include "polyenum128.h"
#include "smooth128.h"
#include "cstd.h"


// Discriminant of Picard curve y^3 = f(x) = a4*x^4 + a3*x^3 + a2*x^2 + a1*x + a0 is 3^9*a4^3*disc(f)^2


#define LOG_CURVES  0   // for debugging purposes, logs all curves enumerated to a smalldisc_log.txt


static inline char *curve_string (char buf[], int f0, int f1, int f2, int f3, int f4)
{
    long f[5];
    
    f[0] = f0; f[1] = f1; f[2] = f2; f[3] = f3; f[4] = f4;
    i_poly_sprint (buf, f, 4);
    return buf;
}


#if LOG_CURVES
FILE *log_fp;

static void log_open (char *filename)
{
    log_fp = fopen(filename, "w");
    if ( ! log_fp ) { printf ("Error creating log file\n"); exit(-1); }
}

static void log_curve (__int128_t disc, int f0, int f1, int f2, int f3, int f4)
{
    char logbuf[256], dbuf[256];
    if ( disc < 0 ) disc = -disc;
    #pragma omp critical(log) 
    { curve_string (logbuf, f0,f1,f2,f3,f4); fprintf (log_fp, "%s:%s\n", itoa128(dbuf,disc), logbuf); }
}

static void log_close (void) { fclose (log_fp); log_fp = 0; }

#else
static void log_open (char *filename) {}
static void log_curve (__int128_t disc, int f0, int f1, int f2, int f3, int f4) {}
static void log_close (void) {}
#endif


static inline void process_curve (mpz_t D, long Dbnd, mpz_t f[5], int f0, int f1, int f2, int f3, int f4, long *cnt, mpz_t w[4*4+2], int tid, long ctr, __int128_t disc)
{
    mpz_set_si (f[4],f4); mpz_set_si (f[3], f3); mpz_set_si (f[2], f2); mpz_set_si (f[1], f1); mpz_set_si (f[0], f0);
    mpz_poly_discriminant (D, f, 4, w);
    if ( !f4 || ! mpz_sgn(D) ) return;
    mpz_abs (D, D);
    if ( ! mpz_is_seven_smooth_or_small (D,Dbnd,w) ) { char buf[64]; gmp_fprintf(stderr, "Skipping misidentified smooth/small discriminant %Zd != %s\n", D, itoa128(buf,disc)); }
    mpz_mul (D, D, D);
    mpz_mul_ui (D,D,19683L*(long)f4*(long)f4*(long)f4);
    char buf[256];
    curve_string (buf, f0, f1, f2, f3, f4);
    #pragma omp critical(smalldisc)
    {
        gmp_fprintf (stdout, "%Zd:%s\n", D, buf);  fflush(stdout);
        (*cnt)++;
    }
}

/*
    General Picard curve C of genus g has the form y^3 = f(x) where deg f = 4 and f has integral discriminant
    The discriminant of C is 3^9*f4^3*disc(f)^2

    We search in a box bounded by |f_i| <= b_i
*/

int main (int argc, char *argv[])
{
    double start;
    int threads, instances, instance_id;
    long c, cnt, tot, dbnd;
    
    if ( argc < 3 ) { puts ("psmoothdisc_box coeff-bound disc-bound [threads instances instance-id]"); return 0; }
    c = atoi (argv[1]);
    assert(c);
    dbnd = atol (argv[2]);
    assert(dbnd >= 0);
    if ( argc > 3 ) threads = atoi(argv[3]); else threads = 0;
    if ( ! threads ) threads = omp_get_max_threads();
    fprintf (stderr, "Using %d threads\n", threads);
    if ( argc > 4 ) {
        if ( argc <= 5 ) { puts ("psmoothdisc_box coeff-bound disc-bound [threads instances instance-id]"); return 0; }
        instances = atoi(argv[4]);
        assert (instances > 0 && instances <= 1000000);
        instance_id = atoi(argv[5]);
        assert (instance_id >= 0 && instance_id <= instances);
        if ( instance_id == instances ) instance_id = 0;
        fprintf (stderr,"Instance %d of %d\n", instance_id, instances);
    } else {
        instances = 1; instance_id = 0;
    }

    tot = c*(c+1)*(2*c+1);
    if ( tot < 10L*threads*instances )
        printf ("Warning: 32*%ld*%ld*%ld = %ld < 10*threads*instances=%ld\n", c, c+1, 2*c+1, tot, 10*(long)threads*(long)instances);
    
    seven_smooth_setup_128 ();

    log_open ("smalldisc.log");
    

    tot *= (2*c+1)*(2*c+1);
    fprintf (stderr, "Beginning scan of %ld of %ld (2^%.2f, 10^%.2f) curve equations with |f_i| <= %ld...\n", tot/instances, tot, log2(tot), log10(tot), c);
    tot /= instances;
    
    cnt = 0;
    start = get_time();
    
    #pragma omp parallel num_threads(threads)
    {
        mpoly128_t mpoly;
        mpz_t D, f[5], w[4*4+2];
        int ltid;
        long ctr, m, d, gtid;
        int b = c, b0 = 2 - b;

        if ( b0 > b ) b0 = b;        
        m = instances*threads;
        ltid = omp_get_thread_num();
        gtid = instance_id*threads+ltid;
        mpoly128_discriminant_init (mpoly, 4);
        for ( int i = 0 ; i <= 4 ; i++ ) mpz_init (f[i]);
        for ( int i = 0 ; i < 4*4+2 ; i++ ) mpz_init (w[i]);
        mpz_init (D);
        
        ctr = 0;
        for ( int f4 = 1 ; f4 <= b ; f4++ ) {
        if ( ! is_seven_smooth(f4) ) continue;
        mpoly128_update (mpoly, 4, f4);
        for ( int f3 = 0 ; f3 <= b ; f3++ ) {
        mpoly128_update (mpoly, 3, f3);
        for ( int f2 = -b ; f2 <= b ; f2++, ctr++ ) {
        d = i_mod(ctr,m);
        if ( d != gtid ) continue;
        mpoly128_update (mpoly, 2, f2);
        for ( register int f1 = -b  ; f1 <= b ; f1++ ) {
            __int128_t Dd[4], dpts[3] = {0,0,0};
            register __int128_t *dptr = dpts, disc, *dpoly;
            register int f0;
            mpoly128_update (mpoly, 1, f1);
            dpoly = mpoly128_upoly_ptr(mpoly);
            for ( f0 = -b ; f0 <= b0 ; f0++ ) {
                disc = *dptr++ = cubic_eval (dpoly, f0);
                if ( disc < 0 ) disc = -disc;
                log_curve (disc,f0,f1,f2,f3,f4);
                if ( disc && is_seven_smooth_or_small(disc,dbnd) ) process_curve (D, dbnd, f, f0, f1, f2, f3, f4, &cnt, w, ltid, ctr, disc);
            }
            if ( f0 <= b ) cubic_enum_setup (Dd, dpoly, -b0, dpts, 6);
            for ( ; f0 <= b ; f0++ ) {
                disc = cubic_enum (Dd);
                if ( disc < 0 ) disc = -disc;
                log_curve (disc,f0,f1,f2,f3,f4);
                if ( disc && is_seven_smooth_or_small(disc,dbnd) ) process_curve (D, dbnd, f, f0, f1, f2, f3, f4, &cnt, w, ltid, ctr, disc);
            }
        }
        }}}
    }
    log_close();
    start = (get_time() - start);
    fprintf (stderr, "Found %ld of %ld curves with smooth discriminant in %.3f secs using %d threads (%.1f ns/curve)\n", cnt, tot, start, threads, (1000000000.0*start*threads)/tot);
}
